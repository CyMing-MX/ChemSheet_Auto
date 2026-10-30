#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""EPA CompTox Chemicals Dashboard API 客户端。

输入英文名、CAS 号、SMILES、InChIKey 或 DTXSID，返回 CompTox 上的
理化性质预测值（ACD/Labs Percepta、OPERA、TEST、EPISUITE…）与实验值，
可选返回环境归趋数据（OPERA）。

API 文档（OpenAPI）：https://comptox.epa.gov/dashboard-api/ccdapp2/v3/api-docs
无需 API key，请控制请求频率。

命令行示例：
    python comptox_client.py aspirin
    python comptox_client.py 50-78-2
    python comptox_client.py "CC(=O)OC1=CC=CC=C1C(=O)O"
    python comptox_client.py -f names.txt -o result.csv
    python comptox_client.py aspirin --prefer opera --with-experimental --envfate

作为库使用：
    from comptox_client import ComptoxClient
    client = ComptoxClient()
    rec = client.lookup("aspirin")
    print(rec["properties"]["logkow"]["value"])
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import os
import re
import sys
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

__version__ = "1.0.0"

BASE = "https://comptox.epa.gov/dashboard-api"
SEARCH_URL = BASE + "/ccdapp1/search/chemical/equal-with-detail/{}"
STRUCTURE_URL = BASE + "/ccdapp1/search/chemical/structure/?smiles={}"
PROPERTY_URL = BASE + "/ccdapp2/chemical-property/search/by-dtxsid?id={}"
ENVFATE_URL = BASE + "/ccdapp2/envfate/search/by-dtxsid?id={}"
DETAIL_URL = BASE + "/ccdapp2/chemical-detail/search/by-dsstoxsid?id={}"
PUBCHEM_INCHIKEY = (
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{}/property/InChIKey/JSON"
)

CAS_RE = re.compile(r"^\d{2,7}-\d{2}-\d$")
INCHIKEY_RE = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")
DTXSID_RE = re.compile(r"^DTXSID\d+$", re.I)
SMILES_HINT = re.compile(r"[=\[\]\(\)#@+\\/]")

# CompTox 性质名 -> 输出列名
PROPERTY_KEYS = {
    "LogKow: Octanol-Water": "logkow",
    "LogKoa: Octanol-Air": "logkoa",
    "LogD7.4": "logd_7_4",
    "LogD5.5": "logd_5_5",
    "LogD": "logd",
    "pKa Acidic Apparent": "pka_acidic",
    "pKa Basic Apparent": "pka_basic",
    "pKa": "pka",
    "Boiling Point": "boiling_point",
    "Melting Point": "melting_point",
    "Vapor Pressure": "vapor_pressure",
    "Water Solubility": "water_solubility",
    "Henry's Law": "henrys_law",
    "Density": "density",
    "Flash Point": "flash_point",
    "Molar Refractivity": "molar_refractivity",
    "Molar Volume": "molar_volume",
    "Surface Tension": "surface_tension",
    "Index of Refraction": "index_of_refraction",
    "Polarizability": "polarizability",
}

ENVFATE_KEYS = {
    "Biodeg. Half-Life": "env_biodeg_half_life",
    "Bioconcentration Factor": "env_bcf",
    "Soil Adsorp. Coeff. (Koc)": "env_koc",
    "Atmos. Hydroxylation Rate": "env_atmos_hydroxylation",
    "Fish Biotrans. Half-Life (Km)": "env_fish_biotrans_half_life",
    "ReadyBiodeg": "env_ready_biodeg",
}

# 输出时优先列出的性质顺序
PROPERTY_ORDER = [
    "logkow", "logd_7_4", "logd_5_5", "logd", "logkoa",
    "pka_acidic", "pka_basic", "pka",
    "water_solubility", "vapor_pressure", "henrys_law",
    "boiling_point", "melting_point", "flash_point",
    "density", "molar_refractivity", "molar_volume",
    "surface_tension", "index_of_refraction", "polarizability",
]

IDENTITY_FIELDS = (
    "query", "query_type", "dtxsid", "dtxcid", "preferred_name",
    "casrn", "formula", "smiles", "inchikey",
)

SOURCE_PREFIXES = {
    "percepta": ("percepta",),
    "opera": ("opera",),
    "test": ("test",),
    "episuite": ("episuite", "epi suite"),
}


class ComptoxError(Exception):
    """CompTox API 调用失败。"""


class NotFoundError(ComptoxError):
    """CompTox 中找不到该物质。"""


def snake(name: str) -> str:
    text = re.sub(r"[^0-9a-zA-Z]+", "_", name).strip("_").lower()
    return text or "unknown"


def property_key(name: str) -> str:
    return PROPERTY_KEYS.get(name) or snake(name)


def envfate_key(name: str) -> str:
    return ENVFATE_KEYS.get(name) or ("env_" + snake(name))


def looks_like_smiles(term: str) -> bool:
    return bool(term) and not re.search(r"\s", term) and bool(SMILES_HINT.search(term))


def _num(value: Any) -> Any:
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        return value
    if isinstance(value, str):
        try:
            return float(value) if "." in value or "e" in value.lower() else int(value)
        except ValueError:
            return value
    return value


class ComptoxClient:
    def __init__(
        self,
        *,
        timeout: float = 40.0,
        min_interval: float = 0.6,
        retries: int = 3,
        cache_dir: Optional[str] = None,
        cache_ttl: float = 7 * 24 * 3600,
        use_cache: bool = True,
        pubchem_fallback: bool = True,
    ) -> None:
        self.timeout = timeout
        self.min_interval = min_interval
        self.retries = max(1, retries)
        self.cache_ttl = cache_ttl
        self.use_cache = use_cache
        self.pubchem_fallback = pubchem_fallback
        self._last_request = 0.0
        default_cache = os.environ.get("COMPTOX_CACHE_DIR") or str(
            Path(tempfile.gettempdir()) / "comptox_client"
        )
        self.cache_dir = Path(cache_dir) if cache_dir else Path(default_cache)
        self.stats = {"requests": 0, "cache_hits": 0, "retries": 0, "cache_error": ""}

    # -- 网络层 -------------------------------------------------------------
    def _throttle(self) -> None:
        wait = self.min_interval - (time.monotonic() - self._last_request)
        if wait > 0:
            time.sleep(wait)
        self._last_request = time.monotonic()

    def _cache_path(self, key: str) -> Path:
        return self.cache_dir / (hashlib.sha1(key.encode("utf-8")).hexdigest() + ".json")

    def _cache_get(self, key: str) -> Tuple[bool, Any]:
        if not self.use_cache:
            return False, None
        path = self._cache_path(key)
        try:
            if time.time() - path.stat().st_mtime > self.cache_ttl:
                return False, None
            with path.open("r", encoding="utf-8") as fh:
                self.stats["cache_hits"] += 1
                return True, json.load(fh)
        except (OSError, ValueError):
            return False, None

    def _cache_put(self, key: str, data: Any) -> None:
        if not self.use_cache:
            return
        try:
            self.cache_dir.mkdir(parents=True, exist_ok=True)
            with self._cache_path(key).open("w", encoding="utf-8") as fh:
                json.dump(data, fh, ensure_ascii=False)
        except OSError as exc:
            self.stats["cache_error"] = str(exc)

    def _get_json(self, url: str, *, refresh: bool = False) -> Any:
        if not refresh:
            hit, data = self._cache_get(url)
            if hit:
                return data
        last_error: Optional[Exception] = None
        for attempt in range(1, self.retries + 1):
            self._throttle()
            self.stats["requests"] += 1
            req = urllib.request.Request(
                url,
                headers={
                    "Accept": "application/json",
                    "User-Agent": f"comptox-client/{__version__}",
                },
            )
            try:
                with urllib.request.urlopen(req, timeout=self.timeout) as resp:
                    data = json.loads(resp.read().decode("utf-8"))
                self._cache_put(url, data)
                return data
            except urllib.error.HTTPError as exc:
                text = exc.read().decode("utf-8", "replace").strip()
                if exc.code == 404:
                    raise NotFoundError(f"HTTP 404: {text or url}") from exc
                if exc.code in (429, 500, 502, 503, 504) and attempt < self.retries:
                    self.stats["retries"] += 1
                    time.sleep(2.0 * attempt)
                    last_error = exc
                    continue
                raise ComptoxError(f"HTTP {exc.code}: {text or exc.reason}") from exc
            except (urllib.error.URLError, TimeoutError, ValueError) as exc:
                last_error = exc
                if attempt < self.retries:
                    self.stats["retries"] += 1
                    time.sleep(2.0 * attempt)
                    continue
                raise ComptoxError(f"请求失败：{exc}") from exc
        raise ComptoxError(f"请求失败：{last_error}")

    # -- 解析标识符 ---------------------------------------------------------
    def _pubchem_inchikey(self, smiles: str, refresh: bool = False) -> str:
        url = PUBCHEM_INCHIKEY.format(urllib.parse.quote(smiles, safe=""))
        data = self._get_json(url, refresh=refresh)
        try:
            return data["PropertyTable"]["Properties"][0]["InChIKey"]
        except (KeyError, IndexError, TypeError) as exc:
            raise NotFoundError(f"无法把 SMILES 转成 InChIKey: {smiles}") from exc

    def resolve(self, term: str, *, refresh: bool = False) -> Dict[str, Any]:
        """把 名称/CAS/SMILES/InChIKey/DTXSID 解析成 CompTox 物质信息。"""
        term = term.strip()
        if not term:
            raise ComptoxError("查询内容为空")

        kind = "name"
        lookup_value = term
        if DTXSID_RE.match(term):
            kind, lookup_value = "dtxsid", term.upper()
        elif CAS_RE.match(term):
            kind = "casrn"
        elif INCHIKEY_RE.match(term):
            kind = "inchikey"
        elif looks_like_smiles(term):
            kind = "smiles"
            try:
                lookup_value = self._pubchem_inchikey(term, refresh=refresh)
            except ComptoxError:
                if not self.pubchem_fallback:
                    raise
                raise

        data = self._get_json(SEARCH_URL.format(urllib.parse.quote(lookup_value, safe="")), refresh=refresh)
        if isinstance(data, dict):
            data = data.get("content") or data.get("_embedded") or []
        if not data:
            raise NotFoundError(f"CompTox 中没有找到 {term!r}")

        best = data[0]
        return {
            "query": term,
            "query_type": kind,
            "lookup_value": lookup_value,
            "dtxsid": best.get("dtxsid"),
            "dtxcid": best.get("dtxcid"),
            "preferred_name": best.get("preferredName"),
            "acd_iupac_name": best.get("acdIupacName"),
            "casrn": best.get("casrn"),
            "formula": best.get("molFormula") or best.get("molecularFormula"),
            "smiles": best.get("smiles"),
            "inchikey": best.get("inchiKey") or best.get("inchikey"),
            "_raw": best,
            "_match_count": len(data),
        }

    # -- 数据 ---------------------------------------------------------------
    def properties(self, dtxsid: str, *, refresh: bool = False) -> List[Dict[str, Any]]:
        data = self._get_json(PROPERTY_URL.format(urllib.parse.quote(dtxsid)), refresh=refresh)
        return list((data or {}).get("data") or [])

    def envfate(self, dtxsid: str, *, refresh: bool = False) -> List[Dict[str, Any]]:
        data = self._get_json(ENVFATE_URL.format(urllib.parse.quote(dtxsid)), refresh=refresh)
        return list((data or {}).get("fateDetails") or [])

    @staticmethod
    def _pick(entries: Sequence[Dict[str, Any]], prefer: str) -> Optional[Dict[str, Any]]:
        """按来源优先级挑一条预测值。"""
        if not entries:
            return None
        order = [prefer] + [k for k in ("percepta", "opera", "test", "episuite") if k != prefer]
        for source_key in order:
            for entry in entries:
                src = str(entry.get("source", "")).lower()
                if any(src.startswith(prefix) for prefix in SOURCE_PREFIXES[source_key]):
                    return entry
        return entries[0]

    def lookup(
        self,
        term: str,
        *,
        prefer: str = "percepta",
        with_experimental: bool = False,
        with_envfate: bool = False,
        refresh: bool = False,
    ) -> Dict[str, Any]:
        info = self.resolve(term, refresh=refresh)
        dtxsid = info.get("dtxsid")
        if not dtxsid:
            raise NotFoundError(f"{term!r} 没有 DTXSID")

        properties: Dict[str, Dict[str, Any]] = {}
        for item in self.properties(dtxsid, refresh=refresh):
            name = item.get("name") or ""
            predicted = ((item.get("predicted") or {}).get("rawData")) or []
            chosen = self._pick(predicted, prefer)
            record: Dict[str, Any] = {
                "name": name,
                "unit": item.get("unit"),
                "value": _num(chosen.get("value")) if chosen else None,
                "source": chosen.get("source") if chosen else None,
                "sources": {
                    str(e.get("source")): _num(e.get("value")) for e in predicted if e.get("value") is not None
                },
            }
            if with_experimental:
                experimental = ((item.get("experimental") or {}).get("rawData")) or []
                if experimental:
                    record["experimental"] = _num(experimental[0].get("value"))
                    record["experimental_source"] = experimental[0].get("source")
            properties[property_key(name)] = record

        record = {k: v for k, v in info.items() if not k.startswith("_")}
        record["properties"] = properties
        record["match_count"] = info.get("_match_count", 1)

        if with_envfate:
            env: Dict[str, Dict[str, Any]] = {}
            for item in self.envfate(dtxsid, refresh=refresh):
                predicted = ((item.get("predicted") or {}).get("rawData")) or []
                chosen = self._pick(predicted, "opera")
                env[envfate_key(item.get("name") or "")] = {
                    "name": item.get("name"),
                    "unit": item.get("unit"),
                    "value": _num(chosen.get("value")) if chosen else None,
                    "source": chosen.get("source") if chosen else None,
                }
            record["envfate"] = env
        return record

    def lookup_many(
        self,
        terms: Iterable[str],
        *,
        prefer: str = "percepta",
        with_experimental: bool = False,
        with_envfate: bool = False,
        refresh: bool = False,
        on_error: str = "raise",
    ) -> List[Dict[str, Any]]:
        records: List[Dict[str, Any]] = []
        for term in terms:
            term = term.strip()
            if not term or term.startswith("#"):
                continue
            try:
                records.append(
                    self.lookup(
                        term,
                        prefer=prefer,
                        with_experimental=with_experimental,
                        with_envfate=with_envfate,
                        refresh=refresh,
                    )
                )
            except ComptoxError as exc:
                if on_error == "raise":
                    raise
                records.append({"query": term, "error": str(exc)})
        return records


# ---------------------------------------------------------------------------
# 输出
# ---------------------------------------------------------------------------
def _cell(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.15g}"
    return str(value)


def _column_plan(records: Sequence[Dict[str, Any]], with_experimental: bool, with_envfate: bool):
    prop_keys: List[str] = []
    env_keys: List[str] = []
    for record in records:
        if record.get("error"):
            continue
        for key in (record.get("properties") or {}):
            if key not in prop_keys:
                prop_keys.append(key)
        for key in (record.get("envfate") or {}):
            if key not in env_keys:
                env_keys.append(key)

    def sort_key(key: str) -> Tuple[int, str]:
        return (PROPERTY_ORDER.index(key) if key in PROPERTY_ORDER else 999, key)

    prop_keys.sort(key=sort_key)
    columns = list(IDENTITY_FIELDS)
    for key in prop_keys:
        columns.append(key)
        columns.append(key + "_src")
        if with_experimental:
            columns.append(key + "_exp")
    if with_envfate:
        columns.extend(env_keys)
        columns.append("env_src")
    return columns, prop_keys, env_keys


def record_to_row(record: Dict[str, Any], prop_keys: Sequence[str], env_keys: Sequence[str],
                  with_experimental: bool, with_envfate: bool) -> Dict[str, Any]:
    row: Dict[str, Any] = {key: record.get(key) for key in IDENTITY_FIELDS}
    properties = record.get("properties") or {}
    for key in prop_keys:
        item = properties.get(key) or {}
        row[key] = item.get("value")
        row[key + "_src"] = item.get("source")
        if with_experimental:
            row[key + "_exp"] = item.get("experimental")
    if with_envfate:
        env = record.get("envfate") or {}
        sources = []
        for key in env_keys:
            item = env.get(key) or {}
            row[key] = item.get("value")
            if item.get("source") and item["source"] not in sources:
                sources.append(item["source"])
        row["env_src"] = "; ".join(sources)
    return row


def records_to_csv(records: Sequence[Dict[str, Any]], with_experimental: bool, with_envfate: bool) -> str:
    columns, prop_keys, env_keys = _column_plan(records, with_experimental, with_envfate)
    buffer = io.StringIO()
    writer = csv.DictWriter(
        buffer, fieldnames=columns + ["error"], extrasaction="ignore", lineterminator="\n"
    )
    writer.writeheader()
    for record in records:
        if record.get("error"):
            row = {key: record.get(key) for key in ("query",)}
            row["error"] = record["error"]
        else:
            row = record_to_row(record, prop_keys, env_keys, with_experimental, with_envfate)
        writer.writerow({key: _cell(row.get(key)) for key in writer.fieldnames})
    return buffer.getvalue()


def records_to_table(records: Sequence[Dict[str, Any]]) -> str:
    blocks: List[str] = []
    for record in records:
        if record.get("error"):
            blocks.append(f"[FAILED] {record.get('query')}\n  {record['error']}")
            continue
        head = (
            f"[{record['dtxsid']}] {record.get('preferred_name') or ''}"
            f"   ({record.get('query')})"
        )
        lines = [head]
        identity = [(k, record[k]) for k in ("casrn", "formula", "smiles", "inchikey") if record.get(k)]
        prop_items = [
            (key, item) for key, item in (record.get("properties") or {}).items()
            if item.get("value") is not None or item.get("experimental") is not None
        ]
        env_items = [
            (key, item) for key, item in (record.get("envfate") or {}).items()
            if item.get("value") is not None
        ]
        width = max(
            [len(k) for k, _ in identity] + [len(k) for k, _ in prop_items] + [len(k) for k, _ in env_items],
            default=16,
        ) + 2
        for key, value in identity:
            lines.append(f"  {key:<{width}}{value}")
        for key, item in prop_items:
            unit = item.get("unit") or ""
            text = f"{_cell(item.get('value'))} {unit}  [{item.get('source')}]"
            if item.get("experimental") is not None:
                text += f"   exp={_cell(item['experimental'])}"
            lines.append(f"  {key:<{width}}{text}")
        for key, item in env_items:
            lines.append(f"  {key:<{width}}{_cell(item['value'])} {item.get('unit') or ''}  [{item.get('source')}]")
        blocks.append("\n".join(lines))
    return "\n\n".join(blocks)


def records_to_long_csv(records: Sequence[Dict[str, Any]]) -> str:
    buffer = io.StringIO()
    fields = ["query", "dtxsid", "property", "unit", "type", "source", "value"]
    writer = csv.DictWriter(buffer, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    for record in records:
        if record.get("error"):
            continue
        base = {"query": record.get("query"), "dtxsid": record.get("dtxsid")}
        for key, item in (record.get("properties") or {}).items():
            for source, value in (item.get("sources") or {}).items():
                writer.writerow({**base, "property": key, "unit": item.get("unit"),
                                 "type": "predicted", "source": source, "value": _cell(value)})
            if item.get("experimental") is not None:
                writer.writerow({**base, "property": key, "unit": item.get("unit"),
                                 "type": "experimental",
                                 "source": item.get("experimental_source"),
                                 "value": _cell(item["experimental"])})
        for key, item in (record.get("envfate") or {}).items():
            if item.get("value") is not None:
                writer.writerow({**base, "property": key, "unit": item.get("unit"),
                                 "type": "predicted", "source": item.get("source"),
                                 "value": _cell(item["value"])})
    return buffer.getvalue()


# ---------------------------------------------------------------------------
# 命令行
# ---------------------------------------------------------------------------
def _read_terms(args: argparse.Namespace) -> List[str]:
    terms: List[str] = list(args.terms or [])
    if args.file:
        text = Path(args.file).read_text(encoding="utf-8-sig")
        terms.extend(line.strip() for line in text.splitlines())
    if not terms and not sys.stdin.isatty():
        terms.extend(line.strip() for line in sys.stdin.read().splitlines())
    return [t for t in terms if t and not t.startswith("#")]


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        prog="comptox_client",
        description="用名称/CAS/SMILES 查 EPA CompTox，返回 Percepta、OPERA 等预测性质。",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "示例：\n"
            "  python comptox_client.py aspirin\n"
            "  python comptox_client.py 50-78-2\n"
            "  python comptox_client.py -f names.txt -o result.csv\n"
            "  python comptox_client.py aspirin --prefer opera --with-experimental\n"
            "  python comptox_client.py aspirin --envfate --long -o long.csv\n"
        ),
    )
    parser.add_argument("terms", nargs="*", help="英文名 / CAS / SMILES / InChIKey / DTXSID")
    parser.add_argument("-f", "--file", help="批量文件，每行一个")
    parser.add_argument("-o", "--output", help="输出文件（.csv/.json 自动识别）")
    parser.add_argument("--format", choices=["auto", "table", "json", "csv"], default="auto")
    parser.add_argument("--prefer", choices=["percepta", "opera", "test", "episuite"], default="percepta",
                        help="主值优先来源（默认 percepta，缺则依次回退）")
    parser.add_argument("--with-experimental", action="store_true", help="额外输出实验值列")
    parser.add_argument("--envfate", action="store_true", help="额外输出环境归趋（OPERA）")
    parser.add_argument("--long", action="store_true", help="长表：每个来源一行")
    parser.add_argument("--refresh", action="store_true", help="忽略缓存")
    parser.add_argument("--cache-dir", help="缓存目录")
    parser.add_argument("--no-cache", action="store_true")
    parser.add_argument("--timeout", type=float, default=40.0)
    parser.add_argument("--ignore-errors", action="store_true")
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("--version", action="version", version=f"comptox_client {__version__}")
    args = parser.parse_args(argv)

    client = ComptoxClient(
        timeout=args.timeout,
        cache_dir=args.cache_dir,
        use_cache=not args.no_cache,
    )

    terms = _read_terms(args)
    if not terms:
        parser.print_help()
        return 2

    try:
        records = client.lookup_many(
            terms,
            prefer=args.prefer,
            with_experimental=args.with_experimental,
            with_envfate=args.envfate,
            refresh=args.refresh,
            on_error="skip" if args.ignore_errors else "raise",
        )
    except ComptoxError as exc:
        print(f"错误：{exc}", file=sys.stderr)
        return 1

    if args.verbose:
        print(
            f"请求 {client.stats['requests']} 次，缓存命中 {client.stats['cache_hits']} 次，"
            f"重试 {client.stats['retries']} 次",
            file=sys.stderr,
        )

    fmt = args.format
    if fmt == "auto":
        suffix = Path(args.output).suffix.lower() if args.output else ""
        fmt = "json" if suffix == ".json" else "csv" if suffix == ".csv" else "table"
    if args.long and fmt != "json":
        fmt = "csv"

    if fmt == "json":
        payload = records[0] if len(records) == 1 else records
        text = json.dumps(payload, ensure_ascii=False, indent=2)
    elif fmt == "csv":
        text = (records_to_long_csv(records) if args.long
                else records_to_csv(records, args.with_experimental, args.envfate)).rstrip("\n")
    else:
        text = records_to_table(records)

    if args.output:
        encoding = "utf-8-sig" if fmt == "csv" else "utf-8"
        Path(args.output).write_text(text + "\n", encoding=encoding)
        print(f"已写入 {args.output}（{len(records)} 条）", file=sys.stderr)
    else:
        try:
            sys.stdout.reconfigure(encoding="utf-8")
        except (AttributeError, ValueError):
            pass
        print(text)

    failures = [r for r in records if r.get("error")]
    if failures:
        print(f"共 {len(failures)} 项失败", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
