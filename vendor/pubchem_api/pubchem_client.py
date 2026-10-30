#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""PubChem PUG REST / PUG View client.

输入英文名、CAS 号、CID、InChIKey 或 SMILES，返回 SMILES、InChI 以及
PubChem 化合物页面上 "Computed Properties" 的全部计算结果。

只依赖 Python 标准库（3.8+），无需 API key，无需安装 requests。

命令行示例:
    python pubchem_client.py aspirin
    python pubchem_client.py 50-78-2 --format json
    python pubchem_client.py -f names.txt -o result.csv
    python pubchem_client.py --serve --port 8765

作为库使用:
    from pubchem_client import PubChemClient
    client = PubChemClient()
    rec = client.lookup("aspirin")[0]
    print(rec["smiles"], rec["xlogp"])
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
from typing import Any, Dict, Iterable, Iterator, List, Optional, Sequence

__version__ = "1.0.0"

PUG_REST = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
PUG_VIEW = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view"
COMPOUND_URL = "https://pubchem.ncbi.nlm.nih.gov/compound/{}"

# ---------------------------------------------------------------------------
# PubChem 属性名（已逐个验证，2026-10 有效）
# ---------------------------------------------------------------------------
DEFAULT_PROPERTIES: tuple[str, ...] = (
    "Title",
    "MolecularFormula",
    "MolecularWeight",
    "SMILES",              # 异构 SMILES（含立体化学）
    "ConnectivitySMILES",  # 连接性 SMILES（不含立体化学），即 Canonical SMILES
    "InChI",
    "InChIKey",
    "IUPACName",
    "XLogP",
    "ExactMass",
    "MonoisotopicMass",
    "TPSA",
    "Complexity",
    "Charge",
    "HBondDonorCount",
    "HBondAcceptorCount",
    "RotatableBondCount",
    "HeavyAtomCount",
    "IsotopeAtomCount",
    "DefinedAtomStereoCount",
    "UndefinedAtomStereoCount",
    "DefinedBondStereoCount",
    "UndefinedBondStereoCount",
    "AtomStereoCount",
    "BondStereoCount",
    "CovalentUnitCount",
)

# PUG REST 属性名 -> 返回字典里的 key
PROPERTY_KEY_MAP: Dict[str, str] = {
    "Title": "title",
    "MolecularFormula": "molecular_formula",
    "MolecularWeight": "molecular_weight",
    "SMILES": "smiles",
    "ConnectivitySMILES": "canonical_smiles",
    "InChI": "inchi",
    "InChIKey": "inchikey",
    "IUPACName": "iupac_name",
    "XLogP": "xlogp",
    "ExactMass": "exact_mass",
    "MonoisotopicMass": "monoisotopic_mass",
    "TPSA": "tpsa",
    "Complexity": "complexity",
    "Charge": "formal_charge",
    "HBondDonorCount": "h_bond_donor_count",
    "HBondAcceptorCount": "h_bond_acceptor_count",
    "RotatableBondCount": "rotatable_bond_count",
    "HeavyAtomCount": "heavy_atom_count",
    "IsotopeAtomCount": "isotope_atom_count",
    "DefinedAtomStereoCount": "defined_atom_stereocenter_count",
    "UndefinedAtomStereoCount": "undefined_atom_stereocenter_count",
    "DefinedBondStereoCount": "defined_bond_stereocenter_count",
    "UndefinedBondStereoCount": "undefined_bond_stereocenter_count",
    "AtomStereoCount": "atom_stereocenter_count",
    "BondStereoCount": "bond_stereocenter_count",
    "CovalentUnitCount": "covalently_bonded_unit_count",
}

# PUG View "Computed Properties" 小节标题 -> 返回字典里的 key
COMPUTED_HEADING_MAP: Dict[str, str] = {
    "Molecular Weight": "molecular_weight",
    "XLogP3": "xlogp",
    "Hydrogen Bond Donor Count": "h_bond_donor_count",
    "Hydrogen Bond Acceptor Count": "h_bond_acceptor_count",
    "Rotatable Bond Count": "rotatable_bond_count",
    "Exact Mass": "exact_mass",
    "Monoisotopic Mass": "monoisotopic_mass",
    "Topological Polar Surface Area": "tpsa",
    "Heavy Atom Count": "heavy_atom_count",
    "Formal Charge": "formal_charge",
    "Complexity": "complexity",
    "Isotope Atom Count": "isotope_atom_count",
    "Defined Atom Stereocenter Count": "defined_atom_stereocenter_count",
    "Undefined Atom Stereocenter Count": "undefined_atom_stereocenter_count",
    "Defined Bond Stereocenter Count": "defined_bond_stereocenter_count",
    "Undefined Bond Stereocenter Count": "undefined_bond_stereocenter_count",
    "Covalently-Bonded Unit Count": "covalently_bonded_unit_count",
    "Compound Is Canonicalized": "compound_is_canonicalized",
}

# 输出列顺序（CSV 表头）
OUTPUT_FIELDS: tuple[str, ...] = (
    "query",
    "query_type",
    "resolved_via",
    "cid",
    "title",
    "molecular_formula",
    "molecular_weight",
    "smiles",
    "canonical_smiles",
    "inchi",
    "inchikey",
    "iupac_name",
    "xlogp",
    "tpsa",
    "exact_mass",
    "monoisotopic_mass",
    "formal_charge",
    "complexity",
    "heavy_atom_count",
    "h_bond_donor_count",
    "h_bond_acceptor_count",
    "rotatable_bond_count",
    "isotope_atom_count",
    "atom_stereocenter_count",
    "defined_atom_stereocenter_count",
    "undefined_atom_stereocenter_count",
    "bond_stereocenter_count",
    "defined_bond_stereocenter_count",
    "undefined_bond_stereocenter_count",
    "covalently_bonded_unit_count",
    "compound_is_canonicalized",
    "url",
)

CAS_RE = re.compile(r"^\d{2,7}-\d{2}-\d$")
INCHIKEY_RE = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")
SMILES_HINT_RE = re.compile(r"[=\[\]\(\)#@+\\]")

# 这些字段 PubChem 有时返回字符串，统一转成数字，方便下游直接计算
NUMERIC_FIELDS = {
    "molecular_weight", "exact_mass", "monoisotopic_mass", "xlogp", "tpsa",
    "complexity", "formal_charge", "heavy_atom_count", "isotope_atom_count",
    "h_bond_donor_count", "h_bond_acceptor_count", "rotatable_bond_count",
    "defined_atom_stereocenter_count", "undefined_atom_stereocenter_count",
    "defined_bond_stereocenter_count", "undefined_bond_stereocenter_count",
    "atom_stereocenter_count", "bond_stereocenter_count", "covalently_bonded_unit_count",
}


def _coerce_number(value: Any) -> Any:
    if not isinstance(value, str):
        return value
    text = value.strip()
    if re.fullmatch(r"-?\d+", text):
        return int(text)
    try:
        return float(text)
    except ValueError:
        return value


# ---------------------------------------------------------------------------
# 异常
# ---------------------------------------------------------------------------
class PubChemError(Exception):
    """PubChem 请求失败。"""


class NotFoundError(PubChemError):
    """PubChem 里找不到该化合物。"""


# ---------------------------------------------------------------------------
# 工具函数
# ---------------------------------------------------------------------------
def is_cas_number(text: str) -> bool:
    """判断是否为格式合法且校验位正确的 CAS 号。"""
    text = text.strip()
    if not CAS_RE.match(text):
        return False
    digits = text.replace("-", "")
    body, check = digits[:-1], int(digits[-1])
    total = sum((i + 1) * int(d) for i, d in enumerate(reversed(body)))
    return total % 10 == check


def classify_identifier(text: str) -> str:
    """粗判查询类型: cid / inchikey / cas / smiles / name。"""
    s = text.strip()
    if re.fullmatch(r"(?i)(cid[:\s]*)?\d{1,9}", s):
        return "cid"
    if INCHIKEY_RE.match(s):
        return "inchikey"
    if is_cas_number(s):
        return "cas"
    if not re.search(r"\s", s) and SMILES_HINT_RE.search(s) and len(s) > 2:
        return "smiles"
    return "name"


def _first_value(info: Dict[str, Any]) -> Any:
    """从 PUG View 的 Information[0].Value 中取出标量值。"""
    value = info.get("Value") or {}
    if "Number" in value and value["Number"]:
        return value["Number"][0]
    if "StringWithMarkup" in value and value["StringWithMarkup"]:
        return "".join(part.get("String", "") for part in value["StringWithMarkup"]).strip()
    if "String" in value:
        return value["String"]
    if "Binary" in value:
        return True
    return None


def _iter_sections(section: Dict[str, Any]) -> Iterator[Dict[str, Any]]:
    yield section
    for child in section.get("Section") or ():
        yield from _iter_sections(child)


class PubChemClient:
    """轻量 PubChem 客户端（标准库实现，内置限速/重试/磁盘缓存）。"""

    def __init__(
        self,
        *,
        timeout: float = 30.0,
        min_interval: float = 0.25,   # PubChem 建议 <= 5 次/秒
        retries: int = 3,
        cache_dir: Optional[str] = None,
        cache_ttl: float = 7 * 24 * 3600,
        use_cache: bool = True,
        user_agent: str = f"pubchem-client/{__version__} (python-urllib)",
    ) -> None:
        self.timeout = timeout
        self.min_interval = min_interval
        self.retries = max(1, retries)
        self.cache_ttl = cache_ttl
        self.use_cache = use_cache
        self.user_agent = user_agent
        self._last_request = 0.0
        default_cache = os.environ.get("PUBCHEM_CACHE_DIR") or str(
            Path(tempfile.gettempdir()) / "pubchem_client"
        )
        self.cache_dir = Path(cache_dir) if cache_dir else Path(default_cache)
        self.stats = {"requests": 0, "cache_hits": 0, "retries": 0, "cache_error": ""}

    # -- 网络层 -------------------------------------------------------------
    def _throttle(self) -> None:
        wait = self.min_interval - (time.monotonic() - self._last_request)
        if wait > 0:
            time.sleep(wait)
        self._last_request = time.monotonic()

    def _cache_path(self, url: str) -> Path:
        digest = hashlib.sha1(url.encode("utf-8")).hexdigest()
        return self.cache_dir / f"{digest}.json"

    def _read_cache(self, url: str, refresh: bool) -> Optional[Dict[str, Any]]:
        if not (self.use_cache and not refresh):
            return None
        path = self._cache_path(url)
        try:
            if time.time() - path.stat().st_mtime > self.cache_ttl:
                return None
            with path.open("r", encoding="utf-8") as fh:
                data = json.load(fh)
            self.stats["cache_hits"] += 1
            return data
        except (OSError, ValueError):
            return None

    def _write_cache(self, url: str, data: Dict[str, Any]) -> None:
        if not self.use_cache:
            return
        try:
            self.cache_dir.mkdir(parents=True, exist_ok=True)
            with self._cache_path(url).open("w", encoding="utf-8") as fh:
                json.dump(data, fh, ensure_ascii=False)
        except OSError as exc:
            self.stats["cache_error"] = str(exc)

    def _request(self, url: str) -> Dict[str, Any]:
        last_error: Optional[Exception] = None
        for attempt in range(1, self.retries + 1):
            self._throttle()
            self.stats["requests"] += 1
            req = urllib.request.Request(
                url, headers={"User-Agent": self.user_agent, "Accept": "application/json"}
            )
            try:
                with urllib.request.urlopen(req, timeout=self.timeout) as resp:
                    return json.loads(resp.read().decode("utf-8"))
            except urllib.error.HTTPError as exc:
                body = exc.read().decode("utf-8", "replace")
                message = body.strip()
                try:
                    fault = json.loads(body).get("Fault", {})
                    message = f"{fault.get('Code', exc.code)}: {fault.get('Message', '')}".strip()
                except ValueError:
                    pass
                if exc.code == 404:
                    raise NotFoundError(message or "not found") from exc
                if exc.code in (429, 500, 502, 503, 504) and attempt < self.retries:
                    self.stats["retries"] += 1
                    time.sleep(0.8 * attempt)
                    last_error = exc
                    continue
                raise PubChemError(f"HTTP {exc.code} {message}") from exc
            except (urllib.error.URLError, TimeoutError, ValueError) as exc:
                last_error = exc
                if attempt < self.retries:
                    self.stats["retries"] += 1
                    time.sleep(0.8 * attempt)
                    continue
                raise PubChemError(f"request failed: {exc}") from exc
        raise PubChemError(f"request failed: {last_error}")

    def _get_json(self, url: str, *, refresh: bool = False) -> Dict[str, Any]:
        cached = self._read_cache(url, refresh)
        if cached is not None:
            return cached
        data = self._request(url)
        self._write_cache(url, data)
        return data

    # -- API ---------------------------------------------------------------
    def resolve(self, identifier: str, *, refresh: bool = False) -> List[Dict[str, Any]]:
        """把任意标识符解析为 CID 列表。

        返回 [{"cid": 2244, "method": "name"}, ...]，method 说明命中的命名空间。
        """
        raw = identifier.strip()
        if not raw:
            raise PubChemError("empty identifier")

        kind = classify_identifier(raw)
        attempts: List[tuple[str, str]] = []

        def add(kind_name: str, value: str) -> None:
            pair = (kind_name, value)
            if pair not in attempts:
                attempts.append(pair)

        if kind == "cid":
            cid_text = re.sub(r"(?i)^cid[:\s]*", "", raw)
            add("cid", cid_text)
            add("name", raw)
        elif kind == "inchikey":
            add("inchikey", raw)
        elif kind == "cas":
            add("name", raw)
            add("cid", raw)
        elif kind == "smiles":
            add("smiles", raw)
            add("name", raw)
        else:
            add("name", raw)

        errors: List[str] = []
        for method, value in attempts:
            quoted = urllib.parse.quote(value, safe="")
            url = f"{PUG_REST}/compound/{method}/{quoted}/cids/JSON"
            try:
                data = self._get_json(url, refresh=refresh)
            except NotFoundError:
                errors.append(f"{method}={value!r}: not found")
                continue
            except PubChemError as exc:
                errors.append(f"{method}={value!r}: {exc}")
                continue
            cids: List[int] = []
            for cid in (data.get("IdentifierList") or {}).get("CID") or []:
                if cid not in cids:
                    cids.append(int(cid))
            if cids:
                return [{"cid": cid, "method": method} for cid in cids]
            errors.append(f"{method}={value!r}: no CID")
        raise NotFoundError(f"无法在 PubChem 中找到 {identifier!r} ({'; '.join(errors)})")

    def properties(
        self,
        cids: Sequence[int],
        names: Sequence[str] = DEFAULT_PROPERTIES,
        *,
        chunk_size: int = 50,
        refresh: bool = False,
    ) -> Dict[int, Dict[str, Any]]:
        """批量取 PUG REST 属性，返回 {cid: {key: value}}。"""
        result: Dict[int, Dict[str, Any]] = {}
        names_param = ",".join(names)
        for start in range(0, len(cids), chunk_size):
            chunk = cids[start:start + chunk_size]
            cid_param = ",".join(str(c) for c in chunk)
            url = f"{PUG_REST}/compound/cid/{cid_param}/property/{names_param}/JSON"
            data = self._get_json(url, refresh=refresh)
            for row in (data.get("PropertyTable") or {}).get("Properties") or []:
                cid = int(row["CID"])
                record: Dict[str, Any] = {}
                for prop, key in PROPERTY_KEY_MAP.items():
                    if prop in row:
                        record[key] = row[prop]
                result[cid] = record
        return result

    def computed_properties(self, cid: int, *, refresh: bool = False) -> Dict[str, Any]:
        """取 PubChem 化合物页 "Computed Properties" 整张表。"""
        heading = urllib.parse.quote_plus("Computed Properties")
        url = f"{PUG_VIEW}/data/compound/{cid}/JSON?heading={heading}"
        data = self._get_json(url, refresh=refresh)
        record = data.get("Record") or {}
        target = None
        for section in _iter_sections(record):
            if section.get("TOCHeading") == "Computed Properties":
                target = section
                break
        if target is None:
            return {}

        output: Dict[str, Any] = {}
        extra: Dict[str, Any] = {}
        for sub in target.get("Section") or ():
            heading_text = sub.get("TOCHeading") or ""
            infos = sub.get("Information") or []
            if not infos:
                continue
            value = _first_value(infos[0])
            if value is None:
                continue
            key = COMPUTED_HEADING_MAP.get(heading_text)
            if key:
                output[key] = value
            else:
                extra[heading_text] = value
        if extra:
            output["other_computed_properties"] = extra
        return output

    def iupac_name_from_view(self, cid: int, *, refresh: bool = False) -> Optional[str]:
        """PUG REST 没给 IUPACName 时，从 PUG View 的 "IUPAC Name" 小节兜底取一个。"""
        heading = urllib.parse.quote_plus("IUPAC Name")
        url = f"{PUG_VIEW}/data/compound/{cid}/JSON?heading={heading}"
        data = self._get_json(url, refresh=refresh)
        target = None
        for section in _iter_sections(data.get("Record") or {}):
            if section.get("TOCHeading") == "IUPAC Name":
                target = section
                break
        if target is None:
            return None
        for info in target.get("Information") or ():
            value = _first_value(info)
            if isinstance(value, str) and value.strip():
                return value.strip()
        return None

    def lookup(
        self,
        identifier: str,
        *,
        include_computed: bool = True,
        refresh: bool = False,
    ) -> List[Dict[str, Any]]:
        """查询单个标识符，返回记录列表（一个名字可能对应多个 CID）。"""
        resolved = self.resolve(identifier, refresh=refresh)
        cids = [item["cid"] for item in resolved]
        props = self.properties(cids, refresh=refresh)

        records: List[Dict[str, Any]] = []
        query_type = classify_identifier(identifier)
        for item in resolved:
            cid, method = item["cid"], item["method"]
            record: Dict[str, Any] = {
                "query": identifier.strip(),
                "query_type": query_type,
                "resolved_via": method,
                "cid": cid,
            }
            record.update(props.get(cid, {}))
            if include_computed:
                try:
                    record.update(self.computed_properties(cid, refresh=refresh))
                except PubChemError as exc:
                    record["computed_properties_error"] = str(exc)
            # IUPAC 名：PUG REST 缺失时用 PUG View 兜底
            if not record.get("iupac_name"):
                try:
                    fallback = self.iupac_name_from_view(cid, refresh=refresh)
                except PubChemError:
                    fallback = None
                if fallback:
                    record["iupac_name"] = fallback
                    record["iupac_name_source"] = "pug_view"
            record["url"] = COMPOUND_URL.format(cid)
            for key in NUMERIC_FIELDS:
                if key in record:
                    record[key] = _coerce_number(record[key])
            records.append(record)
        return records

    def lookup_many(
        self,
        identifiers: Iterable[str],
        *,
        include_computed: bool = True,
        refresh: bool = False,
        on_error: str = "raise",
    ) -> List[Dict[str, Any]]:
        """批量查询；on_error='skip' 时跳过失败项并在结果里保留错误行。"""
        records: List[Dict[str, Any]] = []
        for identifier in identifiers:
            identifier = identifier.strip()
            if not identifier or identifier.startswith("#"):
                continue
            try:
                records.extend(
                    self.lookup(identifier, include_computed=include_computed, refresh=refresh)
                )
            except PubChemError as exc:
                if on_error == "raise":
                    raise
                records.append({"query": identifier, "error": str(exc)})
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


def records_to_csv(records: Sequence[Dict[str, Any]], delimiter: str = ",") -> str:
    buffer = io.StringIO()
    writer = csv.DictWriter(
        buffer,
        fieldnames=list(OUTPUT_FIELDS) + ["error"],
        extrasaction="ignore",
        delimiter=delimiter,
        lineterminator="\n",
    )
    writer.writeheader()
    for record in records:
        writer.writerow({key: _cell(record.get(key)) for key in writer.fieldnames})
    return buffer.getvalue()


def records_to_table(records: Sequence[Dict[str, Any]]) -> str:
    blocks: List[str] = []
    for record in records:
        if record.get("error"):
            blocks.append(f"[FAILED] {record.get('query')}\n  {record['error']}")
            continue
        keys = [
            key for key in OUTPUT_FIELDS
            if key not in ("query", "query_type", "cid", "title") and record.get(key) not in (None, "")
        ]
        if record.get("computed_properties_error"):
            keys.append("computed_properties_error")
        width = max((len(key) for key in keys), default=20) + 2
        lines = [f"[CID {record['cid']}] {record.get('title', '')}   ({record.get('query')})"]
        for key in keys:
            lines.append(f"  {key:<{width}}{_cell(record.get(key))}")
        blocks.append("\n".join(lines))
    return "\n\n".join(blocks)


# ---------------------------------------------------------------------------
# 本地 HTTP 服务（把查询包成 JSON 接口）
# ---------------------------------------------------------------------------
def serve(host: str = "127.0.0.1", port: int = 8765, client: Optional[PubChemClient] = None) -> None:
    from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

    client = client or PubChemClient()

    class Handler(BaseHTTPRequestHandler):
        server_version = f"PubChemClient/{__version__}"

        def _send(self, status: int, payload: Dict[str, Any]) -> None:
            body = json.dumps(payload, ensure_ascii=False, indent=2).encode("utf-8")
            self.send_response(status)
            self.send_header("Content-Type", "application/json; charset=utf-8")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)

        def do_GET(self) -> None:  # noqa: N802
            parsed = urllib.parse.urlparse(self.path)
            query = urllib.parse.parse_qs(parsed.query)
            if parsed.path in ("/health", "/"):
                self._send(200, {"status": "ok", "version": __version__})
                return
            if parsed.path != "/lookup":
                self._send(404, {"error": "unknown endpoint, use /lookup?q=..."})
                return
            term = (query.get("q") or [""])[0].strip()
            if not term:
                self._send(400, {"error": "missing query parameter q"})
                return
            include_computed = (query.get("computed") or ["1"])[0] not in ("0", "false", "no")
            try:
                records = client.lookup(term, include_computed=include_computed)
            except NotFoundError as exc:
                self._send(404, {"error": str(exc), "query": term})
                return
            except PubChemError as exc:
                self._send(502, {"error": str(exc), "query": term})
                return
            if (query.get("first") or ["0"])[0] in ("1", "true", "yes") and records:
                self._send(200, records[0])
            else:
                self._send(200, {"query": term, "count": len(records), "records": records})

        def log_message(self, fmt: str, *args: Any) -> None:  # 安静一点
            sys.stderr.write("[serve] " + fmt % args + "\n")

    httpd = ThreadingHTTPServer((host, port), Handler)
    print(f"PubChem API 已启动: http://{host}:{port}/lookup?q=aspirin", file=sys.stderr)
    print("Ctrl+C 停止", file=sys.stderr)
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print("\n已停止", file=sys.stderr)
    finally:
        httpd.server_close()


# ---------------------------------------------------------------------------
# 交互模式
# ---------------------------------------------------------------------------
def interactive(client: PubChemClient, args: argparse.Namespace) -> int:
    try:
        sys.stdout.reconfigure(errors="replace")  # type: ignore[attr-defined]
    except (AttributeError, ValueError):
        pass
    print("PubChem 查询 —— 输入英文名 / CAS 号 / CID / InChIKey，回车查询；直接回车退出。")
    print("例：aspirin   或   50-78-2\n")
    while True:
        try:
            term = input("查询> ").strip()
        except (EOFError, KeyboardInterrupt):
            print("\n再见。")
            return 0
        if not term:
            print("再见。")
            return 0
        try:
            records = client.lookup(
                term,
                include_computed=not args.no_computed,
                refresh=args.refresh,
            )
        except PubChemError as exc:
            print(f"  ✗ 查询失败：{exc}\n")
            continue
        if args.first:
            records = records[:1]
        print()
        print(records_to_table(records))
        print()


# ---------------------------------------------------------------------------
# 命令行
# ---------------------------------------------------------------------------
def _read_identifiers(args: argparse.Namespace) -> List[str]:
    items: List[str] = list(args.identifiers or [])
    if args.file:
        text = Path(args.file).read_text(encoding="utf-8-sig")
        items.extend(line.strip() for line in text.splitlines())
    if not items and not sys.stdin.isatty():
        items.extend(line.strip() for line in sys.stdin.read().splitlines())
    return [item for item in items if item and not item.startswith("#")]


def _infer_format(args: argparse.Namespace) -> str:
    if args.format != "auto":
        return args.format
    if args.output:
        suffix = Path(args.output).suffix.lower()
        if suffix == ".json":
            return "json"
        if suffix in (".csv", ".tsv"):
            return "tsv" if suffix == ".tsv" else "csv"
    return "table"


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        prog="pubchem_client",
        description="按英文名 / CAS 号 / CID / InChIKey / SMILES 查询 PubChem，返回 SMILES 和全部计算化学属性。",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "示例:\n"
            "  python pubchem_client.py aspirin\n"
            "  python pubchem_client.py 50-78-2 --format json\n"
            "  python pubchem_client.py caffeine ibuprofen 64-17-5 -o out.csv\n"
            "  python pubchem_client.py -f names.txt --no-computed -o out.csv\n"
            "  python pubchem_client.py --serve --port 8765\n"
            "  curl 'http://127.0.0.1:8765/lookup?q=aspirin&first=1'\n"
        ),
    )
    parser.add_argument("identifiers", nargs="*", help="英文名 / CAS 号 / CID / InChIKey / SMILES，可多个")
    parser.add_argument("-i", "--interactive", action="store_true", help="进入交互模式，一条一条查（双击运行时默认进入）")
    parser.add_argument("-f", "--file", help="批量文件，每行一个标识符（# 开头为注释）")
    parser.add_argument("-o", "--output", help="输出文件（.json/.csv/.tsv 自动识别格式）")
    parser.add_argument("--format", choices=["auto", "json", "csv", "tsv", "table"], default="auto")
    parser.add_argument("--no-computed", action="store_true", help="不抓取 Computed Properties 表（省一次请求）")
    parser.add_argument("--first", action="store_true", help="只返回第一个 CID（默认返回全部命中）")
    parser.add_argument("--refresh", action="store_true", help="忽略本地缓存，强制重新请求")
    parser.add_argument("--cache-dir", help="缓存目录（默认 ~/.cache/pubchem_client）")
    parser.add_argument("--no-cache", action="store_true", help="禁用磁盘缓存")
    parser.add_argument("--timeout", type=float, default=30.0, help="单请求超时秒数（默认 30）")
    parser.add_argument("--ignore-errors", action="store_true", help="批量时跳过失败项继续")
    parser.add_argument("-v", "--verbose", action="store_true", help="在 stderr 显示请求/缓存统计")
    parser.add_argument("--serve", action="store_true", help="启动本地 HTTP 服务而不是查询")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8765)
    parser.add_argument("--version", action="version", version=f"pubchem_client {__version__}")

    args = parser.parse_args(argv)

    client = PubChemClient(
        timeout=args.timeout,
        cache_dir=args.cache_dir,
        use_cache=not args.no_cache,
    )

    if args.serve:
        serve(args.host, args.port, client)
        return 0

    stdin_is_tty = False
    try:
        stdin_is_tty = sys.stdin.isatty()
    except (AttributeError, ValueError):
        pass
    if args.interactive or (not args.identifiers and not args.file and stdin_is_tty):
        return interactive(client, args)

    identifiers = _read_identifiers(args)
    if not identifiers:
        parser.print_help()
        return 2

    try:
        records = client.lookup_many(
            identifiers,
            include_computed=not args.no_computed,
            refresh=args.refresh,
            on_error="skip" if args.ignore_errors else "raise",
        )
    except PubChemError as exc:
        print(f"错误: {exc}", file=sys.stderr)
        return 1

    if args.first:
        grouped: Dict[str, Dict[str, Any]] = {}
        for record in records:
            grouped.setdefault(record.get("query", ""), record)
        records = list(grouped.values())

    if args.verbose:
        print(
            f"请求 {client.stats['requests']} 次，缓存命中 {client.stats['cache_hits']} 次，"
            f"重试 {client.stats['retries']} 次",
            file=sys.stderr,
        )

    fmt = _infer_format(args)
    if fmt == "json":
        text = json.dumps(records[0] if len(records) == 1 else records,
                          ensure_ascii=False, indent=2)
    elif fmt in ("csv", "tsv"):
        text = records_to_csv(records, delimiter="\t" if fmt == "tsv" else ",").rstrip("\n")
    else:
        text = records_to_table(records)

    if args.output:
        encoding = "utf-8-sig" if fmt in ("csv", "tsv") else "utf-8"
        Path(args.output).write_text(text + ("\n" if not text.endswith("\n") else ""), encoding=encoding)
        print(f"已写入 {args.output}（{len(records)} 条记录）", file=sys.stderr)
    else:
        try:
            sys.stdout.reconfigure(encoding="utf-8")  # type: ignore[attr-defined]
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
