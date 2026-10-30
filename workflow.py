# -*- coding: utf-8 -*-
"""环境物质信息工作流：名称/CAS -> 合并表（xlsx/csv）。

用法：
    python workflow.py --input input.csv --out result.xlsx [--insecure]
                       [--molgpka auto|on|off] [--molgpka-images]
                       [--optional] [--resume]

输入：CSV/TXT（逗号分隔），可带表头（name/englishname + cas/CASRN 等），
      允许无表头（2 列 = name,cas；1 列 = CAS 或名称），允许只有 CAS 列。
输出：xlsx（results / salt_result / pka_sites / qc / meta）或 csv。
"""
import argparse
import csv
import hashlib
import io
import json
import os
import re
import shutil
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

from plugins._common import RDKIT_LIBS, enable_rdkit

enable_rdkit(RDKIT_LIBS)
from rdkit import Chem, RDLogger  # noqa: E402
from rdkit.Chem import Draw  # noqa: E402
from rdkit.Chem.MolStandardize import rdMolStandardize  # noqa: E402

RDLogger.DisableLog("rdApp.*")

HERE = Path(__file__).resolve().parent
PLUGINS = HERE / "plugins"

NAME_KEYS = {"name", "englishname", "compound", "compoundname", "chemicalname"}
CAS_KEYS = {"cas", "casrn", "casno", "casnumber", "cas_no"}
CAS_RE = re.compile(r"^\d{2,7}-\d{2}-\d$")


def comptox_cert_failed(timeout=8):
    """探测 CompTox 证书是否可用；失败（过期/不受信）时返回 True。"""
    import ssl
    import urllib.error
    import urllib.request

    url = ("https://comptox.epa.gov/dashboard-api/ccdapp1/search/chemical/"
           "equal-with-detail/aspirin")
    try:
        urllib.request.urlopen(url, timeout=timeout)
        return False
    except urllib.error.URLError as exc:
        reason = getattr(exc, "reason", None)
        if isinstance(reason, ssl.SSLCertVerificationError):
            return True
        text = f"{reason or exc}".lower()
        return "certificate" in text and ("expired" in text or "verify" in text)
    except Exception:  # noqa: BLE001
        return False

# ---- 列定义：按性质分类排列 -------------------------------------------------
CORE_COLUMNS = [
    # 标识
    "name", "cas", "iupac_name", "smiles",
    # 分子
    "molecular_formula", "molecular_weight", "exact_mass", "heavy_atom_count",
    "formal_charge",
    # 脂溶性
    "logP_pubchem", "logP_comptox", "logD_7.4", "logD_5.5", "logKoa",
    # 电离（pKa，两个来源放一起）
    "pka_acidic_comptox", "pka_basic_comptox",
    "pka_acid_min_molgpka", "pka_base_max_molgpka",
    # 氢键
    "hba_pubchem", "hbd_pubchem", "lipinski_hba", "lipinski_hbd",
    # 极性与尺寸
    "tpsa", "rotatable_bond_count", "labute_surface_area", "fraction_csp3",
    "atom_count", "heteroatom_count",
    # 环
    "ring_count", "aromatic_ring_count", "aliphatic_ring_count",
    "aromatic_carbocycle_count", "aliphatic_carbocycle_count",
    "aromatic_heterocycle_count", "aliphatic_heterocycle_count", "heterocycle_count",
    # 溶解与挥发
    "water_solubility", "vapor_pressure", "henry_law_constant",
    # 热物性
    "boiling_point", "melting_point", "density", "surface_tension",
    # 环境归趋
    "biodegradation_half_life", "bioconcentration_factor",
    "soil_adsorption_coefficient_koc", "atmospheric_hydroxylation_rate",
    "fish_biotransformation_half_life", "ready_biodegradability",
]

OPTIONAL_COLUMNS = [
    "flash_point", "refractive_index", "molar_refractivity", "molar_volume",
    "polarizability", "thermal_conductivity", "viscosity", "dielectric_constant",
    "bioaccumulation_factor",
]

SALT_COLUMNS = [
    "name", "cas", "iupac_name", "smiles",
    "molecular_formula", "molecular_weight", "exact_mass", "heavy_atom_count",
    "formal_charge", "parent_formula", "parent_mw",
    "logP_pubchem",
    "pka_acid_min_molgpka", "pka_base_max_molgpka",
    "hba_pubchem", "hbd_pubchem", "lipinski_hba", "lipinski_hbd",
    "tpsa", "rotatable_bond_count", "labute_surface_area", "fraction_csp3",
    "atom_count", "heteroatom_count",
    "ring_count", "aromatic_ring_count", "aliphatic_ring_count",
    "aromatic_carbocycle_count", "aliphatic_carbocycle_count",
    "aromatic_heterocycle_count", "aliphatic_heterocycle_count", "heterocycle_count",
]

CT_MAP = {
    "logd_7_4": "logD_7.4", "logd_5_5": "logD_5.5", "logkoa": "logKoa",
    "water_solubility": "water_solubility", "vapor_pressure": "vapor_pressure",
    "henrys_law": "henry_law_constant", "boiling_point": "boiling_point",
    "melting_point": "melting_point", "density": "density",
    "surface_tension": "surface_tension", "pka_acidic": "pka_acidic_comptox",
    "pka_basic": "pka_basic_comptox",
}
ENV_MAP = {
    "env_biodeg_half_life": "biodegradation_half_life",
    "env_bcf": "bioconcentration_factor",
    "env_koc": "soil_adsorption_coefficient_koc",
    "env_atmos_hydroxylation": "atmospheric_hydroxylation_rate",
    "env_fish_biotrans_half_life": "fish_biotransformation_half_life",
    "env_ready_biodeg": "ready_biodegradability",
}
OPT_MAP = {
    "flash_point": "flash_point", "index_of_refraction": "refractive_index",
    "molar_refractivity": "molar_refractivity", "molar_volume": "molar_volume",
    "polarizability": "polarizability", "thermal_conductivity": "thermal_conductivity",
    "viscosity": "viscosity", "dielectric_constant": "dielectric_constant",
    "env_baf": "bioaccumulation_factor",
}

RING_FIELDS = [
    ("num_rings", "ring_count"),
    ("num_aromatic_rings", "aromatic_ring_count"),
    ("num_aliphatic_rings", "aliphatic_ring_count"),
    ("num_aromatic_carbocycles", "aromatic_carbocycle_count"),
    ("num_aliphatic_carbocycles", "aliphatic_carbocycle_count"),
    ("num_aromatic_heterocycles", "aromatic_heterocycle_count"),
    ("num_aliphatic_heterocycles", "aliphatic_heterocycle_count"),
    ("num_heterocycles", "heterocycle_count"),
]


def is_cas(text):
    if not text or not CAS_RE.match(text):
        return False
    digits = text.replace("-", "")
    return sum((i + 1) * int(d) for i, d in enumerate(reversed(digits[:-1]))) % 10 == int(digits[-1])


def parent_inchikey(smiles):
    if not smiles:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
        if len(frags) > 1:
            parent = max(frags, key=lambda m: (m.GetNumHeavyAtoms(),
                                               sum(1 for a in m.GetAtoms() if a.GetAtomicNum() == 6)))
            parent = rdMolStandardize.Uncharger().uncharge(parent)
        else:
            parent = mol
    except Exception:  # noqa: BLE001
        parent = mol
    try:
        return Chem.MolToInchiKey(parent)
    except Exception:  # noqa: BLE001
        return None


def _read_input_lines(path):
    """读取 csv/txt/xlsx，自动判断分隔符与编码，返回二维字符串列表。"""
    p = Path(path)
    suffix = p.suffix.lower()
    if suffix in (".xlsx", ".xlsm"):
        from openpyxl import load_workbook

        wb = load_workbook(p, read_only=True, data_only=True)
        ws = wb.active
        lines = []
        for row in ws.iter_rows(values_only=True):
            cells = ["" if v is None else str(v).strip() for v in row]
            if any(cells):
                lines.append(cells)
        return lines

    text = None
    for encoding in ("utf-8-sig", "gbk", "utf-8"):
        try:
            text = p.read_text(encoding=encoding)
            break
        except UnicodeDecodeError:
            continue
    if text is None:
        text = p.read_text(encoding="utf-8", errors="replace")

    first_line = next((ln for ln in text.splitlines() if ln.strip()), "")
    if "\t" in first_line:
        delimiter = "\t"
    elif ";" in first_line and first_line.count(";") > first_line.count(","):
        delimiter = ";"
    else:
        delimiter = ","

    reader = csv.reader(io.StringIO(text), delimiter=delimiter)
    return [[c.strip() for c in row] for row in reader if any(c.strip() for c in row)]


def load_input(path):
    lines = _read_input_lines(path)
    if not lines:
        raise SystemExit("输入文件为空")

    first = lines[0]
    lower = [c.lower().replace(" ", "").replace("_", "") for c in first]
    name_idx = cas_idx = None
    for i, cell in enumerate(lower):
        if cell in {k.replace("_", "") for k in NAME_KEYS} and name_idx is None:
            name_idx = i
        if cell in {k.replace("_", "") for k in CAS_KEYS} and cas_idx is None:
            cas_idx = i

    has_header = name_idx is not None or cas_idx is not None
    data_lines = lines[1:] if has_header else lines

    rows = []
    for n, cells in enumerate(data_lines, start=1):
        if has_header:
            name = cells[name_idx] if name_idx is not None and name_idx < len(cells) else ""
            cas = cells[cas_idx] if cas_idx is not None and cas_idx < len(cells) else ""
        else:
            if len(cells) >= 2:
                name, cas = cells[0], cells[1]
            elif len(cells) == 1:
                if is_cas(cells[0]):
                    name, cas = "", cells[0]
                else:
                    name, cas = cells[0], ""
            else:
                continue
        if not name and not cas:
            continue
        rows.append({"row_id": n, "name": name, "cas": cas})
    return rows, {"has_header": has_header, "name_idx": name_idx, "cas_idx": cas_idx}


def write_csv(path, rows, columns):
    with open(path, "w", encoding="utf-8-sig", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({c: ("" if row.get(c) is None else row.get(c)) for c in columns})


def save_workbook_with_fallback(wb, out):
    """保存 xlsx；目标被占用或无权限时自动换名/换目录，并返回实际路径。"""
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    alt_name = f"{Path(out).stem}_{stamp}{Path(out).suffix}"
    candidates = [
        Path(out),
        Path(out).with_name(alt_name),
        Path.home() / "Documents" / alt_name,
        Path(os.environ.get("TEMP", ".")) / alt_name,
    ]
    last_error = None
    for path in candidates:
        try:
            path.parent.mkdir(parents=True, exist_ok=True)
            wb.save(path)
            if path != Path(out):
                print(f"提示：{out} 无法写入，已改存为 {path}", file=sys.stderr)
            return path
        except PermissionError as exc:
            last_error = exc
            continue
        except OSError as exc:
            last_error = exc
            continue
    raise SystemExit(
        f"无法写入 {out}：{last_error}\n"
        "常见原因：文件正被 Excel 打开、目标文件夹只读，或该目录有写入限制。\n"
        "请关闭 Excel 后重试，或把输出路径改到 D:\\ 下（例如 D:\\result.xlsx）。"
    )


def run_plugin(name, rows, columns, extra_args, tmpdir, progress=None):
    """运行插件；progress(i, n) 回调用于显示进度。"""
    if not rows:
        return []
    src = Path(tmpdir) / f"{name}_in.csv"
    out = Path(tmpdir) / f"{name}_out.json"
    write_csv(src, rows, columns)
    cmd = [sys.executable, str(PLUGINS / f"{name}_plugin.py"), "--in", str(src), "--out", str(out),
           "--progress"]
    cmd += extra_args
    env = dict(
        os.environ,
        PYTHONIOENCODING="utf-8",
        PUBCHEM_CACHE_DIR=str(HERE / ".cache" / "pubchem"),
        COMPTOX_CACHE_DIR=str(HERE / ".cache" / "comptox"),
    )
    proc = subprocess.Popen(
        cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE,
        text=True, encoding="utf-8", errors="replace", env=env,
    )
    tail = []
    for line in proc.stderr:
        line = line.rstrip("\n")
        if line.startswith("PROGRESS "):
            try:
                _, i, n = line.split()
                if progress:
                    progress(int(i), int(n))
            except ValueError:
                pass
        elif line.strip():
            tail.append(line)
            if len(tail) > 30:
                tail.pop(0)
    proc.wait()
    if proc.returncode != 0:
        detail = "\n".join(tail[-8:]) or "unknown"
        raise SystemExit(f"{name} 插件失败：{detail}")
    return json.loads(out.read_text(encoding="utf-8"))


def pick_comptox(rec, key):
    return ((rec.get("properties") or {}).get(key) or {}).get("value")


def pick_logp(rec):
    item = (rec.get("properties") or {}).get("logkow") or {}
    sources = item.get("sources") or {}
    for src, val in sources.items():
        if src.lower().startswith("percepta"):
            return val
    for src, val in sources.items():
        if src.lower().startswith("opera"):
            return val
    return item.get("experimental")


def pick_env(rec, key):
    return ((rec.get("envfate") or {}).get(key) or {}).get("value")


def to_num(value):
    if value is None or isinstance(value, (int, float)):
        return value
    try:
        return float(value)
    except (TypeError, ValueError):
        return value


def safe_name(text, fallback):
    text = re.sub(r"[^\w\-.]+", "_", (text or "").strip())
    return text or fallback


def render_site_image(mol_block, atom_index, path):
    mol = Chem.MolFromMolBlock(mol_block, removeHs=False) if mol_block else None
    if mol is None:
        return False
    n = mol.GetNumAtoms()
    idx = atom_index if 0 <= atom_index < n else max(0, min(atom_index - 1, n - 1))
    atom = mol.GetAtomWithIdx(idx)
    if atom.GetSymbol() == "H":
        heavy = [nb.GetIdx() for nb in atom.GetNeighbors() if nb.GetAtomicNum() > 1]
        if heavy:
            idx = heavy[0]
    try:
        # 用原子映射号把高亮原子带到去氢后的图上，得到干净的 2D 结构
        mol.GetAtomWithIdx(idx).SetAtomMapNum(999)
        clean = Chem.RemoveHs(mol)
        targets = [a.GetIdx() for a in clean.GetAtoms() if a.GetAtomMapNum() == 999]
        if targets:
            clean.GetAtomWithIdx(targets[0]).SetAtomMapNum(0)
            idx = targets[0]
        img = Draw.MolToImage(clean, size=(320, 320), highlightAtoms=[idx])
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        img.save(path)
        return True
    except Exception:  # noqa: BLE001
        return False


def render_structure_image(smiles, path):
    """用 RDKit 画 2D 结构图（盐按整体画）。"""
    mol = Chem.MolFromSmiles(smiles) if smiles else None
    if mol is None:
        return False
    try:
        img = Draw.MolToImage(mol, size=(300, 220))
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        img.save(path)
        return True
    except Exception:  # noqa: BLE001
        return False


def main():
    ap = argparse.ArgumentParser(
        prog="workflow.py",
        description="环境物质信息采集：输入名称/CAS 表，输出合并后的 Excel。",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "常用示例：\n"
            "  workflow.cmd --input 输入.csv --out 结果.xlsx\n"
            "  workflow.cmd --input 输入.csv --out 结果.xlsx --with-molgpka\n"
            "  workflow.cmd --input 输入.csv --out 结果.xlsx --with-molgpka --no-images\n"
            "\n"
            "参数说明：\n"
            "  --with-molgpka   采集 MolGpKa 位点 pKa（默认不采集，较慢）\n"
            "  --no-images      采集 pKa 时不嵌入结构图\n"
            "  --optional       额外输出 9 个不常用的 CompTox 属性\n"
            "  --insecure       强制跳过 TLS 校验（默认自动检测，一般不用加）\n"
            "  --resume         中断后继续，复用已完成的行\n"
        ),
    )
    ap.add_argument("--input", required=True, help="输入文件（csv/txt，逗号分隔）")
    ap.add_argument("--out", required=True, help="输出文件（.xlsx）")
    ap.add_argument("--format", choices=["xlsx", "csv"], default="xlsx")
    ap.add_argument("--with-molgpka", action="store_true",
                    help="采集 MolGpKa 位点 pKa（默认关闭）")
    ap.add_argument("--no-images", action="store_true",
                    help="采集 pKa 时不把高亮结构图嵌入 Excel")
    ap.add_argument("--optional", action="store_true", help="输出可选 C 组属性")
    ap.add_argument("--insecure", action="store_true",
                    help="强制跳过 TLS 证书校验（默认自动检测，一般不用加）")
    ap.add_argument("--resume", action="store_true", help="复用已完成的行")
    ap.add_argument("--cache-dir", default=str(HERE / ".cache"))
    args = ap.parse_args()
    t0 = time.time()
    do_molgpka = args.with_molgpka
    do_images = args.with_molgpka and not args.no_images

    insecure = args.insecure
    if not insecure:
        if comptox_cert_failed():
            insecure = True
            print("提示：CompTox 证书校验未通过，已自动跳过校验（仅影响公开数据查询）。",
                  file=sys.stderr)

    rows, _ = load_input(args.input)
    if not rows:
        raise SystemExit("没有解析到任何行")

    key = hashlib.sha1(Path(args.input).read_bytes()).hexdigest()[:12]
    cache_dir = Path(args.cache_dir) / key
    cache_dir.mkdir(parents=True, exist_ok=True)

    cached = {}
    if args.resume:
        for f in cache_dir.glob("*.json"):
            try:
                item = json.loads(f.read_text(encoding="utf-8"))
                cached[int(item["_row_id"])] = item
            except Exception:  # noqa: BLE001
                pass
    pending = [r for r in rows if r["row_id"] not in cached]

    tmp_root = HERE / ".work"
    tmpdir = tmp_root / f"run_{os.getpid()}_{int(time.time())}"
    tmpdir.mkdir(parents=True, exist_ok=True)
    extra = ["--insecure"] if insecure else []

    total_rows = len(rows)
    done_cached = total_rows - len(pending)
    print(f"读取到 {total_rows} 个物质；本次查询 {len(pending)} 个"
          + (f"（缓存 {done_cached} 个）" if done_cached else ""))

    stages = ["PubChem", "RDKit", "CompTox"] + (["MolGpKa"] if do_molgpka else [])
    base_pct = (done_cached / total_rows * 100) if total_rows else 0.0
    span_pct = ((len(pending) / total_rows) * 100) if total_rows else 0.0

    last_stage = {"name": None}

    def show_progress(pct, stage):
        if stage != last_stage["name"]:
            if last_stage["name"] is not None:
                sys.stdout.write("\n")
            last_stage["name"] = stage
        sys.stdout.write(f"\r[{pct:5.1f}%] {stage:<10}")
        sys.stdout.flush()

    def stage_cb(index):
        stage_name = stages[index]
        stage_base = base_pct + span_pct * index / len(stages)
        stage_span = span_pct / len(stages)

        def inner(i, n):
            show_progress(stage_base + stage_span * (i / n if n else 1), stage_name)

        return inner

    pub = run_plugin("pubchem", pending, ["row_id", "name", "cas"], extra, tmpdir,
                     progress=stage_cb(0))
    pub_map = {int(r["row_id"]): r for r in pub}

    rd_rows = [
        {"row_id": r["row_id"], "smiles": pub_map[int(r["row_id"])].get("smiles")}
        for r in pending
        if pub_map.get(int(r["row_id"]), {}).get("smiles")
    ]
    rd = run_plugin("rdkit", rd_rows, ["row_id", "smiles"], [], tmpdir, progress=stage_cb(1))
    rd_map = {int(r["row_id"]): r for r in rd}

    ct = run_plugin("comptox", pending, ["row_id", "name", "cas"], extra, tmpdir,
                    progress=stage_cb(2))
    ct_map = {int(r["row_id"]): r for r in ct}

    mol_rows = []
    if do_molgpka:
        for r in pending:
            rid = r["row_id"]
            rk = rd_map.get(rid) or {}
            if not rk.get("parent_smiles"):
                continue
            # 盐也计算，但用最大有机片段（rdkit 插件已给出 parent_smiles）
            if rk.get("is_inorganic"):
                continue
            mol_rows.append({"row_id": rid, "smiles": rk["parent_smiles"]})
        mol = run_plugin("molgpka", mol_rows, ["row_id", "smiles"], [], tmpdir,
                         progress=stage_cb(3))
    else:
        mol = []
    mol_map = {int(r["row_id"]): r for r in mol}
    show_progress(100, "整理结果")
    print()

    results, salts, qc_rows, pka_records = [], [], [], []
    image_dir = Path(str(args.out).rsplit(".", 1)[0] + "_images")
    photo_enabled = args.format == "xlsx"

    for r in rows:
        rid = int(r["row_id"])
        if rid in cached:
            item = cached[rid]
            result = dict(item["result"])
            if photo_enabled and not result.get("photo"):
                fname = f"{safe_name(result.get('name'), r['cas'] or str(rid))}.png"
                photo_path = image_dir / "structures" / fname
                if render_structure_image(result.get("smiles"), photo_path):
                    result["photo"] = str(photo_path)
            qc_rows.append(item["qc"])
            pka_records.append({
                "name": result.get("name"),
                "cas": result.get("cas"),
                "sites": item.get("pka_sites") or [],
            })
            if item["qc"].get("pubchem_status") == "ok":
                if item.get("is_salt"):
                    salts.append(result)
                else:
                    results.append(result)
            continue

        p = pub_map.get(rid) or {}
        rk = rd_map.get(rid) or {}
        c = ct_map.get(rid) or {}
        m = mol_map.get(rid) or {}

        p_ok = p.get("status") == "ok"
        c_ok = c.get("status") == "ok"

        if not p_ok:
            identity = "NO_PUBCHEM"
        elif not c_ok:
            identity = "NO_COMPTOX"
        elif p.get("inchikey") and p.get("inchikey") == c.get("inchikey"):
            identity = "OK"
        else:
            ppk = rk.get("parent_inchikey") or parent_inchikey(p.get("smiles"))
            cpk = parent_inchikey(c.get("smiles"))
            identity = "SALT_FORM_DIFFERS" if (ppk and cpk and ppk == cpk) else "MISMATCH"

        is_salt = bool(rk.get("is_salt"))
        fill_ct = (identity == "OK") and not is_salt

        # CAS：PubChem 同义词优先；没有则 CompTox；再没有用输入值
        cas_list = p.get("cas_found") or []
        if r["cas"] and r["cas"] in cas_list:
            resolved_cas = r["cas"]
        elif cas_list:
            resolved_cas = cas_list[0]
        else:
            resolved_cas = c.get("casrn") or r["cas"]

        row = {
            "name": r["name"] or p.get("title") or c.get("preferred_name") or r["cas"],
            "cas": resolved_cas,
            "iupac_name": p.get("iupac_name") or c.get("iupac_name"),
            "smiles": p.get("smiles"),
            "molecular_formula": p.get("formula"),
            "molecular_weight": to_num(p.get("mw")),
            "exact_mass": to_num(p.get("exact_mass")),
            "heavy_atom_count": p.get("heavy_atoms"),
            "formal_charge": p.get("charge"),
            "logP_pubchem": to_num(p.get("xlogp")),
            "logP_comptox": pick_logp(c) if fill_ct else None,
            "hba_pubchem": p.get("hba"),
            "hbd_pubchem": p.get("hbd"),
            "lipinski_hba": rk.get("lipinski_hba"),
            "lipinski_hbd": rk.get("lipinski_hbd"),
            "tpsa": rk.get("tpsa"),
            "rotatable_bond_count": p.get("rotatable_bonds"),
            "atom_count": rk.get("num_atoms"),
            "heteroatom_count": rk.get("num_heteroatoms"),
            "labute_surface_area": rk.get("labute_asa"),
            "fraction_csp3": rk.get("fraction_csp3"),
            "pka_acid_min_molgpka": m.get("pka_acid_min"),
            "pka_base_max_molgpka": m.get("pka_base_max"),
        }
        for k, out in RING_FIELDS:
            row[out] = rk.get(k)
        for k, out in CT_MAP.items():
            row[out] = pick_comptox(c, k) if fill_ct else None
        for k, out in ENV_MAP.items():
            row[out] = pick_env(c, k) if fill_ct else None
        if args.optional:
            for k, out in OPT_MAP.items():
                if k.startswith("env_"):
                    row[out] = pick_env(c, k) if fill_ct else None
                else:
                    row[out] = pick_comptox(c, k) if fill_ct else None
        if is_salt:
            row["parent_formula"] = rk.get("parent_formula")
            row["parent_mw"] = rk.get("parent_mw")
        if photo_enabled:
            fname = f"{safe_name(row['name'], r['cas'] or str(rid))}.png"
            photo_path = image_dir / "structures" / fname
            if render_structure_image(p.get("smiles"), photo_path):
                row["photo"] = str(photo_path)

        qc = {
            "name": row["name"],
            "input_cas": r["cas"],
            "cas_resolved": resolved_cas,
            "overall_status": "OK" if (p_ok and c_ok) else ("NO_PUBCHEM" if not p_ok else "NO_COMPTOX"),
            "pubchem_status": p.get("status", "not_run"),
            "pubchem_cid": p.get("cid"),
            "pubchem_cas_found": ";".join(p.get("cas_found") or []),
            "comptox_status": c.get("status", "not_run"),
            "comptox_dtxsid": c.get("dtxsid"),
            "identity_status": identity,
            "pubchem_inchikey": p.get("inchikey"),
            "comptox_inchikey": c.get("inchikey"),
            "pubchem_parent_inchikey": rk.get("parent_inchikey"),
            "comptox_parent_inchikey": parent_inchikey(c.get("smiles")),
            "structure_basis": rk.get("structure_basis"),
            "is_salt": rk.get("is_salt"),
            "is_inorganic": rk.get("is_inorganic"),
            "compound_is_canonicalized": p.get("canonicalized"),
            "molgpka_status": m.get("status", "skipped"),
            "molgpka_basis": ("parent_fragment" if is_salt else "molecule")
                              if do_molgpka else None,
            "messages": "; ".join(
                x for x in [p.get("message"), c.get("message"), m.get("message"), rk.get("note")] if x
            ),
        }

        row_pka = []
        for i, site in enumerate(m.get("sites") or [], start=1):
            item = {
                "site_no": i,
                "site_type": site["site_type"],
                "pka": site["pka"],
                "atom_index": site["atom_index"],
            }
            if do_images and m.get("mol_block"):
                fname = f"{safe_name(row['name'], r['cas'] or str(rid))}_{i}.png"
                if render_site_image(m["mol_block"], site["atom_index"], image_dir / fname):
                    item["image_file"] = str(image_dir / fname)
            row_pka.append(item)

        pka_records.append({"name": row["name"], "cas": resolved_cas, "sites": row_pka})

        cache_dir.joinpath(f"{rid}.json").write_text(
            json.dumps({"_row_id": rid, "is_salt": is_salt, "result": row, "qc": qc,
                        "pka_sites": row_pka}, ensure_ascii=False),
            encoding="utf-8",
        )
        if p_ok:
            if is_salt:
                salts.append(row)
            else:
                results.append(row)
        qc_rows.append(qc)

    columns = [c for c in CORE_COLUMNS if do_molgpka or c not in
               ("pka_acid_min_molgpka", "pka_base_max_molgpka")]
    if args.optional:
        columns += OPTIONAL_COLUMNS
    salt_columns = [c for c in SALT_COLUMNS if do_molgpka or c not in
                    ("pka_acid_min_molgpka", "pka_base_max_molgpka",
                     "parent_formula", "parent_mw")]
    if photo_enabled:
        columns.insert(columns.index("smiles") + 1, "photo")
        salt_columns.insert(salt_columns.index("smiles") + 1, "photo")
    pka_records = list(pka_records)
    max_sites = max([len(x["sites"]) for x in pka_records] or [0])
    pka_cols = ["name", "cas", "site_count"]
    for i in range(1, max_sites + 1):
        pka_cols += [f"site{i}_no", f"site{i}_type", f"site{i}_pka", f"site{i}_atom_index"]
        if do_images:
            pka_cols.append(f"site{i}_image")
    pka_wide = []
    for rec in pka_records:
        row = {"name": rec["name"], "cas": rec["cas"], "site_count": len(rec["sites"])}
        for i, site in enumerate(rec["sites"], start=1):
            row[f"site{i}_no"] = site.get("site_no")
            row[f"site{i}_type"] = site.get("site_type")
            row[f"site{i}_pka"] = site.get("pka")
            row[f"site{i}_atom_index"] = site.get("atom_index")
            if do_images:
                row[f"site{i}_image"] = site.get("image_file")
        pka_wide.append(row)
    pka_csv_cols = [c for c in pka_cols if not c.endswith("_image")]
    qc_cols = list(qc_rows[0].keys()) if qc_rows else []
    out = Path(args.out)

    if args.format == "csv":
        write_csv(out.with_name(out.stem + "_results.csv"), results, columns)
        write_csv(out.with_name(out.stem + "_salt_result.csv"), salts, salt_columns)
        written = f"{out.stem}_results.csv / _salt_result.csv"
        if do_molgpka and pka_wide:
            write_csv(out.with_name(out.stem + "_pka_sites.csv"), pka_wide, pka_csv_cols)
            written += " / _pka_sites.csv"
        write_csv(out.with_name(out.stem + "_qc.csv"), qc_rows, qc_cols)
        print(f"已写出 CSV：{written} / _qc.csv", file=sys.stderr)
    else:
        from openpyxl import Workbook
        from openpyxl.drawing.image import Image as XLImage
        from openpyxl.styles import Font
        from openpyxl.utils import get_column_letter

        wb = Workbook()
        ws = wb.active
        ws.title = "results"
        ws.append(columns)
        for rec in results:
            ws.append([rec.get(c) if c != "photo" else None for c in columns])
        for cell in ws[1]:
            cell.font = Font(bold=True)
        ws.freeze_panes = "A2"
        ws.auto_filter.ref = ws.dimensions
        for i, col in enumerate(columns, start=1):
            ws.column_dimensions[get_column_letter(i)].width = min(max(10, len(col) + 2), 34)
        if photo_enabled:
            pcol = columns.index("photo") + 1
            ws.column_dimensions[get_column_letter(pcol)].width = 23
            for r_idx, rec in enumerate(results, start=2):
                p = rec.get("photo")
                if p and Path(p).exists():
                    img = XLImage(p)
                    img.width = 150
                    img.height = 110
                    ws.add_image(img, f"{get_column_letter(pcol)}{r_idx}")
                    ws.row_dimensions[r_idx].height = 88

        ws_salt = wb.create_sheet("salt_result")
        ws_salt.append(salt_columns)
        for rec in salts:
            ws_salt.append([rec.get(c) if c != "photo" else None for c in salt_columns])
        for cell in ws_salt[1]:
            cell.font = Font(bold=True)
        ws_salt.freeze_panes = "A2"
        ws_salt.auto_filter.ref = ws_salt.dimensions
        if photo_enabled:
            pcol = salt_columns.index("photo") + 1
            ws_salt.column_dimensions[get_column_letter(pcol)].width = 23
            for r_idx, rec in enumerate(salts, start=2):
                p = rec.get("photo")
                if p and Path(p).exists():
                    img = XLImage(p)
                    img.width = 150
                    img.height = 110
                    ws_salt.add_image(img, f"{get_column_letter(pcol)}{r_idx}")
                    ws_salt.row_dimensions[r_idx].height = 88

        if do_molgpka and pka_wide:
            ws2 = wb.create_sheet("pka_sites")
            ws2.append(pka_cols)
            for rec in pka_wide:
                ws2.append([rec.get(c) if not c.endswith("_image") else None for c in pka_cols])
            for cell in ws2[1]:
                cell.font = Font(bold=True)
            ws2.freeze_panes = "C2"

            if do_images and max_sites:
                image_cols = {f"site{i}_image": pka_cols.index(f"site{i}_image") + 1
                              for i in range(1, max_sites + 1)}
                for col_name, col_idx in image_cols.items():
                    ws2.column_dimensions[get_column_letter(col_idx)].width = 24
                for r_idx, rec in enumerate(pka_records, start=2):
                    for i, site in enumerate(rec["sites"], start=1):
                        path = site.get("image_file")
                        col_idx = image_cols.get(f"site{i}_image")
                        if path and col_idx and Path(path).exists():
                            img = XLImage(path)
                            img.width = 160
                            img.height = 160
                            ws2.add_image(img, f"{get_column_letter(col_idx)}{r_idx}")
                    ws2.row_dimensions[r_idx].height = 126

        ws3 = wb.create_sheet("qc")
        ws3.append(qc_cols)
        for rec in qc_rows:
            ws3.append([rec.get(c) for c in qc_cols])
        for cell in ws3[1]:
            cell.font = Font(bold=True)
        ws3.freeze_panes = "A2"

        ws4 = wb.create_sheet("meta")
        meta = [
            ["generated_at", datetime.now().isoformat(timespec="seconds")],
            ["input_file", str(Path(args.input).resolve())],
            ["input_sha1", key],
            ["rows", len(rows)],
            ["with_molgpka", do_molgpka],
            ["molgpka_images", do_images],
            ["optional_group", args.optional],
            ["insecure_tls", insecure],
            ["pubchem", "PUG REST, XLogP3 / Cactvs counts"],
            ["comptox", "Dashboard API; Percepta 2023.1.2 > OPERA 2.6 > TEST > EPISUITE; experimental only for logP"],
            ["rdkit", "desalt + neutralized parent; descriptors on parent"],
            ["molgpka", "xundrug.cn port 5001; per-site pKa"],
            ["units", "MW/exact mass g/mol; TPSA/Labute surface area A^2; solubility mol/L; vapor pressure mmHg; Henry atm-m3/mol"],
        ]
        for pair in meta:
            ws4.append(pair)
        ws4.column_dimensions["A"].width = 22
        ws4.column_dimensions["B"].width = 90

        out = save_workbook_with_fallback(wb, out)
        args.out = str(out)
        print(f"已写出 {out}", file=sys.stderr)

    shutil.rmtree(tmpdir, ignore_errors=True)
    elapsed = time.time() - t0
    ok = sum(1 for q in qc_rows if q.get("overall_status") == "OK")
    failed = len(rows) - len(results) - len(salts)
    print()
    print("=" * 52)
    print("任务完成")
    print(f"  输入物质：{len(rows)} 个")
    print(f"  非盐结果：{len(results)} 个")
    print(f"  盐类结果：{len(salts)} 个")
    print(f"  未命中/跳过：{failed} 个（详见 qc）")
    print(f"  PubChem+CompTox 均命中：{ok} 个")
    print(f"  耗时：{int(elapsed // 60)} 分 {elapsed % 60:.1f} 秒")
    print(f"  输出：{args.out}")
    print("=" * 52)


if __name__ == "__main__":
    main()
