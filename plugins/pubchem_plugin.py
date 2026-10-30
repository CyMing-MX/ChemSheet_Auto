# -*- coding: utf-8 -*-
"""PubChem 插件：名称/CAS -> CID、IUPAC 名、SMILES、属性、CAS 同义词。

独立用法：
    python pubchem_plugin.py --in rows.csv --out pubchem.json [--insecure]
输入 CSV 列：row_id,name,cas
"""
import argparse
import re
import sys
import time

from _common import PUBCHEM_DIR, allow_insecure, fail, read_rows, write_json

sys.path.insert(0, str(PUBCHEM_DIR))
from pubchem_client import PubChemClient, PubChemError, is_cas_number  # noqa: E402


def retry(fn, times=3, delay=1.5):
    last = None
    for attempt in range(times):
        try:
            return fn()
        except Exception as exc:  # noqa: BLE001
            last = exc
            if attempt < times - 1:
                time.sleep(delay * (attempt + 1))
    raise last

CAS_RE = re.compile(r"^\d{2,7}-\d{2}-\d$")
PROPS = [
    "Title", "IUPACName", "SMILES", "InChIKey", "MolecularFormula", "MolecularWeight",
    "ExactMass", "HeavyAtomCount", "RotatableBondCount", "HBondDonorCount",
    "HBondAcceptorCount", "Charge", "XLogP",
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", dest="out", required=True)
    ap.add_argument("--insecure", action="store_true")
    ap.add_argument("--progress", action="store_true")
    ap.add_argument("--timeout", type=float, default=45.0)
    ap.add_argument("--retries", type=int, default=4)
    args = ap.parse_args()

    allow_insecure(args.insecure)
    client = PubChemClient(timeout=args.timeout, retries=args.retries, min_interval=0.5)
    rows = read_rows(args.inp)
    results = []
    total = len(rows)

    for index, row in enumerate(rows, start=1):
        if args.progress:
            print(f"PROGRESS {index} {total}", file=sys.stderr, flush=True)
        rid = row.get("row_id", "")
        name = (row.get("name") or "").strip()
        cas = (row.get("cas") or "").strip()
        rec = {"row_id": rid, "input_name": name, "input_cas": cas, "status": ""}

        query, used = None, ""
        if cas and is_cas_number(cas):
            query, used = cas, "cas"
        elif name:
            query, used = name, "name"
        elif cas:
            query, used = cas, "cas_invalid_checksum"
        else:
            rec.update(status="empty")
            results.append(rec)
            continue

        try:
            hits = retry(lambda: client.resolve(query))
            cid = hits[0]["cid"]
        except Exception as exc:  # noqa: BLE001
            rec.update(status="not_found", message=str(exc)[:160])
            results.append(rec)
            continue

        rec.update(query_used=query, query_type=used, cid=cid, method=hits[0].get("method"))
        try:
            props = retry(lambda: client.properties([cid], names=PROPS)).get(cid, {})
        except Exception as exc:  # noqa: BLE001
            rec.update(status="error", message=f"properties: {exc}"[:160])
            results.append(rec)
            continue

        rec.update(
            title=props.get("title"),
            iupac_name=props.get("iupac_name"),
            smiles=props.get("smiles"),
            inchikey=props.get("inchikey"),
            formula=props.get("molecular_formula"),
            mw=props.get("molecular_weight"),
            exact_mass=props.get("exact_mass"),
            heavy_atoms=props.get("heavy_atom_count"),
            rotatable_bonds=props.get("rotatable_bond_count"),
            hbd=props.get("h_bond_donor_count"),
            hba=props.get("h_bond_acceptor_count"),
            charge=props.get("formal_charge"),
            xlogp=props.get("xlogp"),
        )

        # CAS 同义词
        try:
            data = retry(lambda: client._get_json(
                f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/synonyms/JSON"
            ), times=2)
            syns = (data.get("InformationList") or {}).get("Information") or [{}]
            names = syns[0].get("Synonym") or []
            rec["cas_found"] = sorted({s for s in names if CAS_RE.match(s)})
        except Exception as exc:  # noqa: BLE001
            rec["cas_found"] = []
            rec["synonyms_error"] = str(exc)[:80]

        # PubChem 规范化标记
        try:
            pv = retry(lambda: client._get_json(
                "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/"
                f"{cid}/JSON?heading=Computed+Properties"
            ), times=2)
            canon = None
            stack = [pv.get("Record") or {}]
            while stack:
                sec = stack.pop()
                if sec.get("TOCHeading") == "Computed Properties":
                    for sub in sec.get("Section") or ():
                        if sub.get("TOCHeading") == "Compound Is Canonicalized":
                            infos = sub.get("Information") or [{}]
                            val = (infos[0].get("Value") or {}).get("StringWithMarkup") or []
                            canon = "".join(x.get("String", "") for x in val).strip() or None
                    break
                stack.extend(sec.get("Section") or ())
            rec["canonicalized"] = canon
        except Exception:  # noqa: BLE001
            rec["canonicalized"] = None

        rec["status"] = "ok"
        results.append(rec)

    write_json(args.out, results)
    ok = sum(1 for r in results if r.get("status") == "ok")
    print(f"pubchem: {ok}/{len(results)} ok -> {args.out}", file=sys.stderr)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:  # noqa: BLE001
        fail(str(exc))
