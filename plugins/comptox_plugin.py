# -*- coding: utf-8 -*-
"""CompTox 插件：CAS/名称 -> DTXSID、理化性质（预测+实验）、环境归趋。

独立用法：
    python comptox_plugin.py --in rows.csv --out comptox.json [--insecure]
输入 CSV 列：row_id,name,cas
"""
import argparse
import sys

from _common import COMPTOX_DIR, allow_insecure, fail, read_rows, write_json

sys.path.insert(0, str(COMPTOX_DIR))
from comptox_client import ComptoxClient  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", dest="out", required=True)
    ap.add_argument("--insecure", action="store_true")
    ap.add_argument("--timeout", type=float, default=45.0)
    ap.add_argument("--retries", type=int, default=4)
    ap.add_argument("--progress", action="store_true")
    args = ap.parse_args()

    allow_insecure(args.insecure)
    client = ComptoxClient(timeout=args.timeout, retries=args.retries, min_interval=0.8)

    ct_rows = read_rows(args.inp)
    total = len(ct_rows)
    results = []
    for index, row in enumerate(ct_rows, start=1):
        if args.progress:
            print(f"PROGRESS {index} {total}", file=sys.stderr, flush=True)
        rid = row.get("row_id", "")
        name = (row.get("name") or "").strip()
        cas = (row.get("cas") or "").strip()
        rec = {"row_id": rid, "query": cas or name, "status": ""}
        query = cas or name
        if not query:
            rec["status"] = "empty"
            results.append(rec)
            continue
        try:
            data = client.lookup(
                query, prefer="percepta", with_experimental=True, with_envfate=True
            )
        except Exception as exc:  # noqa: BLE001
            rec.update(status="error", message=str(exc)[:160])
            results.append(rec)
            continue

        rec.update(
            status="ok",
            dtxsid=data.get("dtxsid"),
            casrn=data.get("casrn"),
            preferred_name=data.get("preferred_name"),
            iupac_name=data.get("acd_iupac_name"),
            smiles=data.get("smiles"),
            inchikey=data.get("inchikey"),
            formula=data.get("formula"),
            properties=data.get("properties") or {},
            envfate=data.get("envfate") or {},
        )
        results.append(rec)

    write_json(args.out, results)
    ok = sum(1 for r in results if r.get("status") == "ok")
    print(f"comptox: {ok}/{len(results)} ok -> {args.out}", file=sys.stderr)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:  # noqa: BLE001
        fail(str(exc))
