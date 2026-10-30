# -*- coding: utf-8 -*-
"""MolGpKa 插件：SMILES -> 逐位点 pKa（酸/碱）。

独立用法：
    python molgpka_plugin.py --in rows.csv --out molgpka.json
输入 CSV 列：row_id,smiles
"""
import argparse
import json
import sys
import urllib.error
import urllib.request

from _common import fail, read_rows, write_json

URL = "http://xundrug.cn:5001/modules/upload0/"
TOKEN = "O05DriqqQLlry9kmpCwms2IJLC0MuLQ7"


def predict(smiles, token, timeout=40):
    boundary = "----WebKitFormBoundaryMolGpKa"
    body = (
        f"--{boundary}\r\n"
        f'Content-Disposition: form-data; name="Smiles"; filename="tmg"\r\n'
        f"Content-Type: application/octet-stream\r\n\r\n"
        f"{smiles}\r\n"
        f"--{boundary}--\r\n"
    ).encode("utf-8")
    req = urllib.request.Request(
        URL,
        data=body,
        headers={
            "Content-Type": f"multipart/form-data; boundary={boundary}",
            "token": token,
            "User-Agent": "env-pipeline/1.0",
        },
        method="POST",
    )
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        data = json.loads(resp.read().decode("utf-8"))
    gen = data.get("gen_datas") or {}
    sites = []
    for kind, label in (("Acid", "acid"), ("Base", "base")):
        for atom, value in (gen.get(kind) or {}).items():
            try:
                pka = float(value)
            except (TypeError, ValueError):
                continue
            sites.append({"site_type": label, "atom_index": int(atom), "pka": round(pka, 3)})
    return sites, gen.get("MolBlock") or ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", dest="out", required=True)
    ap.add_argument("--token", default=TOKEN)
    ap.add_argument("--timeout", type=float, default=40.0)
    ap.add_argument("--progress", action="store_true")
    args = ap.parse_args()

    mol_rows = read_rows(args.inp)
    total = len(mol_rows)
    results = []
    for index, row in enumerate(mol_rows, start=1):
        if args.progress:
            print(f"PROGRESS {index} {total}", file=sys.stderr, flush=True)
        rid = row.get("row_id", "")
        smi = (row.get("smiles") or "").strip()
        rec = {"row_id": rid, "status": ""}
        if not smi:
            rec["status"] = "empty"
            results.append(rec)
            continue
        try:
            sites, mol_block = predict(smi, args.token, args.timeout)
        except urllib.error.HTTPError as exc:
            rec.update(status="http_error", message=f"HTTP {exc.code}")
            results.append(rec)
            continue
        except Exception as exc:  # noqa: BLE001
            rec.update(status="error", message=str(exc)[:140])
            results.append(rec)
            continue

        acids = [s["pka"] for s in sites if s["site_type"] == "acid"]
        bases = [s["pka"] for s in sites if s["site_type"] == "base"]
        rec.update(
            status="ok",
            sites=sites,
            mol_block=mol_block,
            pka_acid_min=round(min(acids), 3) if acids else None,
            pka_base_max=round(max(bases), 3) if bases else None,
        )
        results.append(rec)

    write_json(args.out, results)
    ok = sum(1 for r in results if r.get("status") == "ok")
    print(f"molgpka: {ok}/{len(results)} ok -> {args.out}", file=sys.stderr)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:  # noqa: BLE001
        fail(str(exc))
