# -*- coding: utf-8 -*-
"""RDKit 插件：去盐/中和、母体信息、17 个描述符。

独立用法：
    python rdkit_plugin.py --in rows.csv --out rdkit.json [--rdkit-path DIR]
输入 CSV 列：row_id,smiles
"""
import argparse
import sys

from _common import RDKIT_LIBS, enable_rdkit, fail, read_rows, write_json

enable_rdkit(RDKIT_LIBS)

from rdkit import Chem, RDLogger  # noqa: E402
from rdkit.Chem import Descriptors, rdMolDescriptors  # noqa: E402
from rdkit.Chem.MolStandardize import rdMolStandardize  # noqa: E402

RDLogger.DisableLog("rdApp.*")


def largest_fragment(mol):
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    return max(frags, key=lambda m: (m.GetNumHeavyAtoms(),
                                     sum(1 for a in m.GetAtoms() if a.GetAtomicNum() == 6)))


def parent_molecule(mol):
    """多片段时取重原子最多的片段，并尽量中和。"""
    frag = largest_fragment(mol)
    try:
        frag = rdMolStandardize.Uncharger().uncharge(frag)
    except Exception:  # noqa: BLE001
        pass
    return frag


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", dest="out", required=True)
    ap.add_argument("--rdkit-path", default=None)
    ap.add_argument("--progress", action="store_true")
    args = ap.parse_args()
    enable_rdkit(args.rdkit_path)

    rd_rows = read_rows(args.inp)
    total = len(rd_rows)
    results = []
    for index, row in enumerate(rd_rows, start=1):
        if args.progress:
            print(f"PROGRESS {index} {total}", file=sys.stderr, flush=True)
        rid = row.get("row_id", "")
        smi = (row.get("smiles") or "").strip()
        rec = {"row_id": rid, "input_smiles": smi}
        mol = Chem.MolFromSmiles(smi) if smi else None
        if mol is None:
            rec.update(status="invalid_smiles")
            results.append(rec)
            continue

        n_frag = len(Chem.GetMolFrags(mol))
        is_multi = n_frag > 1
        parent = mol
        if is_multi:
            parent = parent_molecule(mol)

        carbons = sum(1 for a in parent.GetAtoms() if a.GetAtomicNum() == 6)
        inorganic = carbons == 0

        rec.update(
            status="ok",
            is_salt=is_multi,
            n_fragments=n_frag,
            is_inorganic=inorganic,
            structure_basis="parent_neutral" if is_multi else "as_input",
            parent_smiles=Chem.MolToSmiles(parent),
            parent_formula=rdMolDescriptors.CalcMolFormula(parent),
            parent_mw=round(Descriptors.MolWt(parent), 4),
            parent_inchikey=Chem.MolToInchiKey(parent),
            num_atoms=parent.GetNumAtoms(),
            num_heteroatoms=rdMolDescriptors.CalcNumHeteroatoms(parent),
        )

        if inorganic:
            rec["note"] = "inorganic: prediction descriptors skipped"
            results.append(rec)
            continue

        rec.update(
            tpsa=round(Descriptors.TPSA(parent), 2),
            num_rings=rdMolDescriptors.CalcNumRings(parent),
            num_aromatic_rings=rdMolDescriptors.CalcNumAromaticRings(parent),
            num_aliphatic_rings=rdMolDescriptors.CalcNumAliphaticRings(parent),
            num_aromatic_carbocycles=rdMolDescriptors.CalcNumAromaticCarbocycles(parent),
            num_aliphatic_carbocycles=rdMolDescriptors.CalcNumAliphaticCarbocycles(parent),
            num_aromatic_heterocycles=rdMolDescriptors.CalcNumAromaticHeterocycles(parent),
            num_aliphatic_heterocycles=rdMolDescriptors.CalcNumAliphaticHeterocycles(parent),
            num_heterocycles=rdMolDescriptors.CalcNumHeterocycles(parent),
            lipinski_hbd=rdMolDescriptors.CalcNumLipinskiHBD(parent),
            lipinski_hba=rdMolDescriptors.CalcNumLipinskiHBA(parent),
            labute_asa=round(Descriptors.LabuteASA(parent), 3),
            fraction_csp3=round(rdMolDescriptors.CalcFractionCSP3(parent), 4),
        )
        results.append(rec)

    write_json(args.out, results)
    ok = sum(1 for r in results if r.get("status") == "ok")
    print(f"rdkit: {ok}/{len(results)} ok -> {args.out}", file=sys.stderr)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:  # noqa: BLE001
        fail(str(exc))
