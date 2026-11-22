#!/usr/bin/env python3
"""RDKit validation of stage-1b species: sanitisation, formula, monoisotopic mass vs PubChem, all-cis triol check."""
import argparse, json, sys
from pathlib import Path
import numpy as np
import yaml
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors

HERE = Path(__file__).resolve().parent


def load():
    return yaml.safe_load((HERE / "species.yaml").read_text())


def plain(smiles):
    m = Chem.MolFromSmiles(smiles)
    for a in m.GetAtoms():
        a.SetAtomMapNum(0)
    return m


def formula(m):
    f = rdMolDescriptors.CalcMolFormula(m)
    q = Chem.GetFormalCharge(m)
    return f.rstrip("+-") + ("" if q == 0 else ("+" if q > 0 else "-"))


def triol_all_cis(smiles, seed=7):
    """True if the three O atoms lie on the same side of the mean ring plane (embedded conformer)."""
    m = Chem.AddHs(plain(smiles))
    AllChem.EmbedMolecule(m, randomSeed=seed)
    X = m.GetConformer().GetPositions()
    ring = [a.GetIdx() for a in m.GetAtoms() if a.IsInRing()]
    c = X[ring].mean(0); _, _, vt = np.linalg.svd(X[ring] - c); nrm = vt[2]
    side = []
    for a in m.GetAtoms():
        if a.GetSymbol() == "O":
            cr = [n.GetIdx() for n in a.GetNeighbors() if n.GetSymbol() == "C"][0]
            hs = [n.GetIdx() for n in m.GetAtomWithIdx(cr).GetNeighbors() if n.GetSymbol() == "H"][0]
            side.append(np.sign(np.dot(X[a.GetIdx()] - X[hs], nrm)))   # O vs geminal H: which face
    return len(set(side)) == 1


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--out")
    a = ap.parse_args()
    S = load(); rows, bad = [], 0
    for group, items in S.items():
        for name, d in items.items():
            m = Chem.MolFromSmiles(d["smiles"])
            ok = m is not None
            r = {"group": group, "name": name, "smiles": d["smiles"], "sanitised": ok}
            if ok:
                p = plain(d["smiles"])
                r["formula"] = formula(p); r["charge"] = Chem.GetFormalCharge(p)
                r["mass"] = round(Descriptors.ExactMolWt(p) + (5.4858e-4 * r["charge"]), 6)  # electron mass correction
                r["formula_ok"] = r["formula"] == d["ref_formula"]
                r["mass_ok"] = abs(r["mass"] - d["ref_mass"]) < 2e-3
                r["charge_ok"] = r["charge"] == d["charge"]
                r["maps"] = sorted(at.GetAtomMapNum() for at in m.GetAtoms() if at.GetAtomMapNum())
                if d.get("stereo") == "all-cis":
                    r["all_cis"] = all(triol_all_cis(d["smiles"], s) for s in (1, 7, 42))
            good = ok and r["formula_ok"] and r["mass_ok"] and r["charge_ok"] and r.get("all_cis", True)
            r["ok"] = bool(good); bad += not good; rows.append(r)
            print(f"{'OK ' if good else 'BAD'} {group:9s} {name:9s} {r.get('formula')} {r.get('mass')} maps={r.get('maps')}"
                  + (f" all_cis={r['all_cis']}" if 'all_cis' in r else ""))
    if a.out:
        Path(a.out).write_text(json.dumps(rows, indent=1))
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
