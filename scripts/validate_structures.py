#!/usr/bin/env python3
"""Validate substrates/structures.yaml with RDKit: sanitise, formula, monoisotopic mass, SMARTS matches.

Usage: python scripts/validate_structures.py [--out data/structures_validation.json]
Exit code 1 if any structure fails.
"""
import argparse, json, sys
from pathlib import Path
import yaml
from rdkit import Chem, RDLogger
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem.Descriptors import ExactMolWt

ROOT = Path(__file__).resolve().parents[1]
MASS_TOL = 0.002  # Da; covers electron-mass handling for ions

def check(entry):
    m = Chem.MolFromSmiles(entry["smiles"])
    if m is None:
        return {"ok": False, "error": "RDKit could not parse/sanitise"}
    res = {"canonical_smiles": Chem.MolToSmiles(m), "formula": CalcMolFormula(m),
           "mono_mass": round(ExactMolWt(m), 5), "inchikey": Chem.MolToInchiKey(m)}
    ok = True
    if entry.get("ref_formula"):
        res["formula_ok"] = res["formula"] == entry["ref_formula"]; ok &= res["formula_ok"]
    if entry.get("ref_mass") is not None:
        res["mass_delta"] = round(res["mono_mass"] - float(entry["ref_mass"]), 5)
        res["mass_ok"] = abs(res["mass_delta"]) <= MASS_TOL; ok &= res["mass_ok"]
    res["ok"] = bool(ok)
    return res

def main():
    RDLogger.DisableLog("rdApp.*")
    ap = argparse.ArgumentParser(); ap.add_argument("--out", default=str(ROOT / "data/structures_validation.json"))
    a = ap.parse_args()
    spec = yaml.safe_load(open(ROOT / "substrates/structures.yaml"))
    out, bad = {}, []
    for sect in ("substrates", "fragments"):
        for k, e in spec[sect].items():
            r = check(e); out[k] = r
            if not r["ok"]: bad.append(k)
            print(f"{k:18s} {'OK ' if r['ok'] else 'FAIL'} {r.get('formula','')} {r.get('mono_mass','')} d={r.get('mass_delta','-')}")
    pats = {k: Chem.MolFromSmarts(v) for k, v in spec["smarts"].items() if ">>" not in v}
    for k, p in pats.items():
        if p is None: bad.append("smarts:" + k)
    out["_smarts_matches"] = {k: {s: len(Chem.MolFromSmiles(e["smiles"]).GetSubstructMatches(p))
                                  for s, e in {**spec["substrates"], **spec["fragments"]}.items()}
                              for k, p in pats.items() if p is not None}
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    json.dump(out, open(a.out, "w"), indent=1)
    if bad:
        print("FAILED:", bad); sys.exit(1)
    print("all structures valid")

if __name__ == "__main__":
    main()
