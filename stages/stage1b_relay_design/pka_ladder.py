#!/usr/bin/env python3
"""Relative acidity ladder for the stage-1b donors/shuttles/bases (lever 2: pKa matching), open-source only.

For each acid HA: relax HA and A- (GFN2-xTB/ALPB, tblite) and report the isodesmic proton-exchange energy
  dE(HA) = [E(A-) + E(H2O)] - [E(HA) + E(OH-)]      (kcal/mol, relative to water; electronic, no ZPE/entropy)
Smaller dE = stronger acid. Semi-empirical and implicit-solvent: an ORDERING tool for choosing pKa-matched
shuttles, not a pKa prediction. Literature pKa values are only quoted in the design doc with a checked source.

Usage: python pka_ladder.py --solvent water --out ladder_water.json   (seconds per species; part of the screen)
"""
import argparse, json, os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

PAIRS = {   # name: (acid SMILES, conjugate base SMILES)
    "water": ("O", "[OH-]"), "methanol": ("CO", "C[O-]"), "acac_CH": ("CC(=O)CC(C)=O", "CC(=O)[CH-]C(C)=O"),
    "dienone_CH": ("O=C1CC=CC=C1", "[O-]c1ccccc1"), "acoh": ("CC(=O)O", "CC(=O)[O-]"),
    "phenol": ("Oc1ccccc1", "[O-]c1ccccc1"), "pyridone_NH": ("O=C1C=CC=CN1", "[O-]c1ccccn1"),
    "hfip": ("OC(C(F)(F)F)C(F)(F)F", "[O-]C(C(F)(F)F)C(F)(F)F"), "triol_OH": ("O[C@H]1C[C@@H](O)C[C@@H](O)C1",
                                                                        "[O-][C@H]1C[C@@H](O)C[C@@H](O)C1"),
    "hnme3": ("C[NH+](C)C", "CN(C)C"), "hnmi": ("C[n+]1cc[nH]c1", "Cn1ccnc1"), "hq7_OH": ("Oc1ccc2cccnc2c1",
                                                                                     "[O-]c1ccc2cccnc2c1"),
}


def energy(smiles, solvent, seed=0):
    from ase import Atoms
    from ase.optimize import BFGS
    from tblite.ase import TBLite
    m = Chem.AddHs(Chem.MolFromSmiles(smiles)); q = Chem.GetFormalCharge(m)
    if m.GetNumAtoms() > 1:
        AllChem.EmbedMolecule(m, randomSeed=seed); AllChem.MMFFOptimizeMolecule(m)
        X = m.GetConformer().GetPositions()
    else:
        X = np.zeros((1, 3))
    A = Atoms([a.GetSymbol() for a in m.GetAtoms()], positions=X)
    A.calc = TBLite(method="GFN2-xTB", charge=q, solvation=("alpb", solvent), verbosity=0)
    BFGS(A, logfile=None).run(fmax=0.02, steps=500)
    return float(A.get_potential_energy()) * 23.060548


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--solvent", default="water"); ap.add_argument("--out")
    a = ap.parse_args()
    E = {k: (energy(ha, a.solvent), energy(b, a.solvent)) for k, (ha, b) in PAIRS.items()}
    ref = E["water"][1] - E["water"][0]
    rows = sorted(((k, round((eb - eh) - ref, 2)) for k, (eh, eb) in E.items()), key=lambda r: r[1])
    for k, v in rows:
        print(f"{k:12s} {v:8.2f}")
    if a.out:
        json.dump({"solvent": a.solvent, "level": "GFN2-xTB/ALPB (tblite)", "dE_rel_water_kcal": dict(rows)},
                  open(a.out, "w"), indent=1)


if __name__ == "__main__":
    main()
