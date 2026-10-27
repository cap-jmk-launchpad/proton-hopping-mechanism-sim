#!/usr/bin/env python3
"""Stage-0 smoke test (CPU only, seconds): rigid 1D scan of a proton shuttling
between two water molecules in the Zundel-type cation H5O2+ with GFN2-xTB.

This is a TOY model to validate the toolchain (xtb install, parsing, plotting
data), NOT a result about the proposed mechanism. Rigid scans overestimate
barriers; real work uses constrained optimisation / NEB (see README stages).

Usage: XTB=/path/to/xtb python3 scripts/smoke_proton_scan.py --roo 2.6
"""
import argparse, os, subprocess, tempfile, re, json
import numpy as np

HARTREE_KCAL = 627.509

def geom(roo, x):
    """O1 at 0, O2 at roo on z; shared proton at fraction x along O-O.
    Two terminal H on each O in a fixed, roughly tetrahedral pose."""
    oh, ang = 0.97, np.deg2rad(110.0)
    atoms = [("O", 0, 0, 0), ("O", 0, 0, roo), ("H", 0, 0, x * roo)]
    for o_z, sign, phi in ((0.0, -1, 0.0), (roo, +1, np.pi / 2)):
        for k in (0, 1):
            p = phi + k * np.pi
            atoms.append(("H", oh*np.sin(ang)*np.cos(p), oh*np.sin(ang)*np.sin(p),
                          o_z + sign * oh * abs(np.cos(ang))))
    return atoms

def energy(xtb, atoms):
    with tempfile.TemporaryDirectory() as d:
        with open(os.path.join(d, "m.xyz"), "w") as f:
            f.write(f"{len(atoms)}\n\n" + "".join(f"{a} {x:.5f} {y:.5f} {z:.5f}\n" for a, x, y, z in atoms))
        out = subprocess.run([xtb, "m.xyz", "--gfn", "2", "--chrg", "1", "--sp"],
                             cwd=d, capture_output=True, text=True).stdout
    m = re.search(r"TOTAL ENERGY\s+(-?\d+\.\d+)", out)
    if not m:
        raise RuntimeError("xtb failed:\n" + out[-2000:])
    return float(m.group(1))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--roo", type=float, default=2.6, help="O...O distance / Angstrom")
    ap.add_argument("--n", type=int, default=11)
    ap.add_argument("--out", default="data/smoke_scan.json")
    a = ap.parse_args()
    xtb = os.environ.get("XTB", "xtb")
    xs = np.linspace(0.38, 0.62, a.n)
    es = np.array([energy(xtb, geom(a.roo, x)) for x in xs])
    rel = (es - es.min()) * HARTREE_KCAL
    for x, e in zip(xs, rel):
        print(f"x={x:.3f}  dE={e:7.3f} kcal/mol")
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    json.dump({"model": "H5O2+ rigid scan, GFN2-xTB", "roo": a.roo,
               "x": xs.tolist(), "dE_kcal": rel.tolist(),
               "note": "toolchain smoke test only; not a mechanistic result"},
              open(a.out, "w"), indent=1)

if __name__ == "__main__":
    main()
