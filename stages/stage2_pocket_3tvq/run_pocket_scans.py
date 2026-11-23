#!/usr/bin/env python3
"""Enzyme-pocket (3TVQ TcmN ARO/CYC, 241 atoms) reactant relaxation and 1D relaxed scans.
Method: GFN2-xTB/ALPB(water) via tblite, ASE FixAtoms (24 CA + H-caps frozen) + FixInternals(bondcombos) + BFGS.
Strictly OSI open-source tools (ASE, tblite, NumPy).

Usage:
  python run_pocket_scans.py --dir <pocket_work_dir> relax
  python run_pocket_scans.py --dir <pocket_work_dir> scan <channel> [--npts 15]
  python run_pocket_scans.py --dir <pocket_work_dir> build_seed <seed_int>
  python run_pocket_scans.py --dir <pocket_work_dir> build_hip
  python run_pocket_scans.py --dir <pocket_work_dir> scan_variant <tag> <channel> [--npts 15]
  python run_pocket_scans.py --dir <pocket_work_dir> c10_water_relay [--npts 15]
"""
import argparse, json, time
from pathlib import Path
import numpy as np
from ase import Atoms
from ase.io import read, write
from ase.constraints import FixAtoms, FixInternals
from ase.optimize import BFGS
from tblite.ase import TBLite

HERE = Path(__file__).resolve().parent


def load_cluster_meta(work_dir: Path):
    meta_p = work_dir / "cluster/cluster_3TVQ.json"
    if not meta_p.exists():
        meta_p = HERE / "cluster_3TVQ.json"
    return json.load(open(meta_p))


def load_cluster_xyz(work_dir: Path):
    xyz_p = work_dir / "cluster/cluster_3TVQ.xyz"
    if not xyz_p.exists():
        xyz_p = HERE / "cluster_3TVQ.xyz"
    return read(str(xyz_p))


def cmd_relax(work_dir: Path):
    meta = load_cluster_meta(work_dir)
    atoms = load_cluster_xyz(work_dir)
    frozen_0 = [idx - 1 for idx in meta["frozen_1based"]]
    pos0 = atoms.get_positions().copy()

    atoms.set_constraint(FixAtoms(indices=frozen_0))
    atoms.calc = TBLite(method="GFN2-xTB", charge=meta["charge"], multiplicity=1, solvation=("alpb", "water"), verbosity=0)

    outdir = work_dir / "xtb/ase_opt"
    outdir.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    e0 = float(atoms.get_potential_energy())
    f0 = float(np.max(np.linalg.norm(atoms.get_forces(), axis=1)))
    opt = BFGS(atoms, logfile=str(outdir / "opt.log"), maxstep=0.15)
    converged = opt.run(fmax=0.05, steps=200)
    t1 = time.time()

    e1 = float(atoms.get_potential_energy())
    f1 = float(np.max(np.linalg.norm(atoms.get_forces(), axis=1)))
    pos1 = atoms.get_positions()
    max_frozen_drift = float(np.max(np.linalg.norm(pos1[frozen_0] - pos0[frozen_0], axis=1)))
    rmsd_all = float(np.sqrt(np.mean(np.sum((pos1 - pos0) ** 2, axis=1))))

    write(str(outdir / "reactant_relaxed.xyz"), atoms)
    res = {
        "converged": bool(converged),
        "nsteps": int(opt.get_number_of_steps()),
        "wall_s": round(t1 - t0, 2),
        "E0_eV": round(e0, 6),
        "E1_eV": round(e1, 6),
        "dE_kcal": round((e1 - e0) * 23.060548, 2),
        "fmax0": round(f0, 4),
        "fmax1": round(f1, 4),
        "max_frozen_drift_A": round(max_frozen_drift, 6),
        "rmsd_all_A": round(rmsd_all, 3),
    }
    json.dump(res, open(outdir / "step0_result.json", "w"), indent=2)
    print(json.dumps(res, indent=2), flush=True)


def cmd_scan(work_dir: Path, ch_name: str, NP: int = 15):
    meta = load_cluster_meta(work_dir)
    labels = meta["labels"]
    idx = {lab: i for i, lab in enumerate(labels)}
    frozen_0 = [i - 1 for i in meta["frozen_1based"]]

    CHANNELS = {
        "direct_C10_Asp57": [
            (idx["DQH:H24"], idx["DQH:C10"], idx["ASP57:OD2"]),
        ],
        "relay_O27_W259_W266_His128": [
            (idx["DQH:H31"], idx["DQH:O27"], idx["HOH259:O"]),
            (idx["HOH259:H1"], idx["HOH259:O"], idx["HOH266:O"]),
            (idx["HOH266:H2"], idx["HOH266:O"], idx["HIS128:NE2"]),
        ],
        "relay_O27_W266_His128": [
            (idx["DQH:H31"], idx["DQH:O27"], idx["HOH266:O"]),
            (idx["HOH266:H2"], idx["HOH266:O"], idx["HIS128:NE2"]),
        ],
    }
    mv = CHANNELS[ch_name]
    atoms = read(str(work_dir / "xtb/ase_opt/reactant_relaxed.xyz"))

    if ch_name.startswith("relay_O27"):
        X = atoms.get_positions()
        h0, d0, a0 = mv[0]
        u = (X[a0] - X[d0]) / np.linalg.norm(X[a0] - X[d0])
        X[h0] = X[d0] + 0.98 * u
        h266_1, o266, o259 = idx["HOH266:H1"], idx["HOH266:O"], idx["HOH259:O"]
        u2 = (X[o259] - X[o266]) / np.linalg.norm(X[o259] - X[o266])
        X[h266_1] = X[o266] + 0.97 * u2
        atoms.set_positions(X)
        atoms.set_constraint(FixAtoms(indices=frozen_0))
        atoms.calc = TBLite(method="GFN2-xTB", charge=meta["charge"], multiplicity=1, solvation=("alpb", "water"), verbosity=0, max_iterations=500)
        BFGS(atoms, logfile=None, maxstep=0.1).run(fmax=0.05, steps=80)

    X0 = atoms.get_positions().copy()
    calc = TBLite(method="GFN2-xTB", charge=meta["charge"], multiplicity=1, solvation=("alpb", "water"), verbosity=0, max_iterations=500)
    atoms.calc = calc
    E_ref = float(atoms.get_potential_energy())

    q_fn = lambda pos: float(sum(np.linalg.norm(pos[h] - pos[d]) - np.linalg.norm(pos[h] - pos[a]) for h, d, a in mv))
    q0, q1 = q_fn(X0), 0.85 * len(mv)
    outdir = work_dir / "xtb/scans" / ch_name
    outdir.mkdir(parents=True, exist_ok=True)
    write(str(outdir / "start_reactant.xyz"), atoms)

    targets = np.linspace(q0, q1, NP)
    combo = [[d, h, 1.0] for h, d, a in mv] + [[h, a, -1.0] for h, d, a in mv]
    res = []
    for k, t in enumerate(targets):
        atoms.set_constraint()
        X = atoms.get_positions()
        delta = (t - q_fn(X)) / (2.0 * len(mv))
        for h, d, a in mv:
            u = X[a] - X[d]
            X[h] += delta * u / np.linalg.norm(u)
        atoms.set_positions(X)
        atoms.set_constraint([FixAtoms(indices=frozen_0), FixInternals(bondcombos=[[q_fn(X), combo]])])
        atoms.calc = calc
        t0 = time.time()
        opt = BFGS(atoms, logfile=str(outdir / f"opt_{k:02d}.log"), maxstep=0.1)
        try:
            conv = bool(opt.run(fmax=0.04, steps=120))
            E = float(atoms.get_potential_energy())
            err = None
        except Exception as e:
            conv, E, err = False, None, repr(e)[:200]
            calc = TBLite(method="GFN2-xTB", charge=meta["charge"], multiplicity=1, solvation=("alpb", "water"), verbosity=0, max_iterations=500)
        Xn = atoms.get_positions()
        dE_kcal = round((E - E_ref) * 23.060548, 3) if E is not None else None
        hop_d = [[round(float(np.linalg.norm(Xn[h] - Xn[d])), 3), round(float(np.linalg.norm(Xn[h] - Xn[a])), 3)] for h, d, a in mv]
        pt = {
            "k": k, "target_q": round(float(t), 4), "q": round(q_fn(Xn), 4),
            "E_eV": round(E, 6) if E is not None else None, "dE_kcal": dE_kcal,
            "converged": conv, "steps": int(opt.nsteps), "wall_s": round(time.time() - t0, 2),
            "hop_d_DH_HA": hop_d, "error": err,
        }
        res.append(pt)
        write(str(outdir / f"pt_{k:02d}.xyz"), atoms)
        json.dump({"channel": ch_name, "E_ref_eV": E_ref, "q0": q0, "q1": q1, "points": res}, open(outdir / "scan.json", "w"), indent=2)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default=".", help="Working directory for pocket runs")
    ap.add_argument("cmd", choices=["relax", "scan"])
    ap.add_argument("channel", nargs="?")
    ap.add_argument("--npts", type=int, default=15)
    args = ap.parse_args()
    wdir = Path(args.dir).resolve()
    if args.cmd == "relax":
        cmd_relax(wdir)
    elif args.cmd == "scan":
        cmd_scan(wdir, args.channel, args.npts)
