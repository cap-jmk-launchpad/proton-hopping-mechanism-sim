#!/usr/bin/env python3
"""Stage 1: relay vs direct proton transfer, GFN2-xTB (tblite) + ALPB, endpoint optimisation + CI-NEB.

Each job: build cluster -> optimise reactant R and product P -> classify endpoints -> CI-NEB R->P.
Energies are electronic energies (GFN2-xTB with ALPB implicit solvent); no ZPE, no entropy.
Results: <out>/<job>/result.json (+ R.xyz, P.xyz, band.xyz). Existing results are skipped (checkpointing).

Usage:
  python run_stage1.py --out RUN_DIR [--workers 12] [--seeds 0 1 2] [--only-n 0 1] [--limit K]
"""
import argparse, json, os, signal, sys, time, traceback
from itertools import product
from pathlib import Path
from multiprocessing import Pool

os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from ase import Atoms
from ase.io import write
from ase.optimize import BFGS, FIRE
try:
    from ase.mep import NEB
except ImportError:  # older ASE
    from ase.neb import NEB
from tblite.ase import TBLite

sys.path.insert(0, str(Path(__file__).resolve().parent))
import build  # noqa: E402

EV2KCAL = 23.060548
STEPS = ("deprot", "reprot")
SINKS = ("hydroxide", "methoxide", "nme3")
SOLVENTS = ("water", "methanol")
N_IMAGES = 11
JOB_TIMEOUT_S = 1200      # wall-clock cap per job; a capped job is recorded as status "timeout"


class JobTimeout(BaseException):   # BaseException so NEB retry loops do not swallow it
    pass


def _alarm(signum, frame):
    raise JobTimeout()


def jobs(seeds, only_n=None):
    out = []
    for n in range(5):
        if only_n is not None and n not in only_n:
            continue
        for step, sink, solv, seed in product(STEPS, SINKS, SOLVENTS, seeds):
            tops = ("relay",) if n == 0 else ("relay", "direct")   # n = 0: relay == direct
            for top in tops:
                out.append(dict(step=step, sink=sink, solvent=solv, n=n, topology=top, seed=seed))
    return out


def job_id(j):
    return f"{j['step']}_{j['sink']}_{j['solvent']}_n{j['n']}_{j['topology']}_s{j['seed']}"


class RobustTBLite(TBLite):
    """TBLite that retries a failed SCF with a fresh calculator (no cached wavefunction), then with an EEQ guess.
    The Hamiltonian is unchanged, so converged energies are the same; only the SCF starting point differs."""
    retries = 0

    def calculate(self, atoms=None, properties=None, system_changes=None):
        properties = sorted(set(properties or []) | {"energy", "forces"})   # retries must keep forces
        try:
            return super().calculate(atoms, properties, system_changes or ["positions"])
        except Exception as first:
            RobustTBLite.retries += 1
            params = dict(self.parameters)
            for guess in ("sad", "eeq"):
                fresh = TBLite(**{**params, "guess": guess})
                a = atoms.copy(); a.calc = fresh
                try:
                    e, f = a.get_potential_energy(), a.get_forces()
                except Exception:
                    continue
                self.atoms = atoms.copy()
                self.results = {"energy": e, "free_energy": e, "forces": f}
                return
            raise first


def calc(charge, solvent, etemp=300.0):
    return RobustTBLite(method="GFN2-xTB", charge=charge, solvation=("alpb", solvent), verbosity=0,
                        max_iterations=300, mixer_damping=0.2, electronic_temperature=etemp)


# NEB attempts in order; the first that runs without an SCF failure is kept and recorded in the result
NEB_ATTEMPTS = (("idpp", 300.0), ("linear", 300.0), ("idpp", 1000.0))


def run_neb(R, P, q, solvent, interp, etemp):
    images = [R.copy()] + [R.copy() for _ in range(N_IMAGES - 2)] + [P.copy()]
    for im in images:
        im.calc = calc(q, solvent, etemp)
    neb = NEB(images, climb=False, method="improvedtangent", allow_shared_calculator=False)
    neb.interpolate(interp)
    opt = FIRE(neb, logfile=None, maxstep=0.05)       # plain NEB first, small steps
    opt.run(fmax=0.10, steps=800); n1 = int(opt.nsteps)
    neb.climb = True                                  # then climbing image
    opt = FIRE(neb, logfile=None, maxstep=0.05)
    conv = opt.run(fmax=0.05, steps=1200)
    E = np.array([float(im.get_potential_energy()) for im in images])
    return images, E, bool(conv), n1 + int(opt.nsteps)


def relax(atoms, charge, solvent, fmax=0.03, steps=800):
    atoms.calc = calc(charge, solvent)
    opt = BFGS(atoms, logfile=None, maxstep=0.1)
    conv = opt.run(fmax=fmax, steps=steps)
    return float(atoms.get_potential_energy()), bool(conv), int(opt.nsteps)


def classify(pos, moving):
    d = lambda i, k: float(np.linalg.norm(pos[i] - pos[k]))
    on_donor = all(d(h, a_) > 1.30 and d(h, dn) < 1.25 for h, dn, a_ in moving)
    on_acceptor = all(d(h, a_) < 1.25 for h, dn, a_ in moving)
    h0, d0, a0 = moving[0]                  # primary proton: the C-H (deprot) or B-H (reprot) proton
    return {"all_on_donor": on_donor, "all_on_acceptor": on_acceptor,
            "primary_on_donor": d(h0, d0) < 1.25, "primary_off_donor": d(h0, d0) > 1.60,
            "d_H_donor": [round(d(h, dn), 3) for h, dn, a_ in moving],
            "d_H_acceptor": [round(d(h, a_), 3) for h, dn, a_ in moving]}


def run_job(args):
    j, root = args
    jid = job_id(j); od = Path(root) / jid
    if (od / "result.json").exists():
        return jid, "skip"
    od.mkdir(parents=True, exist_ok=True)
    t0 = time.time(); res = dict(j, job=jid, level="GFN2-xTB/ALPB(%s) via tblite" % j["solvent"])
    RobustTBLite.retries = 0
    signal.signal(signal.SIGALRM, _alarm); signal.alarm(JOB_TIMEOUT_S)
    try:
        s = build.build(j["step"], j["sink"], j["solvent"], j["n"], j["topology"], seed=j["seed"])
        q = s["charge"]; res.update(charge=q, natoms=len(s["symbols"]), moving=s["moving"])
        R = Atoms(s["symbols"], positions=s["positions"])
        P = Atoms(s["symbols"], positions=build.product_guess(s))
        eR, cR, nR = relax(R, q, j["solvent"]); eP, cP, nP = relax(P, q, j["solvent"])
        write(od / "R.xyz", R); write(od / "P.xyz", P)
        kR, kP = classify(R.positions, s["moving"]), classify(P.positions, s["moving"])
        res.update(E_R_eV=eR, E_P_eV=eP, R_converged=cR, P_converged=cP, R_steps=nR, P_steps=nP,
                   R_class=kR, P_class=kP, dE_rxn_kcal=round((eP - eR) * EV2KCAL, 3))
        # wire rearrangements (e.g. an OH- defect hopping along the wire) are recorded, not rejected
        res["wire_rearranged_R"] = not kR["all_on_donor"]; res["wire_rearranged_P"] = not kP["all_on_acceptor"]
        if not kR["primary_on_donor"]:
            res["status"] = "reactant_not_minimum"      # primary proton left its donor during R relaxation
        elif not kP["primary_off_donor"]:
            res["status"] = "product_not_minimum"       # primary proton returned during P relaxation
        else:
            fails = []
            for interp, etemp in NEB_ATTEMPTS:
                try:
                    images, E, conv, nst = run_neb(R, P, q, j["solvent"], interp, etemp)
                    break
                except Exception as e:  # SCF failure on an interpolated image
                    fails.append(f"{interp}/{etemp}: {e!r}")
            else:
                raise RuntimeError("all NEB attempts failed: " + " | ".join(fails))
            write(od / "band.xyz", images)
            imax = int(E.argmax())
            res.update(neb_interp=interp, neb_etemp_K=etemp, neb_failed_attempts=fails)
            res.update(status="ok" if conv else "neb_not_converged", neb_steps=nst,
                       E_band_eV=E.tolist(), i_max=imax,
                       barrier_fwd_kcal=round(float(E[imax] - E[0]) * EV2KCAL, 3),
                       barrier_rev_kcal=round(float(E[imax] - E[-1]) * EV2KCAL, 3),
                       ts_is_endpoint=imax in (0, len(E) - 1))
    except JobTimeout:
        res.update(status="timeout", error=f"job exceeded {JOB_TIMEOUT_S} s")
    except Exception as e:
        res.update(status="error", error=repr(e), tb=traceback.format_exc()[-1500:])
    finally:
        signal.alarm(0)
    res["scf_retries"] = RobustTBLite.retries
    res["wall_s"] = round(time.time() - t0, 1)
    json.dump(res, open(od / "result.json", "w"), indent=1)
    return jid, res["status"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True); ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2])
    ap.add_argument("--only-n", type=int, nargs="+"); ap.add_argument("--limit", type=int)
    ap.add_argument("--jobs", nargs="+", help="run only these job ids (e.g. to retry failed jobs)")
    a = ap.parse_args()
    J = jobs(a.seeds, a.only_n)
    if a.jobs:
        J = [j for j in J if job_id(j) in set(a.jobs)]
    if a.limit:
        J = J[:a.limit]
    Path(a.out).mkdir(parents=True, exist_ok=True)
    print(f"{len(J)} jobs, {a.workers} workers", flush=True)
    t0 = time.time(); done = 0
    with Pool(a.workers) as pool:
        for jid, st in pool.imap_unordered(run_job, [(j, a.out) for j in J]):
            done += 1
            print(f"[{done}/{len(J)}] {time.time() - t0:8.0f}s {jid} {st}", flush=True)
    (Path(a.out) / "DONE").write_text(f"{len(J)} jobs finished {time.strftime('%Y-%m-%d %H:%M:%S %Z')}\n")
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
