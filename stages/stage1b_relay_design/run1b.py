#!/usr/bin/env python3
"""Stage 1b: relay-design screen (GFN2-xTB/ALPB via tblite, ASE CI-NEB), families A, B, C, D and control Q.

Per job: build -> (optional FixBondLength C...base at d) -> relax R and P -> classify -> CI-NEB ->
partial Hessians over the moving protons at R and at the highest image -> H/D KIE (ZPE-only, all-D and C-H-only),
Wigner/Bell/Eckart tunnelling factors. Optional CREST conformer search of the relaxed reactant (--crest).
Energies are GFN2-xTB electronic energies (no ZPE/entropy in the barrier). Screening only; DFT follows for hits.
Results: <out>/<job>/result.json (+ R.xyz, P.xyz, band.xyz). Existing results are skipped; nothing is deleted.

Usage:
  python run1b.py --list [--family A B C Q]              # print job ids and counts, run nothing
  python run1b.py --out DIR --jobs JOB_ID [--crest]       # run selected jobs (the pilot uses exactly one)
  python run1b.py --out DIR --family A B C Q --workers 8   # full screen (ONLY after Julian's go)
"""
import argparse, json, os, shutil, signal, subprocess, sys, time, traceback
from itertools import product
from pathlib import Path
from multiprocessing import Pool

os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent / "stage1_model_relays"))
import build1b  # noqa: E402
import kinetics  # noqa: E402

D_GRID = (3.0, 4.0, 5.0, 6.0, 7.0, None)        # None = natural (unconstrained) geometry
SHUTTLE_SETS = {"none": [], "w1": ["water"], "w2": ["water", "water"], "m1": ["methanol"], "acoh": ["acoh"],
                "pyr": ["pyridone"], "phoh": ["phenol"], "hfip": ["hfip"]}
JOB_TIMEOUT_S = 2400
N_IMAGES = 11
EV2KCAL = 23.060548


def jobs(families=("A", "B", "C", "Q"), seeds_a=(0, 1, 2), seeds=(0, 1, 2)):
    J = []
    if "A" in families:
        for base, solv, n, d, s in product(("nme3", "nmi", "acetate"), ("water", "methanol"), range(4), D_GRID, seeds_a):
            J.append(dict(family="A", base=base, solvent=solv, n=n, d=d, topology="direct" if n == 0 else "relay", seed=s))
    if "B" in families:
        for sub, key, solv, s in product(("acac", "dienone"), SHUTTLE_SETS, ("water", "methanol"), seeds):
            for top in (("relay",) if key == "none" else ("relay", "stepwise")):
                J.append(dict(family="B", substrate=sub, shuttles=key, solvent=solv, topology=top, seed=s))
    if "C" in families:
        for base, top, solv, s in product(("nme3", "acetate"), ("relay", "direct"), ("water", "methanol"), seeds):
            J.append(dict(family="C", base=base, topology=top, solvent=solv, seed=s))
    if "Q" in families:
        for solv, n, s in product(("water", "methanol"), range(4), seeds):
            J.append(dict(family="Q", solvent=solv, n=n, topology="direct" if n == 0 else "relay", seed=s))
    return J


def job_id(j):
    f = j["family"]
    if f == "A":
        d = "nat" if j["d"] is None else f"{j['d']:.1f}"
        core = f"{j['base']}_{j['solvent']}_n{j['n']}_d{d}"
    elif f == "B":
        core = f"{j['substrate']}_{j['shuttles']}_{j['solvent']}"
    elif f == "C":
        core = f"triol_{j['base']}_{j['solvent']}"
    else:
        core = f"hq7_{j['solvent']}_n{j['n']}"
    sh = f"_sh{j['shell']}" if j.get("shell") else ""
    return f"{f}_{core}_{j['topology']}{sh}_s{j['seed']}"


def parse_job_id(jid):
    """Inverse of job_id for the standard grid plus optional _sh<k> (family D = any job id with _sh<k>)."""
    for j in jobs():
        base = job_id(j)
        if jid == base:
            return j
        stem, seed = base.rsplit("_s", 1)
        if jid.startswith(stem + "_sh") and jid.endswith("_s" + seed):
            k = int(jid[len(stem) + 3:].split("_s")[0])
            return dict(j, shell=k)
    raise KeyError(jid)


def build_job(j):
    sh = j.get("shell", 0)
    if j["family"] == "A":
        return build1b.build_tethered(j["base"], j["solvent"], j["n"], d=j["d"], seed=j["seed"], shell=sh)
    if j["family"] == "B":
        return build1b.build_shuttle(j["substrate"], SHUTTLE_SETS[j["shuttles"]], j["topology"], seed=j["seed"],
                                     shell=sh, shell_solvent=j["solvent"])
    if j["family"] == "C":
        return build1b.build_triol(j["base"], j["topology"], seed=j["seed"], solvent=j["solvent"], shell=sh)
    return build1b.build_hq7(j["solvent"], j["n"], seed=j["seed"])


class JobTimeout(BaseException):
    pass


def _alarm(signum, frame):
    raise JobTimeout()


def run_neb(R, P, q, solvent, interp, etemp, calc):
    from ase.mep import NEB
    from ase.optimize import FIRE
    images = [R.copy()] + [R.copy() for _ in range(N_IMAGES - 2)] + [P.copy()]
    for im in images:
        im.calc = calc(q, solvent, etemp)
    neb = NEB(images, climb=False, method="improvedtangent", allow_shared_calculator=False)
    neb.interpolate(interp, apply_constraint=True)
    opt = FIRE(neb, logfile=None, maxstep=0.05); opt.run(fmax=0.10, steps=800); n1 = int(opt.nsteps)
    neb.climb = True
    opt = FIRE(neb, logfile=None, maxstep=0.05); conv = opt.run(fmax=0.05, steps=1200)
    E = np.array([float(im.get_potential_energy()) for im in images])
    return images, E, bool(conv), n1 + int(opt.nsteps)


def partial_hessian(atoms, idx, delta=0.005):
    """Central-difference Hessian block (eV/A^2) over atoms idx, unconstrained forces."""
    m = len(idx); Hh = np.zeros((3 * m, 3 * m)); X0 = atoms.get_positions().copy()
    for a, i in enumerate(idx):
        for c in range(3):
            fs = []
            for sgn in (1, -1):
                X = X0.copy(); X[i, c] += sgn * delta
                atoms.set_positions(X, apply_constraint=False)
                F = atoms.get_forces(apply_constraint=False)
                fs.append(F[idx].reshape(-1))
            Hh[3 * a + c] = -(fs[0] - fs[1]) / (2 * delta)
    atoms.set_positions(X0, apply_constraint=False)
    return 0.5 * (Hh + Hh.T)


def kinetics_block(R, TS, idx, primary, V1, V2):
    HR, HT = partial_hessian(R, idx), partial_hessian(TS, idx)
    mH = [kinetics.M_H] * len(idx)
    mDall = [kinetics.M_D] * len(idx)
    mDp = [kinetics.M_D if i == primary else kinetics.M_H for i in idx]
    f = kinetics.partial_hessian_freqs
    wRH, wTH = f(HR, mH), f(HT, mH)
    out = {"nu_imag_H_cm": float(wTH[0]) if wTH[0] < 0 else None, "ts_partial_n_imag": int((wTH < 0).sum())}
    for tag, mD in (("allD", mDall), ("primaryD", mDp)):
        wRD, wTD = f(HR, mD), f(HT, mD)
        kie = kinetics.kie_semiclassical(wRH, wTH, wRD, wTD)
        out[f"kie_sc_{tag}"] = round(kie, 3)
        if wTH[0] < 0 and wTD[0] < 0:
            nH, nD = abs(wTH[0]), abs(wTD[0])
            out[f"nu_imag_D_{tag}_cm"] = round(float(wTD[0]), 1)
            out[f"kie_wigner_{tag}"] = round(kie * kinetics.wigner(nH) / kinetics.wigner(nD), 3)
            eH, eD = kinetics.eckart(nH, V1, V2), kinetics.eckart(nD, V1, V2)
            out[f"kappa_eckart_H"] = round(eH, 3) if np.isfinite(eH) else None
            out[f"kie_eckart_{tag}"] = round(kie * eH / eD, 3) if np.isfinite(eH) and np.isfinite(eD) else None
    out["kappa_wigner_H"] = round(kinetics.wigner(abs(wTH[0])), 3) if wTH[0] < 0 else None
    return out


def run_crest(od, R, q, solvent, threads, timeout_s, moving):
    """CREST (NCI mode, quick) on the relaxed reactant; counts relay-competent conformers (all links intact)."""
    exe = shutil.which("crest")
    if exe is None:
        return {"crest": "not_found"}
    cd = od / "crest"; cd.mkdir(exist_ok=True)
    from ase.io import write
    write(cd / "R.xyz", R)
    cmd = [exe, "R.xyz", "--gfn2", "--alpb", solvent, "--chrg", str(q), "--nci", "--quick", "-T", str(threads)]
    t0 = time.time()
    try:
        with open(cd / "crest.log", "w") as lg:
            rc = subprocess.run(cmd, cwd=cd, stdout=lg, stderr=subprocess.STDOUT, timeout=timeout_s).returncode
    except subprocess.TimeoutExpired:
        return {"crest": "timeout", "crest_wall_s": round(time.time() - t0, 1)}
    res = {"crest": "ok" if rc == 0 else f"rc{rc}", "crest_wall_s": round(time.time() - t0, 1), "crest_threads": threads}
    f = cd / "crest_conformers.xyz"
    if f.exists():
        lines = f.read_text().splitlines(); confs, E, i = [], [], 0
        while i < len(lines) and lines[i].strip():
            nat = int(lines[i]); E.append(float(lines[i + 1].split()[0]))
            confs.append(np.array([[float(x) for x in l.split()[1:4]] for l in lines[i + 2:i + 2 + nat]]))
            i += nat + 2
        E = np.array(E); rel = (E - E.min()) * 627.509474
        intact = np.array([all(np.linalg.norm(X[h] - X[a]) < 2.6 for h, d, a in moving) for X in confs])
        win = rel <= 3.0
        res.update(crest_n_conf=len(confs), crest_n_within_3kcal=int(win.sum()),
                   crest_relay_competent_within_3kcal=int((intact & win).sum()),
                   crest_lowest_is_relay_competent=bool(intact[int(E.argmin())]))
    return res


def run_job(args):
    j, root, crest_threads, crest_timeout = args
    from ase import Atoms
    from ase.io import write
    from ase.constraints import FixBondLength
    import run_stage1 as s1
    jid = job_id(j); od = Path(root) / jid
    if (od / "result.json").exists():
        return jid, "skip"
    od.mkdir(parents=True, exist_ok=True)
    t0 = time.time(); c0 = time.process_time()
    res = dict(j, job=jid, level=f"GFN2-xTB/ALPB({j['solvent']}) via tblite")
    s1.RobustTBLite.retries = 0
    signal.signal(signal.SIGALRM, _alarm); signal.alarm(JOB_TIMEOUT_S)
    try:
        s = build_job(j)
        q = s["charge"]; res.update(charge=q, natoms=len(s["symbols"]), moving=s["moving"],
                                    d_natural=s.get("d_natural"))
        R = Atoms(s["symbols"], positions=s["positions"])
        P = Atoms(s["symbols"], positions=build1b.product_guess(s))
        if s.get("constraint"):
            for A in (R, P):
                A.set_constraint(FixBondLength(*s["constraint"]))
        eR, cR, nR = s1.relax(R, q, j["solvent"]); eP, cP, nP = s1.relax(P, q, j["solvent"])
        write(od / "R.xyz", R); write(od / "P.xyz", P)
        kR, kP = s1.classify(R.positions, s["moving"]), s1.classify(P.positions, s["moving"])
        res.update(E_R_eV=eR, E_P_eV=eP, R_converged=cR, P_converged=cP, R_steps=nR, P_steps=nP,
                   R_class=kR, P_class=kP, dE_rxn_kcal=round((eP - eR) * EV2KCAL, 3))
        res["wire_rearranged_R"] = not kR["all_on_donor"]; res["wire_rearranged_P"] = not kP["all_on_acceptor"]
        if not kR["primary_on_donor"]:
            res["status"] = "reactant_not_minimum"
        elif not kP["primary_off_donor"]:
            res["status"] = "product_not_minimum"
        else:
            fails = []
            for interp, etemp in s1.NEB_ATTEMPTS:
                try:
                    images, E, conv, nst = run_neb(R, P, q, j["solvent"], interp, etemp, s1.calc)
                    break
                except Exception as e:
                    fails.append(f"{interp}/{etemp}: {e!r}")
            else:
                raise RuntimeError("all NEB attempts failed: " + " | ".join(fails))
            write(od / "band.xyz", images)
            imax = int(E.argmax())
            bf, br = float(E[imax] - E[0]) * EV2KCAL, float(E[imax] - E[-1]) * EV2KCAL
            res.update(neb_interp=interp, neb_etemp_K=etemp, neb_failed_attempts=fails,
                       status="ok" if conv else "neb_not_converged", neb_steps=nst, E_band_eV=E.tolist(),
                       i_max=imax, barrier_fwd_kcal=round(bf, 3), barrier_rev_kcal=round(br, 3),
                       ts_is_endpoint=imax in (0, len(E) - 1))
            if not res["ts_is_endpoint"]:
                idx = sorted({h for h, _, _ in s["moving"]})
                t1 = time.time()
                res["kinetics"] = kinetics_block(images[0], images[imax], idx, s["moving"][0][0], bf, br)
                res["kinetics_wall_s"] = round(time.time() - t1, 1)
        if crest_threads and res.get("status") in ("ok", "neb_not_converged"):
            res.update(run_crest(od, R, q, j["solvent"], crest_threads, crest_timeout, s["moving"]))
    except build1b.InfeasibleGeometry as e:
        res.update(status="geometry_infeasible", error=str(e))
    except JobTimeout:
        res.update(status="timeout", error=f"job exceeded {JOB_TIMEOUT_S} s")
    except Exception as e:
        res.update(status="error", error=repr(e), tb=traceback.format_exc()[-1500:])
    finally:
        signal.alarm(0)
    res["scf_retries"] = s1.RobustTBLite.retries
    res["wall_s"] = round(time.time() - t0, 1); res["cpu_s_main"] = round(time.process_time() - c0, 1)
    json.dump(res, open(od / "result.json", "w"), indent=1)
    return jid, res["status"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out"); ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--family", nargs="+", default=["A", "B", "C", "Q"])
    ap.add_argument("--jobs", nargs="+", help="explicit job ids (family D: append _sh<k> before _s<seed>)")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--crest", type=int, default=0, help="CREST threads (0 = no CREST)")
    ap.add_argument("--crest-timeout", type=int, default=3600)
    a = ap.parse_args()
    J = [parse_job_id(x) for x in a.jobs] if a.jobs else jobs(tuple(a.family))
    if a.list:
        from collections import Counter
        for j in J:
            print(job_id(j))
        print("counts:", dict(Counter(j["family"] for j in J)), "total", len(J), file=sys.stderr)
        return
    Path(a.out).mkdir(parents=True, exist_ok=True)
    print(f"{len(J)} jobs, {a.workers} workers", flush=True)
    t0 = time.time(); done = 0
    args = [(j, a.out, a.crest, a.crest_timeout) for j in J]
    if a.workers == 1 or len(J) == 1:
        it = map(run_job, args)
    else:
        it = Pool(a.workers).imap_unordered(run_job, args)
    for jid, st in it:
        done += 1
        print(f"[{done}/{len(J)}] {time.time() - t0:8.0f}s {jid} {st}", flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
