#!/usr/bin/env python3
"""Stage 2: DFT check of converged stage-1 (GFN2-xTB/ALPB) relay/direct pairs.

Open-source only: PySCF (Apache-2.0) + geomeTRIC (BSD-3) + simple-dftd3 via pyscf-dispersion (LGPL-3).
Per stage-1 job (lowest-barrier converged seed of a condition):
  R  : optimise CI-NEB image 0 (minimum search)
  P  : optimise last CI-NEB image (minimum search)
  TS : geomeTRIC transition-state optimisation from the climbing image, starting from the analytic DFT Hessian
  analytic Hessian at R, TS, P -> harmonic frequencies, ZPE, RRHO thermal corrections (298.15 K, 1 atm)
  optional def2-TZVP single points on the optimised geometries
Default level: wB97X-D3(BJ)/def2-SVP, IEF-PCM (water eps 78.3553, methanol eps 32.613), RI-JK (def2-universal-jkfit).
Every step is checkpointed (one JSON + xyz per step) so a restart skips finished steps. Nothing is ever deleted.
"""
import argparse, json, os, sys, time, traceback
from pathlib import Path
import numpy as np
from pyscf import gto, dft, lib
from pyscf.hessian import thermo
from pyscf.geomopt import geometric_solver

HARTREE2KCAL = 627.509474
EPS = {"water": 78.3553, "methanol": 32.613}


def read_band(path):
    """Read an ASE extxyz band (first 4 columns: symbol x y z)."""
    lines = Path(path).read_text().splitlines()
    imgs, i = [], 0
    while i < len(lines):
        n = int(lines[i]); body = lines[i + 2:i + 2 + n]
        imgs.append(([l.split()[0] for l in body], np.array([[float(x) for x in l.split()[1:4]] for l in body])))
        i += n + 2
    return imgs


def write_xyz(path, sym, xyz, comment=""):
    with open(path, "w") as f:
        f.write(f"{len(sym)}\n{comment}\n")
        for s, (x, y, z) in zip(sym, xyz):
            f.write(f"{s:2s} {x:14.8f} {y:14.8f} {z:14.8f}\n")


def make_mol(sym, xyz, charge, basis):
    return gto.M(atom=[(s, tuple(p)) for s, p in zip(sym, xyz)], basis=basis, charge=charge, spin=0,
                 unit="Angstrom", verbose=int(os.environ.get("PYSCF_VERBOSE", "3")))


def make_mf(mol, solvent, xc, df=True):
    mf = dft.RKS(mol)
    mf.xc = xc
    mf.grids.level = 3
    mf.conv_tol = 1e-9
    mf.max_cycle = 300
    mf.max_memory = int(os.environ.get("PYSCF_MAX_MEMORY", "12000"))
    if df:
        mf = mf.density_fit()
    mf = mf.PCM()
    mf.with_solvent.method = "IEF-PCM"
    mf.with_solvent.eps = EPS[solvent]
    return mf


def hessian_analysis(mf, mol):
    h = mf.Hessian().kernel()
    res = thermo.harmonic_analysis(mol, h, imaginary_freq=False)
    freqs = np.asarray(res["freq_wavenumber"]).real
    th = thermo.thermo(mf, res["freq_au"], 298.15, 101325)
    return h, freqs, {k: float(th[k][0]) for k in ("ZPE", "E_tot", "H_tot", "G_tot", "S_tot") if k in th}


def h_positions(xyz, moving):
    """For each moving proton (h, donor, acceptor): distances to donor and acceptor (A)."""
    return [{"h": h, "d_donor": round(float(np.linalg.norm(xyz[h] - xyz[d])), 3),
             "d_acceptor": round(float(np.linalg.norm(xyz[h] - xyz[a])), 3)} for h, d, a in moving]


def step(out, name, fn):
    js = out / f"{name}.json"
    if js.exists():
        d = json.loads(js.read_text())
        if d.get("status") == "done":
            return d
    t0 = time.time()
    try:
        d = fn()
        d["status"] = "done"
    except Exception as e:
        d = {"status": "error", "error": repr(e), "traceback": traceback.format_exc()}
    d["wall_s"] = round(time.time() - t0, 1)
    js.write_text(json.dumps(d, indent=1))
    print(f"[{time.strftime('%H:%M:%S')}] {out.name}/{name}: {d['status']} in {d['wall_s']} s", flush=True)
    return d


def run(job_dir, outroot, xc, basis, tz_basis, maxsteps):
    job_dir = Path(job_dir)
    meta = json.loads((job_dir / "result.json").read_text())
    out = Path(outroot) / meta["job"]
    out.mkdir(parents=True, exist_ok=True)
    imgs = read_band(job_dir / "band.xyz")
    sym = imgs[0][0]
    i_ts = meta["i_max"]
    charge, solvent, moving = meta["charge"], meta["solvent"], meta["moving"]
    (out / "input.json").write_text(json.dumps({"stage1_job_dir": str(job_dir), "xc": xc, "basis": basis,
                                                "tz_basis": tz_basis, "solvent": solvent, "eps": EPS[solvent],
                                                "charge": charge, "i_ts_guess": i_ts, "moving": moving,
                                                "xtb_barrier_fwd_kcal": meta["barrier_fwd_kcal"],
                                                "omp_threads": lib.num_threads()}, indent=1))
    starts = {"R": imgs[0][1], "TS": imgs[i_ts][1], "P": imgs[-1][1]}
    for k, x in starts.items():
        write_xyz(out / f"{k}_start_xtb.xyz", sym, x, f"{meta['job']} {k} start (GFN2-xTB CI-NEB image)")

    def opt(kind):
        def fn():
            mol = make_mol(sym, starts[kind], charge, basis)
            mf = make_mf(mol, solvent, xc)
            ncyc = []
            kw = dict(maxsteps=maxsteps, callback=lambda env: ncyc.append(float(env["energy"])))
            if kind == "TS":
                kw.update(transition=True, hessian=True)
            conv, mol_eq = geometric_solver.kernel(mf, **kw)
            xyz = mol_eq.atom_coords(unit="Angstrom")
            write_xyz(out / f"{kind}_opt.xyz", sym, xyz, f"{meta['job']} {kind} {xc}/{basis} PCM({solvent}) conv={conv}")
            return {"converged": bool(conv), "n_cycles": len(ncyc), "E_last_opt_Eh": ncyc[-1] if ncyc else None,
                    "h_positions": h_positions(xyz, moving)}
        return fn

    def freq(kind):
        def fn():
            sym2, xyz = read_xyz(out / f"{kind}_opt.xyz")
            mol = make_mol(sym2, xyz, charge, basis)
            mf = make_mf(mol, solvent, xc)
            e = float(mf.kernel())
            if not mf.converged:
                raise RuntimeError("SCF not converged at optimised geometry")
            h, freqs, th = hessian_analysis(mf, mol)
            np.save(out / f"{kind}_hessian.npy", h)
            nimag = int(np.sum(freqs < -50.0))
            return {"E_el_Eh": e, "freqs_cm1": [round(float(f), 2) for f in freqs], "n_imag_lt_-50": nimag,
                    "n_imag_any": int(np.sum(freqs < 0)), "lowest_cm1": round(float(freqs.min()), 2), **th}
        return fn

    def sp_tz(kind):
        def fn():
            sym2, xyz = read_xyz(out / f"{kind}_opt.xyz")
            mf = make_mf(make_mol(sym2, xyz, charge, tz_basis), solvent, xc)
            e = float(mf.kernel())
            if not mf.converged:
                raise RuntimeError("SCF not converged")
            return {"E_el_Eh": e, "basis": tz_basis}
        return fn

    res = {"job": meta["job"]}
    for kind in ("R", "P", "TS"):
        res[f"opt_{kind}"] = step(out, f"opt_{kind}", opt(kind))
    for kind in ("R", "P", "TS"):
        if res[f"opt_{kind}"]["status"] == "done":
            res[f"freq_{kind}"] = step(out, f"freq_{kind}", freq(kind))
    if tz_basis:
        for kind in ("R", "P", "TS"):
            if res[f"opt_{kind}"]["status"] == "done":
                res[f"sp_{kind}"] = step(out, f"sp_tz_{kind}", sp_tz(kind))
    (out / "summary.json").write_text(json.dumps(summarise(res, meta), indent=1))
    return res


def read_xyz(path):
    lines = Path(path).read_text().splitlines()
    n = int(lines[0]); body = lines[2:2 + n]
    return [l.split()[0] for l in body], np.array([[float(x) for x in l.split()[1:4]] for l in body])


def summarise(res, meta):
    s = {"job": meta["job"], "xtb_barrier_fwd_kcal": meta["barrier_fwd_kcal"], "xtb_dE_rxn_kcal": meta["dE_rxn_kcal"]}
    f = {k: res.get(f"freq_{k}", {}) for k in ("R", "TS", "P")}
    ok = all(f[k].get("status") == "done" for k in f)
    s["ts_one_imag"] = bool(ok and f["TS"]["n_imag_lt_-50"] == 1)
    s["minima_ok"] = bool(ok and f["R"]["n_imag_lt_-50"] == 0 and f["P"]["n_imag_lt_-50"] == 0)
    s["opt_converged"] = {k: res[f"opt_{k}"].get("converged") for k in ("R", "TS", "P")}
    if ok:
        E = {k: f[k]["E_el_Eh"] for k in f}
        s["dE_dagger_kcal"] = round((E["TS"] - E["R"]) * HARTREE2KCAL, 2)
        s["dE_rxn_kcal"] = round((E["P"] - E["R"]) * HARTREE2KCAL, 2)
        if all("ZPE" in f[k] for k in f):
            s["dE0_dagger_kcal"] = round((E["TS"] + f["TS"]["ZPE"] - E["R"] - f["R"]["ZPE"]) * HARTREE2KCAL, 2)
        if all("G_tot" in f[k] for k in f):
            s["dG_dagger_kcal"] = round((f["TS"]["G_tot"] - f["R"]["G_tot"]) * HARTREE2KCAL, 2)
            s["dG_rxn_kcal"] = round((f["P"]["G_tot"] - f["R"]["G_tot"]) * HARTREE2KCAL, 2)
        sp = {k: res.get(f"sp_{k}", {}) for k in ("R", "TS", "P")}
        if all(sp[k].get("status") == "done" for k in sp):
            Et = {k: sp[k]["E_el_Eh"] for k in sp}
            s["dE_dagger_tz_kcal"] = round((Et["TS"] - Et["R"]) * HARTREE2KCAL, 2)
            s["dG_dagger_tz_kcal"] = round(((Et["TS"] - E["TS"] + f["TS"]["G_tot"]) -
                                            (Et["R"] - E["R"] + f["R"]["G_tot"])) * HARTREE2KCAL, 2)
    s["h_positions"] = {k: res[f"opt_{k}"].get("h_positions") for k in ("R", "TS", "P")}
    return s


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("job_dirs", nargs="+", help="stage-1 job directories (with result.json and band.xyz)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--xc", default="wb97x-d3bj")
    ap.add_argument("--basis", default="def2-svp")
    ap.add_argument("--tz-basis", default="", help="e.g. def2-tzvp for single points (empty = skip)")
    ap.add_argument("--maxsteps", type=int, default=150)
    a = ap.parse_args()
    for jd in a.job_dirs:
        run(jd, a.out, a.xc, a.basis, a.tz_basis, a.maxsteps)
