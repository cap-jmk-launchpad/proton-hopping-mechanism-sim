#!/usr/bin/env python3
"""List stage-1 job dirs for the DFT check: per (step, sink, solvent, n, topology) the lowest-barrier converged seed
(status ok, 300 K NEB). Priority 1: NMe3 deprotonation, water then methanol, n = 0-4, conditions where both relay and
direct converged plus the n = 0 contact pair. Priority 2: OH-/MeO- relay chains with a non-zero xtb barrier
(barrierless chains, ts_is_endpoint = True, have no TS to refine)."""
import argparse, glob, json
from pathlib import Path

ap = argparse.ArgumentParser()
ap.add_argument("run_dirs", nargs="+")
ap.add_argument("--priority", type=int, choices=[1, 2], default=1)
a = ap.parse_args()
best = {}
for rd in a.run_dirs:
    for f in glob.glob(f"{rd}/*/result.json"):
        if ".failed" in f:
            continue
        r = json.loads(Path(f).read_text())
        if r["status"] != "ok" or r.get("neb_etemp_K", 300.0) != 300.0:
            continue
        k = (r["step"], r["sink"], r["solvent"], r["n"], r["topology"])
        if k not in best or r["barrier_fwd_kcal"] < best[k][0]:
            best[k] = (r["barrier_fwd_kcal"], r["natoms"], str(Path(f).parent), r["ts_is_endpoint"])
out = []
if a.priority == 1:
    for solv in ("water", "methanol"):
        for n in range(5):
            topo = ["relay"] if n == 0 else ["relay", "direct"]
            ks = [("deprot", "nme3", solv, n, t) for t in topo]
            if all(k in best for k in ks):
                out += [best[k] for k in ks]
else:
    for k in sorted(best):
        if k[1] in ("hydroxide", "methoxide") and k[4] == "relay" and k[0] == "deprot" and not best[k][3]:
            out.append(best[k])
for b, nat, d, _ in out:
    print(d)
