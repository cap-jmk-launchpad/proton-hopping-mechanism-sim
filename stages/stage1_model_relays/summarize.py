#!/usr/bin/env python3
"""Summarise stage-1 results: per-job CSV and per-condition relay-vs-direct table (min and spread over seeds).

Usage: python summarize.py RUN_DIR [--out-dir ../../data]
Status meanings: ok = CI-NEB converged; neb_not_converged; reactant_not_minimum = the primary proton left its donor
while the reactant complex was relaxed (no barrier at this level); product_not_minimum = the product complex
relaxed back to the reactant (no stable product at this level); timeout; error.
"""
import argparse, csv, json, statistics as st
from collections import defaultdict
from pathlib import Path

KEYS = ["job", "step", "sink", "solvent", "n", "topology", "seed", "natoms", "charge", "status",
        "dE_rxn_kcal", "barrier_fwd_kcal", "barrier_rev_kcal", "ts_is_endpoint", "wire_rearranged_R",
        "wire_rearranged_P", "neb_interp", "neb_etemp_K", "neb_steps", "scf_retries", "wall_s"]


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("run_dir")
    ap.add_argument("--out-dir", default=str(Path(__file__).resolve().parents[2] / "data"))
    a = ap.parse_args()
    rows = [json.load(open(p)) for p in sorted(Path(a.run_dir).glob("*/result.json"))]
    out = Path(a.out_dir); out.mkdir(parents=True, exist_ok=True)
    with open(out / "stage1_jobs.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=KEYS, extrasaction="ignore"); w.writeheader(); w.writerows(rows)
    groups = defaultdict(list)
    for r in rows:
        groups[(r["step"], r["sink"], r["solvent"], r["n"], r["topology"])].append(r)
    table = []
    for k, rs in sorted(groups.items()):
        ok = [r for r in rs if r["status"] == "ok" and r.get("neb_etemp_K", 300.0) == 300.0]
        bar = [r["barrier_fwd_kcal"] for r in ok]
        status = defaultdict(int)
        for r in rs:
            status[r["status"]] += 1
        table.append(dict(step=k[0], sink=k[1], solvent=k[2], n=k[3], topology=k[4], n_seeds=len(rs),
                          n_ok=len(ok), status_counts=dict(status),
                          barrier_min_kcal=min(bar) if bar else None,
                          barrier_median_kcal=round(st.median(bar), 3) if bar else None,
                          barrier_max_kcal=max(bar) if bar else None,
                          dE_rxn_median_kcal=round(st.median([r["dE_rxn_kcal"] for r in rs if r.get("dE_rxn_kcal") is not None]), 3)
                          if any(r.get("dE_rxn_kcal") is not None for r in rs) else None))
    json.dump({"run_dir": str(a.run_dir), "n_jobs": len(rows), "level": "GFN2-xTB/ALPB (tblite), electronic energies",
               "groups": table}, open(out / "stage1_summary.json", "w"), indent=1)
    print(f"{len(rows)} jobs")
    print(f"{'step':6s} {'sink':9s} {'solv':8s} n {'topology':7s} ok/seeds  barrier min/med/max (kcal/mol)  dE_rxn med  statuses")
    for t in table:
        b = "/".join("-" if t[x] is None else f"{t[x]:.1f}" for x in ("barrier_min_kcal", "barrier_median_kcal", "barrier_max_kcal"))
        print(f"{t['step']:6s} {t['sink']:9s} {t['solvent']:8s} {t['n']} {t['topology']:7s} {t['n_ok']}/{t['n_seeds']}      {b:22s} "
              f"{'-' if t['dE_rxn_median_kcal'] is None else t['dE_rxn_median_kcal']:>8}  {t['status_counts']}")


if __name__ == "__main__":
    main()
