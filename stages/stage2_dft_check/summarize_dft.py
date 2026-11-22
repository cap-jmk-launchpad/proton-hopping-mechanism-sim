#!/usr/bin/env python3
"""Collect stage-2 DFT summaries (summary.json per job) into a CSV and a relay-vs-direct table next to the xtb barriers.
Rejection rule (plan, fixed before runs): the relay is rejected if dG_dagger(relay) is not lower than dG_dagger(direct)
at the same level of theory (or if the sign flips between levels). Only TS with exactly one imaginary frequency
(< -50 cm-1) and minima without one are counted as validated."""
import argparse, csv, glob, json, re
from pathlib import Path

ap = argparse.ArgumentParser()
ap.add_argument("dft_dirs", nargs="+")
ap.add_argument("--csv", default=None)
a = ap.parse_args()
rows = {}
for d in a.dft_dirs:
    for f in glob.glob(f"{d}/*/summary.json"):
        s = json.loads(Path(f).read_text())
        m = re.match(r"(\w+?)_(\w+?)_(water|methanol)_n(\d)_(relay|direct)_s(\d)", s["job"])
        s.update(step=m[1], sink=m[2], solvent=m[3], n=int(m[4]), topology=m[5], seed=int(m[6]), dir=str(Path(f).parent))
        rows[s["job"]] = s
keys = ["job", "step", "sink", "solvent", "n", "topology", "seed", "xtb_barrier_fwd_kcal", "dE_dagger_kcal",
        "dE0_dagger_kcal", "dG_dagger_kcal", "dE_dagger_tz_kcal", "dG_dagger_tz_kcal", "dE_rxn_kcal", "dG_rxn_kcal",
        "ts_one_imag", "minima_ok", "opt_converged"]
if a.csv:
    with open(a.csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, keys, extrasaction="ignore"); w.writeheader()
        for k in sorted(rows): w.writerow(rows[k])
fmt = lambda v: "-" if v is None else (f"{v:.2f}" if isinstance(v, float) else str(v))
print("| solvent | n | topology | xtb dE‡ | DFT dE‡ | dE0‡ (ZPE) | dG‡ 298 K | dG‡ TZVP//SVP | TS 1 imag | minima ok |")
print("|---|---|---|---|---|---|---|---|---|---|")
for k in sorted(rows, key=lambda j: (rows[j]["sink"], rows[j]["solvent"], rows[j]["n"], rows[j]["topology"])):
    s = rows[k]
    print(f"| {s['solvent']} | {s['n']} | {s['topology']} | {fmt(s.get('xtb_barrier_fwd_kcal'))} | {fmt(s.get('dE_dagger_kcal'))} | "
          f"{fmt(s.get('dE0_dagger_kcal'))} | {fmt(s.get('dG_dagger_kcal'))} | {fmt(s.get('dG_dagger_tz_kcal'))} | "
          f"{s.get('ts_one_imag')} | {s.get('minima_ok')} |")
print("\nRejection rule per condition (validated stationary points only):")
conds = {}
for s in rows.values():
    conds.setdefault((s["sink"], s["solvent"], s["n"]), {})[s["topology"]] = s
for (sink, solv, n), t in sorted(conds.items()):
    if "relay" in t and "direct" in t:
        r, d = t["relay"], t["direct"]
        valid = all(x.get("ts_one_imag") and x.get("minima_ok") for x in (r, d))
        if r.get("dG_dagger_kcal") is not None and d.get("dG_dagger_kcal") is not None:
            verdict = "relay rejected" if r["dG_dagger_kcal"] >= d["dG_dagger_kcal"] else "relay favoured"
            print(f"- {sink}/{solv}/n={n}: dG‡ relay {r['dG_dagger_kcal']:.2f} vs direct {d['dG_dagger_kcal']:.2f} kcal/mol"
                  f" -> {verdict}{'' if valid else ' (NOT validated: frequency check failed)'}")
        else:
            print(f"- {sink}/{solv}/n={n}: incomplete")
