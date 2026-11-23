#!/usr/bin/env python3
"""Apply the stage-1b screen rule of docs/relay-system-design.md section 4 (fixed before the screen) to result.json files.

Per job, a barrier is usable if (a) status == "ok" and the band maximum is not an endpoint (barrier_fwd_kcal), or
(b) the sequential-hop path converged (segmented.seg_converged) and its maximum is not the first or last point
(seg_barrier_kcal). If both are usable, the lower one is taken (the lowest path found). This (b)/min clause was added
to section 4 on 22.11.2026 before the screen was launched (the segmented fallback did not exist when section 4 was
written); no other criterion is changed.
Per condition: barrier = median over usable seeds; the condition counts only if >= 2 of 3 seeds are usable.
Usage: python analyze1b.py RUN_DIR [--out summary.json]
"""
import argparse, json, statistics
from collections import defaultdict
from pathlib import Path

DELTA = 2.0          # kcal/mol
A_MAX = 20.0         # kcal/mol, family A absolute ceiling
SEEDS_NEEDED = 2


def job_barrier(r):
    opts = []
    if r.get("status") == "ok" and r.get("ts_is_endpoint") is False and r.get("barrier_fwd_kcal") is not None:
        opts.append(("band", r["barrier_fwd_kcal"]))
    sg = r.get("segmented") or {}
    path = sg.get("seg_E_path_kcal") or []
    if sg.get("seg_converged") and path:
        imax = max(range(len(path)), key=lambda i: path[i])
        if 0 < imax < len(path) - 1:
            opts.append(("segmented", sg["seg_barrier_kcal"]))
    return min(opts, key=lambda x: x[1]) if opts else None


def key(r):
    f = r["family"]
    if f == "A":
        return (f, r["base"], r["solvent"], r["n"], r["d"])
    if f == "B":
        return (f, r["substrate"], r["shuttles"], r["solvent"], r["topology"])
    if f == "C":
        return (f, r["base"], r["solvent"], r["topology"])
    return (f, r["solvent"], r["n"])


def conditions(rows):
    by = defaultdict(list)
    for r in rows:
        by[key(r)].append(r)
    out = {}
    for k, rs in by.items():
        b = [job_barrier(r) for r in rs]
        use = [x[1] for x in b if x]
        out[k] = {"n_jobs": len(rs), "n_usable": len(use), "via": sorted({x[0] for x in b if x}),
                  "statuses": dict(sorted({s: sum(r["status"] == s for r in rs) for s in {r["status"] for r in rs}}.items())),
                  "barrier": round(statistics.median(use), 3) if len(use) >= SEEDS_NEEDED else None,
                  "seed_barriers": sorted(round(u, 2) for u in use)}
    return out


def verdicts(C):
    B = lambda k: C.get(k, {}).get("barrier")
    res = {"A": [], "B": [], "C": [], "Q": []}
    for base in ("nme3", "nmi", "acetate"):
        for solv in ("water", "methanol"):
            for d in (4.0, 5.0, 6.0, 7.0):
                rel = [(B(("A", base, solv, n, d)), n) for n in (1, 2, 3) if B(("A", base, solv, n, d)) is not None]
                best = min(rel) if rel else None
                dr = B(("A", base, solv, 0, d))
                if best is None:
                    v = "no_usable_relay"
                elif dr is None:
                    v = "undetermined_direct_unusable"
                else:
                    v = "HIT" if (best[0] <= dr - DELTA and best[0] <= A_MAX) else "no"
                res["A"].append(dict(base=base, solvent=solv, d=d, relay_best=best[0] if best else None,
                                     relay_best_n=best[1] if best else None, direct=dr, verdict=v))
    for sub in ("acac", "dienone"):
        for sh in ("w1", "w2", "m1", "acoh", "pyr", "phoh", "hfip"):
            for solv in ("water", "methanol"):
                r, st = B(("B", sub, sh, solv, "relay")), B(("B", sub, sh, solv, "stepwise"))
                w1 = B(("B", sub, "w1", solv, "relay")); intra = B(("B", sub, "none", solv, "relay"))
                if r is None or st is None or w1 is None or (sub == "dienone" and intra is None):
                    v = "undetermined"
                else:
                    ok = r <= st - DELTA and r <= w1 - DELTA and (sub != "dienone" or r <= intra - DELTA)
                    v = "HIT" if ok else "no"
                res["B"].append(dict(substrate=sub, shuttle=sh, solvent=solv, relay=r, stepwise=st, relay_w1=w1,
                                     intramolecular=intra, verdict=v))
    for base in ("nme3", "acetate"):
        for solv in ("water", "methanol"):
            r, dr = B(("C", base, solv, "relay")), B(("C", base, solv, "direct"))
            loose = B(("A", base, solv, 3, None))
            v = "undetermined" if None in (r, dr, loose) else ("HIT" if r <= dr - DELTA and r <= loose - DELTA else "no")
            res["C"].append(dict(base=base, solvent=solv, relay=r, direct=dr, loose_w3=loose, verdict=v))
    q_pass = False
    for solv in ("water", "methanol"):
        d0 = B(("Q", solv, 0))
        for n in (2, 3):
            r = B(("Q", solv, n))
            p = r is not None and (d0 is None or d0 >= r + DELTA)
            res["Q"].append(dict(solvent=solv, n=n, relay=r, n0=d0, passes=p))
            q_pass |= p
    res["Q_control_passes"] = q_pass
    res["hits"] = {f: sum(x["verdict"] == "HIT" for x in res[f]) for f in "ABC"}
    return res


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("run_dir"); ap.add_argument("--out")
    a = ap.parse_args()
    rows = [json.load(open(p)) for p in sorted(Path(a.run_dir).glob("*/result.json"))]
    C = conditions(rows); V = verdicts(C)
    st = defaultdict(int)
    for r in rows:
        st[r["status"]] += 1
    summary = {"n_results": len(rows), "statuses": dict(st), "verdicts": V,
               "conditions": {"|".join(map(str, k)): v for k, v in sorted(C.items(), key=lambda x: str(x[0]))}}
    print(json.dumps({"n_results": len(rows), "statuses": dict(st), "hits": V["hits"],
                      "Q_control_passes": V["Q_control_passes"]}, indent=1))
    if a.out:
        Path(a.out).write_text(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
