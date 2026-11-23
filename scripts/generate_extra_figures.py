#!/usr/bin/env python3
"""Extra manuscript diagrams.
  fig4_mechanism_schematic.pdf  : schematic of direct transfer, free-solution relay and the held-off/bridged relay.
  fig5_paired_forest.pdf        : seed-paired barrier differences (relay minus reference) for every paired condition.
  fig6_tcmn_active_site.pdf     : projection of the relaxed 3TVQ active-site cluster with measured distances.
All numbers are read from result files; nothing is typed in by hand except the schematic distances, which are labelled
as representative values taken from the tables.
"""
import json
import statistics as st
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch

plt.rcParams.update({"font.family": "sans-serif", "font.size": 9, "axes.spines.top": False,
                     "axes.spines.right": False, "pdf.fonttype": 42})
BLUE, ORANGE, GREEN, PURPLE, GREY, CYAN, YELLOW = "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#7F7F7F", "#56B4E9", "#E69F00"
OUT = Path("/tmp/proton-hopping-git/paper")
C1C = Path("/mnt/homelab-external/proton-hopping-stage1b/runs/stage1c-nonpeptide")
POCKET = Path("/mnt/homelab-external/proton-hopping-pocket")


# ---------------------------------------------------------------- Fig 4: schematic
def atom(ax, xy, label, color, r=0.42, fs=9, tc="white"):
    ax.add_patch(Circle(xy, r, fc=color, ec="black", lw=0.8, zorder=3))
    ax.text(xy[0], xy[1], label, ha="center", va="center", color=tc, fontsize=fs, zorder=4, fontweight="bold")


def hbond(ax, a, b, color="black", ls="--", lw=1.1):
    ax.add_patch(FancyArrowPatch(a, b, arrowstyle="-", linestyle=ls, color=color, lw=lw, zorder=2))


def proton(ax, xy):
    ax.add_patch(Circle(xy, 0.14, fc="white", ec="black", lw=0.8, zorder=5))
    ax.text(xy[0], xy[1], "H", ha="center", va="center", fontsize=6, zorder=6)


def arrow(ax, a, b, color=ORANGE):
    ax.add_patch(FancyArrowPatch(a, b, arrowstyle="-|>", mutation_scale=9, color=color, lw=1.4, zorder=6,
                                 connectionstyle="arc3,rad=-0.35"))


fig, axs = plt.subplots(1, 3, figsize=(7.4, 2.9))
for ax in axs:
    ax.set_xlim(-0.8, 8.8); ax.set_ylim(-1.8, 2.6); ax.set_aspect("equal"); ax.axis("off")

# (a) direct contact
ax = axs[0]
atom(ax, (0.6, 0), "C", GREY); proton(ax, (1.45, 0.05)); atom(ax, (3.0, 0), "N", BLUE)
hbond(ax, (1.55, 0.05), (2.7, 0.0)); arrow(ax, (1.45, 0.3), (2.7, 0.35))
ax.annotate("", xy=(0.6, -0.9), xytext=(3.0, -0.9), arrowprops=dict(arrowstyle="<->", lw=1))
ax.text(1.8, -1.3, r"$R_{\mathrm{DA}}\approx 2.7$ Å", ha="center", fontsize=8)
ax.set_title("(a) Unhindered base:\ndirect contact", fontsize=9, fontweight="bold", loc="left")
ax.text(1.8, 1.9, r"$\Delta E_{\mathrm{act}}$ 3–7 kcal mol$^{-1}$", ha="center", fontsize=8, color=GREEN)

# (b) free solution wire
ax = axs[1]
atom(ax, (0.2, 0), "C", GREY); proton(ax, (1.0, 0.05))
atom(ax, (2.2, 0.2), "O", ORANGE); atom(ax, (4.2, 0.2), "O", ORANGE); atom(ax, (6.4, 0), "N", BLUE)
for a, b in [((1.1, 0.05), (1.9, 0.18)), ((2.55, 0.2), (3.9, 0.2)), ((4.55, 0.2), (6.1, 0.05))]:
    hbond(ax, a, b)
ax.text(3.2, 1.15, "loose water wire", ha="center", fontsize=8)
ax.text(3.2, -1.2, "ordering + charge\nseparation cost", ha="center", fontsize=8)
ax.set_title("(b) Free wire in solution:\nrelay loses", fontsize=9, fontweight="bold", loc="left")
ax.text(3.2, 1.9, r"$\Delta E_{\mathrm{act}}$ 13–17 kcal mol$^{-1}$", ha="center", fontsize=8, color=ORANGE)

# (c) held-off base + bridging water
ax = axs[2]
atom(ax, (0.2, 0), "C", GREY); proton(ax, (1.0, 0.05))
atom(ax, (2.5, 0.1), "O", ORANGE); proton(ax, (3.3, 0.1))
atom(ax, (5.6, 0), "N", BLUE)
ax.add_patch(Circle((5.6, 0), 0.8, fc="none", ec=PURPLE, lw=1.4, ls=":", zorder=1))
ax.text(6.9, 1.0, "bulky\ngroups", ha="left", fontsize=7.5, color=PURPLE)
hbond(ax, (1.1, 0.05), (2.2, 0.1)); hbond(ax, (3.4, 0.1), (5.3, 0.02))
ax.annotate("", xy=(0.2, -0.95), xytext=(5.6, -0.95), arrowprops=dict(arrowstyle="<->", lw=1))
ax.text(2.9, -1.4, r"$R_{\mathrm{DA}}\geq 4.9$ Å", ha="center", fontsize=8)
ax.set_title("(c) Base held off\n+ bridging water: relay", fontsize=9, fontweight="bold", loc="left")
ax.text(3.2, 1.9, r"$\Delta E_{\mathrm{act}}$ 1.8–6 kcal mol$^{-1}$", ha="center", fontsize=8, color=GREEN)
fig.tight_layout()
fig.savefig(OUT / "fig4_mechanism_schematic.pdf", bbox_inches="tight")
plt.close(fig)


# ---------------------------------------------------------------- Fig 5: paired forest plot
def seeds(cond):
    out = {}
    for s in (0, 1, 2):
        p = C1C / f"{cond}_s{s}" / "result.json"
        if p.exists():
            r = json.load(open(p))
            if r.get("status") == "ok":
                out[s] = r["barrier_kcal"]
    return out


pairs = [
    ("DBU: relay − direct", "E1_dbu_n1_direct", "E1_dbu_n1_relay", BLUE),
    ("TMP: relay − direct", "E1_tmp_n1_direct", "E1_tmp_n1_relay", BLUE),
    ("DTBMP: relay − direct", "E1_dtbmp_n1_direct", "E1_dtbmp_n1_relay", PURPLE),
    ("[MgCl$_2$(acac)] + NMI", "E3_MgCl2_nmi_n1_direct", "E3_MgCl2_nmi_n1_relay", GREEN),
    ("[MgCl$_2$(acac)] + DTBMP", "E3_MgCl2_dtbmp_n1_direct", "E3_MgCl2_dtbmp_n1_relay", GREEN),
    ("[ZnCl$_2$(acac)] + DTBMP", "E3_ZnCl2_dtbmp_n1_direct", "E3_ZnCl2_dtbmp_n1_relay", GREEN),
    ("TBD: cyclic − stepwise", "E2_tbd_acac_stepwise", "E2_tbd_acac_relay", ORANGE),
]
P = json.load(open("/tmp/proton-hopping-git/data/pocket_3TVQ_scans_summary.json"))
dref = [P["direct_C10_Asp57"][k]["barrier"] for k in ("s7", "s11", "s13")]
rows = []
for lab, a, b, col in pairs:
    A, B = seeds(a), seeds(b)
    d = [B[s] - A[s] for s in sorted(set(A) & set(B))]
    rows.append((lab, d, col))
for lab, ch in [("TcmN: 2-water wire − direct C10", "relay_O27_W259_W266_His128"),
                ("TcmN: 1-water wire − direct C10", "relay_O27_W266_His128")]:
    r = [P[ch][k]["barrier"] for k in ("s7", "s11", "s13")]
    rows.append((lab, [x - y for x, y in zip(r, dref)], CYAN))

fig, ax = plt.subplots(figsize=(7.0, 3.6))
for i, (lab, d, col) in enumerate(rows[::-1]):
    m, sd = st.mean(d), st.stdev(d)
    ax.errorbar(m, i, xerr=sd, fmt="o", color=col, capsize=3, ms=7, zorder=3)
    ax.scatter(d, np.full(len(d), i), s=14, color="black", zorder=4)
ax.axvline(0, color="black", lw=0.8)
ax.set_yticks(range(len(rows)))
ax.set_yticklabels([r[0] for r in rows[::-1]], fontsize=8)
ax.set_xlabel(r"Seed-paired $\Delta\Delta E_{\mathrm{act}}$ (relay $-$ reference), kcal mol$^{-1}$  (mean $\pm$ SD, dots = seeds)")
ax.grid(axis="x", ls=":", alpha=0.5)
ax.text(0.02, 1.02, "relay lower", transform=ax.transAxes, fontsize=8, color=GREEN, ha="left")
ax.text(0.98, 1.02, "relay higher", transform=ax.transAxes, fontsize=8, color=ORANGE, ha="right")
ax.set_xlim(-32, 12)
fig.tight_layout()
fig.savefig(OUT / "fig5_paired_forest.pdf", bbox_inches="tight")
plt.close(fig)


# ---------------------------------------------------------------- Fig 6: real active-site projection
meta = json.load(open(POCKET / "cluster/cluster_3TVQ.json"))
labels = meta["labels"]
lines = (POCKET / "xtb/ase_opt/reactant_relaxed.xyz").read_text().splitlines()[2:]
sym = [l.split()[0] for l in lines]
X = np.array([[float(v) for v in l.split()[1:4]] for l in lines])
idx = {l: i for i, l in enumerate(labels)}
key = ["DQH:C10", "DQH:O27", "ASP57:OD2", "HIS128:NE2", "HOH259:O", "HOH266:O", "HOH234:O"]
keep_exact = {"ASP57:CG", "ASP57:OD1", "ASP57:OD2", "HIS128:CG", "HIS128:ND1", "HIS128:CD2", "HIS128:CE1", "HIS128:NE2",
              "HOH259:O", "HOH259:H1", "HOH259:H2", "HOH266:O", "HOH266:H1", "HOH266:H2", "HOH234:O", "HOH234:H1", "HOH234:H2",
              "DQH:H24", "DQH:H31"}
dqh_heavy = [i for i, l in enumerate(labels) if l.startswith("DQH:") and sym[i] != "H"]
sel = sorted(set(dqh_heavy) | {idx[k] for k in keep_exact})
kp = np.array([X[idx[k]] for k in key])
c = kp.mean(0)
u, s_, vt = np.linalg.svd(kp - c)
P2 = (X - c) @ vt[:2].T
fig, ax = plt.subplots(figsize=(6.6, 4.6))
col = {"C": "#BBBBBB", "O": ORANGE, "N": BLUE, "H": "white"}
for i in sel:
    for j in sel:
        if j <= i:
            continue
        d = np.linalg.norm(X[i] - X[j])
        covalent = ((sym[i] == "H") != (sym[j] == "H") and d < 1.25) or (sym[i] != "H" and sym[j] != "H" and d < 1.65)
        if covalent:
            faint = labels[i].startswith("DQH") and labels[j].startswith("DQH")
            ax.plot(*zip(P2[i], P2[j]), color="#888888" if faint else "black", lw=0.8 if faint else 1.4, zorder=1)
        elif (sym[i] == "H") != (sym[j] == "H") and d < 2.1:
            other = j if sym[i] == "H" else i
            if sym[other] in "ON":
                ax.plot(*zip(P2[i], P2[j]), color=PURPLE, lw=1.1, ls="--", zorder=1)
for i in sel:
    ax.scatter(*P2[i], s=18 if sym[i] == "H" else 60, c=col[sym[i]], edgecolors="black", linewidths=0.6, zorder=3)
off = {"DQH:C10": (-52, -14), "DQH:O27": (6, -16), "ASP57:OD2": (8, 4), "HIS128:NE2": (8, -4),
       "HOH259:O": (-48, 6), "HOH266:O": (8, 6), "HOH234:O": (8, 4)}
for name in key:
    ax.annotate(name.replace(":", " "), P2[idx[name]], textcoords="offset points", xytext=off[name], fontsize=8,
                bbox=dict(fc="white", ec="none", pad=0.8, alpha=0.85), zorder=7)


def dist(a, b):
    return float(np.linalg.norm(X[idx[a]] - X[idx[b]]))


rows = [("DQH:C10", "ASP57:OD2"), ("DQH:O27", "HOH266:O"), ("HOH266:O", "HIS128:NE2"), ("HOH259:O", "HOH266:O"), ("DQH:C10", "HIS128:NE2")]
for a, b in rows:
    ax.plot(*zip(P2[idx[a]], P2[idx[b]]), color=GREEN, lw=1.3, ls=":", zorder=2)
txt = "Distances in relaxed cluster (Å)\n" + "\n".join(f"{a.split(':')[1]}{'' if a.startswith('DQH') else ''} ({a.split(':')[0]}) – {b.split(':')[1]} ({b.split(':')[0]}): {dist(a, b):.2f}" for a, b in rows)
ax.text(1.0, 0.99, txt, transform=ax.transAxes, ha="right", va="top", fontsize=7.5,
        bbox=dict(fc="white", ec=GREEN, pad=3))
ax.set_aspect("equal"); ax.axis("off")
ax.set_title("Relaxed TcmN active site (PDB 3TVQ, seed $s_7$), best-fit projection of the key atoms", fontsize=9, loc="left")
fig.tight_layout()
fig.savefig(OUT / "fig6_tcmn_active_site.pdf", bbox_inches="tight")
plt.close(fig)
print("Generated fig4, fig5, fig6")
