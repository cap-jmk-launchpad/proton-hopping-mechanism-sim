#!/usr/bin/env python3
"""Generate publication-quality vector PDF figures for the proton-hopping manuscript:
  - fig1_solution_vs_heldoff.pdf
  - fig2_tcmn_pocket_profiles.pdf
  - fig3_nonpeptide_complexes.pdf
Uses colourblind-safe Okabe-Ito palette, >=8 pt fonts, vector PDF output.
"""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 9,
    "axes.labelsize": 9.5,
    "axes.titlesize": 10,
    "legend.fontsize": 8,
    "xtick.labelsize": 8.5,
    "ytick.labelsize": 8.5,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "pdf.fonttype": 42,
})

C_BLUE = "#0072B2"
C_ORANGE = "#D55E00"
C_GREEN = "#009E73"
C_PURPLE = "#CC79A7"
C_YELLOW = "#E69F00"
C_CYAN = "#56B4E9"
C_GREY = "#7F7F7F"

OUT = Path("/tmp/proton-hopping-git/paper")
OUT.mkdir(parents=True, exist_ok=True)

# ---------------- FIGURE 1: Solution (Stage 1) vs Held-off Base (Stage 1b Family A) ----------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 3.0), dpi=300)

ns = np.array([0, 1, 2, 3])
dir_m1 = np.array([7.87, 6.78, 5.62, 4.25])
dir_sd = np.array([0.03, 0.03, 0.06, 0.16])

ns_rel = np.array([1, 2, 3])
rel_m1 = np.array([15.43, 13.44, 16.72])
rel_sd = np.array([1.27, 0.58, 0.57])
rel_neb = np.array([6.41, 17.82, 15.83])

ax1.errorbar(ns, dir_m1, yerr=dir_sd, fmt="-o", color=C_BLUE, lw=1.8, ms=5.5, capsize=3, label="Direct (1D scan, $n$ spectators)")
ax1.errorbar(ns_rel, rel_m1, yerr=rel_sd, fmt="-s", color=C_ORANGE, lw=1.8, ms=5.5, capsize=3, label="Relay (2D $11\\times 11$ grid)")
ax1.plot(ns_rel, rel_neb, "--^", color=C_YELLOW, lw=1.3, ms=5, label="Relay (static CI-NEB)")
ax1.scatter([1], [6.61], marker="*", s=90, color=C_GREEN, zorder=5, label="DFT $\\Delta G_{\\mathrm{act,TZ}}$ ($n=1$ relay)")

ax1.set_xlabel("Explicit water molecules $n$")
ax1.set_ylabel("Activation barrier $\\Delta E_{\\mathrm{act}}$ (kcal mol$^{-1}$)")
ax1.set_title("(a) Unconstrained base in water ($\\mathrm{NMe}_3$)", loc="left", fontweight="bold")
ax1.set_xticks([0, 1, 2, 3])
ax1.set_ylim(0, 22)
ax1.grid(axis="y", ls=":", alpha=0.5)
ax1.legend(loc="upper left", frameon=False)

d_nmi_meoh = np.array([4.0, 6.0, 7.0])
dir_nmi_meoh = np.array([57.563, 106.067, 111.645])
d_nmi_rel = np.array([5.0, 6.0, 7.0])
rel_nmi_meoh = np.array([9.576, 9.737, 8.301])

d_nme3_meoh = np.array([4.0, 5.0])
dir_nme3_meoh = np.array([13.43, 89.107])
d_nme3_rel = np.array([5.0, 6.0, 7.0])
rel_nme3_meoh = np.array([13.23, 7.652, 13.46])

d_nme3_wat_rel = np.array([5.0, 7.0])
rel_nme3_wat = np.array([6.223, 5.803])

ax2.plot(d_nmi_meoh, dir_nmi_meoh, "--o", color=C_ORANGE, lw=1.5, ms=5, label="Direct: NMI (MeOH)")
ax2.plot(d_nme3_meoh, dir_nme3_meoh, "--s", color=C_PURPLE, lw=1.5, ms=5, label="Direct: $\\mathrm{NMe}_3$ (MeOH)")
ax2.plot(d_nmi_rel, rel_nmi_meoh, "-o", color=C_BLUE, lw=1.8, ms=5.5, label="Relay ($n=1$–$2$): NMI (MeOH)")
ax2.plot(d_nme3_rel, rel_nme3_meoh, "-s", color=C_GREEN, lw=1.8, ms=5.5, label="Relay ($n=1$–$2$): $\\mathrm{NMe}_3$ (MeOH)")
ax2.plot(d_nme3_wat_rel, rel_nme3_wat, "-^", color=C_CYAN, lw=1.8, ms=5.5, label="Relay ($n=1$–$2$): $\\mathrm{NMe}_3$ ($\\mathrm{H}_2\\mathrm{O}$)")

ax2.axhspan(0, 20, color=C_GREEN, alpha=0.08, label="Pre-registered ceiling (20 kcal mol$^{-1}$)")
ax2.set_xlabel("Fixed distance $d(\\mathrm{C3}\\cdots\\mathrm{base})$ (\\AA)")
ax2.set_ylabel("Median barrier $\\Delta E_{\\mathrm{act}}$ (kcal mol$^{-1}$)")
ax2.set_title("(b) Held-off base crossover (Family A)", loc="left", fontweight="bold")
ax2.set_xticks([4.0, 5.0, 6.0, 7.0])
ax2.set_ylim(0, 120)
ax2.grid(axis="y", ls=":", alpha=0.5)
ax2.legend(loc="upper left", frameon=False, fontsize=7.3)

fig.tight_layout()
fig.savefig(OUT / "fig1_solution_vs_heldoff.pdf", bbox_inches="tight")
plt.close(fig)

# ---------------- FIGURE 2: 3TVQ TcmN Enzyme-Pocket Profiles & Multi-Seed Summary ----------------
B = Path("/mnt/homelab-external/proton-hopping-pocket")
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 3.0), dpi=300)

profiles = [
    ("xtb/scans/direct_C10_Asp57/scan.json", "Direct: C10–H24 $\\to$ Asp57", C_ORANGE, "-o"),
    ("xtb/scans/relay_C10_W999_Asp57_tethered/scan.json", "C10 Relay: C10 $\\to$ W999 $\\to$ Asp57", C_PURPLE, "-^"),
    ("xtb/scans/relay_O27_W259_W266_His128/scan.json", "2-Water Wire: O27 $\\to$ W259 $\\to$ W266 $\\to$ His128", C_BLUE, "-s"),
    ("xtb/scans/relay_O27_W266_His128/scan.json", "1-Water Wire: O27 $\\to$ W266 $\\to$ His128", C_GREEN, "-D"),
]

for rel_path, lab, col, fmt in profiles:
    sc = json.load(open(B / rel_path))
    pts = sc["points"]
    emin = min(p["dE_kcal"] for p in pts[:4])
    xs = np.linspace(0, 1, len(pts))
    ys = np.array([p["dE_kcal"] - emin for p in pts])
    ax1.plot(xs, ys, fmt, color=col, lw=1.7, ms=4.2, label=lab)

ax1.set_xlabel("Normalized proton-transfer coordinate $\\xi$")
ax1.set_ylabel("$\\Delta E - \\Delta E_{\\mathrm{R,min}}$ (kcal mol$^{-1}$)")
ax1.set_title("(a) TcmN pocket 1D relaxed scans (PDB 3TVQ, $s_7$)", loc="left", fontweight="bold")
ax1.set_ylim(-2, 40)
ax1.grid(axis="y", ls=":", alpha=0.5)
ax1.legend(loc="upper left", frameon=False, fontsize=7.3)

cats = ["Direct\nC10$\\to$Asp57", "C10 Relay\nvia W999", "2-Water Wire\nO27$\\to$His128", "1-Water Wire\nO27$\\to$His128"]
x = np.arange(len(cats))
w = 0.34

neutral_means = [35.60, 27.96, 12.39, 10.03]
neutral_sds = [1.10, 0.0, 1.88, 4.02]
hip_vals = [33.66, np.nan, 7.08, 6.74]

ax2.bar(x - w/2, neutral_means, width=w, yerr=neutral_sds, capsize=3, color=C_BLUE, alpha=0.88, label="Neutral $\\mathrm{His124}^0$ (3-seed mean $\\pm$ SD)")
ax2.bar(x + w/2, [v if not np.isnan(v) else 0 for v in hip_vals], width=w, color=C_GREEN, alpha=0.88, label="Cationic $\\mathrm{His124}^+$ (HIP)")

seed_data = [
    [35.98, 34.36, 36.46],
    [27.96],
    [14.56, 11.33, 11.29],
    [7.01, 14.60, 8.49],
]
for i, svals in enumerate(seed_data):
    jitter = np.linspace(-0.07, 0.07, len(svals)) if len(svals) > 1 else [0.0]
    ax2.scatter(np.full(len(svals), x[i] - w/2) + jitter, svals, color="black", s=16, zorder=5)

ax2.set_xticks(x)
ax2.set_xticklabels(cats, fontsize=7.8)
ax2.set_ylabel("Activation barrier $\\Delta E_{\\mathrm{act}}$ (kcal mol$^{-1}$)")
ax2.set_title("(b) Pocket barrier robustness across seeds & $\\mathrm{p}K_{\\mathrm{a}}$", loc="left", fontweight="bold")
ax2.set_ylim(0, 42)
ax2.grid(axis="y", ls=":", alpha=0.5)
ax2.legend(loc="upper right", frameon=False, fontsize=7.5)

fig.tight_layout()
fig.savefig(OUT / "fig2_tcmn_pocket_profiles.pdf", bbox_inches="tight")
plt.close(fig)

# ---------------- FIGURE 3: Stage-1c Non-Peptide Bulky Bases & Metal-Chelated Complexes ----------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 3.0), dpi=300)

systems = ["DBU\n(amidine)", "TMP\n(piperidine)", "DTBMP\n(2,6-$t\\mathrm{Bu}_2$Py)", "$[\\mathrm{MgCl}_2(\\mathrm{acac})]$\n+ DTBMP", "$[\\mathrm{ZnCl}_2(\\mathrm{acac})]$\n+ DTBMP"]
x = np.arange(len(systems))
w = 0.35

dir_means = [3.10, 7.00, 12.22, 10.93, 5.70]
dir_sds   = [1.31, 2.40, 0.75,  5.11,  2.52]
rel_means = [9.37, 12.47, 15.35, 5.87, 3.06]
rel_sds   = [0.39, 0.74,  2.48,  3.33, 2.18]

ax1.bar(x - w/2, dir_means, width=w, yerr=dir_sds, capsize=3, color=C_ORANGE, alpha=0.88, label="Direct ($n=1$ spectator $\\mathrm{H}_2\\mathrm{O}$)")
ax1.bar(x + w/2, rel_means, width=w, yerr=rel_sds, capsize=3, color=C_BLUE, alpha=0.88, label="1-Water Relay ($\\mathrm{C3-H}\\to\\mathrm{H}_2\\mathrm{O}\\to\\mathrm{base}$)")

ax1.set_xticks(x)
ax1.set_xticklabels(systems, fontsize=7.5)
ax1.set_ylabel("Activation barrier $\\Delta E_{\\mathrm{act}}$ (kcal mol$^{-1}$)")
ax1.set_title("(a) Unconstrained non-peptide bases & $M\\mathrm{Cl}_2$ complexes", loc="left", fontweight="bold")
ax1.set_ylim(0, 20)
ax1.grid(axis="y", ls=":", alpha=0.5)
ax1.legend(loc="upper right", frameon=False, fontsize=7.5)

C1C = Path("/mnt/homelab-external/proton-hopping-stage1b/runs/stage1c-nonpeptide")
p_list = [
    ("E1_dtbmp_n1_direct_s1/result.json", "Free acac + DTBMP (Direct, $12.22$)", C_GREY, "--o"),
    ("E3_MgCl2_dtbmp_n1_direct_s2/result.json", "$[\\mathrm{MgCl}_2(\\mathrm{acac})]$ + DTBMP (Direct, $10.58$)", C_ORANGE, "--s"),
    ("E3_MgCl2_dtbmp_n1_relay_s2/result.json", "$[\\mathrm{MgCl}_2(\\mathrm{acac})]$ + DTBMP (Relay, $3.54$)", C_BLUE, "-s"),
    ("E3_ZnCl2_dtbmp_n1_relay_s2/result.json", "$[\\mathrm{ZnCl}_2(\\mathrm{acac})]$ + DTBMP (Relay, $1.79$)", C_GREEN, "-D"),
]
for rel_p, lab, col, fmt in p_list:
    r = json.load(open(C1C / rel_p))
    pts = r["points"]
    emin = min(p["dE_kcal"] for p in pts[:4])
    xs = np.linspace(0, 1, len(pts))
    ys = np.array([p["dE_kcal"] - emin for p in pts])
    ax2.plot(xs, ys, fmt, color=col, lw=1.7, ms=4.2, label=lab)

ax2.set_xlabel("Normalized proton-transfer coordinate $\\xi$")
ax2.set_ylabel("$\\Delta E - \\Delta E_{\\mathrm{R,min}}$ (kcal mol$^{-1}$)")
ax2.set_title("(b) Metal-chelated + bulky-base relay rescue", loc="left", fontweight="bold")
ax2.set_ylim(-26, 16)
ax2.axhline(0, color="black", lw=0.6, ls=":")
ax2.grid(axis="y", ls=":", alpha=0.5)
ax2.legend(loc="lower left", frameon=False, fontsize=7.0)

fig.tight_layout()
fig.savefig(OUT / "fig3_nonpeptide_complexes.pdf", bbox_inches="tight")
plt.close(fig)
print("Generated fig1, fig2, and fig3 PDFs with Delta E_act notation.")
