# Handoff for Codebuff: antibiotics paper (proton-relay / proton hopping)

Date: 23 November 2026. Public repo, MIT licence. Verified facts only.

## Hard constraint: compute credits

Julian's compute credits are running out. No heavy or long compute jobs, no new GPU/cloud runs and no large re-runs without Julian's explicit OK. Prefer cheap checks (xtb, small tests) and report expected cost and runtime before starting anything heavy.

## Claim

This work tests whether a Grotthuss-type proton relay over 1–4 explicit shuttle molecules (water or methanol) lowers the barrier of the base-driven deprotonation step in polyketide ring closure, compared with direct transfer. At the GFN2-xTB level, stage-1 NEB and proton-only relaxed scans find the relay is not favoured; once a heavy-atom coordinate is included, the comparison is open. The work is mechanism-level only and uses open-source software only.

## Repo and manuscript state

- Repo: https://github.com/sailblack-research/proton-hopping-mechanism-sim (public, MIT).
- `main` is at `a865feb` (merge of PR #1: MIT licence, derived structures with provenance, stage-1 pipeline). CI on that head is green on the self-hosted Linux runner (runner left unnamed here).
- PR #2 is an open draft (stage 2 DFT check, scripts only, head `1a1796b`).
- PR #3 is an open draft (stage-1b relay-design screen, head `cface5b`).
- No manuscript yet. The repo has the plan (README) plus `docs/literature.md` and `docs/open-questions.md`. paper-review has not been run, so there is no score. A separate "design note" named in an earlier internal draft was not found as a file in this repo.
- Nothing is ready to merge. Any merge needs Julian's yes on a specific green head.

## Model systems (no enzyme so far)

All results below are from small model systems, not from an enzyme:

- model substrate fragment (acetylacetone C3–H in the stage-1 setup);
- a small base (trimethylamine; stage 1 also used hydroxide and methoxide);
- 0–4 explicit water or methanol shuttles;
- implicit solvent (ALPB water);
- GFN2-xTB, 3 seeds.

There is no enzyme, no construct, no mutations or truncations, and no protein structure. Structures were derived from SMILES, with provenance in the repo.

What this can say: whether a relay is intrinsically preferred over direct transfer for a free shuttle chain in implicit solvent. Stage-1 NEB and proton-only scans find it costlier; with a heavy-atom coordinate included, that conclusion no longer holds robustly.

What this does not say about an enzyme active site. The model has no fixed, preorganised water wire, no residue electrostatics or pKa tuning, no enzyme-imposed base–substrate distance, and no protein dynamics. Those factors could favour a relay inside an enzyme, so the present results do not rule out an enzymatic relay.

## Stage 1 (GFN2-xTB, CI-NEB)

324 jobs, 102 converged (`ok`). Status counts already in the repo: 15 NEB not converged, 100 reactant not a minimum, 105 product not a minimum, 2 SCF errors (100 + 105 = 205 endpoint-collapse cases).

Example, trimethylamine in water, mean ± sd over converged seeds, electronic energies, kcal/mol (no ZPE, no entropy):

- n = 3: relay 15.83 ± 3.19 vs direct 1.97 ± 0.53
- n = 1: effectively a tie, 6.41 ± 0.72 vs 6.65 ± 0.11

The full NMe3 table (water and methanol, n = 0–4) is in the README and in `data/stage1_summary.json`. At this screening level the relay is not favoured over direct transfer.

## Stage 1b relay-design screen

660 jobs finished: 253 ok, 110 NEB not converged, 155 geometry infeasible, 82 product not a minimum, 4 reactant not a minimum, 35 error, 21 timeout.

Not yet analysed against the pre-registered pass rule (relay at least 2 kcal/mol lower, same sign across 3 functionals, barrier at most 23 kcal/mol). Julian has not approved evaluating this screen or updating PR #3.

## DFT check (incomplete)

Level: ωB97X-D3(BJ)/def2-SVP, PCM water, def2-TZVP single points, ZPE+RRHO at 298 K, PySCF/gpu4pyscf.

- Trimethylamine / water / n = 1 relay: one valid TS (one imaginary mode, −210 cm⁻¹): ΔE‡ 4.89, ΔE0‡ 6.22, ΔG‡ 8.25; TZVP//SVP ΔE‡ 3.25, ΔG‡ 6.61 kcal/mol.
- n = 0 relay TS is invalid (two imaginary modes; TS ≈ product).
- No DFT barrier for any direct system yet, so “relay not favoured” is not tested at DFT level.

GPU-path TS steps fail from a dtype (object) bug in the TS-start Hessian. Frequency steps fail from a too-small GPU memory cap. The existing DFT queue was left running and, with those bugs, yields reactant and product structures rather than further valid TS or frequency results. Julian has not approved stopping or fixing that queue.

## Relaxed scans (method switch, 23 Nov 2026)

Julian approved (23 Nov 2026) switching from further static TS searches to relaxed scans along the proton coordinate(s) (xtb first, then DFT for selected systems) plus short GFN2-xTB metadynamics.

All scan values below are GFN2-xTB/ALPB water, NMe3/water, seeds s0/s1/s2, kcal/mol, electronic energies only (no entropy).

M1 (coarse 11×11 proton-only grid), direct vs relay:

- n = 1: 6.78 ± 0.03 vs 15.43 ± 1.27
- n = 2: 5.62 ± 0.06 vs 13.44 ± 0.58
- n = 3: 4.25 ± 0.16 vs 16.72 ± 0.57
- n = 0 (contact pair): 7.87 ± 0.03, with no product minimum

M1b (proton-only antisymmetric 1D scan q = Σ[d(D–H) − d(H–A)], 31 points, continuation), direct vs relay:

- n = 1: 6.89 ± 0.04 vs 13.87 ± 1.65
- n = 2: 5.93 ± 0.16 vs 16.27 ± 0.26
- n = 3: 5.35 ± 0.60 vs 14.94 ± 1.02

21×21 continuation grid, relay:

- n = 1: 12.66 ± 2.78
- n = 2: 11.04 ± 3.52
- n = 3: 15.17 ± 1.89

Why these sit above the stage-1 NEB: at the stage-1 NEB maximum (n = 1, seed s1) both protons have already moved; the peak is heavy-atom reorganisation after transfer, so proton-only scans are upper bounds.

M1c (2D: q × first-hop donor–acceptor heavy-atom distance R_DA, 2.55–3.75 Å, both sweeps, lower-of-both bottleneck):

- n = 1: direct 6.53 ± 0.26 vs relay 7.65 ± 3.03 (seeds 9.16 / 9.63 / 4.17)
- n = 2: direct 5.06 ± 0.22 vs relay 4.34 ± 1.16 (seeds 4.66 / 3.06 / 5.31)

Relay depends strongly on sweep direction (3.6–20.9 kcal/mol between sweeps) while direct agrees within about 2 kcal/mol, so hidden wire/H-bond coordinates matter. Several relay bottlenecks sit at the grid edge R_DA = 3.75 Å. The lower-of-both value may underestimate. Energies only, no entropy.

Conclusion: once heavy atoms may move, “relay not favoured” no longer holds robustly (n = 1 is inconclusive; n = 2, the relay is lower in 2 of 3 seeds and in the mean). Static xtb scans cannot settle relay vs direct for these clusters; sampling along q and R_DA (metadynamics or umbrella sampling) is needed. xtb's native RMSD metadynamics can give transfer events and times, not a proper free-energy profile. A short pilot (2 ps, 1 seed) took about 40 s and saw no transfers. The full run (2 systems × 3 seeds × 100 ps) is estimated at about 3.6 CPU-hours and awaits Julian's OK. DFT relaxed scans are estimated at 25–115 CPU-hours for 4 systems and are not proposed under the credit limit.

## Convergence

Failures are mostly not SCF or geometry-optimisation failures (those rarely fail). They are mostly endpoint collapse — the proton moves spontaneously while the endpoints relax (205 cases in stage 1, 86 in stage 1b) — and hydrogen-bond wire rewiring, i.e. a flat multi-well surface. That is why the method was switched to relaxed scans plus short metadynamics.

## New direction: build outward from the enzyme pocket

Julian's direction (23 Nov 2026): start from the enzyme pocket and build outward. Begin with a real enzyme model built from a crystal structure of the antibiotics-paper target, capturing the pocket, the held water wire, and the nearby residues that tune acidity and the base–substrate distance (QM cluster or QM/MM; open-source software only). Choose a dynamic method suited to that pocket-scale system (metadynamics or ab initio MD).

A proposal is pending; no heavy compute until Julian approves it. The proposal must report, before any run: the chosen PDB entry (ID, resolution, organism, ligand state, mutations or missing residues and how they would be handled), the cluster/QM region definition and size, protonation choices, the method and level of theory, collective variables, and estimated runtime and cost. Respect the compute-credit constraint above.

## Where data lives

Raw run data live on the research workstation's external data volume, with an encrypted offsite backup plus a local mirror. Scan results and DFT outputs are not in the repo. The repo holds code, the small stage-1 summaries under `data/`, and docs. Do not commit secrets, tokens, hostnames, bucket names, or internal paths.

## Open approvals (waiting on Julian)

- Stopping or fixing the DFT queue.
- Stage-1b analysis and any update to PR #3.
- The 2 systems × 3 seeds × 100 ps xtb sampling along q and R_DA (about 3.6 CPU-hours).
- Enzyme-model compute (after the proposal).
- Any merge, and only on Julian's yes on a specific green head.
- Any heavy or long compute, new GPU/cloud run, or large re-run (see the credit constraint). DFT relaxed scans (estimated 25–115 CPU-hours for 4 systems) are not proposed under the credit limit.

## Next steps

1. Sample along q and R_DA (metadynamics or umbrella sampling), and widen the R_DA range (several relay bottlenecks sit at the grid edge 3.75 Å). The 2 systems × 3 seeds × 100 ps xtb run (about 3.6 CPU-hours) awaits Julian's OK.
2. Do not start DFT relaxed scans; they are estimated at 25–115 CPU-hours for 4 systems and are not proposed under the credit limit.
3. The enzyme-pocket proposal is pending; wait for approval before any heavy compute.
4. After results warrant it: draft the manuscript, then run paper-review and iterate to a score of at least 8.5.
5. Do not merge PR #2 or PR #3 unless Julian says yes on a named green head.

## Rules

- No merges without Julian's yes on a specific green head.
- Public repo: no secrets, tokens, hostnames, bucket names, or internal paths.
- Mechanism-level work only; open-source software only (no ORCA or Gaussian).
- No invented numbers or citations.
