# Designing a system in which the proton relay beats direct transfer: stage-1b design and screening plan

Researcher, 22.11.2026, for Julian M. Kleber. Repo: `sailblack-research/proton-hopping-mechanism-sim`, branch `stage1b-relay-design`. Mechanism level only: common solvents, reagents and textbook model compounds. No synthesis routes or hazardous structures. Strictly open-source software (tblite/xtb, CREST, ASE, RDKit, PySCF, geomeTRIC). **Nothing beyond the single pilot job has run.** The full screen waits for your go.

## 0. Why stage 1 came out against the relay

Stage 1 (GFN2-xTB/ALPB, tblite, ASE CI-NEB, 324 jobs, 3 seeds; `data/stage1_summary.json`) found the relay not favoured. All numbers below are electronic energies from that run, mean ± sd over converged seeds:

| condition (deprotonation of acac C3–H) | relay | direct (same composition) |
|---|---|---|
| NMe3, water, n = 1 | 6.41 ± 0.72 | 6.65 ± 0.11 |
| NMe3, water, n = 3 | 15.83 ± 3.19 | 1.97 ± 0.53 |
| NMe3, methanol, n = 1 | 12.56 ± 1.93 | 8.71 ± 0.29 |
| NMe3, methanol, n = 3 | 23.78 ± 10.57 | 4.47 ± 2.15 |

With OH⁻ and MeO⁻, the direct reactant is not a minimum: the proton leaves the C–H as soon as the base touches it (status `reactant_not_minimum`). All 90 reprotonation jobs with those bases ended as `product_not_minimum`, because the alkoxide and the solvent anion have similar basicity.

Reading: stage 1 gave direct transfer the easiest case, a free, strong base in contact with the C–H. Every added wire molecule then costs ordering and charge separation and buys nothing. The one hint the other way is the single-water case in water (n = 1), where relay and direct agree within their spread. That matches the single water used as base or proton carrier in the ARO/CYC cyclases TcmN and WhiE (Ames 2008; Lee 2012). So, for the relay to win, the system has to **remove the direct path, make it expensive, or make the relay cheaper**. The five levers below each do one of these.

## 1. The five levers: rationale, falsifiable prediction, minimal model

| # | Lever | Mechanistic rationale (literature) | Falsifiable prediction | Minimal model (family) |
|---|---|---|---|---|
| 1 | **Base held off the C–H** (tether or confinement) | In enzymes, a base that cannot reach the substrate works through water wires. Examples are the rate-limiting proton transfer through water in carbonic anhydrase (Silverman & Lindskog 1988) and long-range transfer through water chains with "proton holes" in solution and enzymes (Riccardi et al. 2006). An intramolecular tether sets the effective molarity and geometry of the base (Kirby 1980). | The direct barrier B_dir(d) rises steeply once the C···base distance d exceeds contact. The best relay barrier min_n B_relay(n, d) stays bounded, so there is a crossover d* above which the relay wins by ≥ 2 kcal/mol. **Falsified** if no d in 4–7 Å gives B_relay ≤ B_dir − 2 kcal/mol with B_relay ≤ 20 kcal/mol. | **A**: acac + base (NMe3, N-methylimidazole, acetate) + (H2O)n or (MeOH)n, n = 0–3. The C3···base distance is fixed at d = 3, 4, 5, 6, 7 Å or left free (ASE `FixBondLength`). |
| 2 | **Weaker, pKa-matched base and bifunctional shuttles** | A carbon acid whose conjugate base is delocalised has a high intrinsic barrier from imbalanced charge delocalisation (Bernasconi 1987). A bifunctional catalyst that accepts and donates a proton in one step avoids the ion pair. The classic example is 2-pyridone in tetramethylglucose mutarotation (Swain & Brown 1952). Water can play the same role in the keto–enol step of acetylacetone (Alagona et al. 2010). HFIP is a strong H-bond-donor solvent (Colomer et al. 2017). | With the bifunctional shuttles (2-pyridone, acetic acid), the concerted cyclic relay (C–H → shuttle → C=O, neutral enol) is lower than the stepwise path (C–H → shuttle only, ion pair). It is also lower than the water-assisted path. The non-bifunctional donors (phenol, HFIP, methanol) give a smaller or no relay advantage. **Falsified** if, for both bifunctional shuttles in both solvents, stepwise ≤ relay − 0 or relay is not below water-assisted. | **B**: substrate (acac; cyclohexa-2,4-dien-1-one) + shuttle ∈ {none, H2O, (H2O)2, MeOH, AcOH, 2-pyridone, phenol, HFIP}, topology relay vs stepwise. The acidity order of all donors comes from `pka_ladder.py`, an isodesmic GFN2/ALPB ordering. |
| 3 | **Preorganised wire** | Most of the relay penalty in stage 1 is the ordering of loose solvent molecules. Armstrong et al. (2014) flag this entropy caveat for computed water relays. Built-in H-bond chains carry multi-proton transfers in bioinspired molecular relays (Odella et al. 2018; Hammes-Schiffer 2015 review). | A covalently preorganised three-OH chain lowers the relay barrier relative to a loose (H2O)3 wire with the same number of hops. The relay-competent arrangement is also the most stable reactant conformer (CREST). **Falsified** if the triol relay barrier is not lower than the loose (H2O)3 relay barrier (stage-1 / family A n = 3 natural) by ≥ 2 kcal/mol, or if CREST finds no relay-competent conformer within 3 kcal/mol. | **C**: acac + all-cis cyclohexane-1,3,5-triol (three 1,3-diaxial OH) + NMe3 or acetate, relay (4 protons) vs direct with the triol as an H-bond spectator (same composition). |
| 4 | **The step after ring closure** | The aromatising keto→enol tautomerisation after aldol closure and dehydration can't happen by a direct intramolecular 1,3-shift; it needs an external acid/base or a shuttle. The phenol ⇌ 2,4-cyclohexadienone equilibrium was measured in water (Capponi et al. 1986). Stage 1 could not test the reprotonation, because the alkoxide endpoints were not minima. | For the dienone → phenol step, a one- or two-molecule shuttle relay beats both the intramolecular shift and the stepwise ion-pair path. It does so by more than for the acac keto→enol step, because aromatisation drives the product side. **Falsified** if the dienone relay advantage is ≤ the acac relay advantage, or if stepwise ≤ relay. | **B** (dienone rows), plus reprotonation with a weak acid moved to stage 2 |
| 5 | **Explicit second shell and tunnelling/KIE as discriminators** | Implicit solvent alone overstates the ordering penalty and leaves out the charge stabilisation of a wire. Primary H/D KIEs with tunnelling, and proton inventories, separate one-proton from multi-proton TSs (Venkatasubban & Schowen 1984; Bigeleisen & Mayer 1947). Tunnelling factors follow Wigner 1932, Bell's parabolic barrier (Skodje & Truhlar 1981) and the unsymmetrical Eckart barrier (Eckart 1930; Johnston & Heicklen 1962). | (i) Adding k = 2 or 4 second-shell waters does not flip the sign of ΔB = B_relay − B_dir for the hits from levers 1–4. (ii) A relay TS has ≥ 2 protons in flight, so KIE(all-D)/KIE(C–H-only D) > 1 and the proton inventory is non-linear. A direct TS gives a ratio of ~1. **Falsified** for a hit if the sign flips with the shell, or if the relay TS is effectively a one-proton TS (ratio ~1). | **D**: hits of A–C re-run with k = 2, 4 explicit second-shell waters. Every converged job gets screening-level KIE and tunnelling factors from partial Hessians over the moving protons. |
| ctrl | **Positive control** | 7-Hydroxyquinoline: the O–H and the ring N are too far apart for a direct transfer. Bridged by an (NH3)3 wire, it shows excited-state H transfer along the wire (Tanner et al. 2003). | The protocol must find a connected relay path for n = 2–3 (ground-state enol → keto tautomer) while n = 0 fails. If it does not, negative results in A–C are **inconclusive** (method failure), not a kill. | **Q**: 7-HQ + (H2O)n or (MeOH)n, n = 0–3 |

## 2. Ranked shortlist (SMILES derived by Researcher, validated with RDKit; PubChem CIDs via PUG-REST on 22.11.2026)

`stages/stage1b_relay_design/validate1b.py` checks each species: RDKit sanitisation, formula, monoisotopic mass against PubChem, and an all-cis embedding check for the triol. Result: all 13 species OK (`species.yaml`). Atom-map numbers mark the reactive atoms.

| Rank | Candidate system | Components (SMILES, PubChem CID) | Why it is ranked here |
|---|---|---|---|
| **1** | **A: held-off base + short wire**: NMe3 / N-methylimidazole / acetate at d = 4–7 Å from acac C3, bridged by 1–3 H2O or MeOH | acac `CC(=O)[CH2:1]C(C)=[O:2]` (31261); NMe3 `C[N:1](C)C` (1146); N-methylimidazole `CN1C=C[N:1]=C1` (1390); acetate `CC(=O)[O-:1]` (175); water (962); methanol (887) | Most direct test of the enzyme-like scenario (base kept away). It reuses the stage-1 deprotonation, and n = 1 in water was already competitive in stage 1. Acetate is the carboxylate type of base of SnoaL (Sultana et al. 2004). N-methylimidazole is a weaker, pKa-closer base than NMe3. |
| **2** | **B-dienone: aromatising tautomerisation with a bifunctional shuttle**: cyclohexa-2,4-dien-1-one → phenol via 2-pyridone or acetic acid (vs H2O, (H2O)2, MeOH, phenol, HFIP) | dienone `[O:2]=C1[CH2:1]C=CC=C1` (141125); 2-pyridone `[O:1]=C1C=CC=C[NH:2]1` (8871); acetic acid `C[C](=[O:1])[OH:2]` (176); phenol (996); HFIP `FC(F)(F)C([OH:1])C(F)(F)F` (13529) | A real step of the polyketide ring-closure sequence (lever 4). The direct intramolecular shift is geometrically poor, and the bifunctional shuttle avoids charge separation (lever 2). Most likely place for the relay to win without any enforced geometry. |
| **3** | **B-acac: keto→enol with a bifunctional shuttle** | acac (as above) + the same shuttles | Benchmark against the published water-assisted acac tautomerisation (Alagona et al. 2010). Separates the shuttle effect from the aromatisation effect of rank 2. |
| **4** | **C: preorganised triol wire**: acac C3–H → O1–H → O2–H → O3–H of all-cis cyclohexane-1,3,5-triol → NMe3 or acetate | triol `[OH:1][C@H]1C[C@@H]([OH:2])C[C@@H]([OH:3])C1` (230351; the all-cis isomer confirmed by embedding) | Tests whether removing the ordering penalty (lever 3) rescues the 3-hop relay that lost by about 14 kcal/mol in stage 1. Ranked lower because it is the most artificial and its 4-proton TS is the hardest to converge. |
| ctrl | **Q: 7-hydroxyquinoline + (H2O)n or (MeOH)n** | 7-HQ `[OH:1]c1ccc2ccc[n:2]c2c1` (135426866) | Positive control of the method, not a candidate. |

Lever 5 (family D) is not a separate candidate. It is the robustness and discrimination layer applied to every hit of ranks 1–4.

## 3. Screening protocol (stage 1b, all on engine, `nice -n 19 ionice -c3`, data disk `/mnt/homelab-external/proton-hopping-stage1b/`)

1. **Species and builds.** `validate1b.py` covers RDKit and PubChem checks. `build1b.py` builds the clusters: H-bonded chain placement as in stage 1, the fixed d through a proportional translation along the C→base axis, a cyclic closure onto the substrate C=O for family B, and a triaxial triol conformer with the O–H torsions set to the 1,3-chain. Compressed wires (n ≥ 1 at small d) are recorded as `geometry_infeasible` and never run.
2. **Acidity ladder (lever 2).** `pka_ladder.py --solvent water|methanol` gives isodesmic GFN2/ALPB proton-exchange energies relative to water for all donors and conjugate acids. It's an ordering only, used to label the pKa match. No pKa numbers are claimed from it.
3. **Endpoints and path** (`run1b.py`, one process per job, OMP 1 thread): GFN2-xTB/ALPB (tblite) BFGS relaxation of R and P (with the `FixBondLength` constraint in family A), then the stage-1 endpoint classification. Then CI-NEB, 11 images, improved tangent, FIRE: plain NEB to fmax 0.10 eV/Å, then climbing image to 0.05 eV/Å, with the stage-1 SCF-fallback order idpp/300 K → linear/300 K → idpp/1000 K. Job cap 2400 s, checkpointed per job.
4. **Kinetics layer (lever 5, screening).** Central-difference partial Hessian over the moving protons at R and at the highest image. From it: the imaginary frequency, the ZPE-only H/D KIE for all moving protons deuterated and for the C–H proton only, the Wigner, Bell and Eckart factors, and the tunnelling-corrected KIE (`kinetics.py`, unit-tested against the Wigner limit and the H/D frequency ratio).
5. **Conformers (CREST 3.0.2, `--nci --quick`, ALPB, 4 threads).** Only for the screen hits (≤ 10 reactant clusters) and the pilot. Output: lowest conformer, number of conformers within 3 kcal/mol, and the fraction that are relay-competent (all H···acceptor links < 2.6 Å). For each hit, the path is re-run from the lowest relay-competent conformer.
6. **Grid** (`run1b.py --list`): A 432 (3 bases × 2 solvents × n 0–3 × 6 d × 3 seeds), B 180, C 24, Q 24 = **660 jobs**. D: ≤ 8 hit conditions × 2 topologies × k ∈ {2, 4} × 3 seeds ≤ 96 jobs.
7. **DFT refinement (stage 2 code, PySCF + geomeTRIC, already in the repo as `stages/stage2_dft_check`).** For the hits only (≤ 8 relay/direct pairs): ωB97X-D3(BJ)/def2-SVP optimisation, TS from the climbing image with the analytic Hessian, RRHO thermochemistry with quasi-RRHO entropies (Grimme 2012), def2-TZVP single points, IEF-PCM. Robustness: M06-2X and B3LYP-D3(BJ) single points, and implicit vs implicit + explicit second shell (family D geometries).

## 4. Success criteria and rejection rule (fixed now, before any screen run)

Barrier of a condition = median of `barrier_fwd_kcal` over seeds with status `ok`. A condition counts only if ≥ 2 of 3 seeds are `ok` and the TS is not an endpoint. Margin δ = 2.0 kcal/mol (the plan's ±2–3 kcal/mol method error; a factor of 29 in rate at 298 K).

**Screen hit (xtb level), per family:**
- A: at some d ∈ {4, 5, 6, 7} Å, min_n B_relay(n, d) ≤ B_dir(d) − δ **and** min_n B_relay(n, d) ≤ 20 kcal/mol (Eyring t½ ≈ 51 s at 298 K, if ΔG‡ ≈ ΔE‡).
- B: for a shuttle and solvent, B_relay ≤ B_stepwise − δ **and** B_relay ≤ B_relay(H2O, n = 1) − δ (bifunctional advantage). For the dienone: also B_relay ≤ B(intramolecular shift) − δ.
- C: B_relay(triol) ≤ B_dir(triol) − δ **and** B_relay(triol) ≤ B_relay(loose (H2O)3, family A d = natural) − δ.
- Q (control): ≥ 2 of 3 seeds `ok` for n = 2 or 3 in at least one solvent, and n = 0 not `ok` or higher by ≥ δ. If Q fails, the method is fixed first, and negative A–C results are inconclusive.

**Relay supported (DFT), per system:** ΔG‡_relay ≤ ΔG‡_direct − δ at ωB97X-D3(BJ)/def2-TZVP//def2-SVP + PCM + qRRHO, **and** the same sign of ΔΔG‡ with M06-2X and B3LYP-D3(BJ), **and** the same sign with the explicit second shell (k = 2 and 4), **and** ΔG‡_relay ≤ 23 kcal/mol (t½ ≈ 2.2 h at 298 K), **and** the relay TS has ≥ 2 protons with a substantial share of the imaginary mode (KIE(all-D)/KIE(C–H-only D) > 1).

**Kill rule:** if Q passes and no system of A–C passes the DFT criterion, the relay hypothesis (H1) is **rejected** for base-driven poly-β-carbonyl ring closure at this level. That negative result gets written up as such. If only A passes, the conclusion is that the relay needs enforced geometry: a base held off the C–H, as in an enzyme pocket. H1 for a free base in solution is then rejected, and H1 is restated as an enzyme-pocket hypothesis (then the WhiE/TcmN cluster model in stage 2). If B-dienone passes, the relay is supported for the aromatising tautomerisation, which is a narrower claim than H1. The paper then says exactly that.

## 5. What would support or kill the idea

- **Support:** a d* ≤ 5 Å crossover in A that survives DFT and the second shell. Or a 2-pyridone/AcOH relay in B-dienone that beats stepwise by ≥ δ at all three functionals, with a two-proton TS (KIE ratio > 1). Testable lab predictions follow: a primary KIE plus a non-linear proton inventory, a rate that depends on shuttle identity at constant base strength, and a stronger effect of 2-pyridone than of phenol at matched acidity.
- **Kill:** Q passes, but every A–C hit vanishes at DFT, or flips sign with the explicit shell. Or the relay wins only at d where its own barrier is > 23 kcal/mol, so it is a relay but too slow to matter.

## 6. Pilot job and cost estimate

**Pilot (the only stage-1b job that has run).** One job, `C_triol_nme3_water_relay_s0`: family C, the largest cluster class (49 atoms, 4 protons in flight). It ran once on engine from 20:30 to 21:10 CET (22.11.2026), with `nice -n 19 ionice -c3`, one OpenMP thread for the path and 4 threads for CREST, on the data disk (`/mnt/homelab-external/proton-hopping-stage1b/runs/pilot-20261122/`). Engine was heavily loaded the whole time (load average about 34–37 on 32 cores, including the separate DFT check), so the timings are conservative.

| Step | Measured |
|---|---|
| Build, R and P relaxation (124 and 107 BFGS steps), CI-NEB (1338 FIRE steps, idpp, first attempt), partial Hessians | **347 s CPU** (one thread, about 6 min wall); the kinetics layer took 1.5 s of that |
| CREST 3.0.2 `--nci --quick`, ALPB water, 4 threads | **34 min wall, not finished.** The 2400 s per-job cap stopped it. Its log shows about 83 CPU-min in the four MD/metadynamics segments alone (8.2 + 33.4 + 8.0 + 33.8 min), with 591 structures in the 10 kcal/mol window still being re-optimised. So CREST on a 49-atom cluster costs **more than 2.3 CPU-h** with `--quick`. |

What the pilot showed about the method (one seed, so not a result):
- R and P were proper endpoints. All four protons sat on their donors in R and on their acceptors in P (ΔE_rxn = +4.8 kcal/mol, GFN2-xTB electronic).
- The single-band CI-NEB for the fully concerted 4-proton path peaked 86.4 kcal/mol above R. At that image the partial Hessian had **4 imaginary modes**, so it is not a first-order saddle. Forcing four protons to move synchronously in one band produced an artefact, not a relay barrier. The partial-Hessian KIEs from this image (42 and 22) are therefore invalid and are not used.
- Consequences, now in the code: (i) a **sequential-hop (segmented) path** for relays. It relaxes the intermediates with the first k protons moved and runs one CI-NEB per hop. It is triggered when a multi-proton band exceeds 40 kcal/mol or does not converge. (ii) KIEs are flagged `kie_valid` only when the TS partial Hessian has exactly one imaginary mode. (iii) CREST no longer runs inside the per-job cap, which now covers the path only, and the path status and timings are recorded before CREST starts. In the pilot, the timeout overwrote the NEB status. The segmented path has **not** been run on engine, since only one pilot job was allowed. It is covered by a mock-calculator unit test.

**Cost estimate for the full stage-1b screen.** Per-job cost comes from measured times. For clusters under 41 atoms, it is the stage-1 90th-percentile wall time per size bin over the 113 stage-1 jobs that reached NEB: 15.5 s (< 24 atoms), 76.6 s (24–30), 44.5 s (31–34) and 150.4 s (35–40). For 41 atoms or more, it is the pilot's 347 s. Atom counts come from building every job locally.

| Family | Jobs | Built / infeasible (skipped in ms) | Atoms | Estimate (CPU-h) | With segmented fallback doubling every multi-proton job |
|---|---|---|---|---|---|
| A | 432 | 277 / 155 | 22–46 | 5.9 | 10.0 |
| B | 180 | 180 / 0 | 13–28 | 2.0 | 3.0 |
| C | 24 | 24 / 0 | 43–49 | 2.3 | 3.5 |
| Q | 24 | 24 / 0 | 18–36 | 0.4 | 0.8 |
| **A–C + Q** | **660** | | | **10.6** | **17.3** |

So the screen (A–C + Q, no CREST) is expected to take **about 11–17 CPU-h, under the 24 h limit**. With 8 workers, that is about 1.5–2.5 h wall if cores are free. Every job is also capped at 2400 s, but that cap is a safety limit, not a forecast. Separate, hit-dependent follow-ups, each needing its own go:
- **D** (second shell, at most 96 jobs): upper bound about 18.5 CPU-h (96 × twice the pilot path cost).
- **CREST** on the hits: more than 2.3 CPU-h per 49-atom cluster (measured lower bound). The proposal is at most 6 hits with a 2 h wall cap at 4 threads (at most 48 CPU-h), or `--squick`.
- **DFT** on at most 8 relay/direct pairs: the separate DFT check now running on engine is a 31-atom stage-1 pair at ωB97X-D3(BJ)/def2-SVP on 16 threads. Its reactant optimisation alone took 48 min wall (20:08–20:56 CET), and the process had used 12.7 CPU-h after 56 min. The cost per full pair is known once that check finishes. DFT will dominate the budget, at days of engine CPU or the Colab A100 route from the paper plan.

## 7. References (each DOI resolved with the Crossref API on 22.11.2026; title/journal/year match)

1. Ames, B. D. et al. *PNAS* 2008, doi:10.1073/pnas.0709223105 (TcmN ARO/CYC structure).
2. Lee, M.-Y.; Ames, B. D.; Tsai, S.-C. *Biochemistry* 2012, doi:10.1021/bi201705q (WhiE-ORFVI).
3. Sultana, A. et al. *EMBO J.* 2004, doi:10.1038/sj.emboj.7600201 (SnoaL).
4. Silverman, D. N.; Lindskog, S. *Acc. Chem. Res.* 1988, doi:10.1021/ar00145a005 (carbonic anhydrase, rate-limiting protolysis of water).
5. Riccardi, D. et al. *J. Am. Chem. Soc.* 2006, 128, 16302, doi:10.1021/ja065451j ("proton holes" in long-range proton transfer).
6. Cukierman, S. *Biochim. Biophys. Acta Bioenerg.* 2006, doi:10.1016/j.bbabio.2005.12.001 (Grotthuss review).
7. Kirby, A. J. *Adv. Phys. Org. Chem.* 1980, 17, doi:10.1016/S0065-3160(08)60129-X (effective molarities).
8. Bernasconi, C. F. *Acc. Chem. Res.* 1987, doi:10.1021/ar00140a006 (intrinsic barriers, nonperfect synchronisation).
9. Swain, C. G.; Brown, J. F. *J. Am. Chem. Soc.* 1952, doi:10.1021/ja01130a024 (Concerted displacement reactions VIII, polyfunctional catalysis).
10. Alagona, G.; Ghio, C.; Nagy, P. I. *Phys. Chem. Chem. Phys.* 2010, doi:10.1039/c003999c (water-assisted keto–enol, acetylacetone).
11. Colomer, I. et al. *Nat. Rev. Chem.* 2017, doi:10.1038/s41570-017-0088 (HFIP).
12. Armstrong, A. et al. *Chem. Sci.* 2014, doi:10.1039/c3sc53416b (Houk–List TSs revisited; water relay and entropy).
13. Odella, E. et al. *J. Am. Chem. Soc.* 2018, doi:10.1021/jacs.8b09724 (bioinspired multi-proton relays).
14. Hammes-Schiffer, S. *J. Am. Chem. Soc.* 2015, doi:10.1021/jacs.5b04087 (PCET perspective).
15. Capponi, M.; Gut, I.; Wirz, J. *Angew. Chem. Int. Ed. Engl.* 1986, doi:10.1002/anie.198603441 (phenol ⇌ 2,4-cyclohexadienone in water).
16. Tanner, C. et al. *Science* 2003, doi:10.1126/science.1091708 (7-hydroxyquinoline·(NH3)3 wire; excited-state H transfer).
17. Venkatasubban, K. S.; Schowen, R. L. *CRC Crit. Rev. Biochem.* 1984, doi:10.3109/10409238409110268 (proton inventory).
18. Bigeleisen, J.; Mayer, M. G. *J. Chem. Phys.* 1947, doi:10.1063/1.1746492 (isotope exchange equilibria).
19. Wigner, E. *Z. Phys. Chem.* 1932, 19B, 203, doi:10.1515/zpch-1932-1920 (Crossref also lists 10.1515/zpch-1932-0120 for the same title; the first one carries the page range).
20. Eckart, C. *Phys. Rev.* 1930, doi:10.1103/PhysRev.35.1303.
21. Johnston, H. S.; Heicklen, J. *J. Phys. Chem.* 1962, 66, 532, doi:10.1021/j100809a040 (unsymmetrical Eckart tunnelling).
22. Skodje, R. T.; Truhlar, D. G. *J. Phys. Chem.* 1981, doi:10.1021/j150606a003 (parabolic tunnelling).
23. Grimme, S. *Chem. Eur. J.* 2012, doi:10.1002/chem.201200497 (quasi-RRHO).
24. Bannwarth, C.; Ehlert, S.; Grimme, S. *J. Chem. Theory Comput.* 2019, doi:10.1021/acs.jctc.8b01176 (GFN2-xTB).
25. Ehlert, S. et al. *J. Chem. Theory Comput.* 2021, doi:10.1021/acs.jctc.1c00471 (ALPB).
26. Pracht, P.; Bohle, F.; Grimme, S. *Phys. Chem. Chem. Phys.* 2020, doi:10.1039/c9cp06869d (CREST).
27. Sun, Q. et al. *J. Chem. Phys.* 2020, doi:10.1063/5.0006074 (PySCF).
28. Wang, L.-P.; Song, C. *J. Chem. Phys.* 2016, doi:10.1063/1.4952956 (geomeTRIC).

Caveats: GFN2-xTB electronic barriers carry no ZPE or entropy. Partial-Hessian KIEs are screening-level only. ALPB has no HFIP parameters, so HFIP runs as an explicit shuttle in ALPB water or methanol. The 7-HQ precedent is an excited-state experiment, and here it serves only as a geometric control.
