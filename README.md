# proton-hopping-mechanism-sim

Simulation research project for Julian M. Kleber's proposed **proton-hopping (Grotthuss-type) mechanism**
for the base-induced polyketide cyclisation cascades behind tetracycline-type ring systems
(origin: *F5-2 API Synthesis Plan*, 2020, DOI 10.13140/RG.2.2.23579.28965; rewritten preprint draft in
`cap-jmk-launchpad/proton-hopping-preprints`).

> Status: structures derived and validated; stage 1 (GFN2-xTB model relays) running/finished on engine, see
> `data/stage1_summary.json`. Stage-1 numbers are semi-empirical electronic energies (no ZPE, no entropy) and are a
> screening step, not a test of H1–H4. Code: MIT (`LICENSE`).

## Hypotheses (Julian, 2026-10-27: mechanism applies to the base-driven ring closure — tentative)
- **H1 (relay deprotonation):** In base-induced cyclisation of poly-β-carbonyl chains, the active-methylene C–H
  proton leaves via a multi-hop Grotthuss relay (protic solvent / base conjugate acid / substrate enol–carbonyl
  groups) rather than by direct transfer to a discrete base.
- **H2 (relay reprotonation + coupling):** The alkoxide formed on ring closure is reprotonated through the same
  wire, and the shuttled proton couples successive ring closures into a cascade.
- **H3 (kinetic signature):** Primary H/D KIE with proton-inventory dependence on hop number; strong dependence
  on protic-solvent structure at constant base strength.
- **H4 (flow relevance):** Local base concentration, water/alcohol content and mixing (controllable in flow) govern
  rate and regioselectivity.
- Possible extension (low priority): analogous relay in acid-initiated epoxy-polyene cyclisation (control system).

Closest prior art (see `docs/literature.md`): ARO/CYC enzymes TcmN/WhiE use a single coordinated water as base
(Ames 2008; Lee 2012); computed water relays in enolisation/organocatalytic aldol (Alagona 2010; Armstrong 2014).

## Substrates (all three are simulated — Julian, 2026-10-27)
| ID | Substrate | Role | Folder |
|---|---|---|---|
| S1 | Platform compound **10** (F5-2) | target platform molecule; base cascade | `substrates/s1_compound10/` |
| S2 | Krohn polyketide **12** (→ 5,6-dideoxypradinone 11) | main base-cascade substrate | `substrates/s2_krohn12/` |
| S3 | Simple epoxy-polyene model | acid-cascade control / extension | `substrates/s3_epoxypolyene/` |
Structures, sources and ambiguities: `substrates/structures.yaml` (Julian, 2026-11-22: derive them from the literature).
S2 is taken from the experimental entry of polyketide 34 in Hayat's Paderborn thesis (Krohn group, 1999), which reports the
same synthesis as Krohn et al. JOC 2000. S1 is a stand-in: the F5-2 scheme for compound 10 could not be retrieved.
Validate with `python scripts/validate_structures.py`.

## Methods (staged)
| Stage | Goal | Method | Compute |
|---|---|---|---|
| 0 | Toolchain smoke test | GFN2-xTB rigid proton scan in H5O2+ (`scripts/smoke_proton_scan.py`) | laptop, seconds |
| 1 | Model relays | Acetylacetone C3–H → base (OH⁻, MeO⁻, NMe3) and reprotonation of an aldol alkoxide by BH, through (H2O)n / (MeOH)n wires, n=0–4, relay vs. direct topology at equal composition, 3 seeds; GFN2-xTB/ALPB (tblite) endpoint relaxation + CI-NEB (ASE) | CPU, hours |
| 2 | Substrate paths | For S1, S2 (base) and S3 (acid control): CREST conformers, then classical vs. relay paths for deprotonation, ring closure, reprotonation; DFT refinement (Psi4 or PySCF/gpu4pyscf, ωB97X-D/def2-TZVP//def2-SVP, implicit + explicit solvent); cut-down fragments first | CPU cluster |
| 3 | Dynamics | Explicit-solvent QM/MM (CP2K) and AIMD/metadynamics on proton-transfer CVs for S2 (then S1, S3) | GPU/HPC — **not before November 2026** |
| 4 | Kinetics | TST rates, Wigner/Bell tunnelling, H/D KIE and proton inventory; relay vs. classical microkinetics | CPU |
| 5 | Experimental proposal | Predictions (KIE, solvent dependence, D-labelling) for wet lab / flow | — |

Success criterion: relay path lower than classical at the same level of theory and robust to functional and
solvation choice; otherwise H1/H2 are rejected and reported. No threshold numbers are fixed until stage 1.

## Layout
```
scripts/smoke_proton_scan.py      stage 0
scripts/archive_to_bucket.sh      encrypted upload
substrates/structures.yaml        all SMILES/SMARTS with provenance
substrates/<id>/substrate.yaml    SMILES, charge, fragments
scripts/validate_structures.py    RDKit validation
stages/stage1_model_relays/       build.py, run_stage1.py, summarize.py
tests/                            pytest (no xtb needed)
stages/stage2_substrate_paths/{s1_compound10,s2_krohn12,s3_epoxypolyene}/
docs/literature.md, docs/open-questions.md
```

## Tools (OSI open-source only — Julian, 2026-11-22: no ORCA)
xtb and tblite (GFN2-xTB, LGPL-3.0), CREST (LGPL-3.0), ASE (LGPL-2.1), RDKit (BSD-3), Psi4 (LGPL-3.0),
PySCF and gpu4pyscf (Apache-2.0), CP2K (GPL-2.0+), PLUMED (LGPL-3.0).

## Environment
```bash
micromamba create -y -n proton-relay -c conda-forge --override-channels python=3.11 xtb crest rdkit ase numpy scipy pyyaml tblite-python pytest
micromamba activate proton-relay
```
Used on engine (2026-11-22): tblite 0.7.0, ASE 3.29.0, RDKit 2026.09.2, NumPy 2.4.6, xtb 6.4.1 (conda-forge resolved
6.4.1 once tblite-python was added), CREST 3.0.2.

## Validate structures and run the tests
```bash
python scripts/validate_structures.py
python -m pytest -q tests
```

## Run stage 1
```bash
cd stages/stage1_model_relays
OMP_NUM_THREADS=1 nice -n 19 ionice -c3 python run_stage1.py --out /tmp/stage1-test --workers 2 --limit 4   # quick check
OMP_NUM_THREADS=1 nohup nice -n 19 ionice -c3 python run_stage1.py --out RUN_DIR --workers 12 > stage1.log 2>&1 &
python summarize.py RUN_DIR          # writes data/stage1_jobs.csv and data/stage1_summary.json
```
Statuses: `ok` (CI-NEB converged), `neb_not_converged`, `reactant_not_minimum` (the C–H/B–H proton left its donor while
the reactant complex relaxed, i.e. no barrier at this level), `product_not_minimum` (the product complex relaxed back),
`timeout` (20 min cap per job), `error`. Finished jobs are skipped on restart.

## Run the smoke test
```bash
python3 -m venv venv && . venv/bin/activate && pip install numpy
# xtb binary: https://github.com/grimme-lab/xtb/releases (6.7.1 tested on 2026-10-27; 6.4.1 from conda-forge re-checked on 2026-11-22)
XTB=/path/to/xtb python3 scripts/smoke_proton_scan.py --roo 2.6   # single-well (shared proton)
XTB=/path/to/xtb python3 scripts/smoke_proton_scan.py --roo 2.9   # double-well
```
Smoke-test sanity check (2026-10-27, xtb 6.7.1): R(O–O)=2.6 Å gives a single minimum at the midpoint,
2.9 Å a symmetric double well with a ~4 kcal/mol *rigid-scan* barrier. That is the expected qualitative
Zundel behaviour; numbers are toolchain checks only and must not be quoted as results.

## Data management
- Raw outputs stay out of git (`data/` is git-ignored except small JSON/CSV summaries; `runs/` is ignored).
- Stage-1 raw runs live on engine under `~/proton-hopping-sim/runs/` and are included in the engine research backup.
- Archive per run: `tar` → encrypt client-side (age or gpg, key from Vault, never in repo) → upload to the
  object-storage bucket under `s3://$S3_BUCKET/research/proton-hopping/<stage>/<date>/` via
  `scripts/archive_to_bucket.sh`. Credentials come from Vault at run time.
- Each archive gets a SHA256 + metadata JSON (code commit, xtb/DFT version, level of theory).

## Not yet decided (questions for Julian)
See `docs/open-questions.md`.
