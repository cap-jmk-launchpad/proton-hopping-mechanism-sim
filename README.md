# proton-hopping-mechanism-sim

Simulation research project for Julian M. Kleber's proposed **proton-hopping (Grotthuss-type) mechanism**
for the base-induced polyketide cyclisation cascades behind tetracycline-type ring systems
(origin: *F5-2 API Synthesis Plan*, 2020, DOI 10.13140/RG.2.2.23579.28965; rewritten preprint draft in
`cap-jmk-launchpad/proton-hopping-preprints`).

> Status: research plan + toolchain smoke test. **No mechanistic result exists yet.** H1–H4 are hypotheses.

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
Structures must be entered by Julian (SMILES in each `substrate.yaml`); cut-down fragments are derived per stage.

## Methods (staged)
| Stage | Goal | Method | Compute |
|---|---|---|---|
| 0 | Toolchain smoke test | GFN2-xTB rigid proton scan in H5O2+ (`scripts/smoke_proton_scan.py`) | laptop, seconds |
| 1 | Model relays | Proton transfer from a β-diketone active methylene through (H2O)n / (MeOH)n, n=0–4, to a base (OH⁻, MeO⁻, amine); relaxed scans + NEB (xtb/ASE) | CPU, hours |
| 2 | Substrate paths | For S1, S2 (base) and S3 (acid control): CREST conformers, then classical vs. relay paths for deprotonation, ring closure, reprotonation; DFT refinement (Psi4/ORCA, ωB97X-D/def2-TZVP//def2-SVP, SMD + explicit solvent); cut-down fragments first | CPU cluster |
| 3 | Dynamics | Explicit-solvent QM/MM (CP2K) and AIMD/metadynamics on proton-transfer CVs for S2 (then S1, S3) | GPU/HPC — **not before November 2026** |
| 4 | Kinetics | TST rates, Wigner/Bell tunnelling, H/D KIE and proton inventory; relay vs. classical microkinetics | CPU |
| 5 | Experimental proposal | Predictions (KIE, solvent dependence, D-labelling) for wet lab / flow | — |

Success criterion: relay path lower than classical at the same level of theory and robust to functional and
solvation choice; otherwise H1/H2 are rejected and reported. No threshold numbers are fixed until stage 1.

## Layout
```
scripts/smoke_proton_scan.py      stage 0
scripts/archive_to_bucket.sh      encrypted upload
substrates/<id>/substrate.yaml    SMILES (TODO Julian), charge, fragments
stages/stage1_model_relays/       (to be filled)
stages/stage2_substrate_paths/{s1_compound10,s2_krohn12,s3_epoxypolyene}/
docs/literature.md, docs/open-questions.md
```

## Tools (open source / free-for-academia)
xtb (GFN2, LGPL), ASE (NEB, constraints), CREST (conformers), Psi4 (LGPL), CP2K (GPL; AIMD, QM/MM),
ORCA (free academic licence – check licence terms for commercial use by sail.black), RDKit, PLUMED.

## Run the smoke test
```bash
python3 -m venv venv && . venv/bin/activate && pip install numpy
# xtb binary: https://github.com/grimme-lab/xtb/releases (6.7.1 tested)
XTB=/path/to/xtb python3 scripts/smoke_proton_scan.py --roo 2.6   # single-well (shared proton)
XTB=/path/to/xtb python3 scripts/smoke_proton_scan.py --roo 2.9   # double-well
```
Smoke-test sanity check (2026-10-27, xtb 6.7.1): R(O–O)=2.6 Å gives a single minimum at the midpoint,
2.9 Å a symmetric double well with a ~4 kcal/mol *rigid-scan* barrier. That is the expected qualitative
Zundel behaviour; numbers are toolchain checks only and must not be quoted as results.

## Data management
- Raw outputs stay out of git (`data/` is git-ignored except small JSON summaries).
- Archive per run: `tar` → encrypt client-side (age or gpg, key from Vault, never in repo) → upload to the
  Hetzner bucket under `s3://jkleberdokstorebackup/research/proton-hopping/<stage>/<date>/` via
  `scripts/archive_to_bucket.sh`. Credentials come from Vault at run time.
- Each archive gets a SHA256 + metadata JSON (code commit, xtb/DFT version, level of theory).

## Not yet decided (questions for Julian)
See `docs/open-questions.md`.
