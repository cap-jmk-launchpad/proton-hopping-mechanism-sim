# proton-hopping-mechanism-sim

Simulation research project for Julian M. Kleber's proposed **proton-hopping (Grotthuss-type) mechanism**
for the acid-/base-mediated cyclisation cascades behind tetracycline-type and terpenoid ring systems
(origin: *F5-2 API Synthesis Plan*, 2020, DOI 10.13140/RG.2.2.23579.28965; rewritten preprint draft in
`cap-jmk-launchpad/proton-hopping-preprints`).

> Status: research plan + toolchain smoke test. **No mechanistic result exists yet.** Everything below
> marked H1..H4 is a hypothesis.

## Hypotheses
- **H1 (relay initiation):** In strongly acidic, H-bonded media (superacid / Lewis-acid-assisted Brønsted acid
  such as [TiCl4OH]⁻/H⁺, or an SHC-like acid site), the initiating proton reaches the alkene/epoxide via a
  Grotthuss relay through an H-bond wire rather than by direct encounter with a free acid.
- **H2 (concerted shuttle in the cascade):** Proton release/re-uptake that terminates or propagates the
  cationic (terpenoid) or anionic/enolate (polyketide, "superbase") cyclisation is shuttled through the same
  wire, lowering the effective barrier vs. the classical stepwise mechanism.
- **H3 (kinetic signature):** If H1/H2 hold, the rate-limiting step shows a primary H/D KIE and a
  pronounced dependence on wire length / solvent H-bond structure; the classical mechanism does not.
- **H4 (flow relevance):** The relay is sensitive to local acid concentration and mixing, i.e. conditions a
  microreactor can control (link to the continuous-flow argument of the original document).

## Methods (staged)
| Stage | Goal | Method | Compute |
|---|---|---|---|
| 0 | Toolchain smoke test | GFN2-xTB rigid proton scan in H5O2+ (`scripts/smoke_proton_scan.py`) | laptop, seconds |
| 1 | Model systems | Proton transfer along (H2O)n / (HF)n / ROH wires, n=1–4; relaxed constrained scans + NEB (xtb via ASE) | CPU, hours |
| 2 | Cascade models | Truncated substrates: epoxy-polyene (Corey/Nishizawa type) and a poly-β-ketone (Harris–Wittek type); classical stepwise vs. relay-assisted path; DFT refinement (ORCA or Psi4, ωB97X-D/def2-TZVP //def2-SVP, CPCM/SMD + explicit waters) | CPU cluster |
| 3 | Dynamics | Explicit-solvent QM/MM (CP2K QM/MM or ORCA/xtb in ASE) and AIMD/metadynamics on the proton coordinate (CV: proton-transfer coordinate δ, coordination numbers) | GPU/HPC — **not before November 2026** |
| 4 | Kinetics | Eyring/TST rates, Bell/Wigner tunnelling corrections, H/D KIE from ZPE; microkinetic comparison relay vs. classical; sensitivity to wire length and acid strength | CPU |
| 5 | Experimental proposal | Predictions (KIE, solvent/acid dependence) handed to wet lab; flow conditions | — |

Success criteria: relay path ≥ ~2–3 kcal/mol below classical at the same level of theory and robust to
functional/solvation choice; otherwise H2 is rejected and reported as such.

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
