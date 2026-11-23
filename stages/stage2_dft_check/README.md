# Stage 2: DFT check of the stage-1 candidates

Open-source only: PySCF 2.14 (Apache-2.0), geomeTRIC 1.1 (BSD-3), pyscf-dispersion / simple-dftd3 (LGPL-3).

Level: wB97X-D3(BJ)/def2-SVP, IEF-PCM (water eps 78.3553, methanol eps 32.613), RI-JK; def2-TZVP single points
on the SVP geometries. Per stage-1 condition the lowest-barrier converged seed is taken (`select_candidates.py`).
For each one, `run_dft.py` re-optimises the reactant (CI-NEB image 0), the product (last image) and the TS (climbing
image, geomeTRIC TS search starting from the analytic DFT Hessian), then computes analytic Hessians at all three
(frequencies, ZPE, RRHO thermal corrections at 298.15 K and 1 atm). Every step writes its own JSON and xyz, so a
restart skips finished steps.

Validation: a TS counts only with exactly one imaginary frequency below -50 cm-1, a minimum only with none.
Rejection rule (fixed in the plan before any run): the relay is rejected if dG‡(relay) is not lower than dG‡(direct)
at the same level of theory, or if the sign flips between levels.

```bash
# environment (on top of environment.yml)
pip install "pyscf==2.14.0" geometric pyscf-dispersion
# candidate list (priority 1 = NMe3 water/methanol n=0-4; 2 = OH-/MeO- relay chains with a barrier)
python stages/stage2_dft_check/select_candidates.py runs/stage1-20261122 --priority 1 > queue.txt
# one job
OMP_NUM_THREADS=16 python stages/stage2_dft_check/run_dft.py runs/stage1-20261122/<job> --out runs/dft-20261122/x --tz-basis def2-tzvp
# a queue at low priority
OMP_NUM_THREADS=10 stages/stage2_dft_check/run_queue.sh queue.txt runs/dft-20261122/nme3 --tz-basis def2-tzvp
# table + rejection rule
python stages/stage2_dft_check/summarize_dft.py runs/dft-20261122/pilot runs/dft-20261122/nme3 --csv data/stage2_dft.csv
```
