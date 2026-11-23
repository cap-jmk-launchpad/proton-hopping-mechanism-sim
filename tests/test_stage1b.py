import sys
from pathlib import Path
import numpy as np
import pytest

D = Path(__file__).resolve().parents[1] / "stages/stage1b_relay_design"
sys.path.insert(0, str(D))
import build1b  # noqa: E402
import kinetics  # noqa: E402


def _ok(s):
    X = s["positions"]
    d = np.linalg.norm(X[:, None] - X[None], axis=-1)[np.triu_indices(len(X), 1)]
    assert d.min() > 0.85
    for h, dn, a in s["moving"]:
        assert np.linalg.norm(X[h] - X[dn]) < 1.2          # every moving proton starts on its donor


@pytest.mark.parametrize("base,q", [("nme3", 0), ("nmi", 0), ("acetate", -1)])
@pytest.mark.parametrize("n", [0, 1, 2, 3])
def test_tethered_natural(base, q, n):
    s = build1b.build_tethered(base, "water", n)
    _ok(s); assert s["charge"] == q and len(s["moving"]) == n + 1 and s["constraint"] is None


@pytest.mark.parametrize("n,d", [(0, 3.0), (0, 6.0), (1, 5.0), (2, 6.0), (3, 7.0)])
def test_tethered_constrained_distance(n, d):
    s = build1b.build_tethered("nme3", "water", n, d=d)
    c, sink = s["constraint"]
    assert abs(np.linalg.norm(s["positions"][sink] - s["positions"][c]) - d) < 1e-6
    _ok(s)


def test_compressed_wire_is_infeasible():
    with pytest.raises(build1b.InfeasibleGeometry):
        build1b.build_tethered("nme3", "water", 3, d=3.0)


@pytest.mark.parametrize("sub", ["acac", "dienone"])
@pytest.mark.parametrize("sh", [[], ["water"], ["water", "water"], ["acoh"], ["pyridone"], ["hfip"]])
def test_shuttle(sub, sh):
    r = build1b.build_shuttle(sub, sh, "relay"); _ok(r)
    assert len(r["moving"]) == len(sh) + 1
    if sh:
        st = build1b.build_shuttle(sub, sh, "stepwise")
        assert len(st["moving"]) == 1 and sorted(st["symbols"]) == sorted(r["symbols"])


@pytest.mark.parametrize("base", ["nme3", "acetate"])
def test_triol_same_composition(base):
    r, d = build1b.build_triol(base, "relay"), build1b.build_triol(base, "direct")
    _ok(r); _ok(d)
    assert len(r["moving"]) == 4 and len(d["moving"]) == 1 and sorted(r["symbols"]) == sorted(d["symbols"])


def test_shell_adds_molecules():
    a, b = build1b.build_tethered("nme3", "water", 2), build1b.build_tethered("nme3", "water", 2, shell=3)
    assert len(b["symbols"]) == len(a["symbols"]) + 9


def test_kinetics_limits():
    assert abs(kinetics.eckart(150, 15, 15) - kinetics.wigner(150)) < 0.01     # shallow-curvature limit
    w = kinetics.partial_hessian_freqs(np.diag([5.0, 5, 5]), [kinetics.M_H])
    assert np.allclose(w, w[0]) and 1100 < w[0] < 1200
    wd = kinetics.partial_hessian_freqs(np.diag([5.0, 5, 5]), [kinetics.M_D])
    assert abs(w[0] / wd[0] - np.sqrt(kinetics.M_D / kinetics.M_H)) < 1e-9


def test_job_grid():
    run1b = pytest.importorskip("run1b")
    J = run1b.jobs()
    ids = [run1b.job_id(j) for j in J]
    assert len(ids) == len(set(ids)) == 660
    assert run1b.parse_job_id("C_triol_nme3_water_relay_sh2_s0")["shell"] == 2


def test_segmented_path_runs_with_mock_calculator():
    """Plumbing test of the sequential-hop fallback with a Morse calculator (no xtb): returns a path and a barrier."""
    run1b = pytest.importorskip("run1b")
    pytest.importorskip("ase")                      # CI installs rdkit/numpy/pyyaml/pytest only
    from types import SimpleNamespace
    from ase import Atoms
    from ase.calculators.morse import MorsePotential
    from ase.optimize import BFGS

    def calc(q, solvent, etemp=300.0):
        return MorsePotential(epsilon=1.0, r0=1.0, rho0=6.0)

    def relax(a, q, solvent, fmax=0.05, steps=200):
        a.calc = calc(q, solvent); opt = BFGS(a, logfile=None); conv = opt.run(fmax=fmax, steps=steps)
        return float(a.get_potential_energy()), bool(conv), int(opt.nsteps)

    def classify(pos, moving):
        d = lambda i, k: float(np.linalg.norm(pos[i] - pos[k]))
        return {"all_on_acceptor": all(d(h, a) < d(h, dn) for h, dn, a in moving)}

    s1 = SimpleNamespace(calc=calc, relax=relax, classify=classify, NEB_ATTEMPTS=(("linear", 300.0),))
    X = np.array([[0, 0, 0], [1.0, 0, 0], [2.6, 0.3, 0], [3.6, 0.3, 0], [5.2, 0, 0]])
    s = {"symbols": ["H"] * 5, "positions": X, "moving": [(1, 0, 2), (3, 2, 4)]}
    R = Atoms("H5", positions=X); relax(R, 0, None)
    out = run1b.run_segmented(R, s, 0, None, s1)
    assert "seg_barrier_kcal" in out and len(out["intermediates"]) == 2 and len(out["segments"]) >= 1
