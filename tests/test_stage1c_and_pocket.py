"""Unit tests for Stage 1c (non-peptide bulky bases, bicyclic guanidine TBD, and [MCl2(acac)] metal complexes)
and Stage 2 (3TVQ TcmN ARO/CYC 241-atom enzyme-pocket QM cluster)."""
import json
import sys
from pathlib import Path
import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "stages/stage1c_nonpeptide"))
import run_stage1c as s1c  # noqa: E402


def _check_cluster(spec):
    X = spec["positions"]
    assert len(X) == len(spec["symbols"])
    d = np.linalg.norm(X[:, None] - X[None], axis=-1)[np.triu_indices(len(X), 1)]
    assert d.min() > 0.85, f"atomic clash: min dist {d.min():.3f} A"
    for h, dn, acc in spec["moving"]:
        assert spec["symbols"][h] == "H"
        assert np.linalg.norm(X[h] - X[dn]) < 1.25, "moving proton must start bonded to its donor"
        assert np.linalg.norm(X[h] - X[acc]) > 1.30, "moving proton must not start bonded to its acceptor"


@pytest.mark.parametrize("base", ["dtbmp", "tmp", "dbu"])
@pytest.mark.parametrize("n,top", [(0, "direct"), (1, "direct"), (1, "relay"), (2, "relay")])
def test_stage1c_bulky_base_clusters(base, n, top):
    spec = s1c.build_bulky_base_cluster(base, n, top, seed=0)
    _check_cluster(spec)
    assert spec["charge"] == 0
    expected_hops = (n + 1) if top == "relay" else 1
    assert len(spec["moving"]) == expected_hops


@pytest.mark.parametrize("base", ["dtbmp", "tmp", "dbu"])
def test_stage1c_bulky_base_n1_composition_parity(base):
    d = s1c.build_bulky_base_cluster(base, 1, "direct", seed=0)
    r = s1c.build_bulky_base_cluster(base, 1, "relay", seed=0)
    assert sorted(d["symbols"]) == sorted(r["symbols"])
    assert len(d["moving"]) == 1 and len(r["moving"]) == 2


@pytest.mark.parametrize("sub", ["acac", "dienone"])
def test_stage1c_tbd_shuttle_parity(sub):
    r = s1c.build_tbd_shuttle(sub, "relay", seed=0)
    st = s1c.build_tbd_shuttle(sub, "stepwise", seed=0)
    _check_cluster(r)
    _check_cluster(st)
    assert sorted(r["symbols"]) == sorted(st["symbols"])
    assert len(r["moving"]) == 2 and len(st["moving"]) == 1


@pytest.mark.parametrize("metal", ["Mg", "Zn"])
@pytest.mark.parametrize("base", ["nmi", "dtbmp"])
def test_stage1c_metal_complex_parity_and_chelation(metal, base):
    d = s1c.build_metal_complex(metal, base, "direct", seed=0)
    r = s1c.build_metal_complex(metal, base, "relay", seed=0)
    _check_cluster(d)
    _check_cluster(r)
    assert sorted(d["symbols"]) == sorted(r["symbols"])
    assert d["symbols"].count(metal) == 1 and d["symbols"].count("Cl") == 2
    assert len(d["moving"]) == 1 and len(r["moving"]) == 2


def test_stage1c_job_list_has_72_unique_jobs():
    jobs = s1c.all_jobs()
    names = [name for name, _ in jobs]
    assert len(names) == len(set(names)) == 72


def test_stage2_3tvq_pocket_cluster_integrity():
    meta = json.load(open(ROOT / "stages/stage2_pocket_3tvq/cluster_3TVQ.json"))
    lines = (ROOT / "stages/stage2_pocket_3tvq/cluster_3TVQ.xyz").read_text().splitlines()
    natoms = int(lines[0].strip())
    assert natoms == meta["natoms"] == 241
    assert meta["charge"] == 0 and meta["uhf"] == 0 and meta["nelectrons"] == 914
    assert len(meta["frozen_1based"]) == 24
    assert len(meta["labels"]) == 241
    # Verify key active-site residues and water wire labels exist
    for required in ("ASP57:OD2", "HIS128:NE2", "HIS124:ND1", "ARG82:NH1", "HOH259:O", "HOH266:O", "DQH:C10", "DQH:H24", "DQH:O27", "DQH:H31"):
        assert required in meta["labels"]
