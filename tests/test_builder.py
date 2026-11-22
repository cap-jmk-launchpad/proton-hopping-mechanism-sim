import sys
from pathlib import Path
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "stages/stage1_model_relays"))
import build  # noqa: E402


@pytest.mark.parametrize("step", ["deprot", "reprot"])
@pytest.mark.parametrize("sink", ["hydroxide", "methoxide", "nme3"])
@pytest.mark.parametrize("n", [0, 2, 4])
@pytest.mark.parametrize("top", ["relay", "direct"])
def test_build_counts(step, sink, n, top):
    s = build.build(step, sink, "water", n, top, seed=0)
    assert len(s["moving"]) == (n + 1 if top == "relay" else 1)
    X = s["positions"]
    d = np.linalg.norm(X[:, None] - X[None], axis=-1)[np.triu_indices(len(X), 1)]
    assert d.min() > 0.85                       # no overlapping atoms
    exp_q = {"deprot": {"hydroxide": -1, "methoxide": -1, "nme3": 0},
             "reprot": {"hydroxide": -1, "methoxide": -1, "nme3": 0}}[step][sink]
    assert s["charge"] == exp_q


def test_same_composition_relay_direct():
    a = build.build("deprot", "hydroxide", "methanol", 3, "relay", seed=1)
    b = build.build("deprot", "hydroxide", "methanol", 3, "direct", seed=1)
    assert sorted(a["symbols"]) == sorted(b["symbols"])


def test_product_guess_moves_protons():
    s = build.build("deprot", "nme3", "water", 2, "relay", seed=0)
    P = build.product_guess(s)
    for h, d, a in s["moving"]:
        assert abs(np.linalg.norm(P[h] - P[a]) - 1.0) < 1e-6
