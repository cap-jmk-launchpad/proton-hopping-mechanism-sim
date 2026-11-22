import subprocess, sys
from pathlib import Path
import yaml
from rdkit import Chem

ROOT = Path(__file__).resolve().parents[1]
SPEC = yaml.safe_load(open(ROOT / "substrates/structures.yaml"))


def test_validator_passes(tmp_path):
    r = subprocess.run([sys.executable, str(ROOT / "scripts/validate_structures.py"), "--out", str(tmp_path / "v.json")],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr


def _m(k):
    e = {**SPEC["substrates"], **SPEC["fragments"]}[k]
    return Chem.MolFromSmiles(e["smiles"])


def test_reaction_centre_smarts():
    p = {k: Chem.MolFromSmarts(v) for k, v in SPEC["smarts"].items() if ">>" not in v}
    assert len(_m("F1_acac").GetSubstructMatches(p["active_methylene_13dicarbonyl"])) == 1
    assert len(_m("S2_krohn12").GetSubstructMatches(p["benzylic_ketone_methylene"])) == 1
    assert len(_m("S3_epoxypolyene").GetSubstructMatches(p["epoxide"])) == 1
    assert len(_m("F1_alkoxide").GetSubstructMatches(p["alkoxide"])) == 1
    assert len(_m("F1_enolate").GetSubstructMatches(p["enolate_carbanion"])) == 1


def test_substrate_files_match_master():
    for d, k in (("s1_compound10", "S1_compound10"), ("s2_krohn12", "S2_krohn12"), ("s3_epoxypolyene", "S3_epoxypolyene")):
        y = yaml.safe_load(open(ROOT / "substrates" / d / "substrate.yaml"))
        a, b = Chem.MolFromSmiles(y["smiles"]), Chem.MolFromSmiles(SPEC["substrates"][k]["smiles"])
        assert Chem.MolToSmiles(a) == Chem.MolToSmiles(b)
