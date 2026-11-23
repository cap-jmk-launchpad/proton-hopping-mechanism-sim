"""Generate crisp 2D chemical structure figure (Figure 7) using RDKit + Matplotlib."""
from __future__ import annotations

import io
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D

OUT_PDF = Path(__file__).resolve().parents[1] / "paper" / "fig7_chemical_structures.pdf"

SPECIES = [
    # Row 1: Carbon acids & preorganised / enzyme substrates
    ("Acetylacetone", r"$\mathrm{C_5H_8O_2}$ · carbon acid", "CC(=O)CC(=O)C"),
    ("Cyclohexa-2,4-dien-1-one", r"$\mathrm{C_6H_6O}$ · tautomerisation", "O=C1CC=CC=C1"),
    (r"$\it{all\text{-}cis}$-Cyclohexane-1,3,5-triol", r"$\mathrm{C_6H_{12}O_3}$ · preorganised wire", "OC1CC(O)CC(O)C1"),
    (r"$\it{trans}$-Dihydroquercetin (DQH)", r"$\mathrm{C_{15}H_{12}O_7}$ · TcmN ligand", "O=C1[C@@H](O)[C@H](c2ccc(O)c(O)c2)Oc2cc(O)cc(O)c21"),
    # Row 2: Bifunctional shuttles, positive control & amidine/guanidine bases
    ("2-Pyridone", r"$\mathrm{C_5H_5NO}$ · bifunctional", "O=c1cccc[nH]1"),
    ("7-Hydroxyquinoline", r"$\mathrm{C_9H_7NO}$ · positive control", "Oc1ccc2cccnc2c1"),
    ("TBD", r"$\mathrm{C_7H_{13}N_3}$ · bicyclic guanidine", "C1CCN2CCCN=C2N1"),
    ("DBU", r"$\mathrm{C_9H_{16}N_2}$ · bicyclic amidine", "C1CCC2=NCCCN2CC1"),
    # Row 3: Sterically hindered bases & Lewis-acid diketone chelates
    ("TMP", r"$\mathrm{C_9H_{19}N}$ · hindered amine", "CC1(C)CCCC(C)(C)N1"),
    ("DTBMP", r"$\mathrm{C_{14}H_{23}N}$ · bulky pyridine", "Cc1cc(C(C)(C)C)nc(C(C)(C)C)c1"),
    (r"$[\mathrm{MgCl_2(acac)}]$", r"Lewis-acid diketone chelate", "CC1=CC(C)=[O]->[Mg](Cl)(Cl)O1"),
    (r"$[\mathrm{ZnCl_2(acac)}]$", r"Lewis-acid diketone chelate", "CC1=CC(C)=[O]->[Zn](Cl)(Cl)O1"),
]


def render_mol(smiles: str, w: int = 640, h: int = 420) -> Image.Image:
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None, f"Invalid SMILES: {smiles}"
    for atom in mol.GetAtoms():
        if atom.GetSymbol() in ("Mg", "Zn"):
            atom.SetNumRadicalElectrons(0)
    AllChem.Compute2DCoords(mol)
    if "preorganised wire" in smiles or smiles == "OC1CC(O)CC(O)C1":
        for b in mol.GetBonds():
            a1, a2 = b.GetBeginAtom().GetSymbol(), b.GetEndAtom().GetSymbol()
            if {a1, a2} == {"C", "O"}:
                if b.GetBeginAtom().GetSymbol() == "O":
                    b.SetBeginAtomIdx if False else None
                b.SetBondDir(Chem.BondDir.BEGINWEDGE)
    drawer = rdMolDraw2D.MolDraw2DCairo(w, h)
    opts = drawer.drawOptions()
    opts.bondLineWidth = 2.4
    opts.padding = 0.14
    opts.minFontSize = 20
    opts.maxFontSize = 32
    opts.useBWAtomPalette()
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    return Image.open(io.BytesIO(drawer.GetDrawingText()))


def main() -> None:
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.size": 8.5,
        "axes.edgecolor": "#D0D7DE",
        "axes.linewidth": 0.8,
    })
    fig, axes = plt.subplots(3, 4, figsize=(10.5, 7.2), dpi=300)
    for ax, (title, subtitle, smi) in zip(axes.flat, SPECIES):
        img = render_mol(smi)
        ax.imshow(img)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(f"{title}\n{subtitle}", fontsize=8.2, pad=4)
    fig.tight_layout(h_pad=1.1, w_pad=0.8)
    fig.savefig(OUT_PDF, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {OUT_PDF}")


if __name__ == "__main__":
    main()
