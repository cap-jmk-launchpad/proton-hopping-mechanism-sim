"""Stage 1 cluster builder: proton source -> H-bond wire of n solvent molecules -> proton sink.

Topologies (same atoms in both, so energies are comparable):
  relay : source-H ... O(wire 1)-H ... O(wire n)-H ... sink   (n+1 protons move)
  direct: source-H ... sink, the n solvent molecules donate H-bonds to spectator sites (1 proton moves)
Deprotonation (step 2): source = acetylacetone C3-H, sink = OH-, MeO- or NMe3.
Reprotonation (step 4): source = conjugate acid BH (H2O, MeOH, HNMe3+), sink = diacetone-alcohol alkoxide.
Geometries are crude starting points; they are relaxed with GFN2-xTB before any path search.
"""
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

R_HB = {"O": 1.80, "N": 1.90}       # H...acceptor distance used for placement (Angstrom)
R_CH_HB = 2.00                       # C-H...acceptor
CLASH = 1.55                         # minimum intermolecular distance accepted at build time

SMILES = {"acac": "CC(=O)CC(C)=O", "alkoxide": "CC(=O)CC(C)(C)[O-]",
          "water": "O", "methanol": "CO", "hydroxide": "[OH-]", "methoxide": "C[O-]",
          "nme3": "CN(C)C", "hnme3": "C[NH+](C)C"}
CHARGE = {"acac": 0, "alkoxide": -1, "water": 0, "methanol": 0, "hydroxide": -1,
          "methoxide": -1, "nme3": 0, "hnme3": 1}
# deprotonation sink -> conjugate acid used as reprotonation source
CONJ_ACID = {"hydroxide": "water", "methoxide": "methanol", "nme3": "hnme3"}


def monomer(name, seed=0):
    m = Chem.AddHs(Chem.MolFromSmiles(SMILES[name]))
    if m.GetNumAtoms() == 1:
        return ["O"] if name else [], np.zeros((1, 3)), m
    AllChem.EmbedMolecule(m, randomSeed=seed)
    try:
        AllChem.MMFFOptimizeMolecule(m)
    except Exception:
        pass
    xyz = m.GetConformer().GetPositions()
    return [a.GetSymbol() for a in m.GetAtoms()], np.array(xyz), m


def _unit(v):
    n = np.linalg.norm(v)
    return v / n if n > 1e-9 else v


def _rot(a, b):
    """Rotation matrix taking unit vector a onto unit vector b."""
    a, b = _unit(a), _unit(b)
    v, c = np.cross(a, b), float(np.dot(a, b))
    if np.linalg.norm(v) < 1e-9:
        if c > 0:
            return np.eye(3)
        p = _unit(np.cross(a, [1, 0, 0] if abs(a[0]) < 0.9 else [0, 1, 0]))
        return 2 * np.outer(p, p) - np.eye(3)
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + vx + vx @ vx * (1 / (1 + c))


def _axis_rot(axis, ang):
    axis = _unit(axis); K = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
    return np.eye(3) + np.sin(ang) * K + (1 - np.cos(ang)) * K @ K


def lone_pair_dirs(sym, xyz, mol, idx):
    """Approximate lone-pair directions of acceptor atom idx."""
    nb = [n.GetIdx() for n in mol.GetAtomWithIdx(idx).GetNeighbors()] if mol.GetNumAtoms() > 1 else []
    if not nb:
        return [np.array([1.0, 0, 0])]
    vs = [_unit(xyz[j] - xyz[idx]) for j in nb]
    back = _unit(-np.sum(vs, axis=0))
    a = mol.GetAtomWithIdx(idx)
    if a.GetSymbol() == "O" and len(nb) == 1 and a.GetFormalCharge() == 0 and \
            mol.GetBondBetweenAtoms(idx, nb[0]).GetBondTypeAsDouble() == 2:      # carbonyl O: sp2, two lone pairs
        c = nb[0]
        others = [n.GetIdx() for n in mol.GetAtomWithIdx(c).GetNeighbors() if n.GetIdx() != idx]
        normal = _unit(np.cross(xyz[others[0]] - xyz[c], xyz[others[1]] - xyz[c]))
        return [_axis_rot(normal, s * np.deg2rad(60)) @ back for s in (1, -1)]
    if len(nb) == 1:                                                           # sp3 O with one neighbour (alkoxide, OH-)
        perp = _unit(np.cross(back, [0, 0, 1] if abs(back[2]) < 0.9 else [0, 1, 0]))
        tilt = _axis_rot(perp, np.deg2rad(70.5)) @ back
        return [_axis_rot(back, k * 2 * np.pi / 3) @ tilt for k in range(3)]
    return [back]                                                              # water, alcohol, amine


class Cluster:
    def __init__(self):
        self.sym, self.xyz, self.charge, self.parts = [], np.zeros((0, 3)), 0, []

    def add(self, name, sym, xyz):
        off = len(self.sym)
        self.sym += sym; self.xyz = np.vstack([self.xyz, xyz]); self.charge += CHARGE[name]
        self.parts.append((name, off, len(sym)))
        return off

    def min_dist(self, xyz):
        if len(self.sym) == 0:
            return 9.9
        d = np.linalg.norm(self.xyz[:, None, :] - xyz[None, :, :], axis=-1)
        return float(d.min())


def _place(cl, name, seed, acc, h_pos, u, donor_h=None, fwd=np.array([1.0, 0, 0]), lp_choice=0, r=None):
    """Place monomer so acceptor atom `acc` accepts an H-bond from the H at h_pos (donor direction u)."""
    sym, xyz, mol = monomer(name, seed)
    lps = lone_pair_dirs(sym, xyz, mol, acc)
    lp = lps[lp_choice % len(lps)]
    R0 = _rot(lp, -u)
    base = (xyz - xyz[acc]) @ R0.T
    r = r if r is not None else R_HB.get(sym[acc], 1.85)
    best, best_s = None, -1e9
    rng = np.random.default_rng(seed)
    for ang in np.linspace(0, 2 * np.pi, 36, endpoint=False) + rng.uniform(0, 0.17):
        X = base @ _axis_rot(u, ang).T + h_pos + r * u
        md = cl.min_dist(np.delete(X, [], axis=0))
        s = (0 if md >= CLASH else -10 * (CLASH - md)) + 0.2 * md
        if donor_h is not None:
            s += float(np.dot(_unit(X[donor_h] - X[acc]), fwd))
        if s > best_s:
            best, best_s = X, s
    return sym, best


def _h_on(mol, heavy):
    return [n.GetIdx() for n in mol.GetAtomWithIdx(heavy).GetNeighbors() if n.GetSymbol() == "H"]


def build(step, sink, solvent, n, topology, seed=0):
    """Return dict(symbols, positions, charge, moving=[(h_index, from_atom, to_atom), ...])."""
    assert step in ("deprot", "reprot") and topology in ("relay", "direct")
    cl = Cluster()
    if step == "deprot":
        src_name, snk_name = "acac", sink
    else:
        src_name, snk_name = CONJ_ACID[sink], "alkoxide"
    sym, xyz, mol = monomer(src_name, seed)
    if src_name == "acac":
        donor_heavy = 3                                                       # central CH2 carbon
    else:
        donor_heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() in ("O", "N")][0]
    dh = _h_on(mol, donor_heavy)[0]
    u = _unit(xyz[dh] - xyz[donor_heavy])
    xyz = (xyz - xyz[donor_heavy]) @ _rot(u, [1, 0, 0]).T
    off = cl.add(src_name, sym, xyz)
    chain = [(off + donor_heavy, off + dh)]                                   # (heavy, H) donors along the wire
    wire_n = n if topology == "relay" else 0
    for k in range(wire_n):
        h_heavy, h = chain[-1]
        u = _unit(cl.xyz[h] - cl.xyz[h_heavy])
        s_sym, s_xyz, s_mol = monomer(solvent, seed + 11 * k + 1)
        o = [a.GetIdx() for a in s_mol.GetAtoms() if a.GetSymbol() == "O"][0]
        hd = _h_on(s_mol, o)[0]
        r = R_CH_HB if cl.sym[h_heavy] == "C" else None
        ps, X = _place(cl, solvent, seed + 11 * k + 1, o, cl.xyz[h], u, donor_h=hd, r=r)
        o2 = cl.add(solvent, ps, X)
        chain.append((o2 + o, o2 + hd))
    # sink
    h_heavy, h = chain[-1]
    u = _unit(cl.xyz[h] - cl.xyz[h_heavy])
    k_sym, k_xyz, k_mol = monomer(snk_name, seed + 97)
    acc = [a.GetIdx() for a in k_mol.GetAtoms() if a.GetSymbol() in ("O", "N") and
           (a.GetFormalCharge() < 0 or a.GetSymbol() == "N" or len(k_sym) <= 2)]
    acc = acc[0] if acc else [a.GetIdx() for a in k_mol.GetAtoms() if a.GetSymbol() == "O"][0]
    r = R_CH_HB if cl.sym[h_heavy] == "C" else None
    ps, X = _place(cl, snk_name, seed + 97, acc, cl.xyz[h], u, r=r)
    sk = cl.add(snk_name, ps, X)
    sink_atom = sk + acc
    # spectators for the direct topology: H-bond donors to carbonyl / alkoxide oxygens
    if topology == "direct" and n > 0:
        sites = []
        tmp_mol = Chem.AddHs(Chem.MolFromSmiles(SMILES[src_name if src_name == "acac" else "alkoxide"]))
        ref_name = "acac" if src_name == "acac" else "alkoxide"
        ref = [p for p in cl.parts if p[0] == ref_name][0]
        ref_xyz = cl.xyz[ref[1]:ref[1] + ref[2]]
        for a in tmp_mol.GetAtoms():
            if a.GetSymbol() == "O":
                for d in lone_pair_dirs(None, ref_xyz, tmp_mol, a.GetIdx()):
                    sites.append((ref[1] + a.GetIdx(), d))
        # fallback sites: the same lone pairs tilted by +-30 degrees (used only if the ideal ones clash)
        extra = []
        for site, d in sites:
            p = _unit(np.cross(d, [0, 0, 1] if abs(d[2]) < 0.9 else [0, 1, 0]))
            extra += [(site, _axis_rot(p, t) @ d) for t in (np.deg2rad(30), -np.deg2rad(30))]
        sites += extra
        used = set()
        placed = 0
        for site, d in sites:
            if placed == n:
                break
            if (site, round(float(d[0]), 1)) in used:
                continue
            s_sym, s_xyz, s_mol = monomer(solvent, seed + 31 * placed + 5)
            o = [a.GetIdx() for a in s_mol.GetAtoms() if a.GetSymbol() == "O"][0]
            hd = _h_on(s_mol, o)[0]
            v = _unit(s_xyz[hd] - s_xyz[o])
            Xs = (s_xyz - s_xyz[hd]) @ _rot(v, -d).T + cl.xyz[site] + 1.85 * d
            best, bs = Xs, -1e9
            for ang in np.linspace(0, 2 * np.pi, 24, endpoint=False):
                Y = (Xs - (cl.xyz[site] + 1.85 * d)) @ _axis_rot(d, ang).T + cl.xyz[site] + 1.85 * d
                md = cl.min_dist(Y)
                if md > bs:
                    best, bs = Y, md
            if bs < 1.4:
                continue
            cl.add(solvent, s_sym, best); placed += 1; used.add((site, round(float(d[0]), 1)))
        if placed < n:   # second shell: donate to the O of an already placed spectator
            for name, o0, ln in [p for p in cl.parts if p[0] == solvent]:
                if placed == n:
                    break
                o_idx = o0 + [i for i in range(ln) if cl.sym[o0 + i] == "O"][0]
                hs = [o0 + i for i in range(ln) if cl.sym[o0 + i] in ("H", "C") and
                      np.linalg.norm(cl.xyz[o0 + i] - cl.xyz[o_idx]) < 1.5]
                d = _unit(-np.sum([_unit(cl.xyz[j] - cl.xyz[o_idx]) for j in hs], axis=0))
                s_sym, s_xyz, s_mol = monomer(solvent, seed + 53 * placed + 7)
                o = [a.GetIdx() for a in s_mol.GetAtoms() if a.GetSymbol() == "O"][0]
                hd = _h_on(s_mol, o)[0]
                v = _unit(s_xyz[hd] - s_xyz[o]); c0 = cl.xyz[o_idx] + 1.85 * d
                Xs = (s_xyz - s_xyz[hd]) @ _rot(v, -d).T + c0
                cands = [(Xs - c0) @ _axis_rot(d, ang).T + c0 for ang in np.linspace(0, 2 * np.pi, 24, endpoint=False)]
                best = max(cands, key=cl.min_dist)
                if cl.min_dist(best) >= 1.4:
                    cl.add(solvent, s_sym, best); placed += 1
        if placed < n:
            raise RuntimeError(f"only {placed}/{n} spectator sites usable")
    # protons that move: each donor H goes to the next acceptor (next wire O or the sink)
    acceptors = [c[0] for c in chain[1:]] + [sink_atom]
    moving = [(h, d, a) for (d, h), a in zip(chain, acceptors)]
    return {"symbols": cl.sym, "positions": cl.xyz, "charge": cl.charge, "moving": moving,
            "parts": cl.parts, "source": src_name, "sink": snk_name}


def product_guess(sys_, r_new=1.00):
    """Move every relaying proton onto its acceptor (along the donor-acceptor line)."""
    P = sys_["positions"].copy()
    for h, d, a in sys_["moving"]:
        P[h] = P[a] + r_new * _unit(P[d] - P[a])
    return P
