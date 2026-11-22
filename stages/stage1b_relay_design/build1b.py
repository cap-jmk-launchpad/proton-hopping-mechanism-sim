"""Stage-1b cluster builder for the relay-design screen (families A-D, see docs/relay-system-design.md).

Family A  tethered base : carbon acid C-H -> (solvent)n wire -> base, with the C...base distance held at d (Angstrom).
Family B  shuttle tautomerisation : substrate alpha C-H -> shuttle acceptor; shuttle donor H -> substrate C=O
          (concerted cyclic relay) vs "stepwise" (only the C-H proton moves to the shuttle).
Family C  preorganised wire : acac C-H -> O1-H -> O2-H -> O3-H of all-cis cyclohexane-1,3,5-triol -> base.
Family D  second shell : k extra solvent molecules donating H-bonds to free lone pairs of the inner cluster.
Same conventions as stage 1: build -> GFN2-xTB relaxation of R and P -> CI-NEB. Geometries are crude starts.
"""
import sys
from pathlib import Path
import numpy as np
import yaml
from rdkit import Chem
from rdkit.Chem import AllChem

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "stage1_model_relays"))
from build import _unit, _rot, _axis_rot, lone_pair_dirs, product_guess  # noqa: E402,F401

HERE = Path(__file__).resolve().parent
SPECIES = yaml.safe_load((HERE / "species.yaml").read_text())
ALL = {k: v for g in SPECIES.values() for k, v in g.items()}
R_HB = {"O": 1.80, "N": 1.90}
R_CH_HB = 2.00
CLASH = 1.55


def monomer(name, seed=0, conf_select=None):
    """Embed species `name`; returns symbols, xyz, mol (with H) and {map number: atom index}."""
    m = Chem.AddHs(Chem.MolFromSmiles(ALL[name]["smiles"]))
    maps = {a.GetAtomMapNum(): a.GetIdx() for a in m.GetAtoms() if a.GetAtomMapNum()}
    if m.GetNumAtoms() > 1:
        if conf_select is None:
            AllChem.EmbedMolecule(m, randomSeed=seed)
            try:
                AllChem.MMFFOptimizeMolecule(m)
            except Exception:
                pass
        else:                                            # pick the best of several conformers by a score
            cids = AllChem.EmbedMultipleConfs(m, numConfs=60, randomSeed=seed)
            AllChem.MMFFOptimizeMoleculeConfs(m)
            best = max(cids, key=lambda c: conf_select(m, m.GetConformer(c).GetPositions(), maps))
            keep = Chem.Conformer(m.GetConformer(best)); m.RemoveAllConformers(); m.AddConformer(keep, assignId=True)
    xyz = m.GetConformer().GetPositions() if m.GetNumAtoms() > 1 else np.zeros((1, 3))
    return [a.GetSymbol() for a in m.GetAtoms()], np.array(xyz), m, maps


def h_on(mol, heavy):
    return [n.GetIdx() for n in mol.GetAtomWithIdx(heavy).GetNeighbors() if n.GetSymbol() == "H"]


class Cluster:
    def __init__(self):
        self.sym, self.xyz, self.charge, self.parts = [], np.zeros((0, 3)), 0, []

    def add(self, name, sym, xyz):
        off = len(self.sym)
        self.sym += list(sym); self.xyz = np.vstack([self.xyz, xyz]); self.charge += ALL[name]["charge"]
        self.parts.append((name, off, len(sym)))
        return off

    def min_dist(self, X, skip=()):
        if not self.sym:
            return 9.9
        keep = [i for i in range(len(self.sym)) if i not in skip]
        d = np.linalg.norm(self.xyz[keep][:, None, :] - X[None, :, :], axis=-1)
        return float(d.min())


def place_acceptor(cl, name, seed, h_pos, u, r, donor_h=None, target=None, lp_all=True):
    """Place monomer so its acceptor (map 1) accepts the H at h_pos along u. If target is given, the monomer's donor H
    (map 2, or map 1 for ROH/H2O) is rotated to point at target (closing a cyclic relay)."""
    sym, xyz, mol, maps = monomer(name, seed)
    acc = maps[1]
    lps = lone_pair_dirs(sym, xyz, mol, acc)
    best, best_s, best_h = None, -1e9, None
    rng = np.random.default_rng(seed)
    dh_heavy = maps.get(2, acc)
    hs = h_on(mol, dh_heavy) if mol.GetNumAtoms() > 1 else []
    for lp in (lps if lp_all else lps[:1]):
        base = (xyz - xyz[acc]) @ _rot(lp, -u).T
        for ang in np.linspace(0, 2 * np.pi, 48, endpoint=False) + rng.uniform(0, 0.13):
            X = base @ _axis_rot(u, ang).T + h_pos + r * u
            md = cl.min_dist(X)
            s = (0 if md >= CLASH else -10 * (CLASH - md)) + 0.2 * md
            hbest = None
            if target is not None and hs:
                dt = [(float(np.linalg.norm(X[h] - target)), h) for h in hs]
                dmin, hbest = min(dt)
                s -= 2.0 * dmin
            elif hs:
                hbest = max(hs, key=lambda h: float(np.dot(_unit(X[h] - X[dh_heavy]), _unit(u))))
                s += float(np.dot(_unit(X[hbest] - X[dh_heavy]), _unit(u)))
            if s > best_s:
                best, best_s, best_h = X, s, hbest
    return sym, best, acc, dh_heavy, best_h


def _start(cl, name, seed):
    sym, xyz, mol, maps = monomer(name, seed)
    c = maps[1]; h = h_on(mol, c)[0]
    u = _unit(xyz[h] - xyz[c])
    xyz = (xyz - xyz[c]) @ _rot(u, [1, 0, 0]).T
    off = cl.add(name, sym, xyz)
    return off, off + c, off + h, {k: off + v for k, v in maps.items()}


def add_spectators(cl, solvent, k, sites, seed):
    """k solvent molecules donate an O-H...acceptor H-bond to lone pairs of the given cluster atoms
    (round-robin over the sites, at most two passes)."""
    placed = 0
    for site in list(sites) * 2:
        if placed == k:
            break
        name, off, ln = [p for p in cl.parts if p[1] <= site < p[1] + p[2]][0]
        sym, xyz, mol, _ = monomer(name, 0)
        cand = lone_pair_dirs(sym, cl.xyz[off:off + ln], mol, site - off) if len(sym) > 1 else [np.array([1.0, 0, 0])]
        nb = [q.GetIdx() + off for q in mol.GetAtomWithIdx(site - off).GetNeighbors()]
        if cl.sym[site] == "O" and len(nb) == 2:            # sp3 O with two neighbours: two tetrahedral lone pairs
            n1, n2 = (_unit(cl.xyz[q] - cl.xyz[site]) for q in nb)
            bk, pp = _unit(-(n1 + n2)), _unit(np.cross(n1, n2)); t = np.deg2rad(54.75)
            cand = [_unit(bk * np.cos(t) + pp * np.sin(t)), _unit(bk * np.cos(t) - pp * np.sin(t))]
        tilted = []
        for d0 in cand:                                   # fallback directions: lone pairs tilted by +-35 degrees
            p = _unit(np.cross(d0, [0, 0, 1] if abs(d0[2]) < 0.9 else [0, 1, 0]))
            tilted += [_axis_rot(p, t) @ d0 for t in (np.deg2rad(35), -np.deg2rad(35))]
            tilted += [_axis_rot(d0, np.pi / 2) @ _axis_rot(p, t) @ d0 for t in (np.deg2rad(35), -np.deg2rad(35))]
        cand = list(cand) + tilted
        for d in cand:
            if placed == k:
                break
            s_sym, s_xyz, s_mol, s_maps = monomer(solvent, seed + 31 * placed + 5)
            o = s_maps[1]; hd = h_on(s_mol, o)[0]
            v = _unit(s_xyz[hd] - s_xyz[o]); c0 = cl.xyz[site] + 1.85 * d
            Xs = (s_xyz - s_xyz[hd]) @ _rot(v, -d).T + c0
            cands = [(Xs - c0) @ _axis_rot(d, a).T + c0 for a in np.linspace(0, 2 * np.pi, 24, endpoint=False)]
            best = max(cands, key=cl.min_dist)
            if cl.min_dist(best) >= 1.45:
                cl.add(solvent, s_sym, best); placed += 1
                break                                     # one spectator per site per pass
    if placed < k:
        raise RuntimeError(f"only {placed}/{k} spectator sites usable")


def acceptor_sites(cl, exclude=()):
    out = []
    for i, s in enumerate(cl.sym):
        if s in ("O", "N") and i not in exclude:
            out.append(i)
    return out


def set_distance(cl, i, j, d, movable_from):
    """Translate every fragment with offset >= movable_from along the i->j axis so that |r_j - r_i| = d.
    Intermediate fragments move in proportion to the projection of their first heavy atom on the axis."""
    ax = cl.xyz[j] - cl.xyz[i]; d0 = float(np.linalg.norm(ax)); ax = ax / d0
    for name, off, ln in cl.parts:
        if off < movable_from:
            continue
        anchor = off + ([q for q in range(ln) if cl.sym[off + q] != "H"] or [0])[0]
        if off <= j < off + ln:
            f = 1.0
        else:
            f = float(np.clip(np.dot(cl.xyz[anchor] - cl.xyz[i], ax) / d0, 0.0, 1.0))
        cl.xyz[off:off + ln] += f * (d - d0) * ax
    return d0


class InfeasibleGeometry(RuntimeError):
    pass


def check_geometry(cl, min_d=0.90, min_inter=1.30):
    """Raise InfeasibleGeometry if atoms overlap (e.g. a wire compressed below its length by the d constraint)."""
    X = cl.xyz; D = np.linalg.norm(X[:, None] - X[None], axis=-1); np.fill_diagonal(D, 9.9)
    frag = np.zeros(len(X), int)
    for k, (_, off, ln) in enumerate(cl.parts):
        frag[off:off + ln] = k
    inter = D[frag[:, None] != frag[None, :]]
    if D.min() < min_d or (inter.size and inter.min() < min_inter):
        raise InfeasibleGeometry(f"overlap: min {D.min():.2f}, intermolecular {inter.min() if inter.size else 9.9:.2f}")


def build_tethered(base, solvent, n, d=None, seed=0, shell=0):
    """Family A. n = 0 is direct transfer; n >= 1 relay through n solvent molecules. d=None: natural (unconstrained)."""
    cl = Cluster()
    off, c, h, maps = _start(cl, "acac", seed)
    chain = [(c, h)]
    for k in range(n):
        hh, hx = chain[-1]
        u = _unit(cl.xyz[hx] - cl.xyz[hh]); r = R_CH_HB if cl.sym[hh] == "C" else R_HB["O"]
        sym, X, acc, dh, hbest = place_acceptor(cl, solvent, seed + 11 * k + 1, cl.xyz[hx], u, r)
        o2 = cl.add(solvent, sym, X)
        chain.append((o2 + acc, o2 + hbest))
    hh, hx = chain[-1]
    u = _unit(cl.xyz[hx] - cl.xyz[hh]); r = R_CH_HB if cl.sym[hh] == "C" else R_HB.get("N", 1.9)
    sym, X, acc, _, _ = place_acceptor(cl, base, seed + 97, cl.xyz[hx], u, r)
    b0 = cl.add(base, sym, X); sink = b0 + acc
    d_nat = float(np.linalg.norm(cl.xyz[sink] - cl.xyz[c]))
    if d is not None:
        set_distance(cl, c, sink, d, movable_from=cl.parts[1][1])
        check_geometry(cl)
    if shell:
        add_spectators(cl, solvent, shell, [maps[2]] + [a for a, _ in chain[1:]], seed)
    acceptors = [a for a, _ in chain[1:]] + [sink]
    moving = [(hx, dn, a) for (dn, hx), a in zip(chain, acceptors)]
    return dict(symbols=cl.sym, positions=cl.xyz, charge=cl.charge, moving=moving, parts=cl.parts,
                constraint=(c, sink) if d is not None else None, d_natural=round(d_nat, 3), d=d)


def build_shuttle(substrate, shuttles, topology, seed=0, shell=0, shell_solvent="water"):
    """Family B. shuttles: list of species names (e.g. ["water", "water"] or ["pyridone"]); [] = intramolecular shift.
    topology 'relay': all protons move (cyclic, the last shuttle H goes to the substrate C=O);
    'stepwise': only the C-H proton moves to the first shuttle (ion pair / protonated shuttle)."""
    cl = Cluster()
    off, c, h, maps = _start(cl, substrate, seed)
    target = cl.xyz[maps[2]]
    donor, hcur, moving = c, h, []
    for k, s in enumerate(shuttles):
        u = _unit(cl.xyz[hcur] - cl.xyz[donor]); r = R_CH_HB if cl.sym[donor] == "C" else R_HB["O"]
        sym, X, acc, dh, hbest = place_acceptor(cl, s, seed + 11 * k + 1, cl.xyz[hcur], u, r, target=target)
        o2 = cl.add(s, sym, X)
        moving.append((hcur, donor, o2 + acc))
        donor, hcur = o2 + dh, o2 + hbest
    moving.append((hcur, donor, maps[2]))
    if shell:
        add_spectators(cl, shell_solvent, shell, [maps[2]] + [m[2] for m in moving[:-1]], seed)
    if shuttles and topology == "stepwise":
        moving = moving[:1]
    return dict(symbols=cl.sym, positions=cl.xyz, charge=cl.charge, moving=moving, parts=cl.parts, constraint=None)


def _triol_chain_score(m, X, maps):
    """Prefer conformers with O1-H -> O2 and O2-H -> O3 H-bonds (triaxial chain)."""
    o1, o2, o3 = maps[1], maps[2], maps[3]
    h1, h2 = h_on(m, o1)[0], h_on(m, o2)[0]
    return -(np.linalg.norm(X[h1] - X[o2]) + np.linalg.norm(X[h2] - X[o3]) + 0.5 * np.linalg.norm(X[o1] - X[o3]))


def rotate_h_toward(X, mol, o, h, target):
    """Rotate hydroxyl H about its C-O bond so that it points at `target` (sets the O-H torsion)."""
    c = [n.GetIdx() for n in mol.GetAtomWithIdx(o).GetNeighbors() if n.GetSymbol() != "H"][0]
    ax = _unit(X[o] - X[c]); best, bd = X[h], 9e9
    for a in np.linspace(0, 2 * np.pi, 72, endpoint=False):
        p = (X[h] - X[o]) @ _axis_rot(ax, a).T + X[o]
        dd = float(np.linalg.norm(p - target))
        if dd < bd:
            best, bd = p, dd
    X = X.copy(); X[h] = best
    return X


def build_triol(base, topology, seed=0, solvent="water", shell=0):
    """Family C. relay: acac C-H -> O1; O1-H -> O2; O2-H -> O3; O3-H -> base (4 protons).
    direct: acac C-H -> base, triol as an H-bond donor spectator on an acac carbonyl O (same composition)."""
    cl = Cluster()
    off, c, h, maps = _start(cl, "acac", seed)
    tsym, tX, tmol, tm = monomer("triol", seed, conf_select=_triol_chain_score)
    o1, o2, o3 = tm[1], tm[2], tm[3]
    h1, h2, h3 = h_on(tmol, o1)[0], h_on(tmol, o2)[0], h_on(tmol, o3)[0]
    tX = rotate_h_toward(tX, tmol, o1, h1, tX[o2])          # O1-H -> O2, O2-H -> O3 (1,3-diaxial chain)
    tX = rotate_h_toward(tX, tmol, o2, h2, tX[o3])
    far = tX[o3] + 2.0 * _unit(tX[o3] - tX[o1])
    tX = rotate_h_toward(tX, tmol, o3, h3, far)             # O3-H points away from the ring, towards the base
    if topology == "relay":
        u = np.array([1.0, 0, 0]); hpos = cl.xyz[h]
        lp = _unit(-(sum(_unit(tX[n.GetIdx()] - tX[o1]) for n in tmol.GetAtomWithIdx(o1).GetNeighbors())))
        best, bs = None, -1e9
        for ang in np.linspace(0, 2 * np.pi, 72, endpoint=False):
            X = ((tX - tX[o1]) @ _rot(lp, -u).T) @ _axis_rot(u, ang).T + hpos + R_CH_HB * u
            s = cl.min_dist(X) - 0.05 * np.linalg.norm(X[h3] - X[o1])
            if s > bs:
                best, bs = X, s
        t0 = cl.add("triol", tsym, best)
        u2 = _unit(cl.xyz[t0 + h3] - cl.xyz[t0 + o3])
        sym, X, acc, _, _ = place_acceptor(cl, base, seed + 97, cl.xyz[t0 + h3], u2, R_HB.get("N", 1.9))
        b0 = cl.add(base, sym, X)
        moving = [(h, c, t0 + o1), (t0 + h1, t0 + o1, t0 + o2), (t0 + h2, t0 + o2, t0 + o3), (t0 + h3, t0 + o3, b0 + acc)]
    else:
        sym, X, acc, _, _ = place_acceptor(cl, base, seed + 97, cl.xyz[h], np.array([1.0, 0, 0]), R_CH_HB)
        b0 = cl.add(base, sym, X)
        site = maps[2]; a_sym, _, a_mol, _ = monomer("acac", 0)
        d = lone_pair_dirs(a_sym, cl.xyz[off:off + len(a_sym)], a_mol, site - off)[0]
        v = _unit(tX[h1] - tX[o1]); c0 = cl.xyz[site] + 1.85 * d
        Xs = (tX - tX[h1]) @ _rot(v, -d).T + c0
        cands = [(Xs - c0) @ _axis_rot(d, a).T + c0 for a in np.linspace(0, 2 * np.pi, 36, endpoint=False)]
        cl.add("triol", tsym, max(cands, key=cl.min_dist))
        moving = [(h, c, b0 + acc)]
    if shell:
        add_spectators(cl, solvent, shell, [maps[2]] + [i for i in acceptor_sites(cl) if i != maps[2]], seed)
    return dict(symbols=cl.sym, positions=cl.xyz, charge=cl.charge, moving=moving, parts=cl.parts, constraint=None)


def build_hq7(solvent, n, seed=0):
    """Positive control: 7-hydroxyquinoline O-H -> (solvent)n -> ring N (enol -> keto tautomer)."""
    cl = Cluster()
    sym, xyz, mol, maps = monomer("hq7", seed)
    o = maps[1]; h = h_on(mol, o)[0]
    u = _unit(xyz[h] - xyz[o]); xyz = (xyz - xyz[o]) @ _rot(u, [1, 0, 0]).T
    off = cl.add("hq7", sym, xyz); target = cl.xyz[off + maps[2]]
    chain = [(off + o, off + h)]
    for k in range(n):
        hh, hx = chain[-1]
        u = _unit(cl.xyz[hx] - cl.xyz[hh])
        s_sym, X, acc, dh, hbest = place_acceptor(cl, solvent, seed + 11 * k + 1, cl.xyz[hx], u, R_HB["O"], target=target)
        o2 = cl.add(solvent, s_sym, X)
        chain.append((o2 + acc, o2 + hbest))
    acceptors = [a for a, _ in chain[1:]] + [off + maps[2]]
    moving = [(hx, dn, a) for (dn, hx), a in zip(chain, acceptors)]
    return dict(symbols=cl.sym, positions=cl.xyz, charge=cl.charge, moving=moving, parts=cl.parts, constraint=None)
