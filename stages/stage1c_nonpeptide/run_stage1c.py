#!/usr/bin/env python3
"""Stage 1c: Non-peptide synthetic proton-relay systems (bulky bases, bicyclic amidines/guanidines,
and Lewis-acid-coordinated metal-acac complexes [Mg(acac)(H2O)2]+ and [Zn(acac)(H2O)2]+).
All runs are unconstrained (no artificial distance constraint on the base) at GFN2-xTB/ALPB(water) via tblite,
using 15-point 1D relaxed scans along the antisymmetric proton coordinate q = sum_i [d(D_i-H_i) - d(H_i-A_i)].

Families tested (3 seeds s0, s1, s2 each):
  E1. Bulky encumbered bases (no distance constraint!):
      - DTBMP (2,6-di-tert-butyl-4-methylpyridine, CID 138316) + acac + (H2O)_n, n=0 (direct), n=1 (relay), n=2 (relay)
      - TMP   (2,2,6,6-tetramethylpiperidine, CID 64963)      + acac + (H2O)_n, n=0 (direct), n=1 (relay), n=2 (relay)
      - DBU   (1,8-diazabicyclo[5.4.0]undec-7-ene, CID 81184) + acac + (H2O)_n, n=0 (direct), n=1 (relay), n=2 (relay)
  E2. Bifunctional guanidine shuttle:
      - TBD   (1,5,7-triazabicyclo[4.4.0]dec-5-ene, CID 79837) on acac (keto->enol) and cyclohexa-2,4-dien-1-one (dienone->phenol),
        comparing concerted cyclic relay (C-H -> N_imine; N_amine-H -> C=O) vs direct ion-pair formation (C-H -> N_imine only).
  E3. Metal-chelated complexes [M(acac)(H2O)2]+ (M = Mg2+, Zn2+, charge +1) + base (NMI or DTBMP):
      - Direct: M2+-chelated acac C3-H -> base (with 1 spectator water + 2 M-coordinated aqua ligands)
      - Outer-sphere water relay: M2+-chelated acac C3-H -> H2O -> base (with 2 M-coordinated aqua ligands)
      - Inner-sphere aqua-ligand relay: M2+-chelated acac C3-H -> base while M-OH2 donates H+ to acac C=O or base
"""
import json, os, sys, time
from pathlib import Path
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "stage1_model_relays"))
from build import _unit, _rot, _axis_rot, lone_pair_dirs  # noqa

OUT = Path(os.environ.get("STAGE1C_OUT", "/tmp/stage1c-nonpeptide"))

SMILES = {
    "acac": "CC(=O)[CH2:1]C(C)=[O:2]",
    "dienone": "[O:2]=C1[CH2:1]C=CC=C1",
    "water": "[OH2:1]",
    "nmi": "CN1C=C[N:1]=C1",
    "dtbmp": "CC(C)(C)c1cc(C)cc(C(C)(C)C)[n:1]1",
    "tmp": "CC1(C)CCCC(C)(C)[NH:1]1",
    "dbu": "C1CCC2=[N:1]CCCN2CC1",
    "tbd": "C1CN2CCC[N:1]=C2[NH:2]1",
}


def embed_monomer(name, seed=0):
    m = Chem.AddHs(Chem.MolFromSmiles(SMILES[name]))
    maps = {a.GetAtomMapNum(): a.GetIdx() for a in m.GetAtoms() if a.GetAtomMapNum()}
    if m.GetNumAtoms() > 1:
        AllChem.EmbedMolecule(m, randomSeed=seed + 101)
        try:
            AllChem.MMFFOptimizeMolecule(m)
        except Exception:
            pass
    xyz = np.array(m.GetConformer().GetPositions()) if m.GetNumAtoms() > 1 else np.zeros((1, 3))
    sym = [a.GetSymbol() for a in m.GetAtoms()]
    return sym, xyz, m, maps


def h_on(mol, heavy):
    return [n.GetIdx() for n in mol.GetAtomWithIdx(heavy).GetNeighbors() if n.GetSymbol() == "H"]


def outward_lp_dir(sym, xyz, mol, acc):
    """Outward direction along the lone pair of acceptor atom `acc` (away from its bonded heavy neighbors)."""
    nbs = [n.GetIdx() for n in mol.GetAtomWithIdx(acc).GetNeighbors()]
    v = sum(_unit(xyz[n] - xyz[acc]) for n in nbs)
    if np.linalg.norm(v) < 1e-4:
        return np.array([1.0, 0.0, 0.0])
    return _unit(-v)


def build_bulky_base_cluster(base_name, n_waters, topology, seed=0):
    """Build acac + (H2O)_n + bulky base (DTBMP, TMP, DBU) in 'relay' or 'direct' topology (equal composition)."""
    sym_a, xyz_a, mol_a, maps_a = embed_monomer("acac", seed)
    c3 = maps_a[1]
    o_carb = maps_a[2]
    h_c3 = h_on(mol_a, c3)[0]
    u = _unit(xyz_a[h_c3] - xyz_a[c3])
    xyz_a = (xyz_a - xyz_a[c3]) @ _rot(u, [1.0, 0.0, 0.0]).T

    sym_all = list(sym_a)
    xyz_all = xyz_a.copy()
    chain = [(c3, h_c3)]

    def min_inter_dist(X_new):
        d = np.linalg.norm(xyz_all[:, None, :] - X_new[None, :, :], axis=-1)
        return float(d.min())

    n_wire = n_waters if topology == "relay" else 0
    n_spec = 0 if topology == "relay" else n_waters

    rng = np.random.default_rng(seed + 17)
    for k in range(n_wire):
        d_idx, h_idx = chain[-1]
        u_dir = _unit(xyz_all[h_idx] - xyz_all[d_idx])
        r_hb = 2.00 if sym_all[d_idx] == "C" else 1.80
        wsym, wxyz, wmol, wmaps = embed_monomer("water", seed + 13 * k + 3)
        o_w = wmaps[1]
        hs_w = h_on(wmol, o_w)
        lp = outward_lp_dir(wsym, wxyz, wmol, o_w)
        base_X = (wxyz - wxyz[o_w]) @ _rot(lp, -u_dir).T
        best_X, best_s, best_h = None, -1e9, hs_w[0]
        for ang in np.linspace(0, 2 * np.pi, 36, endpoint=False) + rng.uniform(0, 0.15):
            X = base_X @ _axis_rot(u_dir, ang).T + xyz_all[h_idx] + r_hb * u_dir
            md = min_inter_dist(X)
            h_pick = max(hs_w, key=lambda h: float(np.dot(_unit(X[h] - X[o_w]), u_dir)))
            s = (0 if md >= 1.55 else -10 * (1.55 - md)) + 0.2 * md + float(np.dot(_unit(X[h_pick] - X[o_w]), u_dir))
            if s > best_s:
                best_X, best_s, best_h = X, s, h_pick
        off = len(sym_all)
        sym_all += wsym
        xyz_all = np.vstack([xyz_all, best_X])
        chain.append((off + o_w, off + best_h))

    # Place bulky base at end of chain
    d_idx, h_idx = chain[-1]
    u_dir = _unit(xyz_all[h_idx] - xyz_all[d_idx])
    r_hb = 2.15 if sym_all[d_idx] == "C" else 1.90
    bsym, bxyz, bmol, bmaps = embed_monomer(base_name, seed + 97)
    n_b = bmaps[1]
    lp = outward_lp_dir(bsym, bxyz, bmol, n_b)
    base_X = (bxyz - bxyz[n_b]) @ _rot(lp, -u_dir).T
    best_X, best_s = None, -1e9
    for r_try in (r_hb, r_hb + 0.25, r_hb + 0.50):
        for ang in np.linspace(0, 2 * np.pi, 48, endpoint=False) + rng.uniform(0, 0.1):
            X = base_X @ _axis_rot(u_dir, ang).T + xyz_all[h_idx] + r_try * u_dir
            md = min_inter_dist(X)
            s = (0 if md >= 1.65 else -20 * (1.65 - md)) + 0.1 * md - 0.5 * r_try
            if s > best_s:
                best_X, best_s = X, s
    off_b = len(sym_all)
    sym_all += bsym
    xyz_all = np.vstack([xyz_all, best_X])
    sink = off_b + n_b

    # Place spectator waters on acac carbonyl O atoms if direct topology
    if n_spec > 0:
        o_sites = [i for i, s in enumerate(sym_a) if s == "O"]
        for k in range(n_spec):
            site = o_sites[k % len(o_sites)]
            c_nb = [n.GetIdx() for n in mol_a.GetAtomWithIdx(site).GetNeighbors()][0]
            u_co = _unit(xyz_all[site] - xyz_all[c_nb])
            wsym, wxyz, wmol, wmaps = embed_monomer("water", seed + 41 * k + 7)
            o_w = wmaps[1]; h_w = h_on(wmol, o_w)[0]
            v_oh = _unit(wxyz[h_w] - wxyz[o_w])
            # try cone around u_co
            perp = _unit(np.cross(u_co, [0, 0, 1] if abs(u_co[2]) < 0.9 else [0, 1, 0]))
            best_w, best_ws = None, -1e9
            for tilt in (0.4, 0.7, -0.5):
                for phi in np.linspace(0, 2 * np.pi, 24, endpoint=False):
                    d_vec = _axis_rot(u_co, phi) @ (_axis_rot(perp, tilt) @ u_co)
                    c0 = xyz_all[site] + 1.85 * d_vec
                    Xw = (wxyz - wxyz[h_w]) @ _rot(v_oh, -d_vec).T + c0
                    md = min_inter_dist(Xw)
                    if md > best_ws:
                        best_w, best_ws = Xw, md
            sym_all += wsym
            xyz_all = np.vstack([xyz_all, best_w])

    acceptors = [a for a, _ in chain[1:]] + [sink]
    moving = [(hx, dn, a) for (dn, hx), a in zip(chain, acceptors)]
    return dict(symbols=sym_all, positions=xyz_all, charge=0, moving=moving, c3=c3, sink=sink)


def build_tbd_shuttle(sub_name, topology, seed=0):
    """Build substrate (acac or dienone) + TBD (bicyclic guanidine, CID 79837) in 'relay' (concerted dual H-bond shuttle)
    vs 'stepwise' (only C-H -> N_imine moves)."""
    sym_s, xyz_s, mol_s, maps_s = embed_monomer(sub_name, seed)
    c_don = maps_s[1]
    o_acc = maps_s[2]
    h_don = h_on(mol_s, c_don)[0]
    u = _unit(xyz_s[h_don] - xyz_s[c_don])
    xyz_s = (xyz_s - xyz_s[c_don]) @ _rot(u, [1.0, 0.0, 0.0]).T

    tsym, txyz, tmol, tmaps = embed_monomer("tbd", seed + 33)
    n_acc = tmaps[1]
    n_don = tmaps[2]
    h_tbd = h_on(tmol, n_don)[0]
    lp = outward_lp_dir(tsym, txyz, tmol, n_acc)
    base_X = (txyz - txyz[n_acc]) @ _rot(lp, -np.array([1.0, 0.0, 0.0])).T

    target = xyz_s[o_acc]
    best_X, best_s = None, -1e9
    for ang in np.linspace(0, 2 * np.pi, 72, endpoint=False):
        X = base_X @ _axis_rot([1.0, 0.0, 0.0], ang).T + xyz_s[h_don] + 2.05 * np.array([1.0, 0.0, 0.0])
        d_mat = np.linalg.norm(xyz_s[:, None, :] - X[None, :, :], axis=-1)
        md = float(d_mat.min())
        d_target = float(np.linalg.norm(X[h_tbd] - target))
        s = (0 if md >= 1.55 else -20 * (1.55 - md)) - 2.0 * d_target
        if s > best_s:
            best_X, best_s = X, s

    off = len(sym_s)
    sym_all = list(sym_s) + list(tsym)
    xyz_all = np.vstack([xyz_s, best_X])
    if topology == "relay":
        moving = [(h_don, c_don, off + n_acc), (off + h_tbd, off + n_don, o_acc)]
    else:
        moving = [(h_don, c_don, off + n_acc)]
    return dict(symbols=sym_all, positions=xyz_all, charge=0, moving=moving, c3=c_don, sink=off + n_acc)


def build_metal_complex(metal_sym, base_name, topology, seed=0):
    """Build [M(acac)(H2O)2]+ (M = Mg or Zn, charge +2 -> complex charge +2 with neutral keto acac, or wait:
    we are testing deprotonation of keto-acac C3-H while its two C=O oxygens are chelated to M2+(H2O)2!
    Wait: in [M(keto-acac)(H2O)2]2+, M2+ chelates the two carbonyl oxygens of keto-acac, polarizing C3-H!
    To keep the complex charge +1 (like a counterion-paired [M(Cl)(keto-acac)(H2O)2]+ or [M(OH)(keto-acac)(H2O)]+),
    let's check both [Mg(Cl)(acac)(H2O)]+ (charge +1, neutral overall with Cl-) or [Mg(acac)(H2O)2]2+ (charge +2)!
    Even cleaner: [Mg(Cl)(acac)(H2O)]+ or [Zn(Cl)(acac)(H2O)]+ with Cl- or AcO- coordinated to M2+ gives a mono-cationic (+1)
    or neutral (with 2 Cl-) complex! Let's use [M(Cl)2(acac)] (neutral, charge 0, tetrahedral MCl2 chelated to the two C=O of keto-acac),
    like SnCl4/TiCl4/ZnCl2/MgCl2 Lewis-acid catalysis in polyketide aldol reactions!
    Wait: MgCl2 and ZnCl2 are classic neutral Lewis-acid salts used with amines (Et3N/DIPEA/NMI) in mild soft-enolisation of beta-dicarbonyls!
    And because MCl2 is neutral (charge 0), [MCl2(acac)] + base + H2O has total charge 0, avoiding any high-charge continuum artifacts!
    """
    sym_a, xyz_a, mol_a, maps_a = embed_monomer("acac", seed)
    c3 = maps_a[1]
    o_indices = [i for i, s in enumerate(sym_a) if s == "O"]
    o1, o2 = o_indices[0], o_indices[1]

    # Force syn-diketone conformation so both carbonyl oxygens chelate MCl2 on the side opposite C3-H
    h_c3 = h_on(mol_a, c3)[0]
    u = _unit(xyz_a[h_c3] - xyz_a[c3])
    xyz_a = (xyz_a - xyz_a[c3]) @ _rot(u, [1.0, 0.0, 0.0]).T

    # Place M (Mg or Zn) chelating o1 and o2, shifted slightly away from the +x (C3-H) direction
    mid_o = 0.5 * (xyz_a[o1] + xyz_a[o2])
    u_away = _unit(mid_o - xyz_a[c3])
    if u_away[0] > 0.2:
        u_away[0] = -abs(u_away[0])
        u_away = _unit(u_away)
    m_pos = mid_o + 1.85 * u_away
    perp = _unit(np.cross(xyz_a[o2] - xyz_a[o1], u_away))
    cl1_pos = m_pos + 1.30 * u_away + 1.85 * perp
    cl2_pos = m_pos + 1.30 * u_away - 1.85 * perp

    sym_all = list(sym_a) + [metal_sym, "Cl", "Cl"]
    xyz_all = np.vstack([xyz_a, [m_pos, cl1_pos, cl2_pos]])

    def min_inter_dist(X_new):
        d = np.linalg.norm(xyz_all[:, None, :] - X_new[None, :, :], axis=-1)
        return float(d.min())

    chain = [(c3, h_c3)]
    rng = np.random.default_rng(seed + 29)

    if topology == "relay":
        # 1 water in wire between C3-H and base
        d_idx, h_idx = chain[-1]
        u_dir = _unit(xyz_all[h_idx] - xyz_all[d_idx])
        wsym, wxyz, wmol, wmaps = embed_monomer("water", seed + 19)
        o_w = wmaps[1]; hs_w = h_on(wmol, o_w)
        lp = outward_lp_dir(wsym, wxyz, wmol, o_w)
        base_X = (wxyz - wxyz[o_w]) @ _rot(lp, -u_dir).T
        best_X, best_s, best_h = None, -1e9, hs_w[0]
        for ang in np.linspace(0, 2 * np.pi, 36, endpoint=False) + rng.uniform(0, 0.15):
            X = base_X @ _axis_rot(u_dir, ang).T + xyz_all[h_idx] + 2.00 * u_dir
            md = min_inter_dist(X)
            h_pick = max(hs_w, key=lambda h: float(np.dot(_unit(X[h] - X[o_w]), u_dir)))
            s = (0 if md >= 1.55 else -10 * (1.55 - md)) + 0.2 * md + float(np.dot(_unit(X[h_pick] - X[o_w]), u_dir))
            if s > best_s:
                best_X, best_s, best_h = X, s, h_pick
        off_w = len(sym_all)
        sym_all += wsym
        xyz_all = np.vstack([xyz_all, best_X])
        chain.append((off_w + o_w, off_w + best_h))

    # Place base at end of chain
    d_idx, h_idx = chain[-1]
    u_dir = _unit(xyz_all[h_idx] - xyz_all[d_idx])
    r_hb = 2.10 if sym_all[d_idx] == "C" else 1.88
    bsym, bxyz, bmol, bmaps = embed_monomer(base_name, seed + 97)
    n_b = bmaps[1]
    lp = outward_lp_dir(bsym, bxyz, bmol, n_b)
    base_X = (bxyz - bxyz[n_b]) @ _rot(lp, -u_dir).T
    best_X, best_s = None, -1e9
    for r_try in (r_hb, r_hb + 0.25, r_hb + 0.50):
        for ang in np.linspace(0, 2 * np.pi, 48, endpoint=False) + rng.uniform(0, 0.1):
            X = base_X @ _axis_rot(u_dir, ang).T + xyz_all[h_idx] + r_try * u_dir
            md = min_inter_dist(X)
            s = (0 if md >= 1.65 else -20 * (1.65 - md)) + 0.1 * md - 0.5 * r_try
            if s > best_s:
                best_X, best_s = X, s
    off_b = len(sym_all)
    sym_all += bsym
    xyz_all = np.vstack([xyz_all, best_X])
    sink = off_b + n_b

    if topology == "direct":
        # Add 1 spectator water H-bonded to Cl or carbonyl O to keep composition identical to n=1 relay
        wsym, wxyz, wmol, wmaps = embed_monomer("water", seed + 53)
        o_w = wmaps[1]; h_w = h_on(wmol, o_w)[0]
        v_oh = _unit(wxyz[h_w] - wxyz[o_w])
        best_w, best_ws = None, -1e9
        for phi in np.linspace(0, 2 * np.pi, 24, endpoint=False):
            d_vec = _unit(np.array([-0.5, np.cos(phi), np.sin(phi)]))
            c0 = xyz_all[o1] + 1.90 * d_vec
            Xw = (wxyz - wxyz[h_w]) @ _rot(v_oh, -d_vec).T + c0
            md = min_inter_dist(Xw)
            if md > best_ws:
                best_w, best_ws = Xw, md
        sym_all += wsym
        xyz_all = np.vstack([xyz_all, best_w])

    acceptors = [a for a, _ in chain[1:]] + [sink]
    moving = [(hx, dn, a) for (dn, hx), a in zip(chain, acceptors)]
    return dict(symbols=sym_all, positions=xyz_all, charge=0, moving=moving, c3=c3, sink=sink)


def run_1d_scan(job_name, spec, NP=15):
    from ase import Atoms
    from ase.io import write
    from ase.constraints import FixInternals
    from ase.optimize import BFGS
    from tblite.ase import TBLite
    jdir = OUT / job_name
    jdir.mkdir(parents=True, exist_ok=True)
    if (jdir / "result.json").exists():
        return json.load(open(jdir / "result.json"))

    atoms = Atoms(spec["symbols"], positions=spec["positions"])
    mv = spec["moving"]
    c3, sink = spec["c3"], spec["sink"]
    calc = TBLite(method="GFN2-xTB", charge=spec["charge"], multiplicity=1, solvation=("alpb", "water"), verbosity=0, max_iterations=500)
    atoms.calc = calc

    t0 = time.time()
    # Step 1: Unconstrained reactant relaxation (check whether proton transfers spontaneously or base stays at natural distance)
    opt_r = BFGS(atoms, logfile=str(jdir / "opt_R.log"), maxstep=0.12)
    try:
        conv_r = bool(opt_r.run(fmax=0.03, steps=200))
        E_ref = float(atoms.get_potential_energy())
    except Exception as e:
        res = {"job": job_name, "status": "scf_error_R", "error": repr(e)[:200]}
        json.dump(res, open(jdir / "result.json", "w"), indent=2)
        return res

    Xr = atoms.get_positions().copy()
    write(str(jdir / "R_relaxed.xyz"), atoms)
    r_da_R = float(np.linalg.norm(Xr[c3] - Xr[sink]))
    h0, d0, a0 = mv[0]
    if np.linalg.norm(Xr[h0] - Xr[d0]) > 1.35:
        res = {"job": job_name, "status": "reactant_not_minimum", "R_DA_R": round(r_da_R, 3), "wall_s": round(time.time()-t0, 2)}
        json.dump(res, open(jdir / "result.json", "w"), indent=2)
        print(f"[{job_name}] reactant_not_minimum (R_DA={r_da_R:.2f} A)", flush=True)
        return res

    # Step 2: 15-point 1D relaxed scan along q = sum_i [d(D_i-H_i) - d(H_i-A_i)]
    q_fn = lambda pos: float(sum(np.linalg.norm(pos[h] - pos[d]) - np.linalg.norm(pos[h] - pos[a]) for h, d, a in mv))
    q0 = q_fn(Xr)
    q1 = 0.85 * len(mv)
    targets = np.linspace(q0, q1, NP)
    combo = [[d, h, 1.0] for h, d, a in mv] + [[h, a, -1.0] for h, d, a in mv]

    pts = []
    for k, t in enumerate(targets):
        atoms.set_constraint()
        pos = atoms.get_positions()
        delta = (t - q_fn(pos)) / (2.0 * len(mv))
        for h, d, a in mv:
            u = (pos[a] - pos[d]) / np.linalg.norm(pos[a] - pos[d])
            pos[h] += delta * u
        atoms.set_positions(pos)
        atoms.set_constraint(FixInternals(bondcombos=[[q_fn(pos), combo]]))
        atoms.calc = calc
        opt = BFGS(atoms, logfile=None, maxstep=0.1)
        try:
            conv = bool(opt.run(fmax=0.04, steps=120))
            E = float(atoms.get_potential_energy())
        except Exception:
            conv, E = False, None
            calc = TBLite(method="GFN2-xTB", charge=spec["charge"], multiplicity=1, solvation=("alpb", "water"), verbosity=0, max_iterations=500)
        Xn = atoms.get_positions()
        dE = round((E - E_ref) * 23.060548, 3) if E is not None else None
        hops = [[round(float(np.linalg.norm(Xn[h] - Xn[d])), 3), round(float(np.linalg.norm(Xn[h] - Xn[a])), 3)] for h, d, a in mv]
        pts.append({"k": k, "q": round(q_fn(Xn), 4), "dE_kcal": dE, "converged": conv, "R_DA": round(float(np.linalg.norm(Xn[c3] - Xn[sink])), 3), "hops": hops})

    valid = [p for p in pts if p["dE_kcal"] is not None]
    if len(valid) < 5:
        res = {"job": job_name, "status": "scan_failed"}
        json.dump(res, open(jdir / "result.json", "w"), indent=2)
        return res

    emin = min(p["dE_kcal"] for p in valid[:4])
    imax = max(range(len(valid)), key=lambda i: valid[i]["dE_kcal"])
    emax = valid[imax]["dE_kcal"]
    eend = valid[-1]["dE_kcal"]
    barrier = round(emax - emin, 2)

    # Step 3: Unconstrained product check
    atoms.set_constraint()
    opt_p = BFGS(atoms, logfile=None, maxstep=0.1)
    try:
        conv_p = bool(opt_p.run(fmax=0.04, steps=120))
        Ep = float(atoms.get_potential_energy())
        dE_p = round((Ep - E_ref) * 23.060548 - emin, 2)
    except Exception:
        conv_p, dE_p = False, None
    Xp = atoms.get_positions()
    write(str(jdir / "P_unconstrained.xyz"), atoms)
    prod_stays = bool(np.linalg.norm(Xp[mv[-1][0]] - Xp[mv[-1][2]]) < 1.25)

    res = {
        "job": job_name,
        "status": "ok",
        "natoms": len(spec["symbols"]),
        "R_DA_R": round(r_da_R, 3),
        "R_DA_TS": valid[imax]["R_DA"],
        "barrier_kcal": barrier,
        "dE_end_constrained_kcal": round(eend - emin, 2),
        "dE_rxn_unconstrained_kcal": dE_p,
        "product_is_minimum": prod_stays,
        "imax": imax,
        "hops_TS": valid[imax]["hops"],
        "wall_s": round(time.time() - t0, 2),
        "points": pts,
    }
    json.dump(res, open(jdir / "result.json", "w"), indent=2)
    print(f"[{job_name:36s}] barrier={barrier:6.2f} kcal/mol | R_DA(R->TS)={r_da_R:.2f}->{valid[imax]['R_DA']:.2f} A | dE_rxn={str(dE_p):6s} (prod_min={prod_stays}) | {res['wall_s']}s", flush=True)
    return res


def all_jobs():
    jobs = []
    # E1: Bulky bases (DTBMP, TMP, DBU) with n=0 (direct), n=1 (direct vs relay), n=2 (relay)
    for base in ("dtbmp", "tmp", "dbu"):
        for s in (0, 1, 2):
            jobs.append((f"E1_{base}_n0_direct_s{s}", lambda b=base, s=s: build_bulky_base_cluster(b, 0, "direct", s)))
            jobs.append((f"E1_{base}_n1_direct_s{s}", lambda b=base, s=s: build_bulky_base_cluster(b, 1, "direct", s)))
            jobs.append((f"E1_{base}_n1_relay_s{s}",  lambda b=base, s=s: build_bulky_base_cluster(b, 1, "relay", s)))
            jobs.append((f"E1_{base}_n2_relay_s{s}",  lambda b=base, s=s: build_bulky_base_cluster(b, 2, "relay", s)))
    # E2: Bifunctional guanidine TBD on acac and dienone
    for sub in ("acac", "dienone"):
        for top in ("relay", "stepwise"):
            for s in (0, 1, 2):
                jobs.append((f"E2_tbd_{sub}_{top}_s{s}", lambda sub=sub, top=top, s=s: build_tbd_shuttle(sub, top, s)))
    # E3: Lewis-acid metal complexes [MgCl2(acac)] and [ZnCl2(acac)] + base (nmi or dtbmp)
    for metal in ("Mg", "Zn"):
        for base in ("nmi", "dtbmp"):
            for top in ("direct", "relay"):
                for s in (0, 1, 2):
                    jobs.append((f"E3_{metal}Cl2_{base}_n1_{top}_s{s}", lambda m=metal, b=base, t=top, s=s: build_metal_complex(m, b, t, s)))
    return jobs


if __name__ == "__main__":
    jobs = all_jobs()
    if len(sys.argv) > 1 and sys.argv[1] == "--list":
        for name, _ in jobs:
            print(name)
        sys.exit(0)
    target = sys.argv[1] if len(sys.argv) > 1 else None
    for name, fn in jobs:
        if target and name != target:
            continue
        run_1d_scan(name, fn())
