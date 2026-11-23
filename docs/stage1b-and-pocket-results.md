# Stage-1b Bioinspired Screen & 3TVQ Enzyme-Pocket Proton-Relay Results (2026-11-23)

Scope: mechanism-level physical organic chemistry only (GFN2-xTB/ALPB via `tblite` and `ASE`, strictly OSI open-source tools).

---

## 1. Stage-1b Bioinspired Lever Screen (660 jobs, `data/stage1b_screen_summary.json`)

The 660-job screen across Families A–C and positive control Q was evaluated with `stages/stage1b_relay_design/analyze1b.py` under the pre-registered criteria in `docs/relay-system-design.md` Section 4 ($\delta = 2.0\text{ kcal/mol}$, $\ge 2/3$ usable seeds per condition, taking $\min(\text{CI-NEB}, \text{segmented sequential-hop})$ when both converge to an interior maximum).

* **Status breakdown (660 jobs):** `253 ok`, `155 geometry_infeasible` (compressed wires at short $d$, skipped in ms), `110 neb_not_converged`, `82 product_not_minimum`, `35 error`, `21 timeout`, `4 reactant_not_minimum`.
* **Positive Control Q (7-hydroxyquinoline + $(\text{H}_2\text{O})_n / (\text{MeOH})_n$):** **PASSES** in water ($n=2$: $20.79\text{ kcal/mol}$; $n=3$: $12.54\text{ kcal/mol}$; $n=0$ unviable).

### 1.1 Family A — Base Held Off the Carbon Acid at Fixed $d(\text{C3}\cdots\text{Base}) \in \{4, 5, 6, 7\}\text{ \AA}$ (4 Formal HITs)

Holding the base off the active-methylene $\text{C3–H}$ of acetylacetone mimics an enzyme pocket or tethered catalyst where direct contact is prevented:

| Base | Solvent | $d(\text{C3}\cdots\text{Base})$ (Å) | Direct ($n=0$) Barrier (kcal/mol) | Best Relay Barrier (kcal/mol) | Best Wire Length $n$ | Verdict |
|---|---|---|---|---|---|---|
| `NMe3` | methanol | 5.0 | 89.11 | **13.23** | $n=1$ | **HIT** ($\Delta\Delta E^\ddagger = -75.88$) |
| `N-methylimidazole` | water | 5.0 | 89.83 | **15.95** | $n=1$ | **HIT** ($\Delta\Delta E^\ddagger = -73.88$) |
| `N-methylimidazole` | methanol | 6.0 | 106.07 | **9.74** | $n=1$ | **HIT** ($\Delta\Delta E^\ddagger = -96.33$) |
| `N-methylimidazole` | methanol | 7.0 | 111.65 | **8.30** | $n=2$ | **HIT** ($\Delta\Delta E^\ddagger = -103.34$) |
| `NMe3` | water | 5.0 | *unusable* (1/3: 11.50) | **6.22** | $n=1$ | `undetermined_direct_unusable` |
| `NMe3` | water | 7.0 | *unusable* (0/3) | **5.80** | $n=2$ | `undetermined_direct_unusable` |
| `NMe3` | methanol | 6.0 | *unusable* (0/3) | **7.65** | $n=2$ | `undetermined_direct_unusable` |
| `NMe3` | methanol | 7.0 | *unusable* (0/3) | **13.46** | $n=2$ | `undetermined_direct_unusable` |
| `N-methylimidazole` | methanol | 5.0 | *unusable* (1/3) | **9.58** | $n=1$ | `undetermined_direct_unusable` |
| `acetate` | water | 5.0 | *unusable* (1/3: 9.85) | **14.24** | $n=1$ | `undetermined_direct_unusable` |
| `acetate` | water | 6.0 | *unusable* (1/3: 108.32) | **12.00** | $n=1$ | `undetermined_direct_unusable` |
| `acetate` | water | 7.0 | *unusable* (0/3) | **18.06** | $n=2$ | `undetermined_direct_unusable` |
| `acetate` | methanol | 6.0 | *unusable* (0/3) | **11.12** | $n=1$ | `undetermined_direct_unusable` |

* **Finding:** Once $d(\text{C3}\cdots\text{Base}) \ge 5.0\text{ \AA}$, direct proton transfer spikes to $89\text{--}112\text{ kcal/mol}$ (or fails to converge), whereas a $1\text{--}2$ molecule water/methanol wire bridges the gap with a flat barrier of **$5.8\text{--}18.1\text{ kcal/mol}$**.

### 1.2 Family B — Bifunctional Shuttles in Aromatising Tautomerisation (`cyclohexa-2,4-dien-1-one` $\to$ `phenol`)

| Substrate | Shuttle | Solvent | Concerted Relay (kcal/mol) | Stepwise Ion-Pair (kcal/mol) | $1\text{ H}_2\text{O}$ Relay (kcal/mol) | Intramolecular Shift (kcal/mol) | Verdict |
|---|---|---|---|---|---|---|---|
| `dienone` | `2-pyridone` | water | **1.00** | 2.20 | 18.69 | 46.73 | `no` ($\Delta = -1.20 < 2.0$ vs. stepwise; $-17.69$ vs. $\text{H}_2\text{O}$) |
| `dienone` | `2-pyridone` | methanol | **0.06** | 0.17 | 0.22 | 47.63 | `no` (near-barrierless for both) |
| `dienone` | `AcOH` | water | **4.44** | 4.35 | 18.69 | 46.73 | `no` (relay $\approx$ stepwise; $-14.25$ vs. $\text{H}_2\text{O}$) |
| `dienone` | `AcOH` | methanol | **4.73** | 4.67 | 0.22 | 47.63 | `no` |
| `acac` | `2-pyridone` | methanol | **16.69** | 16.10 | *unusable* | *unusable* | `undetermined` |
| `acac` | `AcOH` | methanol | **11.00** | *unusable* (1/3: 8.28) | *unusable* | *unusable* | `undetermined` |

* **Finding:** While neither bifunctional shuttle beats the *stepwise* path by the strict $\delta = 2.0\text{ kcal/mol}$ margin (because the dienone aromatisation driving force makes both paths extremely fast), **2-pyridone** ($1.00\text{ kcal/mol}$) and **AcOH** ($4.44\text{ kcal/mol}$) lower the dienone tautomerisation barrier in water by **$14.2\text{--}17.7\text{ kcal/mol}$** relative to a single water molecule ($18.69\text{ kcal/mol}$) and by **$> 42\text{ kcal/mol}$** relative to the unassisted intramolecular shift ($46.73\text{ kcal/mol}$).

### 1.3 Family C — Preorganised *all-cis* Cyclohexane-1,3,5-Triol Wire (3-Hop Relay)

| Base | Solvent | Triol Relay Barrier (kcal/mol) | Loose $(\text{H}_2\text{O})_3$ Relay Barrier (kcal/mol) | Direct with Triol Spectator | Verdict |
|---|---|---|---|---|---|
| `NMe3` | methanol | **12.45** | 19.68 | *unusable* (1/3: 7.92) | `undetermined` (**$-7.23\text{ kcal/mol}$ vs. loose $(\text{H}_2\text{O})_3$**) |
| `acetate` | methanol | **10.52** | 18.86 | *unusable* (1/3: 16.45) | `undetermined` (**$-8.35\text{ kcal/mol}$ vs. loose $(\text{H}_2\text{O})_3$**) |
| `NMe3` | water | **16.74** | 14.54 | *unusable* (1/3: 5.38) | `undetermined` |
| `acetate` | water | **26.30** | 15.61 | *unusable* (0/3) | `undetermined` |

* **Finding:** In methanol, covalently preorganising the three OH groups on *all-cis* cyclohexane-1,3,5-triol lowers the 3-hop relay barrier by **$7.2\text{--}8.4\text{ kcal/mol}$** compared to a loose 3-water wire.

---

## 2. Enzyme-Pocket Proton-Relay Model (`3TVQ` TcmN ARO/CYC, 241–244 atoms, `data/pocket_3TVQ_scans_summary.json`)

### 2.1 Structure Selection & QM Cluster Setup

Four ARO/CYC crystal structures (`3TVQ`, `2RER`, `3TL1`, `3TFZ`) were surveyed:
* **`3TVQ`** (TcmN ARO/CYC from *Streptomyces glaucescens*, **1.67 Å** resolution, wild-type, bound product-mimic *trans*-dihydroquercetin `DQH` at occupancy 1.0) was chosen because it is the highest-resolution WT structure with an explicit crystal water wire bridging the ligand in the interior pocket to a catalytic histidine:
  $$\text{DQH:O27} \xrightarrow{2.87\text{ \AA}} \text{HOH259} \xrightarrow{2.64\text{ \AA}} \text{HOH266} \xrightarrow{2.87\text{ \AA}} \text{His128:NE2}$$
  (with `HOH266` also directly within $3.12\text{ \AA}$ of `DQH:O27` in the crystal, tightening to $2.745\text{ \AA}$ upon relaxation).
* **QM Cluster (`241 atoms`, net charge `0`, singlet, $914\text{ e}^-$):**
  * Contiguous backbone+sidechain segments: `Asp57–Ala58–Asn59`, `His124–Met125`, `His128`.
  * $\text{C}_\alpha$-anchored sidechains: `Trp63`, `Trp65`, `Arg82`, `Phe120`.
  * 10 crystal waters within $4.5\text{ \AA}$ of `DQH`: `HOH176, 225, 233, 234, 255, 259, 264, 266, 280, 298`.
  * Ligand: `DQH` (34 atoms).
  * **Boundary constraints:** **24 frozen atoms** (all terminal $\text{C}_\alpha$ and H-cap atoms) enforced via ASE `FixAtoms` (`0.000000 Å` drift during optimization; heavy/all-atom RMSD = $0.445\text{ \AA}$ from crystal).
  * **Protonation states (PROPKA 3.5.1, pH 7.0):** `Asp57` anionic ($\text{p}K_\text{a}\ 0.55$, $-1$), `Arg82` cationic ($\text{p}K_\text{a}\ 9.41$, $+1$), `His128` neutral `HID` ($\text{p}K_\text{a}\ 5.01$, `ND1–H` protonated so unprotonated `NE2` faces the water wire), `His124` ambiguous ($\text{p}K_\text{a}\ 7.35$, tested both as neutral `HIE` at total charge `0` across 3 independent water/OH orientation seeds `s7, s11, s13` and as doubly-protonated `HIP` at total charge `+1`, 242 atoms).

### 2.2 1D Relaxed Scan Results Across Channels, Seeds, and Protonation States

15-point 1D relaxed scans (`fmax = 0.04 eV/Å`, GFN2-xTB/ALPB(water) via `tblite`, ASE `FixAtoms` + `FixInternals`) along the antisymmetric proton coordinate $q = \sum_i [d(\text{D}_i\text{–H}_i) - d(\text{H}_i\text{–A}_i)]$:

| Channel | Pathway | Seed `s7` $\Delta E^\ddagger$ | Seed `s11` $\Delta E^\ddagger$ | Seed `s13` $\Delta E^\ddagger$ | **3-Seed Mean ± SD (Median)** | `His124` `HIP` (+1) $\Delta E^\ddagger$ | Unconstrained Product $\Delta E_\text{rxn}$ (`s7`) |
|---|---|---|---|---|---|---|---|
| **`direct_C10_Asp57`** | Direct carbon-acid transfer: `DQH:C10–H24` $\to$ `Asp57:OD2` | 35.98 | 34.36 | 36.46 | **35.60 ± 1.10 (35.98)** | **33.66** | **+20.39** (stable min, $d(\text{OD2–H}) = 1.009\text{ \AA}$) |
| **`relay_C10_W999_Asp57_tethered`** | Modelled 1-water carbon-acid relay (244 atoms): `DQH:C10–H24` $\to$ `HOH999` $\to$ `Asp57:OD2` | 27.96 | — | — | **27.96** (1 seed) | — | **+13.58** (stable min, $d(\text{OD2–H}) = 1.006\text{ \AA}$) |
| **`relay_O27_W259_W266_His128`** | 2-water crystal wire (hydroxyl/alkoxide channel): `DQH:O27–H31` $\to$ `HOH259` $\to$ `HOH266` $\to$ `His128:NE2` | 14.56 | 11.33 | 11.29 | **12.39 ± 1.88 (11.33)** | **7.08** | **$-2.64$** (stable min, $d(\text{NE2–H}) = 1.048\text{ \AA}$) |
| **`relay_O27_W266_His128`** | 1-water crystal wire shortcut: `DQH:O27–H31` $\to$ `HOH266` $\to$ `His128:NE2` | 7.01 | 14.60 | 8.49 | **10.03 ± 4.02 (8.49)** | **6.74** | **$-0.48$** (stable min, $d(\text{NE2–H}) = 1.046\text{ \AA}$) |

### 2.3 Key Mechanistic Conclusions from the `3TVQ` Pocket Model

1. **Pocket geometry disfavours direct $\text{C10–H24} \to \text{Asp57}$ transfer (`35.60 ± 1.10 kcal/mol`):**
   Although `Asp57:OD2` is $3.325\text{ \AA}$ from `DQH:C10`, the crystal pocket holds the $\text{C10–H24}$ bond at an acute angle ($\angle(\text{C10–H24}\cdots\text{OD2}) = 55.9^\circ$) while `Asp57` is locked by backbone and `Trp63:HE1` ($1.82\text{ \AA}$) hydrogen bonds.
2. **A bridging water rescues carbon-acid deprotonation by `8.0 kcal/mol` (`27.96` vs. `35.98 kcal/mol`, $\Delta E_\text{rxn}: +13.58$ vs. $+20.39\text{ kcal/mol}$):**
   Inserting a single bridging water (`HOH999`) into the pocket cavity between `DQH:C10–H24` and `Asp57:OD2` ($d(\text{H24}\cdots\text{O}) = 1.89\text{ \AA}$, $d(\text{O}\cdots\text{OD2}) = 2.67\text{ \AA}$) relieves the angular strain via a 2-hop relay ($\text{C10–H24} \to \text{HOH999} \to \text{Asp57:OD2}$) and stabilizes the resulting enolate by $6.8\text{ kcal/mol}$.
3. **The crystallographic water wire (`DQH:O27` $\leftrightarrow$ `HOH259`/`HOH266` $\leftrightarrow$ `His128:NE2`) provides a fast, reversible proton relay (`6.7–14.6 kcal/mol`):**
   * Across all 3 seeds and both `His124` protonation states, the held crystal water wire between `DQH:O27` and `His128:NE2` ($d(\text{O27}\cdots\text{NE2}) = 4.02\text{ \AA}$) operates with a low barrier (**10.03 ± 4.02 kcal/mol** via 1 water; **12.39 ± 1.88 kcal/mol** via 2 waters) and nearly thermoneutral reaction energy ($\Delta E_\text{rxn} = -0.48\text{ to }-2.64\text{ kcal/mol}$).
   * Protonating `His124` (`HIP`, charge `+1`) electrostatically stabilizes the `DQH:O27`$^{-}$ oxyanion and lowers the 2-water and 1-water relay barriers further to **7.08** and **6.74 kcal/mol** ($\Delta E_\text{rxn} = -4.46\text{ kcal/mol}$).
   * This directly supports Hypothesis **H2** (rapid alkoxide reprotonation/deprotonation through a preorganised water wire in a confined pocket) and complements Family A of Stage 1b (Hypothesis **H1** rescued when the base is held $\ge 5.0\text{ \AA}$ from the carbon acid).
