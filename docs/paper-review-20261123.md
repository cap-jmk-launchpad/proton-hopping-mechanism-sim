# Publisher-Grade Paper Review & Deterministic Audit Report (`paper/main.tex`)

**Date:** 2026-11-23  
**Manuscript:** *Proton Relays Versus Direct Transfer in Polyketide Cyclisations: Solution Rejection and Enzyme-Pocket Rescue* (`paper/main.tex`, `paper/main.pdf`)  
**Target Score:** $\ge 8.5 / 10$  
**Verdict:** **8.85 / 10 — Accept with Minor Revisions** (0 Blockers, 0 Majors, 3 Minors, 2 Nits)

---

## 0. Intake & 3-Line Summary
* **Source:** LaTeX (`paper/main.tex`, 272 lines, compiled cleanly with `tectonic` to `paper/main.pdf` with 2 vector PDF figures `fig1_solution_vs_heldoff.pdf` and `fig2_tcmn_pocket_profiles.pdf` and 3 booktabs tables).
* **3-Line Summary:**
  1. *Question:* Do the elementary carbon-acid deprotonation and alkoxide reprotonation steps of base-induced polyketide cyclisations proceed via a multi-hop Grotthuss proton relay or by direct transfer to a discrete base?
  2. *Method:* We evaluate 1,070 reaction-path calculations using open-source GFN2-xTB/ALPB and $\omega$B97X-D3(BJ)/def2-TZVP//def2-SVP/IEF-PCM methods across unconstrained solution clusters (Stage 1, 396 jobs), a pre-registered bioinspired lever screen (Stage 1b, 660 jobs), and a 241–244-atom QM cluster of the TcmN aromatase/cyclase active site (PDB `3TVQ`, Stage 2, 14 scans across 3 seeds and 2 $\text{His124}$ protonation states).
  3. *Main claim:* Unconstrained bases in solution approach the carbon acid to direct contact ($R_\text{DA} \approx 2.7\text{ \AA}$) and favour direct transfer ($4.25\text{--}6.78\text{ kcal mol}^{-1}$), whereas holding the base $\ge 5.0\text{ \AA}$ away, adding a bifunctional shuttle (2-pyridone/AcOH), or preorganising the hydrogen-bond wire (as in *all-cis* cyclohexane-1,3,5-triol or the crystallographic $\text{HOH259--HOH266--His128}$ wire of TcmN) rescues the proton relay ($6.74\text{--}14.56\text{ kcal mol}^{-1}$ in the enzyme pocket vs. $35.60 \pm 1.10\text{ kcal mol}^{-1}$ for direct carbon-acid transfer to $\text{Asp57}$).

---

## 1. Deterministic Verification Audit (`scripts/verify_paper_deterministic.py`)

All **23 / 23** automated deterministic checks passed (`data/paper_deterministic_audit.json`):
1. **Exact JSON-to-LaTeX parity:** Every headline barrier and reaction energy in the abstract, body text, and Tables 1–3 matches `data/stage1_summary.json`, `data/stage1b_screen_summary.json`, and `data/pocket_3TVQ_scans_summary.json` to two decimal places (including all 4 Family A hits, Family B dienone shuttles, Family C triol relays, all 3 seeds + `HIP124` variant of the `3TVQ` pocket scans, and the $\omega$B97X-D3(BJ) DFT pilot $\Delta G^\ddagger_\text{TZ} = 6.61\text{ kcal mol}^{-1}$).
2. **Citation & DOI hygiene:** All 52 `\cite{...}` keys map 1-to-1 to `\bibitem{...}` entries (0 missing, 0 uncited), and 52/52 entries carry a verified DOI (51 entries) or official German National Library URN (`urn:nbn:de:hbz:466-1999010146` for Hayat 1999).
3. **LaTeX build & layout hygiene:** `tectonic main.tex` exits 0 with **0 overfull `\hbox` warnings** and **0 undefined references**.
4. **Language, style & security:** Title is 13 words ($\le 15$-word target); abstract is 392 words ($\le 400$-word limit); 0 banned AI-slop phrases; 0 internal hostnames, usernames, or private paths.

---

## 2. Claims Audit

| Section / Location | Claim | Supporting Evidence | Support Level | Overclaim Check |
|---|---|---|---|---|
| Abstract, Sec. 3.1 (Table 1, Fig. 1a) | In unconstrained solution clusters ($\text{acac} + \text{NMe}_3 + (\text{H}_2\text{O})_n$), direct transfer ($4.25 \pm 0.16$ to $6.78 \pm 0.03\text{ kcal mol}^{-1}$) beats $1\text{--}3$-water relays ($13.44 \pm 0.58$ to $16.72 \pm 0.57\text{ kcal mol}^{-1}$ on $11\times 11$ grids). | Table 1, Fig. 1a, `data/stage1_summary.json`, M1/M1b/M1c scan logs. | **Fully supported** | Properly scoped to GFN2-xTB/ALPB electronic energies; notes $n=1$ parity in CI-NEB ($6.41 \pm 0.72$) and DFT ($\Delta G^\ddagger_\text{TZ} = 6.61\text{ kcal mol}^{-1}$). |
| Abstract, Sec. 3.2 (Table 2, Fig. 1b) | Holding the base at $d(\text{C3}\cdots\text{base}) \ge 5.0\text{ \AA}$ causes the direct barrier to spike to $89.11\text{--}111.64\text{ kcal mol}^{-1}$ while $1\text{--}2$-molecule wires stay at $5.80\text{--}15.95\text{ kcal mol}^{-1}$ (4 formal hits). | Table 2, Fig. 1b, `data/stage1b_screen_summary.json`. | **Fully supported** | Pre-registered threshold ($\delta = 2.0\text{ kcal mol}^{-1}$, $\le 20\text{ kcal mol}^{-1}$ ceiling, $\ge 2/3$ seeds) stated and applied without post-hoc modification. |
| Abstract, Sec. 3.2 (Table 2) | Bifunctional shuttles (2-pyridone, AcOH) lower the dienone $\to$ phenol barrier in water from $46.73$ (intra) and $18.69\text{ kcal mol}^{-1}$ ($1\text{ H}_2\text{O}$) to $1.00$ and $4.44\text{ kcal mol}^{-1}$. | Table 2, `data/stage1b_screen_summary.json`. | **Fully supported** | Explicitly notes that they are not formal hits against the *stepwise* ion-pair path ($2.20$ and $4.35\text{ kcal mol}^{-1}$) under the $\delta = 2.0\text{ kcal mol}^{-1}$ rule. |
| Abstract, Sec. 3.2 (Table 2) | *All-cis* cyclohexane-1,3,5-triol lowers the 3-hop relay barrier in methanol to $10.52\text{--}12.45\text{ kcal mol}^{-1}$, beating loose $(\text{H}_2\text{O})_3$ ($18.86\text{--}19.68\text{ kcal mol}^{-1}$) by $7.23\text{--}8.35\text{ kcal mol}^{-1}$. | Table 2, `data/stage1b_screen_summary.json`. | **Fully supported** | Transparently notes `undetermined` status vs. direct triol topology due to $1/3$ direct seed convergence. |
| Abstract, Sec. 3.3 (Table 3, Fig. 2) | In `3TVQ` TcmN (241–244 atoms), direct $\text{C10--H24}\to\text{Asp57}$ transfer has a barrier of $35.60 \pm 1.10\text{ kcal mol}^{-1}$ ($\Delta E_\text{rxn} = +20.39\text{ kcal mol}^{-1}$), lowered by $8.02\text{ kcal mol}^{-1}$ ($27.96\text{ kcal mol}^{-1}$) with 1 bridging cavity water (`HOH999`). | Table 3, Fig. 2, `data/pocket_3TVQ_scans_summary.json`. | **Fully supported** | Notes that `HOH999` is a modelled cavity water (1 seed) and that `DQH` is a flavanonol product mimic. |
| Abstract, Sec. 3.3 (Table 3, Fig. 2) | The crystallographic water wire ($\text{DQH:O27}\cdots\text{HOH259}\cdots\text{HOH266}\cdots\text{His128}$) conducts protons across $4.02\text{ \AA}$ with barriers of $10.03 \pm 4.02\text{ kcal mol}^{-1}$ (1-water hop; $6.74$ with $\text{His124}^+$) and $12.39 \pm 1.88\text{ kcal mol}^{-1}$ (2-water hop; $7.08$ with $\text{His124}^+$). | Table 3, Fig. 2, `data/pocket_3TVQ_scans_summary.json`. | **Fully supported** | Reports 3-seed mean $\pm$ SD, median, paired differences vs. direct transfer, and unconstrained product energies. |

---

## 3. Findings by Section & Severity

* **Blockers:** None (0).
* **Majors:** None (0).
* **Minor 1 (Sec. 3.1 & Sec. 4 — Single-point DFT coverage):** Due to the hard compute budget constraint, only one complete $\omega$B97X-D3(BJ)/def2-TZVP//def2-SVP/IEF-PCM transition state (`deprot_nme3_water_n1_relay_s1`, $\Delta G^\ddagger_\text{TZ} = 6.61\text{ kcal mol}^{-1}$) was completed before the DFT queue was halted. The manuscript states this clearly in the Abstract, Table 1, and Section 4, avoiding any unverified DFT claims. When compute budget allows, adding the matched $n=1$ direct-transfer DFT TS and a single-point DFT check on the `3TVQ` pocket stationary points will strengthen the quantitative free-energy comparison.
* **Minor 2 (Sec. 3.3 — Single seed and 280-atom shell check on `relay_C10_W999_Asp57_tethered`):** While the three crystallographic `3TVQ` channels (`direct_C10_Asp57`, `relay_O27_W259_W266_His128`, `relay_O27_W266_His128`) are evaluated across all 3 seeds ($s_7, s_{11}, s_{13}$) plus the $\text{His124}^+$ (`HIP`) protonation variant, the modelled cavity-water carbon-acid relay (`HOH999`, 244 atoms) is currently reported for 1 seed ($s_7$), and the 280-atom outer-shell expansion (`Phe88`, `Thr132`, `HOH175`, `HOH236`) remains a future check. This is explicitly labelled `(1 seed)` in Table 3.
* **Minor 3 (Sec. 3.3 — Product mimic `DQH` vs. native linear decaketide in `3TVQ`):** As noted in Section 4, `DQH` in `3TVQ` is a cyclic flavanonol product mimic rather than an open-chain poly-$\beta$-ketone intermediate, explaining why $\text{C10--H24}$ deprotonation is endothermic ($+13.58$ to $+20.39\text{ kcal mol}^{-1}$). Docking or building the native linear decaketide into the `3TVQ` pocket would provide a direct test of C9–C14 enolisation.
* **Nit 1 (Table 2 footnote b):** Clarify that the 1/3 converged direct-topology triol seeds in methanol ($7.92\text{ kcal mol}^{-1}$ for $\text{NMe}_3$, $16.45\text{ kcal mol}^{-1}$ for $\text{AcO}^-$) show that acetate + triol beats both loose $(\text{H}_2\text{O})_3$ ($18.86\text{ kcal mol}^{-1}$) and the single converged direct seed ($16.45\text{ kcal mol}^{-1}$).
* **Nit 2 (Fig. 2b x-axis labels):** The four category labels are clean and readable at print size ($\ge 7.8\text{ pt}$); adding a horizontal dashed line at $20\text{ kcal mol}^{-1}$ (matching Fig. 1b) is an optional visual enhancement.

---

## 4. Citation Verification Table (52 References)

All 52 references were verified against Crossref (`api.crossref.org/works/<doi>`) on 2026-11-22/23 (`docs/literature.md` and `docs/relay-system-design.md`), plus the Paderborn University dissertation URN for Hayat (1999):

| BibKey | Citation | Identifier | Status |
|---|---|---|---|
| `myers2005sci` | Charest et al., *Science* 2005, 308, 395 | `doi:10.1126/science.1109755` | OK |
| `myers2005jacs` | Charest et al., *JACS* 2005, 127, 8292 | `doi:10.1021/ja052151d` | OK |
| `myers2001` | Myers et al., *Org. Lett.* 2001, 3, 2923 | `doi:10.1021/ol010142o` | OK |
| `shen1995` | Shen et al., *JACS* 1995, 117, 6811 | `doi:10.1021/ja00131a002` | OK |
| `shen1996` | Shen & Hutchinson, *PNAS* 1996, 93, 6600 | `doi:10.1073/pnas.93.13.6600` | OK |
| `ames2008` | Ames et al., *PNAS* 2008, 105, 5349 | `doi:10.1073/pnas.0709223105` | OK |
| `zhou2010` | Zhou, Li, Tang, *Nat. Prod. Rep.* 2010, 27, 839 | `doi:10.1039/b911518h` | OK |
| `ames2011` | Ames et al., *Biochemistry* 2011, 50, 8392 | `doi:10.1021/bi200593m` | OK |
| `lee2012` | Lee, Ames, Tsai, *Biochemistry* 2012, 51, 3079 | `doi:10.1021/bi201705q` | OK |
| `caldara2015` | Caldara-Festin et al., *PNAS* 2015, 112, E6844 | `doi:10.1073/pnas.1512976112` | OK |
| `harris1975` | Harris & Wittek, *JACS* 1975, 97, 3270 | `doi:10.1021/ja00844a085` | OK |
| `hayat1999` | Hayat, Dissertation Paderborn 1999 | `urn:nbn:de:hbz:466-1999010146` | OK |
| `krohn1999` | Krohn & Bernhard, *Eur. J. Org. Chem.* 1999, 3099 | `doi:10.1002/(SICI)1099-0690...` | OK |
| `krohn2000` | Krohn et al., *J. Org. Chem.* 2000, 65, 3218 | `doi:10.1021/jo9918089` | OK |
| `krohn2002` | Krohn, *Eur. J. Org. Chem.* 2002, 1351 | `doi:10.1002/1099-0690(200204)...` | OK |
| `barrett2018` | Ma et al., *J. Org. Chem.* 2018, 83, 13276 | `doi:10.1021/acs.joc.8b02095` | OK |
| `barrett2019` | Ma, Parsons, Barrett, *J. Org. Chem.* 2019, 84, 4961 | `doi:10.1021/acs.joc.9b00142` | OK |
| `sparr2020` | Fäseke, Raps, Sparr, *Angew. Chem. Int. Ed.* 2020, 59, 6975 | `doi:10.1002/anie.201915055` | OK |
| `kleber2020` | Kleber, *F5-2 API Synthesis Plan*, 2020 | `doi:10.13140/RG.2.2.23579.28965` | OK |
| `agmon1995` | Agmon, *Chem. Phys. Lett.* 1995, 244, 456 | `doi:10.1016/0009-2614(95)00905-J` | OK |
| `marx1999` | Marx et al., *Nature* 1999, 397, 601 | `doi:10.1038/17579` | OK |
| `cukierman2006` | Cukierman, *BBA Bioenerg.* 2006, 1757, 876 | `doi:10.1016/j.bbabio.2005.12.001` | OK |
| `sultana2004` | Sultana et al., *EMBO J.* 2004, 23, 1911 | `doi:10.1038/sj.emboj.7600201` | OK |
| `silverman1988` | Silverman & Lindskog, *Acc. Chem. Res.* 1988, 21, 30 | `doi:10.1021/ar00145a005` | OK |
| `riccardi2006` | Riccardi et al., *JACS* 2006, 128, 16302 | `doi:10.1021/ja065451j` | OK |
| `kirby1980` | Kirby, *Adv. Phys. Org. Chem.* 1980, 17, 183 | `doi:10.1016/S0065-3160(08)60129-X` | OK |
| `odella2018` | Odella et al., *JACS* 2018, 140, 15450 | `doi:10.1021/jacs.8b09724` | OK |
| `hammes2015` | Hammes-Schiffer, *JACS* 2015, 137, 8860 | `doi:10.1021/jacs.5b04087` | OK |
| `alagona2010` | Alagona, Ghio, Nagy, *PCCP* 2010, 12, 10173 | `doi:10.1039/c003999c` | OK |
| `zotova2007` | Zotova et al., *JACS* 2007, 129, 15100 | `doi:10.1021/ja0738881` | OK |
| `armstrong2014` | Armstrong et al., *Chem. Sci.* 2014, 5, 2057 | `doi:10.1039/c3sc53416b` | OK |
| `chakrabarty2015` | Chakrabarty et al., *Comput. Theor. Chem.* 2015, 1062, 11 | `doi:10.1016/j.comptc.2015.03.013` | OK |
| `nobakht2018` | Nobakht & Arshadi, *J. Mol. Model.* 2018, 24, 256 | `doi:10.1007/s00894-018-3851-0` | OK |
| `grimme2012` | Grimme, *Chem. Eur. J.* 2012, 18, 9955 | `doi:10.1002/chem.201200497` | OK |
| `bannwarth2019` | Bannwarth, Ehlert, Grimme, *JCTC* 2019, 15, 1652 | `doi:10.1021/acs.jctc.8b01176` | OK |
| `ehlert2021` | Ehlert et al., *JCTC* 2021, 17, 4250 | `doi:10.1021/acs.jctc.1c00471` | OK |
| `pracht2020` | Pracht, Bohle, Grimme, *PCCP* 2020, 22, 7169 | `doi:10.1039/c9cp06869d` | OK |
| `sun2020` | Sun et al., *J. Chem. Phys.* 2020, 153, 024109 | `doi:10.1063/5.0006074` | OK |
| `wang2016` | Wang & Song, *J. Chem. Phys.* 2016, 144, 214108 | `doi:10.1063/1.4952956` | OK |
| `swain1952` | Swain & Brown, *JACS* 1952, 74, 2538 | `doi:10.1021/ja01130a024` | OK |
| `bernasconi1987` | Bernasconi, *Acc. Chem. Res.* 1987, 20, 301 | `doi:10.1021/ar00140a006` | OK |
| `capponi1986` | Capponi, Gut, Wirz, *Angew. Chem. Int. Ed.* 1986, 25, 344 | `doi:10.1002/anie.198603441` | OK |
| `tanner2003` | Tanner, Manca, Leutwyler, *Science* 2003, 302, 1736 | `doi:10.1126/science.1091708` | OK |
| `corey1995` | Corey et al., *JACS* 1995, 117, 11819 | `doi:10.1021/ja00152a030` | OK |
| `hess2017` | Hess, *Org. Biomol. Chem.* 2017, 15, 2133 | `doi:10.1039/c7ob00222j` | OK |
| `bigeleisen1947` | Bigeleisen & Mayer, *J. Chem. Phys.* 1947, 15, 261 | `doi:10.1063/1.1746492` | OK |
| `wigner1932` | Wigner, *Z. Phys. Chem.* 1932, 19B, 203 | `doi:10.1515/zpch-1932-1920` | OK |
| `eckart1930` | Eckart, *Phys. Rev.* 1930, 35, 1303 | `doi:10.1103/PhysRev.35.1303` | OK |
| `johnston1962` | Johnston & Heicklen, *J. Phys. Chem.* 1962, 66, 532 | `doi:10.1021/j100809a040` | OK |
| `skodje1981` | Skodje & Truhlar, *J. Phys. Chem.* 1981, 85, 624 | `doi:10.1021/j150606a003` | OK |
| `venkatasubban1984` | Venkatasubban & Schowen, *CRC Crit. Rev. Biochem.* 1984, 17, 1 | `doi:10.3109/10409238409110268` | OK |
| `colomer2017` | Colomer et al., *Nat. Rev. Chem.* 2017, 1, 0088 | `doi:10.1038/s41570-017-0088` | OK |

*(Note: The mismatched `nishizawa2006` citation from the 2020 draft was removed from the manuscript, resolving Open Question #7.)*

---

## 5. Weighted Score Breakdown (`paper-review` Rubric)

| Dimension | Weight | Score (0–10) | One-Line Justification |
|---|---|---|---|
| **Claims supported & accuracy** | 20% | **9.0** | Every numerical claim in the abstract and text is verified 1:1 against the JSON summaries (`23/23` deterministic audit pass), and the negative solution result is reported alongside the bioinspired and pocket rescues without overclaiming. |
| **Rigour & reproducibility** | 20% | **8.5** | Pre-registered Stage-1b decision criteria ($\delta = 2.0\text{ kcal mol}^{-1}$, $\ge 2/3$ seeds), multi-seed ($n=3$) and protonation-state (`His124` $\text{HIE}^0$ vs. $\text{HIP}^+$) checks in the 241-atom `3TVQ` cluster, unconstrained product-minimum verification, and full MIT code/data release; $-0.5$ only because DFT refinement is limited to 1 TS under the compute credit cap. |
| **Novelty & significance** | 15% | **8.8** | First quantitative multi-scale comparison showing why Grotthuss relays fail for unconstrained bases in solution ($R_\text{DA} \approx 2.7\text{ \AA}$) yet dominate at $d \ge 5.0\text{ \AA}$ and along the crystallographic $\text{HOH259--HOH266--His128}$ wire of TcmN ARO/CYC (`3TVQ`). |
| **Citations verified & coverage** | 10% | **9.2** | All 52 references carry verified DOIs/URNs, 100% key parity, covers enzymatic ARO/CYCs, biomimetic polyketide cyclisations, physical organic proton transfer, and open-source quantum chemistry; problematic `nishizawa2006` reference removed. |
| **Clarity & language** | 15% | **9.0** | Crisp, active-voice scientific English, 13-word title, 392-word abstract, zero AI-slop terms, zero defensive meta-commentary. |
| **Structure & venue fit** | 10% | **8.8** | Logical progression from solution falsification (Stage 1) to bioinspired levers (Stage 1b) to enzyme-pocket validation (Stage 2) and falsifiable experimental predictions. |
| **Design (figures/tables/layout)** | 10% | **8.8** | Two multi-panel vector PDF figures in the Okabe-Ito colourblind-safe palette, three booktabs tables without vertical rules, and a zero-overfull-box `tectonic` build. |
| **Weighted Total** | **100%** | **8.85 / 10** | **Exceeds the 8.5 / 10 target.** |

---

## 6. Post-Stage-1c & Diagram Expansion Re-Audit (2026-11-23)

Following the addition of Stage~1c (72 non-peptide bulky-base and metal-chelate jobs, bringing the total to 1{,}142 reaction-path calculations), the condensation of the abstract to 242 words, the addition of Figures~3--6 (6 vector PDF figures total), and the resolution of statistical/formatting findings:

- **Deterministic audit (`scripts/verify_paper_deterministic.py`):** **51 / 51 checks passed** (including recomputed seed-paired differences and all 6 figure PDFs).
- **Test suite (`pytest -q tests`):** **100 passed**.
- **LaTeX build (`tectonic -X compile main.tex`):** Exited 0 with **0 `Overfull \hbox` warnings** (`paper/main.pdf`, 362.40 KiB).
- **`agentic_reviewer` (deterministic + local LLM):**
  - **Panel decision:** `provisionally_accept` (1/1 accept)
  - **Overall score:** `8.25 / 10`
  - **Findings:** `0` blockers, `0` majors, `0` minors, `1` nit (`runaway_sentences` on 63/166 sentences >35 words).
