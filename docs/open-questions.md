# Open chemistry questions for Julian
1. ~~Which step carries the mechanism?~~ **Answered (tentative, 2026-10-27):** the base-driven polyketide ring closure ("glaube ich").
2. Proton wire: solvent (water/HFIP/alcohol), base conjugate acid, or substrate-internal enol/C=O groups? *Assumption (2026-11-22, plan default):* stage 1 tests water and methanol wires; isopropanol (solvent of Krohn's K2CO3 cyclisation) is added in stage 2.
3. ~~Default substrate?~~ **Answered (2026-10-27):** simulate all three. **SMILES (2026-11-22, Julian: derive them yourself):** see `substrates/structures.yaml`; S1 (compound 10) is a documented stand-in because the F5-2 scheme was not retrievable.
4. Concerted relay (single TS) or stepwise hops with shallow intermediates? Open; stage 1 NEB bands show which.
5. Any experimental data (KIE, solvent effects, yields)? None in the original. Krohn/Hayat report 70 % for 34 -> 33 with K2CO3, iPrOH/CH2Cl2, 45 C.
6. Novelty vs. TcmN/WhiE single-water mechanism (Lee 2012): *Assumption:* both the multi-hop relay with ring coupling and the non-enzymatic setting.
7. Nishizawa 2006 (Org. Lett. 8, 5793): *Assumption:* not used as evidence until checked.
8. Which base/solvent for compound 10 and 12 simulations? *Assumption:* K2CO3-type conditions in iPrOH (from Hayat/Krohn) for S2; OH-/MeO- in water/MeOH for stage 1.
9. ~~ORCA licence?~~ **Answered (2026-11-22):** no ORCA; OSI open-source packages only.
10. Repo: ~~public?~~ **Answered (2026-11-22):** public under MIT.
