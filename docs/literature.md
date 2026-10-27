# Literature: proton relays in (polyketide) aldol/Claisen cyclisations

Search date 2026-10-27 (Crossref, web search; Europe PMC/OpenAlex were rate-limited). Every DOI below was
resolved via the Crossref API on that date. Summaries are based on abstracts/excerpts, not full reads unless noted.
**Prior-art verdict:** no paper found that proposes a multi-hop Grotthuss-type relay for base-induced
polyketide cyclisation. Closest prior art: (a) single **active-site water as base** in ARO/CYC enzymes
(TcmN/WhiE), (b) **water/alcohol proton relays** computed for organocatalytic aldol steps and keto–enol
tautomerism. The search was limited; a full-text search (SciFinder/Reaxys) is still needed before claiming novelty.

## Enzymatic polyketide cyclisation
| Ref | Summary | Relevance / support or contradiction |
|---|---|---|
| Ames et al., PNAS 2008, 105, 5349. 10.1073/pnas.0709223105 | TcmN ARO/CYC crystal structure; C9–C14 first-ring and C7–C16 second-ring aldol cyclisation in an interior pocket; Tyr35/Arg69 acid/base for ring 1, Tyr35/Arg82 and an **active-site water** for ring 2. | Strongest enzymatic precedent: water participates in proton transfer. Single water, not a wire → partially supportive; a relay would be an extension. |
| Lee, Ames, Tsai, Biochemistry 2012, 51, 3079. 10.1021/bi201705q | WhiE-ORFVI structure; revised general ARO/CYC mechanism: **Glu34–Tyr35-coordinated water is the base** abstracting H at C14, Arg69 protonates the C13 enolate O, water reprotonates the C9 oxyanion; water "replenished" for ring 2. | Very close to our steps 2–4 (deprotonation → ring closure → reprotonation via water). Key prior art; our novelty must be the *multi-hop relay / coupling of ring closures*, and in a non-enzymatic flow setting. |
| Ames et al., Biochemistry 2011, 50, 8392. 10.1021/bi200593m | ZhuI ARO/CYC (R1128 pathway) structure and characterisation; C7–C12 cyclisation specificity. | Pocket-defined folding; mechanistic proton details not checked in full text. |
| Caldara-Festin et al., PNAS 2015, 112. 10.1073/pnas.1512976112 | Structures of two di-domain ARO/CYCs (C7–C12 first ring). | Background for regio-control; not a proton-relay paper. |
| Sultana et al., EMBO J 2004, 23, 1911. 10.1038/sj.emboj.7600201 | SnoaL: intramolecular aldol by **direct Asp121 acid/base**, no Schiff base/metal; enolate stabilised by delocalisation. | Counter-example: direct (classical) proton transfer suffices in this cyclase. |
| Shen & Hutchinson, PNAS 1996, 93, 6600. 10.1073/pnas.93.13.6600; Shen et al., JACS 1995, 117, 6811. 10.1021/ja00131a002 | Without TcmN the decaketide cyclises spontaneously to shunt products; TcmN dictates regiospecificity. | Shows uncatalysed cyclisation occurs but unselectively — the selectivity question a relay could address. |
| Zhou, Li, Tang, Nat. Prod. Rep. 2010, 27, 839. 10.1039/b911518h | Review of aromatic polyketide cyclisation (bacteria, fungi). | Background. |

## Biomimetic (non-enzymatic) polyketide cyclisation
| Harris & Wittek, JACS 1975, 97, 3270. 10.1021/ja00844a085 | Biogenetic-type syntheses of polycyclic polyketides from partially protected β-hexa/heptaketones (base-mediated). | Experimental model reaction for the base cascade; no mechanistic proton data. |

## Proton relays in organocatalysis / enolisation (computational)
| Armstrong et al., Chem. Sci. 2014, 5, 2057. 10.1039/c3sc53416b | Houk–List TS revisited; includes **a water molecule as active proton-transfer relay** in proline aldol; reviews Patil/Sunoj methanol relays; notes entropy penalties. | Methodological template; warns that free energies (not total energies) must be compared. |
| Alagona, Ghio, Nagy, PCCP 2010, 12, 10173. 10.1039/c003999c | Keto–enol tautomerism of pyruvate and acetylacetone: **water-assisted (1–2 H2O) paths more favourable than direct**. | Supports relay-assisted enolisation of β-dicarbonyls — directly relevant to step 2. |
| Zotova et al., JACS 2007, 129, 15100. 10.1021/ja0738881 | Experimental role of water in proline aldol (kinetics). | Water affects rates/catalyst speciation; not specifically a relay. |
| Chakrabarty et al., Comput. Theor. Chem. 2015, 1062, 11. 10.1016/j.comptc.2015.03.013 | DFT: water-assisted H-bond network stabilises aldol TSs. | Supportive (organocatalysis). |
| Nobakht & Arshadi, J. Mol. Model. 2018, 24. 10.1007/s00894-018-3851-0 | DFT on proline aldol and water effect. | Supportive/context. |

## Grotthuss mechanism
| Agmon, Chem. Phys. Lett. 1995, 244, 456. 10.1016/0009-2614(95)00905-J | Classic description of Grotthuss proton hopping. | Conceptual basis. |
| Marx et al., Nature 1999, 397, 601. 10.1038/17579 | AIMD of the hydrated excess proton (Zundel/Eigen). | Method basis for stage 3. |

## Cationic cascades (epoxy-polyene model)
| Hess, Org. Biomol. Chem. 2017, 15, 2133. 10.1039/c7ob00222j | Computational studies of squalene cyclisation to steroids/hopenes. | Reference for epoxy-polyene model. |
| Corey et al., JACS 1995, 117, 11819. 10.1021/ja00152a030 | Cyclisation pathway of (S)-2,3-oxidosqualene. | Epoxy-polyene model background. |
| Nishizawa, Asai, Imagawa, Org. Lett. 2006, 8, 5793. 10.1021/ol062337x | Title per Crossref: "TiCl4 Induced Anti-Markovnikov Rearrangement". | ⚠ Original F5-2 cites this for ring formation with [TiCl4OH]⁻; title does not match that claim — Julian to check. |

## Gaps / next searches
SciFinder/Reaxys: "proton relay" AND (polyketide OR "poly-β-ketone") AND cyclization; QM/MM on TcmN/WhiE
(none found here); KIE data on ARO/CYC or Harris-type cyclisations (none found).
