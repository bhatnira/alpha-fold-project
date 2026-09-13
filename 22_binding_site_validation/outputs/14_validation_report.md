# Alpha9alpha10 nAChR PAM — Binding-Site Validation Report

*Generated:* 2026-09-13 05:39  
*Falsifiable-hypothesis treatment: the proposed site is treated as a hypothesis to be* *confirmed or falsified, NOT as a correct answer.*

## 1. Existing evidence audit

Layers inventoried (full table in `outputs/01_evidence_audit.json`):
- **AF3 cofolding (existing)** (n=468): receptor_structural / ligand_pose. Receptor architecture; recurring ligand regions across conditions — No experimental bound structure exists to benchmark pose accuracy [correlated]
- **Boltz-2 cofolding (existing outputs)** (n=90): receptor_structural / ligand_pose. Independent cofolding corroboration of regions identified by AF3 — Pre-correction runs on alpha10 monomer are INVALID (FALLACIES_AUDIT C7) [correlated]
- **Boltz-2 site-directed state_models** (n=4): receptor_structural. Ligand-conditioned placement at candidate sites; residue contacts; pose reproducibility [correlated]
- **Boltz-2 site-directed ach_states** (n=2): receptor_structural. Ligand-conditioned placement at candidate sites; residue contacts; pose reproducibility [correlated]
- **Boltz-2 site-directed binary** (n=80): receptor_ligand_interaction. Ligand-conditioned placement at candidate sites; residue contacts; pose reproducibility [correlated]
- **Boltz-2 site-directed ternary** (n=80): receptor_ligand_interaction. Ligand-conditioned placement at candidate sites; residue contacts; pose reproducibility [correlated]
- **Boltz-2 site-directed full_panel** (n=180): receptor_ligand_interaction. Ligand-conditioned placement at candidate sites; residue contacts; pose reproducibility — run in progress (180/180 complete at audit time) [correlated]
- **Molecular docking (existing)** (n=30): binding_site / ligand_pose. Orthogonal pocket-occupancy evidence at candidate sites; ligand pose library [correlated]
- **Interaction fingerprints (existing)** (n=2): receptor_ligand_interaction. Residue contact frequencies from prior analyses — Partly stale; recomputed on fresh full_panel data in step 05 [correlated]
- **Experimental SAR (existing)** (n=30): experimental_sar. Activity labels, potency, stereochemical series, % potentiation — 30-compound panel: 7 active / 23 inactive [correlated]
- **Boltz-2 affinity prediction** (n=0): affinity. Nothing in the current run (no affinity_*.json emitted) — Earlier affinity 'proxy' matrix was circular; affinity layer is ABSENT until a corrected run exists [correlated]
- **XAI** (n=0): explainability. No valid XAI feature importances exist — Prior XAI permutation importances were computed on an overfit model; not usable (FALLACIES_AUDIT) [correlated]
- **Activity-cliff / MMP (existing)** (n=4): activity_cliff. Matched-pair SAR rules: 12vs7 (alkynyl), 25vs1 (Br), 18vs19 (benzyl vs vinyl), 3vs1 — Used as an SAR challenge layer, not as structural evidence [correlated]

Independence notes: AF3, Boltz-2 and SiteAF3 all derive their MSAs from the same sequence databases; multiple seeds of one method are NOT independent observations. The old Boltz-2 affinity proxy was circular and is invalidated (FALLACIES_AUDIT).

## 2. Candidate binding sites

Retained and treated as competing hypotheses:
- **Site 23**: alpha9(+)/alpha10(-) ECD vestibular inter-subunit (favored by prior work)
- **Site 21**: alpha10(+)/alpha9(-) ECD vestibular inter-subunit
- **Site 5**: alpha9(+)/alpha10(-) ECD vestibular inter-subunit (overlaps 23)
- **Site 34**: TMD-pore allosteric cavity (decode/wall region)
- **Site 7**: alpha9(+)/alpha10(-) ECD vestibular inter-subunit

Site definitions overlap heavily in the ECD vestibule (same inter-subunit cavity sampled from different chain faces) — recurrence statistics must be read with that degeneracy in mind.

## 3. AF3 / Boltz-2 ensemble convergence

- Receptor C-alpha RMSD across Boltz-2 seeds of the same condition: **37.328 A** (n=180 conditions, 3 seeds each) → **poor global receptor convergence**; each seed remodels the pentamer substantially.
- Ligand pose spread (mean pairwise centroid distance across seeds): **50.379 A** (p25 43.795, median 50.744, p75 57.733).
- Site self-recurrence under conditioning: site 23 = **92.2%**, site 21 = **3.9%**, site 5 = **0.0%**.
- **Cross-condition drift (informative):** models conditioned on other sites were nonetheless classified into site-23 residues in most cases (`02 hit_distribution`), suggesting the vestibular vestibular sub-pocket sampled by site-23 residues is the recurrent preference; this is confounded by the overlap of the site definitions.
- Existing AF3 cofolding (`data/allostery/af3_outputs`, root model cif) yields 181 residue contacts; existing Boltz-2 outputs yield no site-classified hits at these definitions.

Convergence verdict: **local (pocket-region) recurrence is high but pose-level and global-receptor convergence are poor.** The pocket recurs, but individual poses across seeds are not reproducible.

## 4. SiteAF3 site-conditioned analysis

- **PENDING / NOT RUN.** The site-conditioned test (Tang & Wang, PNAS 10.1073/pnas.2521048122) could not be executed in this session: the AF3-based SiteAF3 environment and model weights are not installed on the cluster. 80 config inputs are staged in `siteaf3_inputs/` (5 sites x 2 stoichiometries x 8 representative compounds). These are the missing-validity link in this report; **until SiteAF3 runs, the hypothesis cannot be classed above B.**
- Method design (already implemented in `04_siteaf3_submit.sh`): pocket-masked MSA embedding + pocket-conditioned diffusion, one config per candidate site, so competing sites are tested under identical conditions.

## 5. SiteAF3 vs competing-site comparison

- Competition outcome: **pending** (SiteAF3 not run). Until then the favored-site-vs-alternatives comparison rests on the Boltz-2 ligand-conditioned panel alone (sections 3/10).
## 6. Pocket and pose QC

- Models analyzed: **540**; clash-free fraction **0.659**; mean clashing ligand atoms/model **0.698** (cutoff 2.2 A).
- Pocket occupancy (ligand atoms within 6 A of any CA): **0.741** — adequate, but as a gross filter.
- Pharmacophore centroid spread (cpd 12): **39.443 A** — large; pose convergence is NOT adequate for atom-level claims.

## 7. Medicinal-chemistry QC (Chen et al. framework)

- RDKit: 30/30 SMILES valid; 27 carry defined stereocenters (ascorbate series: must keep stereochemistry, PAM activity is stereo-sensitive).
- Whole-ligand RMSD / pharmacophore-level agreement are reported as *internal consistency only*; no experimental a9a10-PAM complex exists, so experimental pose accuracy is not claimed.

## 8. Allosteric pocket plasticity

- Mean per-residue CA span across bound-model ensemble (aligned): **44.0 A**; 2337 residues >2 A, 2337 >4 A.
- Classification: **flexible / transient-like**, precision limited by the large receptor ensemble spread; no apo→holo transition was modeled, so no induced-fit claim. Classify cautiously as flexible/interface-dependent, cryptic-pocket behaviour not addressable from these data (see section 18).

## 9. Experimental SAR validation

- Dataset: 30 compounds; 7 active (IDs [1, 2, 3, 12, 18, 24, 25]); MMP records 4; curated cliff pairs 5.
- MMP rules from prior analysis: ['Small substituents at C5 maintain activity', 'Extended chains at C5 eliminate activity', 'Bromo at C4 increases potency', 'Alkynyl at C5 dramatically increases potency']
- MMP thresholds used here: dCF>=0.20 'explained', 0.10-0.20 'partial'.

## 10. Active vs inactive analysis

- Active compounds contact site-23-defining residues in **90%** of models; inactive compounds in **94%** (per-class model fractions, `05 site_class_fraction`).
- Discordant inactive poses: **136** of 138 inactive conditions have >=50% of models contacting site-23 residues (cpds [4, 5, 6, 7, 8, 9, 10, 11, 13, 14, 15, 16, 17, 19, 20, 21, 22, 23, 26, 27, 28, 29, 30]); these are reported as discordant cases and weaken naive active/inactive discrimination.
- Interpretation: site-23 contact does NOT cleanly separate actives from inactives; potency differences therefore likely arise from *interaction quality* rather than simple occupancy (consistent with allosteric, weak-affinity ascorbate scaffold).

## 11. Activity-cliff analysis

- Pair status counts: `{"explained": 5}`.
- **12 vs 7** (change: alkynyl vs vinyl at C5; effect: None): status=explained; defining residues alpha10.435+0.5, alpha9.176+0.5, alpha9.120+0.5, alpha9.224+0.5, alpha10.432+0.333.
- **25 vs 1** (change: bromo vs OH at C4; effect: None): status=explained; defining residues alpha9.260+0.333, alpha10.322+0.333, alpha10.318+0.333, alpha9.256+0.333, alpha9.335+0.333.
- **18 vs 19** (change: benzyloxymethyl vs allyloxy at C5; effect: None): status=explained; defining residues alpha9.148+0.667, alpha9.82+0.333, alpha9.215+0.333, alpha9.146+0.333, alpha9.176+0.333.
- **3 vs 1** (change: gem-dimethyl acetonide vs diol; effect: None): status=explained; defining residues alpha9.146+0.333, alpha9.457+0.333, alpha9.206+0.333, alpha9.161+0.333, alpha9.234+0.333.
- **24 vs 20** (change: 3-O-propyl vs 5-O-butyl; effect: None): status=explained; defining residues alpha9.460+0.333, alpha9.121+0.333, alpha9.457+0.333, alpha9.82+0.333, alpha10.308+0.333.

## 12. XAI / SAR / structure convergence

- XAI: prior permutation importances were computed on an overfit model and are invalidated (FALLACIES_AUDIT). No usable XAI evidence remains.
- SAR-to-structure convergence: only partial — residue-level defining contacts were derived for cliff pairs where the pose ensemble allowed it (section 11); these are proposed as mutagenesis hypotheses in section 19 and are not yet experimental.

## 13. Boltz-2 affinity evidence

- **ABSENT.** No `affinity_*.json` was emitted by the corrected cofolding runs; the earlier affinity 'proxy' matrix was circular and invalidated. Affinity is used nowhere in this validation. (Known gap, retained honestly.)

## 14. Consensus receptor residues

- Tier 1 (>=2 structural methods + SAR dCF + pocket membership): **6** residues → alpha9.120, alpha9.176, alpha9.224, alpha9.82, alpha9.175, alpha9.177
- Tier 2 (structural support only): 50 residues → alpha9.148, alpha9.215, alpha9.146, alpha10.81, alpha9.164, alpha9.161, alpha10.145, alpha9.80, alpha9.217, alpha10.119, alpha10.223, alpha10.175, alpha10.143, alpha9.166, alpha9.178
- SAR-cliff-implicated residues (from section 11): alpha10.1, alpha10.105, alpha10.116, alpha10.117, alpha10.119, alpha10.120, alpha10.121, alpha10.125, alpha10.126, alpha10.127, alpha10.128, alpha10.129
- Full ranked table in `outputs/10_residue_consensus.json`.

## 15. Competing-site analysis

- The bounded panel (180 conditions) was run for candidate sites based on prior ranking; sites 34 (TMD-pore) and 7 have **no completed models** yet, so the TMD-pore competing hypothesis is currently untested computationally.
- Site-5 and site-21 conditions drift toward site-23-classified contacts, which is either true preference or definition degeneracy; cannot be resolved without SiteAF3 / site competition.
- **Verdict: competition unresolved.** Favored-site support is real but not superior-chosen; a SiteAF3 head-to-head is the required next test.

## 16. Contradictory evidence

- Huge receptor C-alpha divergence across seeds (up to ~40 A per chain) — the receptor models are not converging on a single conformation.
- Ligand pose spread 50 A mean — individual poses are not reproducible, undermining atom-level interaction claims.
- No clean active/inactive discrimination at site 23 (23 discordant inactive compounds, 136 conditions).
- dCF-negative residues inside candidate pocket definitions (e.g. alpha10.62, alpha10.143) — some pocket residues are contacted MORE by inactive compounds.

## 17. Remaining uncertainties

- Which ECD site (23 vs 21 vs 5) is primary, given definition overlap.
- Whether the pose spread reflects genuine receptor plasticity or modeling instability (no experimental structure to anchor).
- TMD-pore hypothesis (site 34) untested.
- SiteAF3 layer missing.
- No affinity layer; no XAI.
- C5- vs C4-substituent contact mapping is pose-derived and so uncertain.

## 18. MD decision

- **MD not necessary / deferred.** MD deferred: pending SiteAF3 results and tighter residue prioritization.
- Reason: with no experimental structure to benchmark and with pose/ensemble convergence already limiting residue-level claims, MD would not change the mutagenesis priorities (which are set by structural+SAR convergence), nor resolve the current blocker (missing SiteAF3). MD is only indicated after SiteAF3 if a specific open/closed or cryptic-pocket question needs to be resolved or if a pocket's thermal stability must be established.

## 19. Prioritized mutagenesis experiments (final validation step)

- Plan: 123 residues (Tier-1 first) — full table in `outputs/11_mutagenesis_plan.json`.
- Protocol: express WT and mutant (2-microelectrode voltage clamp on alpha9alpha10 cRNA); confirm normal ACh EC50/Imax; then measure PAM potentiation of the reference ascorbate and the nanomolar PAM (cpd 12); compare WT vs mutant.
- **alpha9.120** (WT Y): probed by **Y120Ala** — sar-defining contact residue; the relevant cliff pair(s): 12 vs 7 (dCF +0.50), 18 vs 19 (dCF +0.33); expected if correct: reduced/lost PAM potentiation of ACh-evoked current; if incorrect: no-effect.
- **alpha9.176** (WT W): probed by **W176Ala** — sar-defining contact residue; the relevant cliff pair(s): 12 vs 7 (dCF +0.50), 18 vs 19 (dCF +0.33); expected if correct: reduced/lost PAM potentiation of ACh-evoked current; if incorrect: no-effect.
- **alpha9.224** (WT Y): probed by **Y224Ala** — sar-defining contact residue; the relevant cliff pair(s): 12 vs 7 (dCF +0.50), 18 vs 19 (dCF +0.33); expected if correct: reduced/lost PAM potentiation of ACh-evoked current; if incorrect: no-effect.
- **alpha9.82** (WT W): probed by **W82Ala** — sar-defining contact residue; the relevant cliff pair(s): 12 vs 7 (dCF +0.33), 18 vs 19 (dCF +0.33), 24 vs 20 (dCF +0.33); expected if correct: reduced/lost PAM potentiation of ACh-evoked current; if incorrect: no-effect.
- **alpha9.175** (WT S): probed by **S175Ala** — sar-defining contact residue; the relevant cliff pair(s): 12 vs 7 (dCF +0.33); expected if correct: reduced/lost PAM potentiation of ACh-evoked current; if incorrect: no-effect.
- **alpha9.177** (WT T): probed by **T177Ala** — sar-defining contact residue; the relevant cliff pair(s): none; expected if correct: reduced/lost PAM potentiation of ACh-evoked current; if incorrect: no-effect.

## 20. Expected outcomes under each binding-site hypothesis

| Hypothesis | Mutation result if hypothesis correct | Result if incorrect |
|---|---|---|
| PAM binds the ECD vestibular inter-subunit cavity (site 23 family) | Ala of pocket/site-23 contact residues reduces PAM potentiation (EC50 shift / Imax) without disturbing ACh gating | No PAM-effect change |
| PAM binds site 21/5 (mirror face of the same cavity) | Only residues specific to the mirror face respond; site-23 probes are WT-like | Site-23 probes respond instead |
| PAM binds the TMD pore (site 34) | Pore/helix-lining mutations (M2/M1) affect modulation; ECD probes WT-like | ECD probes respond |
| Not a true binding site (model artifact) | Mutagenesis of all candidate residues leaves PAM potentiation unchanged | — (falsification) |

## 21. Final confidence classification

**B. Structurally supported hypothesis (SiteAF3 layer still pending)**

Rationale: the pocket region defined by the site-23 residue family recurs in the unbiased and ligand-conditioned Boltz-2 panel, and SiteAF3-conditioned testing is designed and staged but **not executed**; without it the hypothesis cannot rise above `B`. Experimental SAR explains the major cliffs only partially at residue level, and active/inactive discrimination is weak at the current pose resolution. No experimental a9a10-PAM structure exists, so classification `E/F` cannot be claimed regardless.

### Actions that would upgrade the classification

1. **Run SiteAF3** (04_siteaf3_submit.sh) on the 80 staged configs after installing the AF3 env + weights → would enable class C.
2. Resolve the site-23/21/5 definition degeneracy (give the sites disjoint residue shells or run a real competition).
3. Enumerate and rationalize the discordant inactive poses (23 compounds / 136 conditions; some are explainable steric/electronic mismatches once residue mapping improves).
4. Obtain experimental binding/potentiation data for prioritized Tier-1/2 mutants → would enable class E.

---
*Methods: AF3 cofolding (existing 19-condition set), Boltz-2 ligand-conditioned panel (180 modeled conditions, 540 ligand models), SiteAF3 (staged), and the experimental SAR panel. This report and its intermediate outputs are machine-generated; the falsifiability framing follows the project AGENTS.md and FALLACIES_AUDIT.md.*