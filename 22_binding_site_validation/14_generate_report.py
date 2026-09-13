#!/usr/bin/env python3
"""14. Final Binding-Site Validation Report (prompt sections 1-21).

Aggregates every module output into the full report and assigns the
prompt-section-20 confidence classification:
  A weak / B structurally supported / C computationally validated /
  D SAR-supported / E experimentally supported / F experimentally validated.

Outputs: outputs/14_validation_report.md (+ .json)
"""

import json
from datetime import datetime
from pathlib import Path

from config import RESULTS


def load(o):
    p = RESULTS / o
    if p.exists():
        with open(p) as f:
            return json.load(f)
    return {}


def fmt_residue_list(rows, key="residue", limit=12):
    return ", ".join(str(r[key]) for r in rows[:limit])


def main():
    audit = load("01_evidence_audit.json")
    ens = load("02_ensemble_convergence.json")
    b2 = load("05_boltz2_consensus.json")
    conv = load("06_three_method_convergence.json")
    qc = load("07_medchem_qc.json")
    plast = load("08_pocket_plasticity.json")
    sar = load("09_sar_activity_cliff_validation.json")
    consensus = load("10_residue_consensus.json")
    mplan = load("11_mutagenesis_plan.json")
    md = load("12_md_decision.json")
    matrix = load("13_evidence_matrix.json")

    lc = ens.get("local_convergence") or {}
    rc = ens.get("receptor_convergence") or {}
    pc = ens.get("pose_convergence") or {}
    site23_rec = (lc.get("23") or {}).get("self_recurrence")
    site21_rec = (lc.get("21") or {}).get("self_recurrence")
    site5_rec = (lc.get("5") or {}).get("self_recurrence")
    disc = b2.get("discordant_summary", {})
    site_frac = b2.get("site_class_fraction", {})
    active_sf = (site_frac.get("active") or {}).get("23")
    inactive_sf = (site_frac.get("inactive") or {}).get("23")
    tier1 = consensus.get("tier_1", [])
    tier2 = consensus.get("tier_2", [])
    cliff_statuses = {p.get("status"): 0 for p in sar.get("pairs", [])}
    for p in sar.get("pairs", []):
        cliff_statuses[p.get("status")] = cliff_statuses.get(p.get("status"), 0) + 1

    # ---- confidence classification (section 20) ----
    siteaf3_done = not (conv.get("siteaf3_status") or "").startswith("NOT RUN")
    n_explained = cliff_statuses.get("explained", 0)
    if not siteaf3_done and not tier1 and n_explained == 0:
        level = "A. Weak computational hypothesis (no Tier-1 residue, no SAR cliff, SiteAF3 pending)"
    elif not siteaf3_done:
        level = "B. Structurally supported hypothesis (SiteAF3 layer still pending)"
    elif n_explained >= 1:
        level = "D. SAR-supported binding-site hypothesis (computational + SAR convergence)"
    else:
        level = "C. Computationally validated hypothesis (AF3/Boltz-2/SiteAF3 convergence)"
    classification = level

    L = []
    A = L.append

    # ============ header ============
    A("# Alpha9alpha10 nAChR PAM — Binding-Site Validation Report")
    A("")
    A(f"*Generated:* {datetime.now().strftime('%Y-%m-%d %H:%M')}  ")
    A("*Falsifiable-hypothesis treatment: the proposed site is treated as a hypothesis to be* "
      "*confirmed or falsified, NOT as a correct answer.*")
    A("")

    # ============ 1. evidence audit ============
    A("## 1. Existing evidence audit")
    A("")
    A("Layers inventoried (full table in `outputs/01_evidence_audit.json`):")
    for r in audit.get("layers", []):
        corr = " [correlated]" if r.get("correlated") else ""
        caveat = f" — {r['caveat']}" if r.get("caveat") else ""
        A(f"- **{r['layer']}** (n={r.get('n', 0)}): {r.get('evidence_type', '')}. "
          f"{r.get('supports', '')}{caveat}{corr}")
    A("")
    A("Independence notes: AF3, Boltz-2 and SiteAF3 all derive their MSAs from the "
      "same sequence databases; multiple seeds of one method are NOT independent "
      "observations. The old Boltz-2 affinity proxy was circular and is invalidated "
      "(FALLACIES_AUDIT).")
    A("")

    # ============ 2. candidate sites ============
    A("## 2. Candidate binding sites")
    A("")
    A("Retained and treated as competing hypotheses:")
    for s, info in [("23", "alpha9(+)/alpha10(-) ECD vestibular inter-subunit (favored by prior work)"),
                    ("21", "alpha10(+)/alpha9(-) ECD vestibular inter-subunit"),
                    ("5", "alpha9(+)/alpha10(-) ECD vestibular inter-subunit (overlaps 23)"),
                    ("34", "TMD-pore allosteric cavity (decode/wall region)"),
                    ("7", "alpha9(+)/alpha10(-) ECD vestibular inter-subunit")]:
        A(f"- **Site {s}**: {info}")
    A("")
    A("Site definitions overlap heavily in the ECD vestibule (same inter-subunit "
      "cavity sampled from different chain faces) — recurrence statistics must be read "
      "with that degeneracy in mind.")
    A("")

    # ============ 3. AF3 / Boltz-2 convergence ============
    A("## 3. AF3 / Boltz-2 ensemble convergence")
    A("")
    rc_mean = rc.get("chain_ca_rmsd_mean")
    A(f"- Receptor C-alpha RMSD across Boltz-2 seeds of the same condition: "
      f"**{rc_mean} A** (n={rc.get('n_conditions')} conditions, 3 seeds each) → "
      f"**poor global receptor convergence**; each seed remodels the pentamer substantially.")
    A(f"- Ligand pose spread (mean pairwise centroid distance across seeds): "
      f"**{pc.get('ligand_centroid_spread_mean_ang')} A** (p25 {pc.get('spread_percentiles', {}).get('p25')}, "
      f"median {pc.get('spread_percentiles', {}).get('median')}, p75 {pc.get('spread_percentiles', {}).get('p75')}).")
    pct23 = f"{site23_rec*100:.1f}%" if site23_rec is not None else "?"
    pct21 = f"{site21_rec*100:.1f}%" if site21_rec is not None else "?"
    pct5 = f"{site5_rec*100:.1f}%" if site5_rec is not None else "?"
    A(f"- Site self-recurrence under conditioning: site 23 = **{pct23}**, "
      f"site 21 = **{pct21}**, site 5 = **{pct5}**.")
    A(f"- **Cross-condition drift (informative):** models conditioned on other sites "
      f"were nonetheless classified into site-23 residues in most cases "
      f"(`02 hit_distribution`), suggesting the vestibular vestibular sub-pocket sampled "
      f"by site-23 residues is the recurrent preference; this is confounded by the "
      f"overlap of the site definitions.")
    A(f"- Existing AF3 cofolding (`data/allostery/af3_outputs`, root model cif) yields "
      f"{conv.get('af3_n_contacts', '?')} residue contacts; existing Boltz-2 outputs "
      f"yield no site-classified hits at these definitions.")
    A("")
    A("Convergence verdict: **local (pocket-region) recurrence is high but pose-level "
      "and global-receptor convergence are poor.** The pocket recurs, but individual "
      "poses across seeds are not reproducible.")
    A("")

    # ============ 4. SiteAF3 site-conditioned analysis ============
    A("## 4. SiteAF3 site-conditioned analysis")
    A("")
    sa = conv.get("siteaf3_status", "")
    if "NOT RUN" in str(sa):
        A("- **PENDING / NOT RUN.** The site-conditioned test (Tang & Wang, PNAS "
          "10.1073/pnas.2521048122) could not be executed in this session: the AF3-based "
          "SiteAF3 environment and model weights are not installed on the cluster. "
          "80 config inputs are staged in `siteaf3_inputs/` (5 sites x 2 stoichiometries "
          "x 8 representative compounds). These are the missing-validity link in this report; "
          "**until SiteAF3 runs, the hypothesis cannot be classed above B.**")
        A("- Method design (already implemented in `04_siteaf3_submit.sh`): pocket-masked "
          "MSA embedding + pocket-conditioned diffusion, one config per candidate site, "
          "so competing sites are tested under identical conditions.")
    else:
        A(f"- {sa}")
    A("")

    # ============ 5. SiteAF3 vs competing sites ============
    A("## 5. SiteAF3 vs competing-site comparison")
    A("")
    A("- Competition outcome: **pending** (SiteAF3 not run). Until then the "
      "favored-site-vs-alternatives comparison rests on the Boltz-2 ligand-conditioned "
      "panel alone (sections 3/10).")

    # ============ 6. pocket & pose QC ============
    A("## 6. Pocket and pose QC")
    A("")
    A(f"- Models analyzed: **{qc.get('n_models_analyzed')}**; clash-free fraction "
      f"**{qc.get('fraction_clash_free')}**; mean clashing ligand atoms/model "
      f"**{qc.get('mean_clashing_ligand_atoms_per_model')}** (cutoff "
      f"{qc.get('clash_cutoff_A')} A).")
    A(f"- Pocket occupancy (ligand atoms within 6 A of any CA): "
      f"**{qc.get('mean_pocket_occupancy_fraction')}** — adequate, but as a gross filter.")
    A(f"- Pharmacophore centroid spread (cpd 12): "
      f"**{qc.get('pharmacophore_centroid_spread_A_per_cpds', {}).get('12')} A** — "
      f"large; pose convergence is NOT adequate for atom-level claims.")
    A("")

    # ============ 7. Chen et al. med-chem QC ============
    A("## 7. Medicinal-chemistry QC (Chen et al. framework)")
    A("")
    rd = qc.get("rdkit", {})
    A(f"- RDKit: {rd.get('n_valid')}/{rd.get('n_smiles')} SMILES valid; "
      f"{rd.get('n_with_stereo')} carry defined stereocenters (ascorbate series: "
      f"must keep stereochemistry, PAM activity is stereo-sensitive).")
    A(f"- Whole-ligand RMSD / pharmacophore-level agreement are reported as "
      f"*internal consistency only*; no experimental a9a10-PAM complex exists, so "
      f"experimental pose accuracy is not claimed.")
    A("")

    # ============ 8. allosteric plasticity ============
    A("## 8. Allosteric pocket plasticity")
    A("")
    A(f"- Mean per-residue CA span across bound-model ensemble (aligned): "
      f"**{plast.get('mean_pocket_ca_std_A')} A**; {plast.get('n_residues_gt_2A_span')} "
      f"residues >2 A, {plast.get('n_residues_gt_4A_span')} >4 A.")
    A(f"- Classification: **flexible / transient-like**, precision limited by the "
      f"large receptor ensemble spread; no apo→holo transition was modeled, so no "
      f"induced-fit claim. Classify cautiously as flexible/interface-dependent, "
      f"cryptic-pocket behaviour not addressable from these data (see section 18).")
    A("")

    # ============ 9. experimental SAR validation ============
    A("## 9. Experimental SAR validation")
    A("")
    A(f"- Dataset: {sar.get('n_sar')} compounds; {len(sar.get('active_ids', []))} active "
      f"(IDs {sar.get('active_ids')}); MMP records {sar.get('n_mmp_records')}; "
      f"curated cliff pairs {len(sar.get('pairs', []))}.")
    A(f"- MMP rules from prior analysis: {sar.get('mmp_rules', [])}")
    A(f"- MMP thresholds used here: dCF>=0.20 'explained', 0.10-0.20 'partial'.")
    A("")

    # ============ 10. active/inactive analysis ============
    A("## 10. Active vs inactive analysis")
    A("")
    A(f"- Active compounds contact site-23-defining residues in "
      f"**{active_sf*100:.0f}%** of models; inactive compounds in "
      f"**{inactive_sf*100:.0f}%** (per-class model fractions, `05 site_class_fraction`).")
    A(f"- Discordant inactive poses: **{disc.get('n_discordant_conditions', 0)}** of "
      f"{disc.get('n_inactive_conditions', '?')} inactive conditions have >=50% of models "
      f"contacting site-23 residues (cpds {disc.get('discordant_cpds')}); these are "
      f"reported as discordant cases and weaken naive active/inactive discrimination.")
    A(f"- Interpretation: site-23 contact does NOT cleanly separate actives from "
      f"inactives; potency differences therefore likely arise from *interaction quality* "
      f"rather than simple occupancy (consistent with allosteric, weak-affinity "
      f"ascorbate scaffold).")
    A("")

    # ============ 11. activity cliffs ============
    A("## 11. Activity-cliff analysis")
    A("")
    A(f"- Pair status counts: `{json.dumps(cliff_statuses)}`.")
    for p in sar.get("pairs", []):
        st = p.get("status")
        resid = p.get("defining_residues_sorted", [])
        resid_s = ", ".join(f"{r}+{v}" for r, v in resid[:5]) if resid else "none called"
        A(f"- **{p['pair']}** (change: {p.get('change')}; effect: {p.get('effect')}): "
          f"status={st}; defining residues {resid_s}.")
    A("")

    # ============ 12. XAI/SAR/structure convergence ============
    A("## 12. XAI / SAR / structure convergence")
    A("")
    A("- XAI: prior permutation importances were computed on an overfit model and are "
      "invalidated (FALLACIES_AUDIT). No usable XAI evidence remains.")
    A("- SAR-to-structure convergence: only partial — residue-level defining contacts "
      "were derived for cliff pairs where the pose ensemble allowed it (section 11); "
      "these are proposed as mutagenesis hypotheses in section 19 and are not yet "
      "experimental.")
    A("")

    # ============ 13. Boltz-2 affinity ============
    A("## 13. Boltz-2 affinity evidence")
    A("")
    A("- **ABSENT.** No `affinity_*.json` was emitted by the corrected cofolding runs; "
      "the earlier affinity 'proxy' matrix was circular and invalidated. Affinity is "
      "used nowhere in this validation. (Known gap, retained honestly.)")
    A("")

    # ============ 14. consensus receptor residues ============
    A("## 14. Consensus receptor residues")
    A("")
    A(f"- Tier 1 (>=2 structural methods + SAR dCF + pocket membership): "
      f"**{len(tier1)}** residues → {fmt_residue_list([{'residue': r} for r in tier1], limit=15)}")
    A(f"- Tier 2 (structural support only): {len(tier2)} residues → "
      f"{fmt_residue_list([{'residue': r} for r in tier2], limit=15)}")
    A(f"- SAR-cliff-implicated residues (from section 11): "
      f"{', '.join(consensus.get('cliff_residues', [])[:12]) or 'none'}")
    A("- Full ranked table in `outputs/10_residue_consensus.json`.")
    A("")

    # ============ 15. competing site analysis ============
    A("## 15. Competing-site analysis")
    A("")
    A("- The bounded panel (180 conditions) was run for candidate sites based on prior "
      "ranking; sites 34 (TMD-pore) and 7 have **no completed models** yet, so the "
      "TMD-pore competing hypothesis is currently untested computationally.")
    A("- Site-5 and site-21 conditions drift toward site-23-classified contacts, which "
      "is either true preference or definition degeneracy; cannot be resolved without "
      "SiteAF3 / site competition.")
    A("- **Verdict: competition unresolved.** Favored-site support is real but not "
      "superior-chosen; a SiteAF3 head-to-head is the required next test.")
    A("")

    # ============ 16. contradictory evidence ============
    A("## 16. Contradictory evidence")
    A("")
    A("- Huge receptor C-alpha divergence across seeds (up to ~40 A per chain) — the "
      "receptor models are not converging on a single conformation.")
    A("- Ligand pose spread 50 A mean — individual poses are not reproducible, "
      "undermining atom-level interaction claims.")
    n_disc_cond = disc.get("n_discordant_conditions", 0)
    n_disc_cpd = len(disc.get("discordant_cpds", []))
    A("- No clean active/inactive discrimination at site 23 "
      f"({n_disc_cpd} discordant inactive compounds, {n_disc_cond} conditions).")
    A("- dCF-negative residues inside candidate pocket definitions (e.g. alpha10.62, "
      "alpha10.143) — some pocket residues are contacted MORE by inactive compounds.")
    A("")

    # ============ 17. remaining uncertainties ============
    A("## 17. Remaining uncertainties")
    A("")
    A("- Which ECD site (23 vs 21 vs 5) is primary, given definition overlap.")
    A("- Whether the pose spread reflects genuine receptor plasticity or modeling "
      "instability (no experimental structure to anchor).")
    A("- TMD-pore hypothesis (site 34) untested.")
    A("- SiteAF3 layer missing.")
    A("- No affinity layer; no XAI.")
    A("- C5- vs C4-substituent contact mapping is pose-derived and so uncertain.")
    A("")

    # ============ 18. MD decision ============
    A("## 18. MD decision")
    A("")
    A(f"- **{'MD necessary' if md.get('necessary') else 'MD not necessary / deferred'}.** {md.get('decision')}")
    A("- Reason: with no experimental structure to benchmark and with pose/ensemble "
      "convergence already limiting residue-level claims, MD would not change the "
      "mutagenesis priorities (which are set by structural+SAR convergence), nor "
      "resolve the current blocker (missing SiteAF3). MD is only indicated after SiteAF3 "
      "if a specific open/closed or cryptic-pocket question needs to be resolved or if "
      "a pocket's thermal stability must be established.")
    A("")

    # ============ 19. prioritized mutagenesis ============
    A("## 19. Prioritized mutagenesis experiments (final validation step)")
    A("")
    A(f"- Plan: {len(mplan.get('plan', []))} residues (Tier-1 first) — full table in "
      f"`outputs/11_mutagenesis_plan.json`.")
    A("- Protocol: express WT and mutant (2-microelectrode voltage clamp on "
      "alpha9alpha10 cRNA); confirm normal ACh EC50/Imax; then measure PAM "
      "potentiation of the reference ascorbate and the nanomolar PAM (cpd 12); "
      "compare WT vs mutant.")
    for p in [x for x in mplan.get("plan", []) if x.get("tier") == 1][:12]:
        A(f"- **{p['residue']}** (WT {p.get('wt_aa')}): probed by **{p.get('probe_mutation')}** — "
          f"{p.get('predicted_role', '').lower()}; the relevant cliff pair(s): "
          f"{', '.join(p.get('cliff_pairs', [])) or 'none'}; expected if correct: "
          f"{p.get('expected_if_site_correct', '').split(' with ')[0]}; if incorrect: no-effect.")
    if not tier1:
        A("- *No Tier-1 residue passed the honest threshold yet; the plan below is "
          "Tier-2 only and lower priority.*")
    A("")

    # ============ 20. expected outcomes under each hypothesis ============
    A("## 20. Expected outcomes under each binding-site hypothesis")
    A("")
    A("| Hypothesis | Mutation result if hypothesis correct | Result if incorrect |")
    A("|---|---|---|")
    A("| PAM binds the ECD vestibular inter-subunit cavity (site 23 family) | Ala of "
      "pocket/site-23 contact residues reduces PAM potentiation (EC50 shift / Imax) "
      "without disturbing ACh gating | No PAM-effect change |")
    A("| PAM binds site 21/5 (mirror face of the same cavity) | Only residues specific "
      "to the mirror face respond; site-23 probes are WT-like | Site-23 probes respond "
      "instead |")
    A("| PAM binds the TMD pore (site 34) | Pore/helix-lining mutations (M2/M1) affect "
      "modulation; ECD probes WT-like | ECD probes respond |")
    A("| Not a true binding site (model artifact) | Mutagenesis of all candidate "
      "residues leaves PAM potentiation unchanged | — (falsification) |")
    A("")

    # ============ 21. final classification ============
    A("## 21. Final confidence classification")
    A("")
    A(f"**{classification}**")
    A("")
    A("Rationale: the pocket region defined by the site-23 residue family recurs in the "
      "unbiased and ligand-conditioned Boltz-2 panel, and SiteAF3-conditioned testing is "
      "designed and staged but **not executed**; without it the hypothesis cannot rise "
      "above `B`. Experimental SAR explains the major cliffs only partially at residue "
      "level, and active/inactive discrimination is weak at the current pose resolution. "
      "No experimental a9a10-PAM structure exists, so classification `E/F` cannot be "
      "claimed regardless.")
    A("")
    A("### Actions that would upgrade the classification")
    A("")
    A("1. **Run SiteAF3** (04_siteaf3_submit.sh) on the 80 staged configs after "
      "installing the AF3 env + weights → would enable class C.")
    A("2. Resolve the site-23/21/5 definition degeneracy (give the sites disjoint "
      "residue shells or run a real competition).")
    A("3. Enumerate and rationalize the discordant inactive poses "
      f"({n_disc_cpd} compounds / {n_disc_cond} conditions; some are explainable "
      "steric/electronic mismatches once residue mapping improves).")
    A("4. Obtain experimental binding/potentiation data for prioritized Tier-1/2 "
      "mutants → would enable class E.")
    A("")
    A("---")
    A(f"*Methods: AF3 cofolding (existing 19-condition set), Boltz-2 ligand-conditioned "
      f"panel ({ens.get('full_panel', {}).get('n_conditions', '?')} modeled conditions, "
      f"{b2.get('n_valid_ligand_models', '?')} ligand models), SiteAF3 (staged), and the "
      f"experimental SAR panel. This report and its intermediate outputs are "
      f"machine-generated; the falsifiability framing follows the project AGENTS.md and "
      f"FALLACIES_AUDIT.md.*")

    md_path = RESULTS / "14_validation_report.md"
    md_path.write_text("\n".join(L))
    with open(RESULTS / "14_validation_report.json", "w") as f:
        json.dump({"classification": classification, "generated": str(datetime.now()),
                   "siteaf3_run": siteaf3_done, "n_tier1": len(tier1),
                   "n_explained_cliffs": n_explained}, f, indent=2)

    print("\n".join(L[:40]))
    print(f"\nwrote {md_path}")


if __name__ == "__main__":
    main()