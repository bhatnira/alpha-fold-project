#!/usr/bin/env python3
"""13. Final evidence matrix (falsifiability tally).

Assembles the prompt-section-19 matrix from the outputs of scripts 01, 02,
05-10. Each layer is classified strictly:

  supports   = positive, reasonably-assessed signal for the proposed site
  contradicts = definitive negative signal for the proposed site
  neutral     = pending, absent, not-directly-informative, or too weak

Notes distinguish correlated evidence (multiple seeds of one method) from
genuinely orthogonal evidence (different methods / experiment).

Output: outputs/13_evidence_matrix.json
"""

import json
from pathlib import Path

from config import RESULTS, POCKET_RESIDUES


def load(o):
    p = RESULTS / o
    if p.exists():
        with open(p) as f:
            return json.load(f)
    return {}


def main():
    audit = load("01_evidence_audit.json")
    ens = load("02_ensemble_convergence.json")
    b2 = load("05_boltz2_consensus.json")
    conv = load("06_three_method_convergence.json")
    qc = load("07_medchem_qc.json")
    plast = load("08_pocket_plasticity.json")
    sar = load("09_sar_activity_cliff_validation.json")
    consensus = load("10_residue_consensus.json")
    md = load("12_md_decision.json")

    lc = ens.get("local_convergence") or {}
    site23_rec = (lc.get("23") or {}).get("self_recurrence")
    rec_ok = (site23_rec or 0) > 0.5

    # AF3 -> residues contacting + mapping onto site-23 pocket definition
    n_af3_res = 0
    af3_site23 = 0
    for row in (conv.get("residue_agreement") or []):
        if row.get("af3_n_contacts", 0) > 0:
            n_af3_res += 1
            if row.get("site_candidate") == 23:
                af3_site23 += 1

    # active / inactive site-23 contact fractions, discordant count
    site_frac = b2.get("site_class_fraction", {})
    active_site23 = (site_frac.get("active", {}) or {}).get("23")
    inactive_site23 = (site_frac.get("inactive", {}) or {}).get("23")
    disc = b2.get("discordant_summary", {})
    n_disc = disc.get("n_discordant_conditions", 0)
    n_inact = max(1, disc.get("n_inactive_conditions", 0) or 1)
    disc_frac = n_disc / n_inact

    active_coh = bool(active_site23 and active_site23 > 0.5)
    inactive_disc = bool(inactive_site23 and inactive_site23 > 0.5 and disc_frac > 0.5)
    tier1 = consensus.get("tier_1", [])

    cliff_status = [p.get("status") for p in sar.get("pairs", [])]
    n_explained = cliff_status.count("explained")
    n_partial = cliff_status.count("partial")

    pharm_spread_med = None
    sp12 = ((qc.get("pharmacophore_centroid_spread_A_per_cpds") or {}).get("12"))
    pharm_spread_med = sp12

    layers = [
        {"layer": "AF3 cofolding",
         "result": f"{n_af3_res} contacting residues; {af3_site23} map to site-23 pocket",
         "supports_site": af3_site23 >= 3, "strength": "low-medium",
         "independent": False,
         "note": "Unbiased cofolding; weak direct support for site-23 residues."},
        {"layer": "Boltz-2 cofolding",
         "result": f"site 23 recurrence {site23_rec*100:.1f}% of modeled conditions",
         "supports_site": rec_ok, "strength": "medium",
         "independent": False,
         "note": "Ligand-conditioned; 92% self-recurrence but site definitions "
                "overlap in the ECD vestibule -> partly circular (see 02 hit_distribution)."},
        {"layer": "Boltz-2 affinity",
         "result": "ABSENT (0 affinity_*.json)",
         "supports_site": None, "strength": "none",
         "independent": True,
         "note": "Affinity layer never produced; prior proxy invalidated (FALLACIES_AUDIT)."},
        {"layer": "SiteAF3",
         "result": "NOT RUN (results/siteaf3 empty)",
         "supports_site": None, "strength": "pending",
         "independent": True,
         "note": "80 configs prepared; run gated on AF3 env/weights + sdcofold completion."},
        {"layer": "SiteAF3 vs AF3/Boltz convergence",
         "result": "pending",
         "supports_site": None, "strength": "pending",
         "independent": True, "note": ""},
        {"layer": "Stoichiometry comparison",
         "result": "both 2to3 and 3to2 modeled; site-23 recurrence stable across both",
         "supports_site": rec_ok, "strength": "medium",
         "independent": False,
         "note": "Site persists in both stoichiometries (see 02 stoichiometry_convergence)."},
        {"layer": "Pocket quality",
         "result": f"pocket occupancy {qc.get('mean_pocket_occupancy_fraction')}",
         "supports_site": bool(qc.get("mean_pocket_occupancy_fraction") and
                              qc["mean_pocket_occupancy_fraction"] > 0.7),
         "strength": "medium", "independent": False,
         "note": "Occupancy = ligand atoms within 6 A of any CA; gross filter only."},
        {"layer": "Pose QC",
         "result": f"clash-free {qc.get('fraction_clash_free')}",
         "supports_site": bool((qc.get("fraction_clash_free") or 0) >= 0.8),
         "strength": "medium", "independent": False,
         "note": "Intra-model quality only; says nothing about the site itself."},
        {"layer": "Pharmacophore QC",
         "result": f"cpd-12 pharmacophore centroid spread {pharm_spread_med} A",
         "supports_site": bool(pharm_spread_med is not None and pharm_spread_med <= 10.0),
         "strength": "low", "independent": False,
         "note": "Large spread => weak pose convergence; limits residue-level claims."},
        {"layer": "Allosteric plasticity",
         "result": f"mean CA span {plast.get('mean_pocket_ca_std_A')} A",
         "supports_site": None, "strength": "low", "independent": False,
         "note": "High ensemble variance; classification flexible/transient-like, "
                "precision limited."},
        {"layer": "Experimental SAR",
         "result": f"{sar.get('n_sar', 30)} compounds; "
                   f"{len(sar.get('active_ids', []))} active",
         "supports_site": None, "strength": "high", "independent": True,
         "note": "Ground-truth challenge layer; supports the site only via the "
                "SAR-to-structure mapping below."},
        {"layer": "Active-compound explanation",
         "result": (f"active compounds contact site-23 residues in "
                    f"{active_site23*100:.0f}% of models" if active_site23 is not None else "n/a"),
         "supports_site": active_coh, "strength": "medium",
         "independent": False,
         "note": "Actives do form site-23-classified poses."},
        {"layer": "Inactive-compound challenge",
         "result": (f"{n_disc}/{n_inact} inactive-model conditions are discordant "
                    f"(>=50% of models contact site-23)"),
         "supports_site": not inactive_disc, "strength": "low",
         "independent": False,
         "note": "Discordant inactive poses weaken active/inactive discrimination "
                "unless the pocket can explain them."},
        {"layer": "Activity cliffs",
         "result": f"{n_explained} explained / {n_partial} partial / "
                   f"{len(cliff_status)} total pairs",
         "supports_site": n_explained >= 1, "strength": "medium",
         "independent": True,
         "note": "See outputs/09 for defining residues per cliff pair."},
        {"layer": "XAI",
         "result": "existing XAI invalidated (FALLACIES_AUDIT)",
         "supports_site": None, "strength": "none", "independent": True,
         "note": "Not usable as evidence."},
        {"layer": "Residue convergence",
         "result": f"Tier1={len(tier1)}, Tier2={len(consensus.get('tier_2', []))}",
         "supports_site": bool(tier1), "strength": "medium",
         "independent": False,
         "note": "Tier-1 residues now require 2 structural methods + SAR dCF + "
                "pocket membership (see outputs/10)."},
        {"layer": "Competing-site analysis",
         "result": "only the favored-site panel is complete (sites 34/7 not modeled)",
         "supports_site": None, "strength": "pending",
         "independent": True,
         "note": "Head-to-head site competition (SiteAF3) pending."},
    ]

    n_support = sum(1 for L in layers if L["supports_site"] is True)
    n_contra = sum(1 for L in layers if L["supports_site"] is False)
    n_pend = len(layers) - n_support - n_contra

    with open(RESULTS / "13_evidence_matrix.json", "w") as f:
        json.dump({"layers": layers, "n_support": n_support,
                   "n_contradict": n_contra, "n_neutral_pending": n_pend,
                   "md_decision": md,
                   "note": ("Layers are NOT summed as independent observations; "
                            "AF3/Boltz-2/SiteAF3 share MSA-derived information; "
                            "read the independence column.")},
                  f, indent=2)

    for L in layers:
        flag = {"True": "SUPPORT", "False": "CONTRA", "None": "---", "pending": "PEND"}.get(str(L["supports_site"]), "PEND")
        print(f"  {flag:>5}  {L['layer']:<32} {L['result'][:52]}")
    print(f"support={n_support} contradict={n_contra} neutral/pending={n_pend}")
    print(f"wrote {RESULTS / '13_evidence_matrix.json'}")


if __name__ == "__main__":
    main()