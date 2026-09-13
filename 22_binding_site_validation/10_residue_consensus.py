#!/usr/bin/env python3
"""10. Residue-level consensus table (Tier 1/2/3).

Tier 1: structural support from >=2 independent methods (AF3 + Boltz-2, or
        SiteAF3 + one more) AND SAR support (positive dCF for active
        compounds) AND membership in a candidate pocket.
Tier 2: strong structural support but no SAR enrichment.
Tier 3: prediction-only evidence.

Key change vs the first pass: "SAR support" is now evidence-based -- a
residue must be contacted more often by active than inactive compounds
(dCF > SAR_DCF_MIN). It is no longer implied by merely being in an
alpha9/alpha10 subunit (which covers every residue).

Inputs: outputs/06_three_method_convergence.json, outputs/05_boltz2_consensus.json,
outputs/02_ensemble_convergence.json, outputs/09_sar_activity_cliff_validation.json.
Output: outputs/10_residue_consensus.json
"""

import json
from pathlib import Path

from config import RESULTS, POCKET_RESIDUES, ACTIVE_CPDS

SAR_DCF_MIN = 0.05   # minimum active-vs-inactive contact-frequency difference
MAX_TIER_LIST = 50


def load(o):
    p = RESULTS / o
    if p.exists():
        with open(p) as f:
            return json.load(f)
    return None


def main():
    b2 = load("05_boltz2_consensus.json") or {}
    conv = load("06_three_method_convergence.json") or {}
    ens = load("02_ensemble_convergence.json") or {}
    sar = load("09_sar_activity_cliff_validation.json") or {}

    cf = {r["residue"]: r for r in b2.get("residue_table", [])}
    rec = {}
    for site, info in (ens.get("local_convergence") or {}).items():
        rec[int(site)] = info.get("self_recurrence")

    site_res = {}
    for site in reversed(list(POCKET_RESIDUES)):
        for sub, reslist in POCKET_RESIDUES[site].items():
            for r in reslist:
                site_res[(sub, r)] = site

    agre = conv.get("residue_agreement", [])
    table = []
    for row in agre:
        sub = row["subunit"]
        rn = row["resnum"]
        key = f"{sub}.{rn}"
        cfr = cf.get(key, {})
        cfv = cfr.get("cf_all", 0.0)
        dcf = cfr.get("dcf", 0.0)
        site = row.get("site_candidate")
        af3_n = row.get("af3_n_contacts", 0)
        sa_n = row.get("siteaf3_n_contacts", 0)
        n_methods = row.get("n_methods_agree", 0)
        n_struct = (1 if (af3_n > 0) else 0) + (1 if (cfv > 0) else 0) + (1 if (sa_n > 0) else 0)

        sar_support = dcf >= SAR_DCF_MIN
        tier = 3
        if n_struct >= 2:
            tier = 2
        if n_struct >= 2 and sar_support and (site is not None):
            tier = 1
        table.append({
            "residue": key, "site": site, "site_recurrence_pct": rec.get(site),
            "af3_n": af3_n, "boltz2_cf": cfv, "siteaf3_n": sa_n,
            "n_structural_methods": n_struct,
            "dcf": round(dcf, 4), "sar_support": sar_support, "tier": tier,
        })
    table.sort(key=lambda r: (-r["tier"], -r["boltz2_cf"], -r["dcf"]))

    tiers = {}
    for t in (1, 2, 3):
        tiers[t] = [r["residue"] for r in table if r["tier"] == t][:MAX_TIER_LIST]

    cliff_residues = set()
    for p in (sar.get("pairs") or []):
        for sr in p.get("site_results", []):
            cliff_residues.update(sr.get("defining_residues", {}).keys())
            for r_res in sr.get("partial_residues", {}).keys():
                cliff_residues.add(r_res)

    with open(RESULTS / "10_residue_consensus.json", "w") as f:
        json.dump({"thresholds": {"sar_dcf_min": SAR_DCF_MIN},
                   "rows": table, "tier_1": tiers[1], "tier_2": tiers[2],
                   "tier_3": tiers[3],
                   "cliff_residues": sorted(cliff_residues),
                   "note": ("Tier 1 requires >=2 structural methods + dCF>=%.2f + "
                            "candidate-pocket membership." % SAR_DCF_MIN),
                   "active_cpds": ACTIVE_CPDS}, f, indent=2)

    for t in (1, 2, 3):
        print(f"Tier {t}: {len(tiers[t])} residues -> {tiers[t][:8]}...")
    print(f"SAR-cliff residues (context): {sorted(cliff_residues)[:12]}")
    print(f"wrote {RESULTS / '10_residue_consensus.json'}")


if __name__ == "__main__":
    main()