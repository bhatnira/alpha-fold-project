#!/usr/bin/env python3
"""5. Boltz-2 ligand-conditioned consensus (CF / Delta-CF) across the full panel.

Contact frequency per residue:

    CF(r)  = #valid predictions contacting r  /  #valid predictions
    CFa(r) = CF over active compounds
    CFi(r) = CF over inactive compounds
    dCF(r) = CFa(r) - CFi(r)

Reads the parsed models from script 02 output (or re-parses). Reports
recurring regions, pose consensus, and active/inactive differences.

Also implements the prompt's inactive-compound challenge (section 11):
*discordant poses* are enumerated for inactive compounds that nevertheless
produce apparently favorable poses at the candidate site (high site-23
contact rate / high pocket occupancy despite being experimentally inactive).

Affinity layer is absent (no affinity_*.json) -> recorded as not available.

Output: outputs/05_boltz2_consensus.json
"""

import json
from collections import defaultdict
from pathlib import Path

import numpy as np

from config import (RESULTS, ACTIVE_CPDS, INACTIVE_CPDS, CHAIN_SUBUNIT,
                    COFOLDING_RESULTS, CONTACT_CUTOFF)

CA = CONTACT_CUTOFF


def parse_name(name):
    parts = name.split("_")
    if len(parts) >= 4:
        return {"stoich": parts[1], "site": int(parts[2].replace("site", "")),
                "cpd": int(parts[3].replace("cpd", ""))}
    return None


def main():
    data = None
    prev = RESULTS / "02_ensemble_convergence.json"
    if prev.exists():
        with open(prev) as f:
            data = json.load(f)

    active_den = {"all": 0, "active": 0, "inactive": 0}
    contact_counts = defaultdict(lambda: {"active": 0, "inactive": 0, "all": 0})
    per_condition = []
    discordant = []

    # Build compound class sets once.
    act = set(ACTIVE_CPDS)
    inact = set(INACTIVE_CPDS)

    if data and data.get("full_panel"):
        for c in data["full_panel"]["conditions"]:
            meta = parse_name(c["name"])
            if not meta:
                continue
            cpd = meta["cpd"]
            if cpd in act:
                cls = "active"
            elif cpd in inact:
                cls = "inactive"
            else:
                cls = "unknown"
            stoich = meta["stoich"]
            sub = CHAIN_SUBUNIT.get(stoich, {})
            n_valid = 0
            site23_models = 0
            occupancy_scores = []
            for m in c.get("models", []):
                if not m.get("contacts"):
                    continue
                n_valid += 1
                active_den["all"] += 1
                active_den[cls] += 1
                for key in m["contacts"]:
                    # contacts are stored as "<chain>.<resnum>" strings
                    if key.count(".") != 1:
                        continue
                    chain, rn = key.rsplit(".", 1)
                    try:
                        rn = int(rn)
                    except ValueError:
                        continue
                    res = f"{sub.get(chain, chain)}.{rn}"
                    contact_counts[res]["all"] += 1
                    contact_counts[res][cls] += 1
                if m.get("site_hit") == 23:
                    site23_models += 1
            condition_row = {
                "condition": c["name"], "cpd": cpd, "class": cls,
                "site_conditioned": meta["site"], "stoich": stoich,
                "n_valid_models": n_valid,
                "n_models": c.get("n_models", 0),
                "n_site23_contacting_models": site23_models,
                "site23_contact_fraction": round(site23_models / max(1, n_valid), 3),
            }
            per_condition.append(condition_row)

            if cls == "inactive" and n_valid and site23_models / n_valid >= 0.5:
                discordant.append({
                    "cpd": cpd,
                    "condition": c["name"],
                    "site23_contact_fraction": round(site23_models / n_valid, 3),
                    "n_models": n_valid,
                })

    rows = []
    for res, c in sorted(contact_counts.items()):
        rows.append({
            "residue": res,
            "cf_all": c["all"] / max(1, active_den["all"]),
            "cf_active": c["active"] / max(1, active_den["active"]),
            "cf_inactive": c["inactive"] / max(1, active_den["inactive"]),
            "dcf": c["active"] / max(1, active_den["active"]) - c["inactive"] / max(1, active_den["inactive"]),
            "n_contacts_all": c["all"],
        })
    rows.sort(key=lambda r: r["dcf"], reverse=True)

    discordant_by_cpd = defaultdict(list)
    for d in discordant:
        discordant_by_cpd[d["cpd"]].append(d)

    # Per-compound-class site statistics (which candidate site do ligands contact most).
    class_site_tally = defaultdict(lambda: defaultdict(int))
    for cond in (data.get("full_panel") or {}).get("conditions", []):
        meta = parse_name(cond["name"])
        if not meta:
            continue
        cls = "active" if meta["cpd"] in act else ("inactive" if meta["cpd"] in inact else "unknown")
        for m in cond.get("models", []):
            if m.get("site_hit") is not None:
                class_site_tally[cls][m["site_hit"]] += 1
    site_fraction = {}
    for cls, tally in class_site_tally.items():
        tot = sum(tally.values())
        if tot:
            site_fraction[cls] = {str(s): round(n / tot, 3) for s, n in
                                  sorted(tally.items(), key=lambda kv: -kv[1])}

    summary = {
        "n_conditions_with_models": len(per_condition),
        "n_valid_ligand_models": active_den["all"],
        "contact_cutoff_A": CA,
        "active_cpds": sorted(act),
        "inactive_cpds": sorted(inact),
        "affinity_layer": "ABSENT (no affinity_*.json in results)",
        "residue_table": rows,
        "top_dcf_pos": [r["residue"] for r in rows[:20]],
        "top_dcf_neg": [r["residue"] for r in rows if r["dcf"] < 0][:10],
        "site_class_fraction": site_fraction,
        "discordant_inactive_conditions": sorted(discordant,
                                                 key=lambda d: -d["site23_contact_fraction"]),
        "discordant_summary": {
            "n_inactive_conditions": sum(1 for c in per_condition if c["class"] == "inactive"),
            "n_discordant_conditions": len(discordant),
            "discordant_cpds": sorted(discordant_by_cpd.keys()),
            "note": "Discordant = experimentally inactive compound with >=50% of models "
                    "contacting site-23-defining residues (favorable computational pose).",
        },
    }

    with open(RESULTS / "05_boltz2_consensus.json", "w") as f:
        json.dump(summary, f, indent=2)

    print(f"[1/4] panel conditions parsed: {summary['n_conditions_with_models']}; "
          f"valid ligand models: {summary['n_valid_ligand_models']}")
    print(f"[2/4] residues with contacts: {len(rows)}")
    print("\nTop 15 residues by Delta-CF (active - inactive):")
    for r in rows[:15]:
        print(f"  {r['residue']:<8} CF(all)={r['cf_all']:.3f} dCF={r['dcf']:+.3f} "
              f"(a={r['cf_active']:.2f} i={r['cf_inactive']:.2f})")
    print("\nSite contact fraction by compound class:")
    for cls, fracs in summary["site_class_fraction"].items():
        print(f"  {cls}: {fracs}")
    print(f"\n[3/4] discordant inactive conditions: {len(discordant)} "
          f"(cpds: {sorted(discordant_by_cpd.keys())})")
    print("[4/4] affinity layer: ABSENT")
    print(f"wrote {RESULTS / '05_boltz2_consensus.json'}")


if __name__ == "__main__":
    main()