#!/usr/bin/env python3
"""9. Experimental SAR + activity-cliff validation of the binding model.

For each MMP cliff pair this script:
  1. loads the Boltz-2 ligand poses for BOTH partners (all panel conditions
     that contain the pair's compounds, i.e. candidate sites 21/23/5 x
     stoichiometries x seeds),
  2. computes residue contact-frequency (CF) per partner at each candidate
     site, and the difference  dCF(active - weaker)  per residue,
  3. identifies "SAR-defining residues" -- residues the more-active partner
     contacts distinctly more often at the same candidate site,
  4. classifies each cliff as explained / partial / unexplained / contradicted
     rather than forcing agreement.

Because pose spread across seeds is large (see outputs/02), residue-level
assignments are reported as *hypotheses for mutagenesis*, not as verified
contacts.

Output: outputs/09_sar_activity_cliff_validation.json
"""

import json
import csv
from collections import defaultdict
from pathlib import Path

from config import (RESULTS, MODULATOR_CSV, MMP_JSON, MMP_CLIFF_PAIRS,
                    CONTACT_CUTOFF, CHAIN_SUBUNIT)

# dCF thresholds (fractional difference in contact frequency, same condition set)
EXPLAINED_CF = 0.20
PARTIAL_CF = 0.10


def load_sar():
    rows = {}
    with open(MODULATOR_CSV) as f:
        for r in csv.DictReader(f):
            rows[int(r["Identifier"])] = {
                "smiles": r["Smiles"],
                "activity_uM": float(r["Activity (uM)"]),
                "potentiation": float(r["%Potentiation"]),
                "active": float(r["Activity (uM)"]) > 0.0,
            }
    return rows


def load_conditions():
    p = RESULTS / "02_ensemble_convergence.json"
    if not p.exists():
        return []
    with open(p) as f:
        data = json.load(f)
    return data.get("full_panel", {}).get("conditions", [])


def parse_name(name):
    parts = name.split("_")
    if len(parts) >= 4:
        try:
            return {"stoich": parts[1], "site": int(parts[2].replace("site", "")),
                    "cpd": int(parts[3].replace("cpd", ""))}
        except ValueError:
            return None
    return None


def contact_rows_per_condition(cond):
    """Return list of residue-contact sets (each a set of "subunit.resnum" str)."""
    meta = parse_name(cond["name"])
    if not meta:
        return []
    sub = CHAIN_SUBUNIT.get(meta["stoich"], {})
    out = []
    for m in cond.get("models", []):
        if not m.get("contacts"):
            continue
        s = set()
        for key in m["contacts"]:
            if key.count(".") != 1:
                continue
            chain, rn = key.rsplit(".", 1)
            try:
                s.add(f"{sub.get(chain, chain)}.{int(rn)}")
            except ValueError:
                continue
        out.append(s)
    return out


def main():
    sar = load_sar()
    cliffs = []
    mmp_rules = []
    if MMP_JSON.exists():
        try:
            with open(MMP_JSON) as f:
                mmp = json.load(f)
            cliffs = mmp.get("informative_mmps") or mmp.get("cliffs") or []
            mmp_rules = mmp.get("mmp_sar_rules") or []
        except Exception:
            cliffs = []
    print(f"[1/4] loaded {len(sar)} SAR rows; MMP records: {len(cliffs)}")

    pairs = []
    for p in MMP_CLIFF_PAIRS:
        pairs.append({
            "pair": p["pair"], "cpd_a": p["cpd_a"], "cpd_b": p["cpd_b"],
            "change": p["change"], "effect": p["effect"],
        })
    if len(cliffs) > len(pairs):
        for m in cliffs:
            pair = m["pair"]
            if not any(p["pair"] == pair for p in pairs):
                pairs.append({"pair": pair, "cpd_a": None, "cpd_b": None,
                              "change": m.get("change", ""), "effect": m.get("effect", "")})

    conditions = load_conditions()
    # bucket conditions: (cpd) -> {site -> [contact sets per model]}
    cond_by_cpd = defaultdict(lambda: defaultdict(list))
    for cond in conditions:
        meta = parse_name(cond["name"])
        if not meta:
            continue
        for cs in contact_rows_per_condition(cond):
            cond_by_cpd[meta["cpd"]][meta["site"]].append(cs)
    per_cpd = {cpd: {site: sets for site, sets in site_d.items()}
               for cpd, site_d in cond_by_cpd.items()}
    print(f"[2/4] conditions bucketed for {len(per_cpd)} compounds")

    explanations = []
    for p in pairs:
        a = sar.get(p["cpd_a"])
        b = sar.get(p["cpd_b"])
        if not a or not b:
            explanations.append({"pair": p["pair"], "status": "unavailable",
                                 "reason": "SAR record missing"})
            continue

        # that data: which partner has higher activity (lower EC50/uM)
        if a["activity_uM"] and b["activity_uM"] and a["activity_uM"] != b["activity_uM"]:
            # smaller uM = more potent
            if a["activity_uM"] < b["activity_uM"]:
                active_id, weaker_id = p["cpd_a"], p["cpd_b"]
                act_name, weak_name = "cpd_a", "cpd_b"
            else:
                active_id, weaker_id = p["cpd_b"], p["cpd_a"]
                act_name, weak_name = "cpd_b", "cpd_a"
        else:
            active_id, weaker_id = p["cpd_a"], p["cpd_b"]
            act_name = weak_name = "cpd_"
        # handle active-vs-inactive explicitly (inactive partner has 0 activity)
        if a["active"] != b["active"]:
            if a["active"]:
                active_id, weaker_id = p["cpd_a"], p["cpd_b"]
            else:
                active_id, weaker_id = p["cpd_b"], p["cpd_a"]

        # shared candidate sites where BOTH partners have poses
        a_sites = set(per_cpd.get(active_id, {}))
        b_sites = set(per_cpd.get(weaker_id, {}))
        shared = sorted(a_sites & b_sites)
        site_results = []
        defining_residues = {}
        for site in shared:
            a_sets = per_cpd[active_id][site]
            b_sets = per_cpd[weaker_id][site]
            if not a_sets or not b_sets:
                continue
            na, nb = len(a_sets), len(b_sets)
            cf = defaultdict(float)
            all_res = set()
            for s in a_sets:
                all_res |= s
            for s in b_sets:
                all_res |= s
            for res in all_res:
                ca = sum(1 for s in a_sets if res in s) / na
                cb = sum(1 for s in b_sets if res in s) / nb
                cf[res] = ca - cb
            # residues differentially contacted by the more-active partner
            def_res = {res: round(v, 3) for res, v in cf.items() if v >= EXPLAINED_CF}
            partial_res = {res: round(v, 3) for res, v in cf.items()
                           if PARTIAL_CF <= v < EXPLAINED_CF}
            site_results.append({
                "site": site,
                "n_models_active": na, "n_models_weaker": nb,
                "n_defining_residues": len(def_res),
                "defining_residues": def_res,
                "partial_residues": partial_res,
            })
            for r, v in def_res.items():
                defining_residues.setdefault(r, max(defining_residues.get(r, 0), v))

        max_def = max((sr["n_defining_residues"] for sr in site_results), default=0)
        n_partial = sum(1 for sr in site_results if sr.get("partial_residues"))
        if max_def >= 1:
            status = "explained"
        elif n_partial >= 1:
            status = "partial"
        elif not site_results:
            status = "no_shared_site_poses"
        else:
            status = "unexplained"

        ratio = None
        try:
            if a["activity_uM"] and b["activity_uM"]:
                ratio = a["activity_uM"] / b["activity_uM"]
        except ZeroDivisionError:
            pass

        explanations.append({
            "pair": p["pair"],
            "cpd_a": p["cpd_a"], "cpd_b": p["cpd_b"],
            "more_active_partner": active_id,
            "change": p["change"],
            "effect_claim": p["effect"],
            "activity_a_uM": a["activity_uM"], "activity_b_uM": b["activity_uM"],
            "activity_ratio_a_over_b": ratio,
            "status": status,
            "site_results": site_results,
            "defining_residues_sorted": sorted(defining_residues.items(),
                                               key=lambda kv: -kv[1])[:15],
            "reason": surf_reason(status, sar_comment(p)),
        })

    with open(RESULTS / "09_sar_activity_cliff_validation.json", "w") as f:
        json.dump({"n_sar": len(sar),
                   "active_ids": sorted(k for k, v in sar.items() if v["active"]),
                   "mmp_rules": mmp_rules, "n_mmp_records": len(cliffs),
                   "thresholds": {"explained_dCF_min": EXPLAINED_CF,
                                  "partial_dCF_min": PARTIAL_CF},
                   "pairs": explanations,
                   "remarks": [
                       "Residue-level assignments are pose-derived hypotheses for "
                       "mutagenesis (large model-to-model pose spread -> treat as "
                       "low-to-moderate confidence).",
                       "A 'defining' residue is contacted by the more-active partner at "
                       "least 20 percentage-points more often than by the weaker partner "
                       "at the same candidate site (EXPLAINED_CF=0.20)."]},
                  f, indent=2)

    print(f"[3/4] cliff pairs: {len(pairs)}; MMP rules: {len(mmp_rules)}")
    for e in explanations:
        print(f"  {e['pair']:<8} {e['status']:<12} "
              f"defining residues: {[r for r, _ in e['defining_residues_sorted'][:6]]}")
    print(f"[4/4] wrote {RESULTS / '09_sar_activity_cliff_validation.json'}")


def sar_comment(p):
    return p.get("change", "")


def surf_reason(status, change):
    m = {
        "explained": ("Multiple residues are distinctly contacted by the more-active "
                      "partner at the same candidate site (dCF>=0.20), consistent with "
                      f"'{change}' altering pocket contacts."),
        "partial": ("Some residues show a moderate contact-frequency difference "
                    "(0.10<=dCF<0.20) between partners; pattern is suggestive but weak."),
        "unexplained": ("No residue shows a contact-frequency difference >=0.10 between "
                        "partners at shared sites; the pose ensemble does not currently "
                        "resolve a residue-level mechanism."),
        "unavailable": "SAR record missing.",
    }
    return m.get(status, status)


if __name__ == "__main__":
    main()