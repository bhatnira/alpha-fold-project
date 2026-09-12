#!/usr/bin/env python3
"""
Phase 12: Ensemble Docking Analysis

Since Vina CLI is not available, perform ensemble analysis using existing
docking results across both stoichiometries (2to3 and 3to2).

This provides robustness testing: does the site prediction hold across
different receptor conformations?
"""
import json, csv, os
from pathlib import Path
from collections import defaultdict

ROOT = Path(__file__).resolve().parent.parent
OUTDIR = ROOT / "09_docking" / "flexible_ensemble"
OUTDIR.mkdir(parents=True, exist_ok=True)

DOCKING_CSV = ROOT / "09_docking" / "full_docking_results.csv"
SAR_CSV = ROOT / "10_sar_validation" / "full_30compound_sar_matrix.csv"


def load_docking_data():
    results = []
    with open(DOCKING_CSV) as f:
        for row in csv.DictReader(f):
            row["binding_energy"] = float(row["binding_energy"])
            row["is_active"] = row["is_active"] == "True"
            results.append(row)
    return results


def ensemble_analysis(docking_data):
    """Compute consensus scores across stoichiometries for each compound x site."""
    compound_site = defaultdict(list)
    for row in docking_data:
        key = (row["compound_id"], row["site_id"])
        compound_site[key].append(row)

    ensemble_results = []
    for (compound_id, site_id), rows in compound_site.items():
        energies = [r["binding_energy"] for r in rows]
        stoichiometries = [r["stoichiometry"] for r in rows]
        is_active = rows[0]["is_active"]

        ensemble_results.append({
            "compound_id": compound_id,
            "site_id": site_id,
            "is_active": is_active,
            "mean_energy": sum(energies) / len(energies),
            "min_energy": min(energies),
            "max_energy": max(energies),
            "std_energy": (sum((e - sum(energies)/len(energies))**2 for e in energies) / len(energies)) ** 0.5 if len(energies) > 1 else 0,
            "n_stoichiometries": len(set(stoichiometries)),
            "consensus": len(set(stoichiometries)) > 1,
        })

    return ensemble_results


def discrimination_analysis(ensemble_results):
    """Evaluate active vs inactive discrimination at each site."""
    site_analysis = defaultdict(lambda: {"active": [], "inactive": []})

    for r in ensemble_results:
        site = r["site_id"]
        if r["is_active"]:
            site_analysis[site]["active"].append(r["mean_energy"])
        else:
            site_analysis[site]["inactive"].append(r["mean_energy"])

    discrimination = {}
    for site, data in site_analysis.items():
        active_energies = data["active"]
        inactive_energies = data["inactive"]

        if active_energies and inactive_energies:
            active_mean = sum(active_energies) / len(active_energies)
            inactive_mean = sum(inactive_energies) / len(inactive_energies)
            discrimination[site] = {
                "active_mean": active_mean,
                "inactive_mean": inactive_mean,
                "delta": inactive_mean - active_mean,
                "direction": "active_binds_stronger" if active_mean < inactive_mean else "inactive_binds_stronger",
                "n_active": len(active_energies),
                "n_inactive": len(inactive_energies),
                "separation_achieved": active_mean < inactive_mean,
            }

    return discrimination


def main():
    print("Phase 12: Ensemble Docking Analysis")
    print("=" * 60)

    docking_data = load_docking_data()
    print(f"Loaded {len(docking_data)} docking results")

    ensemble_results = ensemble_analysis(docking_data)
    print(f"Ensemble results: {len(ensemble_results)} compound-site combinations")

    discrimination = discrimination_analysis(ensemble_results)

    print("\nActive vs Inactive Discrimination by Site:")
    for site in sorted(discrimination.keys()):
        d = discrimination[site]
        status = "PASS" if d["separation_achieved"] else "FAIL"
        print(f"  Site {site}: active={d['active_mean']:.2f} vs inactive={d['inactive_mean']:.2f} "
              f"(delta={d['delta']:.2f}) [{status}]")

    csv_path = OUTDIR / "ensemble_docking_results.csv"
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "compound_id", "site_id", "is_active", "mean_energy", "min_energy",
            "max_energy", "std_energy", "n_stoichiometries", "consensus"
        ])
        writer.writeheader()
        for r in sorted(ensemble_results, key=lambda x: (x["site_id"], x["mean_energy"])):
            writer.writerow(r)

    report = {
        "method": "ensemble_docking_across_stoichiometries",
        "description": "Docking robustness analysis using both 2to3 and 3to2 receptor models",
        "n_compound_sites": len(ensemble_results),
        "n_sites": len(discrimination),
        "discrimination": discrimination,
        "sites_with_separation": [s for s, d in discrimination.items() if d["separation_achieved"]],
        "sites_without_separation": [s for s, d in discrimination.items() if not d["separation_achieved"]],
        "note": "Since Vina CLI is unavailable, ensemble analysis uses existing docking data across both stoichiometries",
    }

    report_path = OUTDIR / "ensemble_docking_report.json"
    report_path.write_text(json.dumps(report, indent=2))
    print(f"\nResults: {csv_path}")
    print(f"Report: {report_path}")


if __name__ == "__main__":
    main()
