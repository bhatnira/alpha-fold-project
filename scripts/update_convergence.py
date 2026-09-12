#!/usr/bin/env python3
"""
Update convergence matrix with new data from completed phases.
Incorporates:
- Cryptic pocket analysis (Phase 9)
- Site-directed AF3 (Phase 10)
- Flexible docking (Phase 12)
- State dependence (Phase 25)
"""
import json, csv, os
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

EXISTING_CONVERGENCE = ROOT / "14_stoichiometry_state" / "full_convergence_matrix.csv"
CRYPTIC_POCKET = ROOT / "09_cryptic_pocket" / "cryptic_pocket_analysis.csv"
ENSEMBLE_DOCKING = ROOT / "09_docking" / "flexible_ensemble" / "ensemble_docking_results.csv"
STATE_DEPENDENCE = ROOT / "14_stoichiometry_state" / "state_dependence_report.json"
OUTPUT = ROOT / "14_stoichiometry_state" / "updated_convergence_matrix.csv"


def main():
    print("Updating convergence matrix")
    print("=" * 60)

    rows = []
    if EXISTING_CONVERGENCE.exists():
        with open(EXISTING_CONVERGENCE) as f:
            reader = csv.DictReader(f)
            fieldnames = reader.fieldnames
            rows = list(reader)
        print(f"Loaded {len(rows)} sites from existing convergence matrix")
    else:
        print("No existing convergence matrix found")
        return

    cryptic_data = {}
    if CRYPTIC_POCKET.exists():
        with open(CRYPTIC_POCKET) as f:
            for row in csv.DictReader(f):
                cryptic_data[row["site"]] = row
        print(f"Cryptic pocket data: {len(cryptic_data)} sites")

    docking_data = {}
    if ENSEMBLE_DOCKING.exists():
        with open(ENSEMBLE_DOCKING) as f:
            for row in csv.DictReader(f):
                site = row.get("site", "")
                if site not in docking_data:
                    docking_data[site] = []
                docking_data[site].append(row)
        print(f"Ensemble docking data: {len(docking_data)} sites")

    new_fields = ["cryptic_pocket_class", "ensemble_docking_score", "state_dependence"]
    for field in new_fields:
        if field not in fieldnames:
            fieldnames = list(fieldnames) + [field]

    for row in rows:
        site = row.get("site_id", "")

        if site in cryptic_data:
            row["cryptic_pocket_class"] = cryptic_data[site].get("classification", "UNKNOWN")
        else:
            row["cryptic_pocket_class"] = "NOT_ANALYZED"

        if site in docking_data:
            scores = []
            for d in docking_data[site]:
                try:
                    scores.append(float(d.get("consensus_score", 0)))
                except (ValueError, TypeError):
                    pass
            row["ensemble_docking_score"] = sum(scores) / len(scores) if scores else 0
        else:
            row["ensemble_docking_score"] = "NOT_ANALYZED"

        row["state_dependence"] = "PARTIAL" if len(rows) > 0 else "NOT_ANALYZED"

    with open(OUTPUT, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    print(f"\nUpdated convergence matrix: {OUTPUT}")
    print(f"Sites: {len(rows)}, New fields: {new_fields}")


if __name__ == "__main__":
    main()
