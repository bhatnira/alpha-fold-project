#!/usr/bin/env python
"""
PHASE 1+2: Data QC + Experimental SAR Analysis
Master Prompt Compliance: Sections 5-6

Runs BEFORE any structural modeling.
Analyzes the 30-compound dataset independently.
"""

import csv
import json
import numpy as np
from pathlib import Path
from collections import defaultdict, Counter
from itertools import combinations

PIPELINE_DIR = Path("/cluster/home/nbhatt04/lean_pipeline")
CSV_PATH = PIPELINE_DIR / "modulator-dataset-a9a10.csv"
OUTPUT_DIR = PIPELINE_DIR / "pam_project/01_qc"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ============================================================
# PHASE 1: DATA QC
# ============================================================

def load_raw_data():
    """Load and QC the raw dataset."""
    compounds = []
    with open(CSV_PATH, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            compounds.append({
                "id": int(row["Identifier"]),
                "smiles_raw": row["Smiles"],
                "activity_uM": float(row["Activity (uM)"]),
                "potentiation_pct": float(row["%Potentiation"]),
            })
    return compounds


def qc_dataset(compounds):
    """Comprehensive QC of the dataset."""
    qc = {
        "total_compounds": len(compounds),
        "n_active": sum(1 for c in compounds if c["activity_uM"] > 0),
        "n_inactive": sum(1 for c in compounds if c["activity_uM"] == 0),
        "checks": {},
    }
    
    # Check: duplicate IDs
    ids = [c["id"] for c in compounds]
    qc["checks"]["duplicate_ids"] = len(ids) != len(set(ids))
    
    # Check: duplicate SMILES
    smiles = [c["smiles_raw"] for c in compounds]
    qc["checks"]["duplicate_smiles"] = len(smiles) != len(set(smiles))
    
    # Check: invalid SMILES (basic check)
    qc["checks"]["empty_smiles"] = any(not s.strip() for s in smiles)
    
    # Check: missing values
    qc["checks"]["missing_activity"] = any(c["activity_uM"] is None for c in compounds)
    qc["checks"]["missing_potentiation"] = any(c["potentiation_pct"] is None for c in compounds)
    
    # Check: consistency of active/inactive labels
    inconsistent = []
    for c in compounds:
        if c["activity_uM"] == 0 and c["potentiation_pct"] > 0:
            inconsistent.append(c["id"])
        elif c["activity_uM"] > 0 and c["potentiation_pct"] == 0:
            inconsistent.append(c["id"])
    qc["checks"]["inconsistent_labels"] = inconsistent
    
    # Check: stereochemistry
    has_stereo = sum(1 for s in smiles if "@" in s)
    qc["checks"]["compounds_with_stereochemistry"] = has_stereo
    
    # Check: protonation states (anions)
    has_anion = sum(1 for s in smiles if "[O-]" in s or "[N+]" in s)
    qc["checks"]["compounds_with_charges"] = has_anion
    
    # Activity distribution
    active = sorted([c for c in compounds if c["activity_uM"] > 0], 
                    key=lambda x: x["activity_uM"])
    qc["activity_distribution"] = {
        "active_compounds": [{"id": c["id"], "activity_uM": c["activity_uM"], 
                              "potentiation_pct": c["potentiation_pct"]} for c in active],
        "potency_range": {
            "min_uM": active[0]["activity_uM"] if active else None,
            "max_uM": active[-1]["activity_uM"] if active else None,
            "range_fold": active[-1]["activity_uM"] / active[0]["activity_uM"] if active and active[0]["activity_uM"] > 0 else None,
        },
        "n_inactive": sum(1 for c in compounds if c["activity_uM"] == 0),
    }
    
    return qc


# ============================================================
# PHASE 2: EXPERIMENTAL SAR (NO STRUCTURAL MODELING)
# ============================================================

def classify_compounds(compounds):
    """Classify by activity level."""
    classified = []
    for c in compounds:
        act = c["activity_uM"]
        pot = c["potentiation_pct"]
        
        if act == 0 and pot == 0:
            cls = "inactive"
        elif act < 1:
            cls = "highly_potent"
        elif act < 10:
            cls = "potent"
        elif act < 2000:
            cls = "moderate"
        else:
            cls = "weak"
        
        classified.append({**c, "activity_class": cls})
    return classified


def experimental_sar_analysis(compounds):
    """Purely experimental SAR - no structural modeling."""
    
    # Rank by potency
    active = sorted([c for c in compounds if c["activity_uM"] > 0], 
                    key=lambda x: x["activity_uM"])
    
    # Rank by potentiation
    by_potentiation = sorted([c for c in compounds if c["potentiation_pct"] > 0],
                            key=lambda x: -x["potentiation_pct"])
    
    sar = {
        "potency_ranking": [
            {"rank": i+1, "id": c["id"], "activity_uM": c["activity_uM"], 
             "potentiation_pct": c["potentiation_pct"]}
            for i, c in enumerate(active)
        ],
        "potentiation_ranking": [
            {"rank": i+1, "id": c["id"], "activity_uM": c["activity_uM"],
             "potentiation_pct": c["potentiation_pct"]}
            for i, c in enumerate(by_potentiation)
        ],
    }
    
    # SAR anchors
    sar["sar_anchors"] = {
        "strongest_active": {"id": active[0]["id"], "activity_uM": active[0]["activity_uM"],
                             "smiles": active[0]["smiles_raw"],
                             "potentiation": active[0]["potentiation_pct"]},
        "second_strongest": {"id": active[1]["id"], "activity_uM": active[1]["activity_uM"],
                             "smiles": active[1]["smiles_raw"],
                             "potentiation": active[1]["potentiation_pct"]} if len(active) > 1 else None,
        "weakest_active": {"id": active[-1]["id"], "activity_uM": active[-1]["activity_uM"],
                           "smiles": active[-1]["smiles_raw"],
                           "potentiation": active[-1]["potentiation_pct"]},
        "highly_potent": [c for c in active if c["activity_uM"] < 1],
        "potent": [c for c in active if 1 <= c["activity_uM"] < 10],
        "moderate": [c for c in active if 10 <= c["activity_uM"] < 2000],
        "weak": [c for c in active if c["activity_uM"] >= 2000],
    }
    
    # Closest active/inactive pairs (structurally similar)
    inactive = [c for c in compounds if c["activity_uM"] == 0]
    
    # Structural similarity based on SMILES length as proxy
    sar["structural_observations"] = {
        "active_smiles_lengths": [len(c["smiles_raw"]) for c in active],
        "inactive_smiles_lengths": [len(c["smiles_raw"]) for c in inactive],
        "avg_active_smiles_len": np.mean([len(c["smiles_raw"]) for c in active]),
        "avg_inactive_smiles_len": np.mean([len(c["smiles_raw"]) for c in inactive]),
    }
    
    # Potentiation vs potency correlation
    active_pots = [c["potentiation_pct"] for c in active]
    active_acts = [c["activity_uM"] for c in active]
    if len(active) > 2:
        from scipy import stats
        rho, pval = stats.spearmanr(active_acts, active_pots)
        sar["potency_potentiation_correlation"] = {
            "spearman_rho": float(rho),
            "p_value": float(pval),
            "interpretation": "positive" if rho > 0 else "negative",
        }
    
    # Activity class distribution
    classes = Counter(c["activity_class"] for c in compounds)
    sar["activity_classes"] = dict(classes)
    
    # Experimental SAR hypothesis (without structural modeling)
    sar["experimental_sar_hypothesis"] = {
        "observation_1": "Potency spans 3+ orders of magnitude (0.2 to 6077 uM)",
        "observation_2": "Potentiation ranges 150-500% and does NOT correlate with potency",
        "observation_3": "The most potent compound (Cpd 12, 0.2 uM) has LOWEST potentiation (150%)",
        "observation_4": "The highest potentiation compound (Cpd 18, 500%) has MODERATE potency (1288 uM)",
        "observation_5": "23 of 30 compounds are completely inactive (0 activity, 0 potentiation)",
        "observation_6": "Potency and potentiation appear INDEPENDENT - suggesting different structural determinants",
        "key_prediction": "The binding site must accommodate compounds spanning ~30,000-fold potency range while allowing independent potentiation mechanisms",
    }
    
    return sar


def generate_compound_table(compounds):
    """Generate formatted compound table."""
    table = []
    for c in sorted(compounds, key=lambda x: x["id"]):
        act = c["activity_uM"]
        pot = c["potentiation_pct"]
        
        if act == 0 and pot == 0:
            cls = "INACTIVE"
        elif act < 1:
            cls = "HIGHLY_POTENT"
        elif act < 10:
            cls = "POTENT"
        elif act < 2000:
            cls = "MODERATE"
        else:
            cls = "WEAK"
        
        table.append({
            "ID": c["id"],
            "SMILES": c["smiles_raw"],
            "Activity_uM": act,
            "Potentiation_pct": pot,
            "Class": cls,
            "Log_Activity": f"{np.log10(act):.2f}" if act > 0 else "N/A",
        })
    return table


def main():
    print("=" * 70)
    print("PHASE 1+2: DATA QC + EXPERIMENTAL SAR")
    print("=" * 70)
    
    # Load raw data
    compounds = load_raw_data()
    print(f"\nLoaded {len(compounds)} compounds")
    
    # QC
    qc = qc_dataset(compounds)
    print(f"\n--- QC RESULTS ---")
    print(f"Total compounds: {qc['total_compounds']}")
    print(f"Active: {qc['n_active']}")
    print(f"Inactive: {qc['n_inactive']}")
    print(f"Duplicate IDs: {qc['checks']['duplicate_ids']}")
    print(f"Duplicate SMILES: {qc['checks']['duplicate_smiles']}")
    print(f"Compounds with stereochemistry: {qc['checks']['compounds_with_stereochemistry']}")
    print(f"Compounds with charges: {qc['checks']['compounds_with_charges']}")
    print(f"Inconsistent labels: {qc['checks']['inconsistent_labels']}")
    
    # Master prompt says 28 compounds (7 active, 21 inactive)
    # Dataset has 30 compounds (7 active, 23 inactive)
    print(f"\n*** DISCREPANCY: Master prompt specifies 28 compounds (7 active, 21 inactive)")
    print(f"*** Actual dataset: {qc['total_compounds']} compounds ({qc['n_active']} active, {qc['n_inactive']} inactive)")
    print(f"*** Two extra inactive compounds: IDs 29, 30")
    
    # Classify compounds
    classified = classify_compounds(compounds)
    
    # Experimental SAR
    sar = experimental_sar_analysis(classified)
    
    print(f"\n--- EXPERIMENTAL SAR ---")
    print(f"\nPotency Ranking:")
    for r in sar["potency_ranking"]:
        print(f"  #{r['rank']}: Cpd {r['id']} = {r['activity_uM']} uM ({r['potentiation_pct']}% pot)")
    
    print(f"\nPotentiation Ranking:")
    for r in sar["potentiation_ranking"][:5]:
        print(f"  #{r['rank']}: Cpd {r['id']} = {r['potentiation_pct']}% pot ({r['activity_uM']} uM)")
    
    print(f"\nActivity Classes:")
    for cls, count in sar["activity_classes"].items():
        print(f"  {cls}: {count}")
    
    print(f"\n--- KEY SAR OBSERVATIONS ---")
    for key, val in sar["experimental_sar_hypothesis"].items():
        if key.startswith("observation"):
            print(f"  {val}")
    
    # Generate compound table
    table = generate_compound_table(compounds)
    
    # Save outputs
    with open(OUTPUT_DIR / "qc_report.json", "w") as f:
        json.dump(qc, f, indent=2)
    
    with open(OUTPUT_DIR / "experimental_sar.json", "w") as f:
        json.dump(sar, f, indent=2, default=str)
    
    with open(OUTPUT_DIR / "compound_table.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["ID", "SMILES", "Activity_uM", "Potentiation_pct", "Class", "Log_Activity"])
        writer.writeheader()
        writer.writerows(table)
    
    # Copy raw data
    import shutil
    shutil.copy(CSV_PATH, OUTPUT_DIR / "raw_dataset.csv")
    
    print(f"\nOutputs saved to: {OUTPUT_DIR}")
    
    return qc, sar


if __name__ == "__main__":
    main()
