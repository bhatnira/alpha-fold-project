#!/usr/bin/env python3
"""
Analyze Boltz-2 cofolding results from job 3332533 (site21, 30 compounds).
Compare active vs inactive binding.
"""

import json
import csv
from pathlib import Path
from collections import defaultdict

PIPELINE = Path("/cluster/home/nbhatt04/lean_pipeline")
COFOLDING = PIPELINE / "site_directed_cofolding" / "results"
OUTPUT = PIPELINE / "phase2_siteaf3_validation" / "analysis"

COMPOUNDS = {
    1: {"smiles": "C1([C@@H]([C@H](O)CO)OC(=O)C=1O)O", "activity": "active", "uM": 1797},
    2: {"smiles": "C(O)C(C1OC(=O)C(O)=C1O)O", "activity": "active", "uM": 1316},
    3: {"smiles": "CC1(O[C@H]([C@H]2OC(=O)C(O)=C2O)CO1)C", "activity": "active", "uM": 6077},
    12: {"smiles": "C(#C)COC1[C@@H]([C@@H]2OC(C)(C)OC2)OC(=O)C=1O", "activity": "active", "uM": 0.198},
    18: {"smiles": "C1C=CC(COC2[C@@H]([C@H](O)CO)OC(=O)C=2O)=CC=1", "activity": "active", "uM": 1288},
    24: {"smiles": "CCCOC1C(=O)O[C@H]([C@H](O)CO)C=1O", "activity": "active", "uM": 1202},
    25: {"smiles": "C(Br)[C@H]([C@H]1OC(=O)C(O)=C1O)O", "activity": "active", "uM": 2.63},
}


def parse_confidence(result_dir):
    for p in Path(result_dir).rglob("confidence_*.json"):
        try:
            with open(p) as f:
                return json.load(f)
        except Exception:
            continue
    return {}


def parse_affinity(result_dir):
    for p in Path(result_dir).rglob("affinity_*.json"):
        try:
            with open(p) as f:
                return json.load(f)
        except Exception:
            continue
    return {}


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    
    print("=" * 70)
    print("BOLTZ-2 COFOLDING ANALYSIS — Site21, 30 Compounds")
    print("=" * 70)
    
    # Scan binary and ternary results
    results = {"binary": [], "ternary": []}
    
    for complex_type in ["binary", "ternary"]:
        site_dir = COFOLDING / complex_type
        if not site_dir.exists():
            print(f"  {complex_type}: NOT FOUND")
            continue
        
        for result_dir in sorted(site_dir.iterdir()):
            if not result_dir.is_dir():
                continue
            name = result_dir.name
            if "site21" not in name:
                continue
            
            # Find boltz_results subdirectory
            boltz_dirs = list(result_dir.glob("boltz_results_*"))
            if boltz_dirs:
                actual = boltz_dirs[0] / "predictions"
            else:
                actual = result_dir
            
            conf = parse_confidence(actual)
            aff = parse_affinity(actual)
            
            if conf or aff:
                record = {
                    "name": name,
                    "complex_type": complex_type,
                    "confidence": conf.get("confidence_score"),
                    "iptm": conf.get("iptm"),
                    "ligand_iptm": conf.get("ligand_iptm"),
                    "complex_plddt": conf.get("complex_plddt"),
                    "affinity_pred": aff.get("affinity_pred_value"),
                    "affinity_prob": aff.get("affinity_probability_binary"),
                }
                results[complex_type].append(record)
    
    # Print summary
    for ct in ["binary", "ternary"]:
        recs = results[ct]
        print(f"\n{ct.upper()}: {len(recs)} predictions")
        if not recs:
            continue
        
        actives = [r for r in recs if any(f"cpd{c}_" in r["name"] or f"cpd{c}." in r["name"] for c in COMPOUNDS if COMPOUNDS[c]["activity"] == "active")]
        inactives = [r for r in recs if any(f"cpd{c}_" in r["name"] or f"cpd{c}." in r["name"] for c in COMPOUNDS if COMPOUNDS[c]["activity"] == "inactive")]
        
        active_confs = [r["confidence"] for r in actives if r.get("confidence")]
        inactive_confs = [r["confidence"] for r in inactives if r.get("confidence")]
        
        active_iptms = [r["ligand_iptm"] for r in actives if r.get("ligand_iptm")]
        inactive_iptms = [r["ligand_iptm"] for r in inactives if r.get("ligand_iptm")]
        
        if active_confs:
            print(f"  Actives (n={len(actives)}): avg_conf={sum(active_confs)/len(active_confs):.3f}")
        if inactive_confs:
            print(f"  Inactives (n={len(inactives)}): avg_conf={sum(inactive_confs)/len(inactive_confs):.3f}")
        if active_iptms:
            print(f"  Actives iptm: {sum(active_iptms)/len(active_iptms):.3f}")
        if inactive_iptms:
            print(f"  Inactives iptm: {sum(inactive_iptms)/len(inactive_iptms):.3f}")
    
    # Write CSV
    all_records = results["binary"] + results["ternary"]
    if all_records:
        csv_path = OUTPUT / "cofolding_site21_results.csv"
        with open(csv_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=all_records[0].keys())
            writer.writeheader()
            writer.writerows(all_records)
        print(f"\nCSV: {csv_path}")
    
    print("\n" + "=" * 70)


if __name__ == "__main__":
    main()
