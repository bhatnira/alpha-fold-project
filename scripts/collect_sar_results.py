#!/usr/bin/env python3
"""
Collect and analyze Boltz-2 SAR affinity results for all 7 active compounds.
"""
import json, csv, os, glob
from pathlib import Path

ROOT = Path("/cluster/home/nbhatt04/lean_pipeline")
OUT_DIR = ROOT / "phase08_boltz_affinity" / "boltz2_sar_outputs"
RESULTS_CSV = ROOT / "phase08_boltz_affinity" / "boltz2_sar_results.csv"

COMPOUNDS = {
    "1": {"smiles": "C1([C@@H]([C@H](O)CO)OC(=O)C=1O)O", "activity_uM": 1797, "potentiation": 286},
    "2": {"smiles": "C(O)C(C1OC(=O)C(O)=C1O)O", "activity_uM": 1316, "potentiation": 380},
    "3": {"smiles": "CC1(O[C@H]([C@H]2OC(=O)C(O)=C2O)CO1)C", "activity_uM": 6077, "potentiation": 154},
    "12": {"smiles": "C(#C)COC1[C@@H]([C@@H]2OC(C)(C)OC2)OC(=O)C=1O", "activity_uM": 0.1981, "potentiation": 4820},
    "18": {"smiles": "C1C=CC(COC2[C@@H]([C@H](O)CO)OC(=O)C=2O)=CC=1", "activity_uM": 1288, "potentiation": 243},
    "24": {"smiles": "CCCOC1C(=O)O[C@H]([C@H](O)CO)C=1O", "activity_uM": 1202, "potentiation": 271},
    "25": {"smiles": "C(Br)[C@H]([C@H]1OC(=O)C(O)=C1O)O", "activity_uM": 2.63, "potentiation": 3700},
}


def collect_results():
    """Collect affinity results from all predictions."""
    results = []
    
    for pred_dir in sorted(OUT_DIR.rglob("sar_cpd*")):
        if not pred_dir.is_dir():
            continue
        
        name = pred_dir.name
        cpd_id = name.replace("sar_cpd", "").replace("_predictions", "")
        
        # Find confidence JSON
        conf_files = list(pred_dir.rglob("*confidence*"))
        if conf_files:
            with open(conf_files[0]) as f:
                conf = json.load(f)
        else:
            conf = {}
        
        # Find affinity JSON
        aff_files = list(pred_dir.rglob("*affinity*"))
        if aff_files:
            with open(aff_files[0]) as f:
                aff = json.load(f)
        else:
            aff = {}
        
        results.append({
            "compound_id": cpd_id,
            "smiles": COMPOUNDS.get(cpd_id, {}).get("smiles", ""),
            "activity_uM": COMPOUNDS.get(cpd_id, {}).get("activity_uM", 0),
            "potentiation": COMPOUNDS.get(cpd_id, {}).get("potentiation", 0),
            "affinity_value": aff.get("affinity_pred_value", None),
            "affinity_probability": aff.get("affinity_probability_binary", None),
            "confidence_score": conf.get("confidence_score", None),
            "iptm": conf.get("iptm", None),
            "ptm": conf.get("ptm", None),
        })
    
    return results


def analyze_sar(results):
    """Analyze SAR consistency."""
    print("\n" + "=" * 70)
    print("BOLTZ-2 SAR ANALYSIS - ALL 7 ACTIVE COMPOUNDS")
    print("=" * 70)
    
    # Sort by potency (most potent first)
    results.sort(key=lambda x: x["activity_uM"])
    
    print(f"\n{'Cpd':>4} {'Activity':>10} {'Potent%':>8} {'Affinity':>10} {'Prob':>8} {'Conf':>6} {'iPTM':>6}")
    print("-" * 70)
    
    for r in results:
        aff = f"{r['affinity_value']:.3f}" if r['affinity_value'] else "N/A"
        prob = f"{r['affinity_probability']:.3f}" if r['affinity_probability'] else "N/A"
        conf = f"{r['confidence_score']:.3f}" if r['confidence_score'] else "N/A"
        iptm = f"{r['iptm']:.3f}" if r['iptm'] else "N/A"
        print(f"{r['compound_id']:>4} {r['activity_uM']:>10.2f} {r['potentiation']:>8} {aff:>10} {prob:>8} {conf:>6} {iptm:>6}")
    
    # Check correlation
    valid = [r for r in results if r['affinity_value'] is not None]
    if len(valid) >= 3:
        # Lower activity_uM = more potent
        activities = [r['activity_uM'] for r in valid]
        affinities = [r['affinity_value'] for r in valid]
        
        # Calculate rank correlation
        from scipy.stats import spearmanr
        corr, pval = spearmanr(activities, affinities)
        print(f"\nSpearman correlation (activity vs affinity): r={corr:.3f}, p={pval:.4f}")
        
        if corr < 0:
            print("  → NEGATIVE correlation (more potent = higher affinity) ✓")
        else:
            print("  → POSITIVE correlation (unexpected)")
    
    return results


def main():
    print("Collecting Boltz-2 SAR results...")
    results = collect_results()
    
    if not results:
        print("No results found. Check if predictions completed.")
        return
    
    # Save CSV
    with open(RESULTS_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=results[0].keys())
        writer.writeheader()
        writer.writerows(results)
    print(f"Results saved: {RESULTS_CSV}")
    
    # Analyze
    analyze_sar(results)


if __name__ == "__main__":
    main()
