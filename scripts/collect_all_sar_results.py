#!/usr/bin/env python3
"""
Collect ALL Boltz-2 SAR results (active + inactive) and compute SAR metrics.
"""
import json, csv, os, glob
from pathlib import Path

ROOT = Path("/cluster/home/nbhatt04/lean_pipeline")
ACTIVE_DIR = ROOT / "phase08_boltz_affinity" / "boltz2_sar_outputs"
INACTIVE_DIR = ROOT / "phase08_boltz_affinity" / "boltz2_inactive_outputs"
RESULTS_CSV = ROOT / "phase08_boltz_affinity" / "boltz2_full_sar_results.csv"

# All compounds with experimental data
COMPOUNDS = {}
with open(ROOT / "09_docking" / "full_docking_results.csv") as f:
    for row in csv.DictReader(f):
        if row['stoichiometry'] == '2to3' and row['site_id'] == '23':
            COMPOUNDS[row['compound_id']] = {
                "smiles": row['smiles'],
                "is_active": row['is_active'] == 'True',
                "activity_uM": float(row.get('activity_uM', 0)),
                "potentiation": float(row.get('potentiation_pct', 0)),
            }


def collect_from_dir(base_dir, prefix):
    """Collect affinity results from prediction directories."""
    results = []
    for pred_dir in sorted(base_dir.rglob(f"{prefix}*")):
        if not pred_dir.is_dir():
            continue
        
        name = pred_dir.name
        cpd_id = name.replace(f"{prefix}_cpd", "").replace("_predictions", "")
        
        conf_files = list(pred_dir.rglob("*confidence*"))
        aff_files = list(pred_dir.rglob("*affinity*"))
        
        conf = {}
        if conf_files:
            with open(conf_files[0]) as f:
                conf = json.load(f)
        
        aff = {}
        if aff_files:
            with open(aff_files[0]) as f:
                aff = json.load(f)
        
        meta = COMPOUNDS.get(cpd_id, {})
        results.append({
            "compound_id": cpd_id,
            "is_active": meta.get("is_active", False),
            "activity_uM": meta.get("activity_uM", 0),
            "potentiation": meta.get("potentiation", 0),
            "affinity_value": aff.get("affinity_pred_value", None),
            "affinity_probability": aff.get("affinity_probability_binary", None),
            "confidence_score": conf.get("confidence_score", None),
            "iptm": conf.get("iptm", None),
        })
    
    return results


def compute_sar_metrics(results):
    """Compute SAR discrimination metrics."""
    active = [r for r in results if r['is_active'] and r['affinity_value'] is not None]
    inactive = [r for r in results if not r['is_active'] and r['affinity_value'] is not None]
    
    if not active or not inactive:
        print("Insufficient data for SAR analysis")
        return
    
    active_aff = [r['affinity_value'] for r in active]
    inactive_aff = [r['affinity_value'] for r in inactive]
    
    active_mean = sum(active_aff) / len(active_aff)
    inactive_mean = sum(inactive_aff) / len(inactive_aff)
    
    active_prob = [r['affinity_probability'] for r in active if r['affinity_probability'] is not None]
    inactive_prob = [r['affinity_probability'] for r in inactive if r['affinity_probability'] is not None]
    
    print("\n" + "=" * 70)
    print("FULL SAR ANALYSIS - BOLTZ-2 AFFINITY PREDICTIONS")
    print("=" * 70)
    
    print(f"\nActive compounds: {len(active)}")
    print(f"Inactive compounds: {len(inactive)}")
    
    print(f"\n{'Metric':<25} {'Active':>12} {'Inactive':>12} {'Delta':>12} {'Direction':>15}")
    print("-" * 70)
    print(f"{'Affinity (mean)':<25} {active_mean:>12.3f} {inactive_mean:>12.3f} {active_mean-inactive_mean:>12.3f} {'HIGHER=better' if active_mean > inactive_mean else 'LOWER=better':>15}")
    
    if active_prob and inactive_prob:
        active_prob_mean = sum(active_prob) / len(active_prob)
        inactive_prob_mean = sum(inactive_prob) / len(inactive_prob)
        print(f"{'Probability (mean)':<25} {active_prob_mean:>12.3f} {inactive_prob_mean:>12.3f} {active_prob_mean-inactive_prob_mean:>12.3f} {'HIGHER=better' if active_prob_mean > inactive_prob_mean else 'LOWER=better':>15}")
    
    # Rank correlation
    if len(active) + len(inactive) >= 5:
        from scipy.stats import spearmanr
        all_compounds = active + inactive
        activities = [r['activity_uM'] for r in all_compounds]
        affinities = [r['affinity_value'] for r in all_compounds]
        corr, pval = spearmanr(activities, affinities)
        print(f"\nSpearman correlation: r={corr:.3f}, p={pval:.4f}")
        if corr < 0:
            print("  ✓ Negative correlation (more potent = higher affinity)")
        else:
            print("  ✗ Positive correlation (unexpected)")
    
    # Per-compound ranking
    print(f"\n{'Cpd':>4} {'Active':>7} {'Activity':>10} {'Affinity':>10} {'Prob':>8}")
    print("-" * 50)
    for r in sorted(active + inactive, key=lambda x: x['activity_uM'] if x['activity_uM'] > 0 else 99999):
        aff = f"{r['affinity_value']:.3f}" if r['affinity_value'] else "N/A"
        prob = f"{r['affinity_probability']:.3f}" if r['affinity_probability'] else "N/A"
        act = "YES" if r['is_active'] else "no"
        print(f"{r['compound_id']:>4} {act:>7} {r['activity_uM']:>10.2f} {aff:>10} {prob:>8}")


def main():
    print("Collecting Boltz-2 SAR results...")
    
    active_results = collect_from_dir(ACTIVE_DIR, "sar")
    inactive_results = collect_from_dir(INACTIVE_DIR, "inactive")
    
    all_results = active_results + inactive_results
    
    if not all_results:
        print("No results found.")
        return
    
    # Save CSV
    with open(RESULTS_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=all_results[0].keys())
        writer.writeheader()
        writer.writerows(all_results)
    print(f"Results saved: {RESULTS_CSV}")
    
    # Analyze
    compute_sar_metrics(all_results)


if __name__ == "__main__":
    main()
