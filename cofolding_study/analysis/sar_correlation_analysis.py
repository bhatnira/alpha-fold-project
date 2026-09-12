#!/usr/bin/env python
"""
SAR Correlation Analysis

Tests two main hypotheses:
1. Do all analogs share the same/near binding site?
2. Does interaction correlate with experimental potency?

Outputs:
- Binding site sharing analysis
- Potency vs interaction correlation
- Active/inactive discrimination
- Matched molecular pair analysis
- Design rule extraction
"""

import csv
import json
import numpy as np
from pathlib import Path
from collections import defaultdict
from scipy import stats

# Configuration
PIPELINE_DIR = Path("/cluster/home/nbhatt04/lean_pipeline")
COFOLDING_DIR = PIPELINE_DIR / "cofolding_study"
ANALYSIS_DIR = COFOLDING_DIR / "analysis"
ANALYSIS_DIR.mkdir(parents=True, exist_ok=True)

CSV_PATH = PIPELINE_DIR / "modulator-dataset-a9a10.csv"

# Key residues for pharmacophore
PHARMACOPHORE = {
    "HBD1": {"residues": ["W176", "S175"], "type": "H-bond donor", "essential": True},
    "HBA1": {"residues": ["Y120", "W81"], "type": "H-bond acceptor", "essential": True},
    "NEG1": {"residues": ["D145", "R83"], "type": "Negative charge", "essential": True},
    "HYD1": {"residues": ["Y224", "Y217"], "type": "Hydrophobic", "essential": False},
    "POL1": {"residues": ["W55"], "type": "Polar", "essential": False},
}


def load_compounds():
    compounds = []
    with open(CSV_PATH, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            activity = float(row["Activity (uM)"])
            potentiation = float(row["%Potentiation"])
            compounds.append({
                "id": int(row["Identifier"]),
                "smiles": row["Smiles"],
                "activity_uM": activity,
                "potentiation_pct": potentiation,
                "activity_class": classify(activity, potentiation),
                "is_active": activity > 0,
                "log_activity": np.log10(activity) if activity > 0 else None,
            })
    return compounds


def classify(activity, potentiation):
    if activity == 0 and potentiation == 0:
        return "inactive"
    elif activity < 1:
        return "highly_potent"
    elif activity < 10:
        return "potent"
    elif activity < 2000:
        return "moderate"
    else:
        return "weak"


def test_hypothesis_1_binding_site_sharing(compounds, fingerprints_by_site):
    """
    H1: Do all analogs share the same/near binding site?
    
    Test: If compounds bind the same site, they should share:
    1. Similar contact residue profiles
    2. Conserved anchor interactions
    3. Overlapping pocket occupancy
    """
    results = {}
    
    for site_id, fps in fingerprints_by_site.items():
        print(f"\n=== Site {site_id}: Binding Site Sharing Analysis ===")
        
        # Build contact matrix: compounds × residues
        all_residues = set()
        for fp in fps.values():
            all_residues.update(fp.keys())
        all_residues = sorted(all_residues)
        
        # Create binary contact matrix
        contact_matrix = []
        compound_ids = []
        for cid in sorted(fps.keys()):
            fp = fps[cid]
            row = [1 if fp.get(r, {}).get("contact", False) else 0 for r in all_residues]
            contact_matrix.append(row)
            compound_ids.append(cid)
        
        contact_matrix = np.array(contact_matrix)
        
        if len(compound_ids) < 2:
            print("  Insufficient compounds with data")
            continue
        
        # Compute pairwise Jaccard similarity
        jaccard_scores = []
        for i in range(len(compound_ids)):
            for j in range(i + 1, len(compound_ids)):
                intersection = np.sum(contact_matrix[i] & contact_matrix[j])
                union = np.sum(contact_matrix[i] | contact_matrix[j])
                jaccard = intersection / union if union > 0 else 0
                jaccard_scores.append(jaccard)
        
        avg_jaccard = np.mean(jaccard_scores)
        min_jaccard = np.min(jaccard_scores)
        max_jaccard = np.max(jaccard_scores)
        
        # Compute conserved anchor interactions
        anchor_residues = ["W176", "S175", "Y120", "W81", "D145"]
        anchor_contact_rates = {}
        for res in anchor_residues:
            if res in all_residues:
                idx = all_residues.index(res)
                contact_rates = contact_matrix[:, idx]
                anchor_contact_rates[res] = {
                    "contact_rate": float(np.mean(contact_rates)),
                    "active_rate": float(np.mean(contact_rates[:7])) if len(contact_rates) >= 7 else None,
                    "inactive_rate": float(np.mean(contact_rates[7:])) if len(contact_rates) > 7 else None,
                }
        
        # Statistical test: Are contacts more similar within active vs between active/inactive?
        active_ids = [c["id"] for c in compounds if c["is_active"]][:7]
        inactive_ids = [c["id"] for c in compounds if not c["is_active"]][:23]
        
        active_jaccard = []
        inactive_jaccard = []
        cross_jaccard = []
        
        for i, ci in enumerate(compound_ids):
            for j, cj in enumerate(compound_ids):
                if i >= j:
                    continue
                intersection = np.sum(contact_matrix[i] & contact_matrix[j])
                union = np.sum(contact_matrix[i] | contact_matrix[j])
                jaccard = intersection / union if union > 0 else 0
                
                if ci in active_ids and cj in active_ids:
                    active_jaccard.append(jaccard)
                elif ci in inactive_ids and cj in inactive_ids:
                    inactive_jaccard.append(jaccard)
                else:
                    cross_jaccard.append(jaccard)
        
        results[site_id] = {
            "n_compounds": len(compound_ids),
            "n_residues": len(all_residues),
            "avg_jaccard": float(avg_jaccard),
            "min_jaccard": float(min_jaccard),
            "max_jaccard": float(max_jaccard),
            "anchor_contact_rates": anchor_contact_rates,
            "active_internal_similarity": float(np.mean(active_jaccard)) if active_jaccard else None,
            "inactive_internal_similarity": float(np.mean(inactive_jaccard)) if inactive_jaccard else None,
            "cross_class_similarity": float(np.mean(cross_jaccard)) if cross_jaccard else None,
            "hypothesis_supported": avg_jaccard > 0.3,  # Threshold for site sharing
        }
        
        print(f"  Compounds analyzed: {len(compound_ids)}")
        print(f"  Residues tracked: {len(all_residues)}")
        print(f"  Average Jaccard similarity: {avg_jaccard:.3f}")
        print(f"  Active internal similarity: {results[site_id]['active_internal_similarity']}")
        print(f"  Inactive internal similarity: {results[site_id]['inactive_internal_similarity']}")
        print(f"  Cross-class similarity: {results[site_id]['cross_class_similarity']}")
        print(f"  H1 supported: {results[site_id]['hypothesis_supported']}")
        
        print(f"\n  Anchor residue contact rates:")
        for res, info in anchor_contact_rates.items():
            print(f"    {res}: {info['contact_rate']:.2f}")
    
    return results


def test_hypothesis_2_sar_correlation(compounds, fingerprints_by_site):
    """
    H2: Does interaction correlate with experimental potency?
    
    Test: 
    1. Spearman correlation between contact count and potency
    2. Active/inactive discrimination
    3. Potency continuum correlation
    """
    results = {}
    
    for site_id, fps in fingerprints_by_site.items():
        print(f"\n=== Site {site_id}: SAR Correlation Analysis ===")
        
        # Extract features for each compound
        compound_features = []
        for c in compounds:
            if c["id"] in fps:
                fp = fps[c["id"]]
                n_contacts = sum(1 for v in fp.values() if v.get("contact", False))
                min_distances = [v["min_distance"] for v in fp.values() if v.get("min_distance") is not None]
                avg_distance = np.mean(min_distances) if min_distances else None
                
                # Pharmacophore compliance
                pharmacophore_hits = 0
                for feature, info in PHARMACOPHORE.items():
                    if any(fp.get(r, {}).get("contact", False) for r in info["residues"]):
                        pharmacophore_hits += 1
                
                compound_features.append({
                    "id": c["id"],
                    "activity_uM": c["activity_uM"],
                    "log_activity": c["log_activity"],
                    "is_active": c["is_active"],
                    "activity_class": c["activity_class"],
                    "n_contacts": n_contacts,
                    "avg_distance": avg_distance,
                    "pharmacophore_hits": pharmacophore_hits,
                })
        
        if len(compound_features) < 3:
            print("  Insufficient data")
            continue
        
        # 1. Spearman correlation: contact count vs potency
        active_features = [f for f in compound_features if f["is_active"]]
        
        if len(active_features) >= 3:
            log_activities = [f["log_activity"] for f in active_features if f["log_activity"] is not None]
            contact_counts = [f["n_contacts"] for f in active_features if f["log_activity"] is not None]
            
            if len(log_activities) >= 3:
                rho, p_rho = stats.spearmanr(contact_counts, log_activities)
            else:
                rho, p_rho = None, None
        else:
            rho, p_rho = None, None
        
        # 2. Active vs inactive discrimination
        active_contacts = [f["n_contacts"] for f in compound_features if f["is_active"]]
        inactive_contacts = [f["n_contacts"] for f in compound_features if not f["is_active"]]
        
        if active_contacts and inactive_contacts:
            t_stat, p_ttest = stats.ttest_ind(active_contacts, inactive_contacts)
            effect_d = (np.mean(active_contacts) - np.mean(inactive_contacts)) / \
                       np.sqrt((np.var(active_contacts) + np.var(inactive_contacts)) / 2)
        else:
            t_stat, p_ttest, effect_d = None, None, None
        
        # 3. Pharmacophore compliance by activity class
        pharmacophore_by_class = defaultdict(list)
        for f in compound_features:
            pharmacophore_by_class[f["activity_class"]].append(f["pharmacophore_hits"])
        
        pharmacophore_summary = {}
        for cls, hits in pharmacophore_by_class.items():
            pharmacophore_summary[cls] = {
                "mean_hits": float(np.mean(hits)),
                "std_hits": float(np.std(hits)),
                "n": len(hits),
            }
        
        # 4. ANOVA across activity classes
        class_groups = []
        for cls in ["inactive", "weak", "moderate", "potent", "highly_potent"]:
            if cls in pharmacophore_by_class and len(pharmacophore_by_class[cls]) > 1:
                class_groups.append(pharmacophore_by_class[cls])
        
        if len(class_groups) >= 2:
            f_stat, p_anova = stats.f_oneway(*class_groups)
        else:
            f_stat, p_anova = None, None
        
        results[site_id] = {
            "n_compounds": len(compound_features),
            "spearman_rho": float(rho) if rho is not None else None,
            "spearman_p": float(p_rho) if p_rho is not None else None,
            "active_mean_contacts": float(np.mean(active_contacts)) if active_contacts else None,
            "inactive_mean_contacts": float(np.mean(inactive_contacts)) if inactive_contacts else None,
            "cohens_d": float(effect_d) if effect_d is not None else None,
            "ttest_p": float(p_ttest) if p_ttest is not None else None,
            "pharmacophore_by_class": pharmacophore_summary,
            "anova_f": float(f_stat) if f_stat is not None else None,
            "anova_p": float(p_anova) if p_anova is not None else None,
            "hypothesis_supported": (
                (rho is not None and p_rho is not None and p_rho < 0.05) or
                (effect_d is not None and abs(effect_d) > 0.8)
            ),
        }
        
        print(f"\n  Spearman correlation (contacts vs potency):")
        print(f"    rho = {rho}, p = {p_rho}")
        print(f"\n  Active vs inactive discrimination:")
        print(f"    Active contacts: {np.mean(active_contacts):.1f} ± {np.std(active_contacts):.1f}")
        print(f"    Inactive contacts: {np.mean(inactive_contacts):.1f} ± {np.std(inactive_contacts):.1f}")
        print(f"    Cohen's d: {effect_d:.3f}")
        print(f"    t-test p: {p_ttest:.4f}")
        print(f"\n  Pharmacophore compliance:")
        for cls, info in pharmacophore_summary.items():
            print(f"    {cls}: {info['mean_hits']:.1f} ± {info['std_hits']:.1f} hits (n={info['n']})")
        print(f"\n  ANOVA across classes: F={f_stat}, p={p_anova}")
        print(f"  H2 supported: {results[site_id]['hypothesis_supported']}")
    
    return results


def generate_design_rules(compounds, fingerprints_by_site, h1_results, h2_results):
    """Generate design rules from validated interactions."""
    rules = []
    
    for site_id in fingerprints_by_site:
        if site_id not in h1_results or site_id not in h2_results:
            continue
        
        if not h1_results[site_id].get("hypothesis_supported", False):
            continue
        
        fps = fingerprints_by_site[site_id]
        
        # Find most conserved interactions among active compounds
        active_fps = {cid: fp for cid, fp in fps.items() 
                      if any(c["id"] == cid and c["is_active"] for c in compounds)}
        
        if not active_fps:
            continue
        
        # Count contact frequency across active compounds
        contact_freq = defaultdict(int)
        for fp in active_fps.values():
            for res, info in fp.items():
                if info.get("contact", False):
                    contact_freq[res] += 1
        
        # Design rules for high-frequency contacts
        for res, freq in sorted(contact_freq.items(), key=lambda x: -x[1]):
            if freq >= len(active_fps) * 0.5:  # >50% of active compounds
                rule = {
                    "site_id": site_id,
                    "residue": res,
                    "contact_frequency": freq / len(active_fps),
                    "rule_type": "preserve",
                    "description": f"Maintain interaction with {res} (conserved in {freq}/{len(active_fps)} active compounds)",
                }
                rules.append(rule)
    
    return rules


def main():
    print("=" * 60)
    print("SAR Correlation Analysis")
    print("=" * 60)
    
    compounds = load_compounds()
    print(f"Loaded {len(compounds)} compounds")
    print(f"Active: {sum(1 for c in compounds if c['is_active'])}")
    print(f"Inactive: {sum(1 for c in compounds if not c['is_active'])}")
    
    # Load fingerprints from previous analysis
    fp_path = ANALYSIS_DIR / "interaction_fingerprints.csv"
    fingerprints_by_site = {}
    
    if fp_path.exists():
        with open(fp_path, "r") as f:
            reader = csv.DictReader(f)
            for row in reader:
                site_id = int(row["site_id"])
                compound_id = int(row["compound_id"])
                
                if site_id not in fingerprints_by_site:
                    fingerprints_by_site[site_id] = {}
                
                # Parse key residue contacts
                key_contacts = row.get("key_residue_contacts", "").split(";")
                key_contacts = [c for c in key_contacts if c]
                
                fp = {}
                for res in ["W176", "S175", "Y120", "Y224", "Y217", "W55",
                           "D145", "R83", "W81", "R143"]:
                    fp[res] = {
                        "contact": res in key_contacts,
                        "min_distance": 3.0 if res in key_contacts else None,
                        "interaction_type": "unknown",
                        "residue_type": "key",
                        "subunit": "alpha9" if res.startswith(("W1", "S1", "Y")) else "alpha10",
                    }
                
                fingerprints_by_site[site_id][compound_id] = fp
    
    print(f"\nFingerprints loaded for sites: {list(fingerprints_by_site.keys())}")
    
    # Test Hypothesis 1
    print("\n" + "=" * 60)
    print("TESTING HYPOTHESIS 1: Binding Site Sharing")
    print("=" * 60)
    h1_results = test_hypothesis_1_binding_site_sharing(compounds, fingerprints_by_site)
    
    # Test Hypothesis 2
    print("\n" + "=" * 60)
    print("TESTING HYPOTHESIS 2: SAR Correlation")
    print("=" * 60)
    h2_results = test_hypothesis_2_sar_correlation(compounds, fingerprints_by_site)
    
    # Generate design rules
    print("\n" + "=" * 60)
    print("DESIGN RULES")
    print("=" * 60)
    design_rules = generate_design_rules(compounds, fingerprints_by_site, h1_results, h2_results)
    for rule in design_rules:
        print(f"  Site {rule['site_id']}: {rule['description']}")
    
    # Save comprehensive results
    output = {
        "hypothesis_1": h1_results,
        "hypothesis_2": h2_results,
        "design_rules": design_rules,
        "n_compounds": len(compounds),
        "compounds": compounds,
    }
    
    output_path = ANALYSIS_DIR / "sar_correlation_results.json"
    with open(output_path, "w") as f:
        json.dump(output, f, indent=2, default=str)
    
    print(f"\nResults saved to: {output_path}")
    
    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    for site_id in h1_results:
        h1 = h1_results[site_id].get("hypothesis_supported", False)
        h2 = h2_results.get(site_id, {}).get("hypothesis_supported", False)
        print(f"  Site {site_id}: H1={'SUPPORTED' if h1 else 'NOT SUPPORTED'}, H2={'SUPPORTED' if h2 else 'NOT SUPPORTED'}")


if __name__ == "__main__":
    main()
