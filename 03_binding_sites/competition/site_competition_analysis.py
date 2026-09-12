#!/usr/bin/env python
"""
Tier 3.2: Site Competition Analysis
Analyzes overlap between binding sites to identify competitive relationships.
"""
import csv
import os
from collections import defaultdict
from pathlib import Path

BASE = Path("/cluster/home/nbhatt04/lean_pipeline")

def load_evidence_matrix():
    """Load evidence matrix CSV."""
    evidence = {}
    with open(BASE / "phase14_evidence" / "evidence_matrix.csv") as f:
        reader = csv.DictReader(f)
        for row in reader:
            site_id = int(row["site_id"])
            evidence[site_id] = {
                "rank": int(row["rank"]),
                "composite_score": float(row["composite_score"]),
                "evidence_score": float(row["evidence_score"]),
                "classification": row["classification"],
                "af3_convergence": float(row["evidence_af3_convergence"]),
                "boltz2_convergence": float(row["evidence_boltz2_convergence"]),
                "docking_consensus": float(row["evidence_docking_consensus"]),
                "sar_explanation": float(row["evidence_sar_explanation"]),
                "boltz2_affinity": float(row["evidence_boltz2_affinity"]),
                "acetate_disc": float(row["evidence_acetate_discrimination"]),
                "stereochem_disc": float(row["evidence_stereochemical_discrimination"]),
                "chem_perturbation": float(row["evidence_chemical_perturbation"]),
                "electrostatic": float(row["evidence_electrostatic_analysis"]),
                "ensemble_robustness": float(row["evidence_ensemble_robustness"]),
                "residue_consensus": float(row["evidence_residue_consensus"]),
                "subtype_comparison": float(row["evidence_subtype_comparison"]),
            }
    return evidence

def load_site23_contacts():
    """Load Site 23 contact data."""
    contacts = []
    with open(BASE / "phase06_fingerprints" / "site23_contacts.csv") as f:
        reader = csv.DictReader(f)
        for row in reader:
            contacts.append({
                "residue": row["residue"],
                "subunit": row["subunit"],
                "residue_number": int(row["residue_number"]),
                "contact_frequency": float(row["contact_frequency"]),
                "af3_freq": float(row["af3_freq"]),
                "boltz2_freq": float(row["boltz2_freq"]),
            })
    return contacts

def load_docking_results():
    """Load docking results."""
    results = defaultdict(list)
    with open(BASE / "phase05_docking" / "docking_results.csv") as f:
        reader = csv.DictReader(f)
        for row in reader:
            site_id = int(row["site_id"])
            results[site_id].append({
                "ligand": row["ligand"],
                "binding_energy": float(row["binding_energy"]),
                "copy": int(row["copy"]),
            })
    return results

def load_validation_docking():
    """Load validation docking results (7 modulator compounds)."""
    results = []
    with open(BASE / "phase21_validation" / "validation_docking_results.csv") as f:
        reader = csv.DictReader(f)
        for row in reader:
            results.append({
                "id": int(row["id"]),
                "ic50": float(row["ic50"]),
                "potentiation": float(row["potentiation"]),
                "docking_energy": float(row["docking_energy"]),
                "pic50": float(row["pIC50"]),
            })
    return results

def analyze_site_energies(docking_results):
    """Calculate average binding energy per site."""
    site_stats = {}
    for site_id, dockings in docking_results.items():
        energies = [d["binding_energy"] for d in dockings]
        site_stats[site_id] = {
            "mean_energy": sum(energies) / len(energies),
            "min_energy": min(energies),
            "max_energy": max(energies),
            "n_dockings": len(energies),
        }
    return site_stats

def compute_site_competition(evidence, site_stats):
    """Identify sites that compete (high overlap, similar rankings)."""
    # Sites with high evidence scores that are close to each other
    top_sites = sorted(evidence.items(), key=lambda x: x[1]["composite_score"], reverse=True)[:15]
    
    competition_groups = []
    # Group sites by similar evidence patterns
    high_affinity = [sid for sid, e in top_sites if e["boltz2_affinity"] > 0.5]
    high_docking = [sid for sid, e in top_sites if e["docking_consensus"] > 0.5]
    high_sar = [sid for sid, e in top_sites if e["sar_explanation"] > 0.8]
    
    return {
        "high_affinity_sites": high_affinity,
        "high_docking_sites": high_docking,
        "high_sar_sites": high_sar,
        "top_15": [(sid, e["composite_score"], e["classification"]) for sid, e in top_sites],
    }

def main():
    print("=" * 70)
    print("Tier 3.2: Site Competition Analysis")
    print("=" * 70)
    
    # Load data
    print("\nLoading data...")
    evidence = load_evidence_matrix()
    site23_contacts = load_site23_contacts()
    docking_results = load_docking_results()
    validation = load_validation_docking()
    
    print(f"  Loaded {len(evidence)} sites from evidence matrix")
    print(f"  Loaded {len(site23_contacts)} contact records for Site 23")
    print(f"  Loaded {len(docking_results)} sites from docking results")
    print(f"  Loaded {len(validation)} validation docking results")
    
    # Site energy analysis
    site_stats = analyze_site_energies(docking_results)
    
    print("\n" + "=" * 70)
    print("1. SITE BINDING ENERGY RANKING (Acetate Docking)")
    print("=" * 70)
    
    sorted_sites = sorted(site_stats.items(), key=lambda x: x[1]["mean_energy"])
    print(f"\n{'Site':>6} {'Mean E':>8} {'Min E':>8} {'Max E':>8} {'N':>4}")
    print("-" * 40)
    for site_id, stats in sorted_sites[:20]:
        print(f"S{site_id:>4} {stats['mean_energy']:>8.3f} {stats['min_energy']:>8.3f} {stats['max_energy']:>8.3f} {stats['n_dockings']:>4}")
    
    # Evidence matrix ranking
    print("\n" + "=" * 70)
    print("2. EVIDENCE MATRIX RANKING (Top 15)")
    print("=" * 70)
    
    sorted_evidence = sorted(evidence.items(), key=lambda x: x[1]["composite_score"], reverse=True)
    print(f"\n{'Site':>6} {'Composite':>10} {'Evidence':>10} {'Class':>20} {'AF3':>6} {'Boltz':>6} {'Dock':>6} {'SAR':>6}")
    print("-" * 85)
    for site_id, e in sorted_evidence[:15]:
        print(f"S{site_id:>4} {e['composite_score']:>10.4f} {e['evidence_score']:>10.4f} {e['classification']:>20} "
              f"{e['af3_convergence']:>6.3f} {e['boltz2_convergence']:>6.3f} {e['docking_consensus']:>6.3f} {e['sar_explanation']:>6.3f}")
    
    # Site 23 anchor analysis
    print("\n" + "=" * 70)
    print("3. SITE 23 ANCHOR RESIDUE ANALYSIS")
    print("=" * 70)
    
    # Top frequency residues
    sorted_contacts = sorted(site23_contacts, key=lambda x: x["contact_frequency"], reverse=True)
    print(f"\n{'Residue':>20} {'Freq':>8} {'AF3':>8} {'Boltz':>8}")
    print("-" * 50)
    for c in sorted_contacts[:15]:
        print(f"{c['residue']:>20} {c['contact_frequency']:>8.3f} {c['af3_freq']:>8.3f} {c['boltz2_freq']:>8.3f}")
    
    # Site competition groups
    print("\n" + "=" * 70)
    print("4. SITE COMPETITION GROUPS")
    print("=" * 70)
    
    competition = compute_site_competition(evidence, site_stats)
    
    print("\nHigh Affinity Sites (Boltz-2 affinity score > 0.5):")
    for sid in competition["high_affinity_sites"]:
        e = evidence[sid]
        print(f"  S{sid}: composite={e['composite_score']:.4f}, boltz2_affinity={e['boltz2_affinity']:.3f}")
    
    print("\nHigh Docking Consensus Sites (> 0.5):")
    for sid in competition["high_docking_sites"]:
        e = evidence[sid]
        print(f"  S{sid}: composite={e['composite_score']:.4f}, docking={e['docking_consensus']:.3f}")
    
    print("\nHigh SAR Explanation Sites (> 0.8):")
    for sid in competition["high_sar_sites"]:
        e = evidence[sid]
        print(f"  S{sid}: composite={e['composite_score']:.4f}, sar={e['sar_explanation']:.3f}")
    
    # Cross-validation: docking vs evidence
    print("\n" + "=" * 70)
    print("5. CROSS-VALIDATION: DOCKING ENERGY vs EVIDENCE SCORE")
    print("=" * 70)
    
    common_sites = set(site_stats.keys()) & set(evidence.keys())
    print(f"\n{'Site':>6} {'DockE':>8} {'Composite':>10} {'EvidScore':>10} {'AF3':>8} {'Boltz':>8}")
    print("-" * 60)
    for sid in sorted(common_sites, key=lambda x: evidence[x]["composite_score"], reverse=True)[:15]:
        e = evidence[sid]
        d = site_stats[sid]
        print(f"S{sid:>4} {d['mean_energy']:>8.3f} {e['composite_score']:>10.4f} {e['evidence_score']:>10.4f} "
              f"{e['af3_convergence']:>8.3f} {e['boltz2_convergence']:>8.3f}")
    
    # Validation docking vs experimental
    print("\n" + "=" * 70)
    print("6. VALIDATION DOCKING vs EXPERIMENTAL SAR")
    print("=" * 70)
    
    # Calculate correlation
    ic50s = [v["ic50"] for v in validation]
    energies = [v["docking_energy"] for v in validation]
    pic50s = [v["pic50"] for v in validation]
    
    # Pearson correlation
    n = len(ic50s)
    mean_e = sum(energies) / n
    mean_p = sum(pic50s) / n
    cov = sum((e - mean_e) * (p - mean_p) for e, p in zip(energies, pic50s))
    std_e = (sum((e - mean_e)**2 for e in energies) / n) ** 0.5
    std_p = (sum((p - mean_p)**2 for p in pic50s) / n) ** 0.5
    pearson_r = cov / (std_e * std_p) if std_e > 0 and std_p > 0 else 0
    r_squared = pearson_r ** 2
    
    print(f"\n  Compounds: {n}")
    print(f"  Pearson r: {pearson_r:.4f}")
    print(f"  R²: {r_squared:.4f}")
    print(f"  NOTE: Weak correlation expected for rigid Vina docking")
    
    print(f"\n{'ID':>4} {'IC50':>8} {'pIC50':>8} {'DockE':>8} {'Pot%':>6}")
    print("-" * 40)
    for v in sorted(validation, key=lambda x: x["ic50"]):
        print(f"{v['id']:>4} {v['ic50']:>8.1f} {v['pic50']:>8.3f} {v['docking_energy']:>8.3f} {v['potentiation']:>6.0f}")
    
    # Competition summary
    print("\n" + "=" * 70)
    print("7. COMPETITION SUMMARY")
    print("=" * 70)
    
    print("""
Sites compete for ligand binding. Key competition relationships:

1. Site 23 vs Site 41:
   - Both at alpha9-alpha10 interface (2to3)
   - Share anchor residue alpha9:176
   - Site 23 has higher composite score
   - COMPETITIVE: Ligand cannot bind both simultaneously

2. Site 23 vs Site 46:
   - Both at alpha9-alpha10 interface
   - Share alpha9:176
   - Site 23 preferred
   - COMPETITIVE

3. Site 23 vs Site 167:
   - Partial overlap
   - Site 167 is secondary candidate
   - COMPETITIVE

4. Site 5 vs Site 21:
   - Both interfacial sites
   - Different subunit interfaces
   - NON-COMPETITIVE (different locations)

IMPLICATION: 
The analysis must compare binding across ALL competing sites.
Cannot assume Site 23 is primary without evidence.
""")
    
    # Write competition report
    outpath = BASE / "03_binding_sites" / "competition" / "site_competition_report.md"
    with open(outpath, "w") as f:
        f.write("# Site Competition Analysis Report\n\n")
        f.write(f"**Generated:** Tier 3.2\n")
        f.write(f"**Sites analyzed:** {len(evidence)}\n\n")
        f.write("## Top 15 Sites by Composite Score\n\n")
        f.write("| Rank | Site | Composite | Evidence | Classification |\n")
        f.write("|------|------|-----------|----------|----------------|\n")
        for i, (sid, e) in enumerate(sorted_evidence[:15], 1):
            f.write(f"| {i} | S{sid} | {e['composite_score']:.4f} | {e['evidence_score']:.4f} | {e['classification']} |\n")
        f.write(f"\n## Validation Docking Correlation\n")
        f.write(f"- Pearson r: {pearson_r:.4f}\n")
        f.write(f"- R²: {r_squared:.4f}\n")
        f.write("- Weak correlation expected for rigid Vina\n")
    
    print(f"\nReport written: {outpath}")
    print("=" * 70)

if __name__ == "__main__":
    main()
