#!/usr/bin/env python3
"""
Phase K-L: Boltz-2 Analysis and Multi-Method Convergence
Independent validation layer and convergence matrix.
"""
import csv
import json
import os
from pathlib import Path
from collections import defaultdict

PIPELINE_DIR = Path("/cluster/home/nbhatt04/lean_pipeline")
ALLOSTERY_DIR = Path("/cluster/scratch/nbhatt04/allostery")
BOLTZ2_DIR = ALLOSTERY_DIR / "boltz2" / "outputs"
PUB_TABLES = Path("/cluster/home/nbhatt04/alpha9alpha10_publication/11_tables")

OUTPUT_DIR = PIPELINE_DIR / "11_boltz2_analysis"
OUTPUT_DIR.mkdir(exist_ok=True)

def load_boltz2_confidences():
    """Load Boltz-2 confidence metrics."""
    results = []
    
    for result_dir in sorted(BOLTZ2_DIR.iterdir()):
        if not result_dir.is_dir():
            continue
        
        name = result_dir.name
        
        # Find confidence JSON files
        for json_file in result_dir.rglob('*confidence*.json'):
            try:
                with open(json_file, 'r') as f:
                    conf = json.load(f)
                
                results.append({
                    'name': name,
                    'confidence_score': conf.get('confidence_score', 0),
                    'ptm': conf.get('ptm', 0),
                    'iptm': conf.get('iptm', 0),
                    'ligand_iptm': conf.get('ligand_iptm', 0),
                    'protein_iptm': conf.get('protein_iptm', 0),
                    'complex_plddt': conf.get('complex_plddt', 0),
                })
            except:
                pass
    
    return results

def load_convergence_matrix():
    """Load the convergence matrix from publication."""
    conv_file = PUB_TABLES / "CONVERGENCE_MATRIX.csv"
    if conv_file.exists():
        with open(conv_file, 'r') as f:
            reader = csv.DictReader(f)
            return list(reader)
    return []

def load_af3_results():
    """Load AF3 results."""
    af3_file = PUB_TABLES / "Table4_af3_results.csv"
    if af3_file.exists():
        with open(af3_file, 'r') as f:
            reader = csv.DictReader(f)
            return list(reader)
    return []

def analyze_boltz2():
    """Analyze Boltz-2 results."""
    print("\n--- BOLTZ-2 ANALYSIS ---")
    
    confidences = load_boltz2_confidences()
    print(f"  Total Boltz-2 models with confidence: {len(confidences)}")
    
    if confidences:
        # Statistics
        conf_scores = [c['confidence_score'] for c in confidences]
        iptm_scores = [c['iptm'] for c in confidences]
        lig_iptm = [c['ligand_iptm'] for c in confidences]
        
        print(f"\n  Confidence score statistics:")
        print(f"    Mean: {sum(conf_scores)/len(conf_scores):.4f}")
        print(f"    Min: {min(conf_scores):.4f}")
        print(f"    Max: {max(conf_scores):.4f}")
        
        print(f"\n  iPTM statistics:")
        print(f"    Mean: {sum(iptm_scores)/len(iptm_scores):.4f}")
        print(f"    Min: {min(iptm_scores):.4f}")
        print(f"    Max: {max(iptm_scores):.4f}")
        
        print(f"\n  Ligand iPTM statistics:")
        print(f"    Mean: {sum(lig_iptm)/len(lig_iptm):.4f}")
        print(f"    Min: {min(lig_iptm):.4f}")
        print(f"    Max: {max(lig_iptm):.4f}")
    
    # Group by ligand type
    by_ligand = defaultdict(list)
    for c in confidences:
        # Parse ligand from name
        name = c['name']
        parts = name.replace('boltz_results_', '').split('_')
        ligand = '_'.join(parts[2:-1])  # Remove stoichiometry and seed
        by_ligand[ligand].append(c)
    
    print(f"\n  By ligand type:")
    for ligand, models in by_ligand.items():
        mean_conf = sum(m['confidence_score'] for m in models) / len(models)
        print(f"    {ligand}: n={len(models)}, mean_conf={mean_conf:.4f}")
    
    return confidences

def analyze_convergence():
    """Analyze cross-method convergence."""
    print("\n--- MULTI-METHOD CONVERGENCE ---")
    
    # Load convergence data
    convergence_file = PIPELINE_DIR / "phase04_convergence" / "convergence_map.csv"
    convergence = []
    if convergence_file.exists():
        with open(convergence_file, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                convergence.append(row)
    
    print(f"  Convergence entries: {len(convergence)}")
    
    # Classify convergence
    strong = [c for c in convergence if c.get('convergence_class') == 'STRONG']
    moderate = [c for c in convergence if c.get('convergence_class') == 'MODERATE']
    single = [c for c in convergence if c.get('convergence_class') == 'SINGLE_METHOD']
    
    print(f"  Strong convergence: {len(strong)}")
    print(f"  Moderate convergence: {len(moderate)}")
    print(f"  Single method: {len(single)}")
    
    # Load cross-method report
    cross_file = PUB_TABLES / "CROSS_METHOD_REPORT.txt"
    if cross_file.exists():
        with open(cross_file, 'r') as f:
            content = f.read()
        # Parse site clusters
        site_clusters = []
        for line in content.split('\n'):
            if line.startswith('--- site'):
                parts = line.split()
                site_id = parts[2].rstrip(':')
                n_models = int(parts[4])
                site_clusters.append({'site_id': site_id, 'n_models': n_models})
        
        print(f"\n  Cross-method site clusters: {len(site_clusters)}")
        for sc in site_clusters[:10]:
            print(f"    Site {sc['site_id']}: {sc['n_models']} models")
    
    return convergence

def build_convergence_matrix():
    """Build the multi-method convergence matrix."""
    print("\n--- BUILDING CONVERGENCE MATRIX ---")
    
    # Load all data sources
    site_ranking = []
    ranking_file = PUB_TABLES / "FINAL_SITE_RANKING.csv"
    if ranking_file.exists():
        with open(ranking_file, 'r') as f:
            reader = csv.DictReader(f)
            site_ranking = list(reader)
    
    # Build matrix for top sites
    matrix = []
    for site in site_ranking[:15]:
        site_id = site.get('site_id', '')
        
        # AF3 support
        n_af3 = int(site.get('n_af3', 0))
        af3_support = 'STRONG' if n_af3 >= 20 else ('MODERATE' if n_af3 >= 5 else 'WEAK')
        
        # Boltz-2 support
        n_boltz2 = int(site.get('n_boltz2', 0))
        boltz2_support = 'STRONG' if n_boltz2 >= 10 else ('MODERATE' if n_boltz2 >= 3 else 'WEAK')
        
        # Docking support (from existing analysis)
        # Would need to cross-reference with docking results
        
        # SAR support
        stereo_test = site.get('stereo_test', '')
        acetate_test = site.get('acetate_test', '')
        analogue_test = site.get('analogue_test', '')
        
        sar_support = 'STRONG' if ('both' in str(stereo_test) and 'both' in str(analogue_test)) else \
                      ('MODERATE' if 'both' in str(stereo_test) or 'both' in str(analogue_test) else 'WEAK')
        
        # Overall convergence
        methods = [af3_support, boltz2_support, sar_support]
        strong_count = methods.count('STRONG')
        overall = 'STRONG' if strong_count >= 2 else ('MODERATE' if strong_count >= 1 else 'WEAK')
        
        matrix.append({
            'site_id': site_id,
            'af3_support': af3_support,
            'af3_n_models': n_af3,
            'boltz2_support': boltz2_support,
            'boltz2_n_models': n_boltz2,
            'sar_support': sar_support,
            'overall_convergence': overall,
            'site_score': float(site.get('SiteScore', 0)),
            'evidence_level': site.get('evidence_level', ''),
        })
    
    print(f"  Convergence matrix built for {len(matrix)} sites")
    
    # Summary
    strong_sites = [m for m in matrix if m['overall_convergence'] == 'STRONG']
    moderate_sites = [m for m in matrix if m['overall_convergence'] == 'MODERATE']
    weak_sites = [m for m in matrix if m['overall_convergence'] == 'WEAK']
    
    print(f"\n  Overall convergence:")
    print(f"    STRONG: {len(strong_sites)} sites")
    print(f"    MODERATE: {len(moderate_sites)} sites")
    print(f"    WEAK: {len(weak_sites)} sites")
    
    return matrix

def main():
    print("=" * 70)
    print("PHASE K-L: BOLTZ-2 & CONVERGENCE ANALYSIS")
    print("=" * 70)
    
    boltz2_results = analyze_boltz2()
    convergence = analyze_convergence()
    matrix = build_convergence_matrix()
    
    # Save results
    report = {
        'boltz2_models': len(boltz2_results),
        'convergence_entries': len(convergence),
        'convergence_matrix': matrix,
    }
    
    with open(OUTPUT_DIR / "boltz2_convergence_report.json", 'w') as f:
        json.dump(report, f, indent=2, default=str)
    
    # Save convergence matrix
    if matrix:
        with open(OUTPUT_DIR / "convergence_matrix.csv", 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=matrix[0].keys())
            writer.writeheader()
            writer.writerows(matrix)
    
    print("\n" + "=" * 70)
    print("PHASE K-L COMPLETE")
    print("=" * 70)

if __name__ == "__main__":
    main()
