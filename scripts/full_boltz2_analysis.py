#!/usr/bin/env python3
"""
Full Boltz-2 Analysis: Process all existing Boltz-2 outputs.
Confidence metrics, convergence with AF3, SAR consistency.
"""
import csv
import json
import os
from pathlib import Path
from collections import defaultdict

PIPELINE_DIR = Path("/cluster/home/nbhatt04/lean_pipeline")
SCRATCH = Path("/cluster/scratch/nbhatt04/allostery")
BOLTZ2_DIR = SCRATCH / "boltz2" / "outputs"
OUTPUT_DIR = PIPELINE_DIR / "11_boltz2_analysis"
OUTPUT_DIR.mkdir(exist_ok=True)

def load_all_boltz2_confidences():
    """Load confidence metrics from all Boltz-2 outputs."""
    results = []
    
    for result_dir in sorted(BOLTZ2_DIR.iterdir()):
        if not result_dir.is_dir():
            continue
        
        name = result_dir.name
        
        # Find all confidence JSON files
        for json_file in result_dir.rglob('*confidence*.json'):
            try:
                with open(json_file, 'r') as f:
                    conf = json.load(f)
                
                # Parse name components
                parts = name.replace('boltz_results_', '').split('_')
                stoichiometry = parts[0] + '_' + parts[1] if len(parts) >= 2 else 'unknown'
                ligand = '_'.join(parts[2:-1]) if len(parts) >= 4 else 'unknown'
                seed = parts[-1] if len(parts) >= 3 else 'unknown'
                
                results.append({
                    'name': name,
                    'stoichiometry': stoichiometry,
                    'ligand': ligand,
                    'seed': seed,
                    'confidence_score': conf.get('confidence_score', 0),
                    'ptm': conf.get('ptm', 0),
                    'iptm': conf.get('iptm', 0),
                    'ligand_iptm': conf.get('ligand_iptm', 0),
                    'protein_iptm': conf.get('protein_iptm', 0),
                    'complex_plddt': conf.get('complex_plddt', 0),
                    'n_protein_chains': len([c for c in conf.get('chains', []) if c.get('entity_type') == 'protein']),
                    'n_ligand_chains': len([c for c in conf.get('chains', []) if c.get('entity_type') == 'ligand']),
                })
            except Exception as e:
                pass
    
    return results

def analyze_confidence_distributions(results):
    """Analyze confidence score distributions."""
    metrics = ['confidence_score', 'ptm', 'iptm', 'ligand_iptm', 'protein_iptm', 'complex_plddt']
    
    distributions = {}
    for metric in metrics:
        values = [r[metric] for r in results if r[metric] > 0]
        if values:
            distributions[metric] = {
                'n': len(values),
                'mean': sum(values)/len(values),
                'min': min(values),
                'max': max(values),
                'median': sorted(values)[len(values)//2],
            }
    
    return distributions

def analyze_by_ligand(results):
    """Analyze results grouped by ligand type."""
    by_ligand = defaultdict(list)
    for r in results:
        by_ligand[r['ligand']].append(r)
    
    summary = {}
    for ligand, models in by_ligand.items():
        confs = [m['confidence_score'] for m in models if m['confidence_score'] > 0]
        iptms = [m['ligand_iptm'] for m in models if m['ligand_iptm'] > 0]
        
        summary[ligand] = {
            'n_models': len(models),
            'mean_confidence': sum(confs)/len(confs) if confs else 0,
            'mean_ligand_iptm': sum(iptms)/len(iptms) if iptms else 0,
            'seeds': list(set(m['seed'] for m in models)),
        }
    
    return summary

def analyze_convergence_with_af3(results):
    """Analyze convergence between Boltz-2 and AF3 predictions."""
    # Load AF3 results
    af3_file = PIPELINE_DIR / "phase04_convergence" / "convergence_map.csv"
    convergence = []
    if af3_file.exists():
        with open(af3_file, 'r') as f:
            reader = csv.DictReader(f)
            convergence = list(reader)
    
    # Cross-reference
    convergent = [c for c in convergence if c.get('convergence_class') == 'STRONG']
    
    return {
        'total_convergence_entries': len(convergence),
        'strong_convergence': len(convergent),
        'convergent_sites': [c.get('site', '') for c in convergent],
    }

def main():
    print("="*70)
    print("FULL BOLTZ-2 ANALYSIS")
    print("="*70)
    
    # Load all Boltz-2 confidences
    results = load_all_boltz2_confidences()
    print(f"Loaded {len(results)} Boltz-2 models with confidence data")
    
    # Confidence distributions
    distributions = analyze_confidence_distributions(results)
    print(f"\nConfidence distributions:")
    for metric, stats in distributions.items():
        print(f"  {metric}: n={stats['n']}, mean={stats['mean']:.4f}, "
              f"range=[{stats['min']:.4f}, {stats['max']:.4f}]")
    
    # By ligand
    by_ligand = analyze_by_ligand(results)
    print(f"\nBy ligand type:")
    for ligand, stats in by_ligand.items():
        print(f"  {ligand}: n={stats['n_models']}, "
              f"mean_conf={stats['mean_confidence']:.4f}, "
              f"mean_lig_iptm={stats['mean_ligand_iptm']:.4f}")
    
    # Convergence with AF3
    convergence = analyze_convergence_with_af3(results)
    print(f"\nConvergence with AF3:")
    print(f"  Total entries: {convergence['total_convergence_entries']}")
    print(f"  Strong convergence: {convergence['strong_convergence']}")
    
    # SAR consistency
    active_ligands = ['LASC_ANION', 'OETHYL_USER_SUPPLIED', 'OETHYL_CORRECTED_3O']
    inactive_ligands = ['ACETATE_ANION', 'DASC_ANION']
    
    active_confs = [r['confidence_score'] for r in results 
                   if r['ligand'] in active_ligands and r['confidence_score'] > 0]
    inactive_confs = [r['confidence_score'] for r in results 
                     if r['ligand'] in inactive_ligands and r['confidence_score'] > 0]
    
    sar_consistency = {
        'active_mean_confidence': sum(active_confs)/len(active_confs) if active_confs else 0,
        'inactive_mean_confidence': sum(inactive_confs)/len(inactive_confs) if inactive_confs else 0,
        'separation': (sum(active_confs)/len(active_confs) if active_confs else 0) > 
                      (sum(inactive_confs)/len(inactive_confs) if inactive_confs else 0),
    }
    
    print(f"\nSAR consistency:")
    print(f"  Active mean confidence: {sar_consistency['active_mean_confidence']:.4f}")
    print(f"  Inactive mean confidence: {sar_consistency['inactive_mean_confidence']:.4f}")
    print(f"  Separation: {sar_consistency['separation']}")
    
    # Save results
    report = {
        'total_models': len(results),
        'confidence_distributions': distributions,
        'by_ligand': by_ligand,
        'convergence_with_af3': convergence,
        'sar_consistency': sar_consistency,
    }
    
    with open(OUTPUT_DIR / "full_boltz2_analysis.json", 'w') as f:
        json.dump(report, f, indent=2, default=str)
    
    # Save detailed results
    if results:
        fieldnames = list(results[0].keys())
        with open(OUTPUT_DIR / "boltz2_all_confidences.csv", 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(results)
    
    print(f"\nBoltz-2 analysis complete: {OUTPUT_DIR / 'full_boltz2_analysis.json'}")

if __name__ == "__main__":
    main()
