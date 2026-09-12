#!/usr/bin/env python3
"""
Phase I-J: Docking Analysis and SAR Validation
Analyzes docking results, interaction fingerprints,
and validates against experimental SAR.
"""
import csv
import json
import os
import glob
from pathlib import Path
from collections import defaultdict

PIPELINE_DIR = Path("/cluster/home/nbhatt04/lean_pipeline")
OUTPUT_DIR = PIPELINE_DIR / "09_docking"
OUTPUT_DIR.mkdir(exist_ok=True)

def load_docking_results():
    """Load existing docking results."""
    docking_file = PIPELINE_DIR / "phase05_docking" / "docking_results.csv"
    results = []
    if docking_file.exists():
        with open(docking_file, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                results.append(row)
    return results

def load_ifp_data():
    """Load interaction fingerprint data."""
    ifp_file = PIPELINE_DIR / "phase06_fingerprints" / "ifp_summary.csv"
    ifps = []
    if ifp_file.exists():
        with open(ifp_file, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                ifps.append(row)
    return ifps

def load_sar_rules():
    """Load SAR rules."""
    sar_file = PIPELINE_DIR / "phase07_sar_model" / "sar_rules.csv"
    rules = []
    if sar_file.exists():
        with open(sar_file, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                rules.append(row)
    return rules

def analyze_docking():
    """Analyze docking results across all compounds and sites."""
    print("\n--- DOCKING ANALYSIS ---")
    
    results = load_docking_results()
    print(f"  Total docking results: {len(results)}")
    
    if not results:
        print("  No docking results found")
        return {}
    
    # Group by ligand
    by_ligand = defaultdict(list)
    for r in results:
        ligand = r.get('ligand', 'unknown')
        by_ligand[ligand].append(r)
    
    print(f"  By ligand: { {k: len(v) for k, v in by_ligand.items()} }")
    
    # Analyze binding energies
    energies = []
    for r in results:
        try:
            energy = float(r.get('binding_energy', 0))
            energies.append(energy)
        except:
            pass
    
    if energies:
        print(f"\n  Binding energy statistics:")
        print(f"    Mean: {sum(energies)/len(energies):.3f} kcal/mol")
        print(f"    Min: {min(energies):.3f} kcal/mol")
        print(f"    Max: {max(energies):.3f} kcal/mol")
    
    # Active vs inactive comparison
    print("\n  Active vs Inactive Docking:")
    active_ligands = ['ascorbate', 'O-ethyl_ascorbate']  # Known active
    inactive_ligands = ['acetate']  # Known inactive
    
    active_energies = []
    inactive_energies = []
    
    for r in results:
        ligand = r.get('ligand', '')
        try:
            energy = float(r.get('binding_energy', 0))
            if ligand in active_ligands:
                active_energies.append(energy)
            elif ligand in inactive_ligands:
                inactive_energies.append(energy)
        except:
            pass
    
    if active_energies:
        print(f"    Active: mean={sum(active_energies)/len(active_energies):.3f}, "
              f"n={len(active_energies)}")
    if inactive_energies:
        print(f"    Inactive: mean={sum(inactive_energies)/len(inactive_energies):.3f}, "
              f"n={len(inactive_energies)}")
    
    return {
        'total_results': len(results),
        'by_ligand': {k: len(v) for k, v in by_ligand.items()},
        'energy_stats': {
            'mean': sum(energies)/len(energies) if energies else 0,
            'min': min(energies) if energies else 0,
            'max': max(energies) if energies else 0,
        }
    }

def analyze_ifp():
    """Analyze interaction fingerprints."""
    print("\n--- INTERACTION FINGERPRINT ANALYSIS ---")
    
    ifps = load_ifp_data()
    print(f"  Total IFP entries: {len(ifps)}")
    
    # Load individual site contact files
    contact_files = sorted(glob.glob(str(PIPELINE_DIR / "phase06_fingerprints" / "site*_contacts.csv")))
    print(f"  Site contact files: {len(contact_files)}")
    
    # Analyze contact patterns
    all_contacts = defaultdict(int)
    for cf in contact_files:
        try:
            with open(cf, 'r') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    residue = row.get('residue', '')
                    if residue:
                        all_contacts[residue] += 1
        except:
            pass
    
    # Top contact residues
    top_contacts = sorted(all_contacts.items(), key=lambda x: x[1], reverse=True)[:20]
    print(f"\n  Top 20 contact residues:")
    for residue, count in top_contacts:
        print(f"    {residue}: {count} contacts")
    
    return {
        'total_ifp_entries': len(ifps),
        'contact_files': len(contact_files),
        'top_contacts': top_contacts,
    }

def analyze_sar_validation():
    """Validate SAR against docking results."""
    print("\n--- SAR VALIDATION ---")
    
    # Load SAR rules
    rules = load_sar_rules()
    print(f"  SAR rules: {len(rules)}")
    
    # Load active compounds
    active_csv = PIPELINE_DIR / "01_qc" / "active_compounds.csv"
    active = []
    if active_csv.exists():
        with open(active_csv, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                active.append(row)
    
    print(f"  Active compounds: {len(active)}")
    
    # Load docking results
    docking = load_docking_results()
    
    # Cross-reference
    print("\n  Active compound docking summary:")
    for a in active:
        comp_id = a.get('Identifier', '')
        smiles = a.get('Smiles', '')
        activity = a.get('Activity (uM)', '')
        
        # Find docking results for this compound
        # Note: docking uses ligand names, not IDs
        print(f"    Compound {comp_id}: {activity} uM")
    
    return {
        'sar_rules': len(rules),
        'active_compounds': len(active),
    }

def main():
    print("=" * 70)
    print("PHASE I-J: DOCKING & SAR VALIDATION")
    print("=" * 70)
    
    docking_stats = analyze_docking()
    ifp_stats = analyze_ifp()
    sar_stats = analyze_sar_validation()
    
    # Save report
    report = {
        'docking': docking_stats,
        'interaction_fingerprints': ifp_stats,
        'sar_validation': sar_stats,
    }
    
    with open(OUTPUT_DIR / "docking_sar_report.json", 'w') as f:
        json.dump(report, f, indent=2, default=str)
    
    print("\n" + "=" * 70)
    print("PHASE I-J COMPLETE")
    print("=" * 70)

if __name__ == "__main__":
    main()
