#!/usr/bin/env python3
"""
Phase A: Dataset Quality Control
Verifies compound counts, validates SMILES, normalizes structures,
and generates clean computational dataset.
"""
import csv
import json
import os
import sys
from pathlib import Path
from collections import Counter

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, AllChem, rdMolDescriptors
    from rdkit.Chem import Draw
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False
    print("WARNING: RDKit not available. Using basic QC only.")

PIPELINE_DIR = Path("/cluster/home/nbhatt04/lean_pipeline")
INPUT_CSV = PIPELINE_DIR / "modulator-dataset-a9a10.csv"
OUTPUT_DIR = PIPELINE_DIR / "01_qc"
OUTPUT_DIR.mkdir(exist_ok=True)

def load_dataset(csv_path):
    """Load the raw dataset."""
    compounds = []
    with open(csv_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            compounds.append(row)
    return compounds

def validate_smiles(smiles):
    """Validate SMILES string using RDKit."""
    if not HAS_RDKIT:
        return True, "RDKit not available"
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False, "Invalid SMILES"
        return True, "Valid"
    except Exception as e:
        return False, str(e)

def normalize_smiles(smiles):
    """Normalize SMILES to canonical form."""
    if not HAS_RDKIT:
        return smiles
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return smiles
        return Chem.MolToSmiles(mol, canonical=True)
    except:
        return smiles

def compute_descriptors(smiles):
    """Compute molecular descriptors."""
    if not HAS_RDKIT:
        return {}
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {}
        return {
            'mw': round(Descriptors.MolWt(mol), 2),
            'logp': round(Descriptors.MolLogP(mol), 2),
            'hbd': rdMolDescriptors.CalcNumHBD(mol),
            'hba': rdMolDescriptors.CalcNumHBA(mol),
            'rotatable_bonds': rdMolDescriptors.CalcNumRotatableBonds(mol),
            'tpsa': round(Descriptors.TPSA(mol), 2),
            'heavy_atoms': mol.GetNumHeavyAtoms(),
            'formal_charge': Chem.GetFormalCharge(mol),
            'num_rings': rdMolDescriptors.CalcNumRings(mol),
            'num_aromatic_rings': rdMolDescriptors.CalcNumAromaticRings(mol),
        }
    except:
        return {}

def check_duplicates(compounds):
    """Check for duplicate structures and IDs."""
    smiles_list = [c['Smiles'] for c in compounds]
    id_list = [c['Identifier'] for c in compounds]
    
    smiles_counts = Counter(smiles_list)
    id_counts = Counter(id_list)
    
    dup_smiles = {s: cnt for s, cnt in smiles_counts.items() if cnt > 1}
    dup_ids = {i: cnt for i, cnt in id_counts.items() if cnt > 1}
    
    return dup_smiles, dup_ids

def classify_activity(compounds):
    """Classify compounds as active/inactive based on experimental data."""
    active = []
    inactive = []
    
    for c in compounds:
        activity = float(c.get('Activity (uM)', 0) or 0)
        potentiation = float(c.get('%Potentiation', 0) or 0)
        
        if activity > 0 and potentiation > 0:
            active.append(c)
        else:
            inactive.append(c)
    
    return active, inactive

def main():
    print("=" * 70)
    print("PHASE A: DATASET QUALITY CONTROL")
    print("=" * 70)
    
    # Load raw dataset
    print("\n1. Loading raw dataset...")
    compounds = load_dataset(INPUT_CSV)
    print(f"   Total compounds loaded: {len(compounds)}")
    
    # Check duplicates
    print("\n2. Checking for duplicates...")
    dup_smiles, dup_ids = check_duplicates(compounds)
    if dup_smiles:
        print(f"   WARNING: Duplicate SMILES found: {dup_smiles}")
    else:
        print("   No duplicate SMILES found")
    if dup_ids:
        print(f"   WARNING: Duplicate IDs found: {dup_ids}")
    else:
        print("   No duplicate IDs found")
    
    # Validate SMILES
    print("\n3. Validating SMILES...")
    invalid_smiles = []
    for c in compounds:
        valid, msg = validate_smiles(c['Smiles'])
        if not valid:
            invalid_smiles.append((c['Identifier'], c['Smiles'], msg))
            print(f"   INVALID: Compound {c['Identifier']}: {msg}")
    if not invalid_smiles:
        print("   All SMILES valid")
    
    # Classify activity
    print("\n4. Classifying activity...")
    active, inactive = classify_activity(compounds)
    print(f"   Active compounds: {len(active)}")
    print(f"   Inactive compounds: {len(inactive)}")
    print(f"   Total: {len(active) + len(inactive)}")
    
    # Print active compounds
    print("\n   Active compounds:")
    for c in sorted(active, key=lambda x: float(x['Activity (uM)'])):
        print(f"     ID {c['Identifier']}: {c['Activity (uM)']} uM, "
              f"{c['%Potentiation']}% potentiation")
    
    # Compute descriptors
    print("\n5. Computing molecular descriptors...")
    normalized_data = []
    for c in compounds:
        norm_smiles = normalize_smiles(c['Smiles'])
        desc = compute_descriptors(norm_smiles)
        
        activity = float(c.get('Activity (uM)', 0) or 0)
        potentiation = float(c.get('%Potentiation', 0) or 0)
        is_active = 1 if (activity > 0 and potentiation > 0) else 0
        
        normalized_data.append({
            'compound_id': c['Identifier'],
            'raw_smiles': c['Smiles'],
            'canonical_smiles': norm_smiles,
            'activity_uM': activity,
            'potentiation_pct': potentiation,
            'is_active': is_active,
            'activity_class': 'ACTIVE' if is_active else 'INACTIVE',
            **desc
        })
    
    # Save normalized dataset
    output_csv = OUTPUT_DIR / "normalized_dataset.csv"
    if normalized_data:
        fieldnames = list(normalized_data[0].keys())
        with open(output_csv, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(normalized_data)
        print(f"\n6. Saved normalized dataset to {output_csv}")
    
    # Save QC report
    report = {
        'total_compounds': len(compounds),
        'active_compounds': len(active),
        'inactive_compounds': len(inactive),
        'invalid_smiles': len(invalid_smiles),
        'duplicate_smiles': len(dup_smiles),
        'duplicate_ids': len(dup_ids),
        'active_ids': [c['Identifier'] for c in active],
        'inactive_ids': [c['Identifier'] for c in inactive],
    }
    
    report_file = OUTPUT_DIR / "qc_report.json"
    with open(report_file, 'w') as f:
        json.dump(report, f, indent=2)
    print(f"   QC report saved to {report_file}")
    
    # Save active/inactive split
    active_csv = OUTPUT_DIR / "active_compounds.csv"
    inactive_csv = OUTPUT_DIR / "inactive_compounds.csv"
    
    if active:
        with open(active_csv, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=['Identifier', 'Smiles', 'Activity (uM)', '%Potentiation'])
            writer.writeheader()
            writer.writerows(active)
    
    if inactive:
        with open(inactive_csv, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=['Identifier', 'Smiles', 'Activity (uM)', '%Potentiation'])
            writer.writeheader()
            writer.writerows(inactive)
    
    print(f"\n   Active compounds CSV: {active_csv}")
    print(f"   Inactive compounds CSV: {inactive_csv}")
    
    # Summary statistics
    if active:
        activities = [float(c['Activity (uM)']) for c in active]
        potentiations = [float(c['%Potentiation']) for c in active]
        print(f"\n7. Active compound statistics:")
        print(f"   Activity range: {min(activities):.3f} - {max(activities):.3f} uM")
        print(f"   Potentiation range: {min(potentiations):.0f} - {max(potentiations):.0f}%")
    
    print("\n" + "=" * 70)
    print("PHASE A COMPLETE")
    print("=" * 70)
    
    return report

if __name__ == "__main__":
    main()
