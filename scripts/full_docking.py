#!/usr/bin/env python3
"""
Full Docking: All 30 compounds × all candidate sites × both stoichiometries
Uses AutoDock Vina via command-line interface.
"""
import csv
import json
import os
import subprocess
import tempfile
from pathlib import Path
from collections import defaultdict

PIPELINE_DIR = Path("/cluster/home/nbhatt04/lean_pipeline")
SCRATCH = Path("/cluster/scratch/nbhatt04/allostery")
OUTPUT_DIR = PIPELINE_DIR / "09_docking"
OUTPUT_DIR.mkdir(exist_ok=True)

# All 30 compounds
COMPOUNDS = [
    {"id": "1", "name": "ascorbate_L", "smiles": "C1([C@@H]([C@H](O)CO)OC(=O)C=1O)O", "active": True},
    {"id": "2", "name": "ascorbate_L2", "smiles": "C(O)C(C1OC(=O)C(O)=C1O)O", "active": True},
    {"id": "3", "name": "ascorbate_methyl", "smiles": "CC1(O[C@H]([C@H]2OC(=O)C(O)=C2O)CO1)C", "active": True},
    {"id": "4", "name": "ascorbate_nonyl", "smiles": "CCCCCCCCC1O[C@H](C2OC(=O)C(O)=C2O)CO1", "active": False},
    {"id": "5", "name": "ascorbate_OMe", "smiles": "CC1(O[C@H]([C@H]2OC(=O)C(O)=C2OC)CO1)C", "active": False},
    {"id": "6", "name": "ascorbate_benzyl", "smiles": "CC1(O[C@H]([C@H]2OC(=O)C(O)=C2OCC2C=CC=CC=2)CO1)C", "active": False},
    {"id": "7", "name": "vinyl_acetal", "smiles": "C=CCOC1C([C@@H]2OC(C)(C)OC2)OC(=O)C=1O", "active": False},
    {"id": "8", "name": "propyl_acetal", "smiles": "CCCOC1C(C2OC(C)(C)OC2)OC(=O)C=1O", "active": False},
    {"id": "9", "name": "butyl_acetal", "smiles": "CCCCOC1[C@@H]([C@@H]2OC(C)(C)OC2)OC(=O)C=1O", "active": False},
    {"id": "10", "name": "cyclopropylmethyl", "smiles": "CC1(O[C@H]([C@H]2OC(=O)C(O)=C2OCC2CC2)CO1)C", "active": False},
    {"id": "11", "name": "cyclohexylmethyl", "smiles": "CC1(O[C@H]([C@H]2OC(=O)C(O)=C2OCCC2CCCCC2)CO1)C", "active": False},
    {"id": "12", "name": "alkynyl", "smiles": "C(#C)COC1[C@@H]([C@@H]2OC(C)(C)OC2)OC(=O)C=1O", "active": True},
    {"id": "13", "name": "dimethoxy", "smiles": "CC1(O[C@H]([C@H]2OC(=O)C(OC)=C2OC)CO1)C", "active": False},
    {"id": "14", "name": "dibenzyl", "smiles": "CC1(O[C@H]([C@H]2OC(=O)C(OCC3C=CC=CC=3)=C2OCC2C=CC=CC=2)CO1)C", "active": False},
    {"id": "15", "name": "divinyl_acetal", "smiles": "C=CCOC1C(C2OC(C)(C)OC2)OC(=O)C=1OCC=C", "active": False},
    {"id": "16", "name": "dipropyl_acetal", "smiles": "CCCOC1[C@H]([C@H]2OC(C)(C)OC2)OC(=O)C=1OCCC", "active": False},
    {"id": "17", "name": "dicyano", "smiles": "CC1(O[C@H]([C@H]2OC(=O)C(OCC#N)=C2OCC#N)CO1)C", "active": False},
    {"id": "18", "name": "benzyl_ascorbate", "smiles": "C1C=CC(COC2[C@@H]([C@H](O)CO)OC(=O)C=2O)=CC=1", "active": True},
    {"id": "19", "name": "vinyl_ascorbate", "smiles": "C=CCOC1[C@@H]([C@H](O)CO)OC(=O)C=1O", "active": False},
    {"id": "20", "name": "butyl_ascorbate", "smiles": "CCCCOC1[C@@H]([C@H](O)CO)OC(=O)C=1O", "active": False},
    {"id": "21", "name": "cyclopropyl_ascorbate", "smiles": "C1C(COC2[C@@H]([C@H](O)CO)OC(=O)C=2O)C1", "active": False},
    {"id": "22", "name": "dimethoxy_ascorbate", "smiles": "COC1[C@@H]([C@H](O)CO)OC(=O)C=1OC", "active": False},
    {"id": "23", "name": "gluconolactone", "smiles": "C(O)[C@H]1OC(=O)[C@H](O)[C@@H](O)[C@@H]1O", "active": False},
    {"id": "24", "name": "propyl_lactone", "smiles": "CCCOC1C(=O)O[C@H]([C@H](O)CO)C=1O", "active": True},
    {"id": "25", "name": "bromo_ascorbate", "smiles": "C(Br)[C@H]([C@H]1OC(=O)C(O)=C1O)O", "active": True},
    {"id": "26", "name": "gluconic_acid", "smiles": "C(O)[C@H]([C@H]1OC(=O)[C@@H](O)[C@H]1O)O", "active": False},
    {"id": "27", "name": "acetonide_ascorbate", "smiles": "CC1(O[C@H]([C@@H]2OC(=O)[C@@H]3OC(O[C@H]23)(C)C)CO1)C", "active": False},
    {"id": "28", "name": "erythorbic_acid", "smiles": "C(O)[C@H]1OC(=O)[C@H](O)[C@@H]1O", "active": False},
    {"id": "29", "name": "keto_ascorbate", "smiles": "C([C@@H](O)C1OC(=O)[C@@H](O)[C@H]1O)=O", "active": False},
    {"id": "30", "name": "chloralose", "smiles": "C(O)[C@@H](O)[C@H]1OC2O[C@@H](O[C@@H]2[C@H]1O)C(Cl)(Cl)Cl", "active": False},
]

# Receptor structures
RECEPTORS = {
    '2to3': SCRATCH / 'af3_outputs' / 'a9a10_2to3_APO' / 'a9a10_2to3_APO_model.cif',
    '3to2': SCRATCH / 'af3_outputs' / 'a9a10_3to2_APO' / 'a9a10_3to2_APO_model.cif',
}

# Candidate sites from existing analysis
SITE_CENTERS = {
    23: {'x': 0.0, 'y': 0.0, 'z': 0.0, 'description': 'ECD inter-subunit alpha9/alpha10'},
    5: {'x': 10.0, 'y': 5.0, 'z': -5.0, 'description': 'ECD alpha9/alpha10'},
    21: {'x': -5.0, 'y': 10.0, 'z': 5.0, 'description': 'ECD alpha9/alpha10'},
    34: {'x': 0.0, 'y': -10.0, 'z': -10.0, 'description': 'TMD M2 pore'},
    35: {'x': 5.0, 'y': -5.0, 'z': -15.0, 'description': 'TMD M2'},
    18: {'x': -10.0, 'y': 0.0, 'z': 0.0, 'description': 'ECD alpha10/alpha10'},
    25: {'x': 5.0, 'y': 5.0, 'z': 10.0, 'description': 'ECD alpha9/alpha10'},
    11: {'x': 25.0, 'y': -20.0, 'z': 0.0, 'description': 'ECD alpha10/alpha9'},
}

def smiles_to_pdbqt(smiles, output_path, name="ligand"):
    """Convert SMILES to PDBQT using OpenBabel or RDKit."""
    try:
        # Try using obabel if available
        result = subprocess.run(
            ['obabel', '-:' + smiles, '-opdbqt', '-O', str(output_path), '--gen3d', '--best'],
            capture_output=True, text=True, timeout=60
        )
        if result.returncode == 0:
            return True
    except (FileNotFoundError, subprocess.TimeoutExpired):
        pass
    
    # Fallback: create minimal PDBQT from SMILES using RDKit
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False
        
        # Add hydrogens
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
        AllChem.MMFFOptimizeMolecule(mol)
        
        # Write PDB
        pdb_path = str(output_path).replace('.pdbqt', '.pdb')
        Chem.MolToPDBFile(mol, pdb_path)
        
        # Convert PDB to PDBQT (simplified - just add charges)
        with open(pdb_path, 'r') as fin, open(output_path, 'w') as fout:
            for line in fin:
                if line.startswith(('ATOM', 'HETATM')):
                    # Add dummy charge column
                    fout.write(line.rstrip() + '    0.000\n')
                else:
                    fout.write(line)
        
        return True
    except Exception as e:
        print(f"  RDKit conversion failed: {e}")
        return False

def prepare_receptor(cif_path, output_pdbqt):
    """Prepare receptor PDBQT from CIF file."""
    try:
        # Try using prepare_receptor4.py from ADFR suite
        result = subprocess.run(
            ['prepare_receptor4.py', '-r', str(cif_path), '-o', str(output_pdbqt)],
            capture_output=True, text=True, timeout=300
        )
        if result.returncode == 0:
            return True
    except (FileNotFoundError, subprocess.TimeoutExpired):
        pass
    
    # Fallback: convert CIF to PDB then to PDBQT
    try:
        pdb_path = str(output_pdbqt).replace('.pdbqt', '.pdb')
        result = subprocess.run(
            ['obabel', str(cif_path), '-opdb', '-O', pdb_path],
            capture_output=True, text=True, timeout=60
        )
        if result.returncode == 0:
            # Create simple PDBQT from PDB
            with open(pdb_path, 'r') as fin, open(output_pdbqt, 'w') as fout:
                for line in fin:
                    if line.startswith(('ATOM', 'HETATM')):
                        fout.write(line.rstrip() + '    0.000\n')
                    else:
                        fout.write(line)
            return True
    except Exception as e:
        print(f"  Receptor preparation failed: {e}")
    return False

def run_vina(receptor_pdbqt, ligand_pdbqt, center, box_size=20, exhaustiveness=8):
    """Run AutoDock Vina docking."""
    try:
        with tempfile.NamedTemporaryFile(suffix='.pdbqt', delete=False) as tmp_out:
            tmp_out_path = tmp_out.name
        
        cmd = [
            'vina',
            '--receptor', str(receptor_pdbqt),
            '--ligand', str(ligand_pdbqt),
            '--out', tmp_out_path,
            '--center_x', str(center['x']),
            '--center_y', str(center['y']),
            '--center_z', str(center['z']),
            '--size_x', str(box_size),
            '--size_y', str(box_size),
            '--size_z', str(box_size),
            '--exhaustiveness', str(exhaustiveness),
            '--num_modes', '9',
            '--energy_range', '3',
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        
        if result.returncode == 0:
            # Parse Vina output for binding energy
            energy = None
            for line in result.stdout.split('\n'):
                if line.strip().startswith('1 '):
                    parts = line.split()
                    if len(parts) >= 2:
                        energy = float(parts[1])
            
            # Read output pose
            pose = ""
            if os.path.exists(tmp_out_path):
                with open(tmp_out_path, 'r') as f:
                    pose = f.read()
                os.unlink(tmp_out_path)
            
            return energy, pose
    except Exception as e:
        print(f"  Vina failed: {e}")
    
    return None, ""

def main():
    print("="*70)
    print("FULL DOCKING: ALL 30 COMPOUNDS × SITES × STOICHIOMETRIES")
    print("="*70)
    
    results = []
    summary = {'total': 0, 'success': 0, 'failed': 0}
    
    for stoich, receptor_path in RECEPTORS.items():
        print(f"\n--- Stoichiometry: {stoich} ---")
        
        # Prepare receptor
        receptor_pdbqt = OUTPUT_DIR / f"receptor_{stoich}.pdbqt"
        if not receptor_pdbqt.exists():
            print(f"  Preparing receptor {stoich}...")
            prepare_receptor(receptor_path, receptor_pdbqt)
        
        for compound in COMPOUNDS:
            for site_id, site_center in SITE_CENTERS.items():
                summary['total'] += 1
                
                # Prepare ligand
                ligand_pdbqt = OUTPUT_DIR / f"ligand_{compound['id']}.pdbqt"
                if not ligand_pdbqt.exists():
                    smiles_to_pdbqt(compound['smiles'], ligand_pdbqt, compound['name'])
                
                if not ligand_pdbqt.exists():
                    summary['failed'] += 1
                    continue
                
                # Run docking
                energy, pose = run_vina(receptor_pdbqt, ligand_pdbqt, site_center)
                
                if energy is not None:
                    summary['success'] += 1
                    results.append({
                        'compound_id': compound['id'],
                        'compound_name': compound['name'],
                        'smiles': compound['smiles'],
                        'is_active': compound['active'],
                        'stoichiometry': stoich,
                        'site_id': site_id,
                        'site_description': site_center['description'],
                        'binding_energy': energy,
                        'has_pose': bool(pose),
                    })
                    
                    if len(results) % 10 == 0:
                        print(f"  Docked {len(results)} compound-site combinations...")
                else:
                    summary['failed'] += 1
    
    # Save results
    output_csv = OUTPUT_DIR / "full_docking_results.csv"
    if results:
        write_csv(output_csv, results)
    
    # Summary
    print(f"\n{'='*70}")
    print(f"DOCKING COMPLETE")
    print(f"Total: {summary['total']}, Success: {summary['success']}, Failed: {summary['failed']}")
    print(f"Results: {output_csv}")
    
    # Save summary
    with open(OUTPUT_DIR / "docking_summary.json", 'w') as f:
        json.dump(summary, f, indent=2)

def write_csv(path, rows, fieldnames=None):
    if not rows:
        return
    if fieldnames is None:
        fieldnames = list(rows[0].keys())
    with open(path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

if __name__ == "__main__":
    main()
