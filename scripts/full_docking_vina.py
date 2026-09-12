#!/usr/bin/env python3
"""
Full Docking: All 30 compounds × candidate sites × both stoichiometries
Uses AutoDock Vina Python API.
"""
import csv
import json
import math
import os
import sys
from pathlib import Path
from collections import defaultdict

sys.path.insert(0, os.path.expanduser('~/.local/lib/python3.9/site-packages'))

try:
    from vina import Vina
    HAS_VINA = True
except ImportError:
    HAS_VINA = False
    print("ERROR: vina not importable")

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdMolDescriptors
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False

PIPELINE_DIR = Path("/cluster/home/nbhatt04/lean_pipeline")
SCRATCH = Path("/cluster/scratch/nbhatt04/allostery")
OUTPUT_DIR = PIPELINE_DIR / "09_docking"
OUTPUT_DIR.mkdir(exist_ok=True)

# All 30 compounds
COMPOUNDS = [
    {"id": "1", "smiles": "C1([C@@H]([C@H](O)CO)OC(=O)C=1O)O", "active": True, "activity_uM": 1797, "potentiation": 286},
    {"id": "2", "smiles": "C(O)C(C1OC(=O)C(O)=C1O)O", "active": True, "activity_uM": 1316, "potentiation": 300},
    {"id": "3", "smiles": "CC1(O[C@H]([C@H]2OC(=O)C(O)=C2O)CO1)C", "active": True, "activity_uM": 6077, "potentiation": 293},
    {"id": "4", "smiles": "CCCCCCCCC1O[C@H](C2OC(=O)C(O)=C2O)CO1", "active": False, "activity_uM": 0, "potentiation": 0},
    {"id": "5", "smiles": "CC1(O[C@H]([C@H]2OC(=O)C(O)=C2OC)CO1)C", "active": False, "activity_uM": 0, "potentiation": 0},
    {"id": "6", "smiles": "CC1(O[C@H]([C@H]2OC(=O)C(O)=C2OCC2C=CC=CC=2)CO1)C", "active": False, "activity_uM": 0, "potentiation": 0},
    {"id": "7", "smiles": "C=CCOC1C([C@@H]2OC(C)(C)OC2)OC(=O)C=1O", "active": False, "activity_uM": 0, "potentiation": 0},
    {"id": "8", "smiles": "CCCOC1C(C2OC(C)(C)OC2)OC(=O)C=1O", "active": False, "activity_uM": 0, "potentiation": 0},
    {"id": "9", "smiles": "CCCCOC1[C@@H]([C@@H]2OC(C)(C)OC2)OC(=O)C=1O", "active": False, "activity_uM": 0, "potentiation": 0},
    {"id": "10", "smiles": "CC1(O[C@H]([C@H]2OC(=O)C(O)=C2OCC2CC2)CO1)C", "active": False, "activity_uM": 0, "potentiation": 0},
    {"id": "11", "smiles": "CC1(O[C@H]([C@H]2OC(=O)C(O)=C2OCCC2CCCCC2)CO1)C", "active": False, "activity_uM": 0, "potentiation": 0},
    {"id": "12", "smiles": "C(#C)COC1[C@@H]([C@@H]2OC(C)(C)OC2)OC(=O)C=1O", "active": True, "activity_uM": 0.1981, "potentiation": 150},
    {"id": "13", "smiles": "CC1(O[C@H]([C@H]2OC(=O)C(OC)=C2OC)CO1)C", "active": False, "activity_uM": 0, "potentiation": 0},
    {"id": "14", "smiles": "CC1(O[C@H]([C@H]2OC(=O)C(OCC3C=CC=CC=3)=C2OCC2C=CC=CC=2)CO1)C", "active": False, "activity_uM": 0, "potentiation": 0},
    {"id": "15", "smiles": "C=CCOC1C(C2OC(C)(C)OC2)OC(=O)C=1OCC=C", "active": False, "activity_uM": 0, "potentiation": 0},
    {"id": "16", "smiles": "CCCOC1[C@H]([C@H]2OC(C)(C)OC2)OC(=O)C=1OCCC", "active": False, "activity_uM": 0, "potentiation": 0},
    {"id": "17", "smiles": "CC1(O[C@H]([C@H]2OC(=O)C(OCC#N)=C2OCC#N)CO1)C", "active": False, "activity_uM": 0, "potentiation": 0},
    {"id": "18", "smiles": "C1C=CC(COC2[C@@H]([C@H](O)CO)OC(=O)C=2O)=CC=1", "active": True, "activity_uM": 1288, "potentiation": 500},
    {"id": "19", "smiles": "C=CCOC1[C@@H]([C@H](O)CO)OC(=O)C=1O", "active": False, "activity_uM": 0, "potentiation": 0},
    {"id": "20", "smiles": "CCCCOC1[C@@H]([C@H](O)CO)OC(=O)C=1O", "active": False, "activity_uM": 0, "potentiation": 0},
    {"id": "21", "smiles": "C1C(COC2[C@@H]([C@H](O)CO)OC(=O)C=2O)C1", "active": False, "activity_uM": 0, "potentiation": 0},
    {"id": "22", "smiles": "COC1[C@@H]([C@H](O)CO)OC(=O)C=1OC", "active": False, "activity_uM": 0, "potentiation": 0},
    {"id": "23", "smiles": "C(O)[C@H]1OC(=O)[C@H](O)[C@@H](O)[C@@H]1O", "active": False, "activity_uM": 0, "potentiation": 0},
    {"id": "24", "smiles": "CCCOC1C(=O)O[C@H]([C@H](O)CO)C=1O", "active": True, "activity_uM": 1202, "potentiation": 190},
    {"id": "25", "smiles": "C(Br)[C@H]([C@H]1OC(=O)C(O)=C1O)O", "active": True, "activity_uM": 2.63, "potentiation": 180},
    {"id": "26", "smiles": "C(O)[C@H]([C@H]1OC(=O)[C@@H](O)[C@H]1O)O", "active": False, "activity_uM": 0, "potentiation": 0},
    {"id": "27", "smiles": "CC1(O[C@H]([C@@H]2OC(=O)[C@@H]3OC(O[C@H]23)(C)C)CO1)C", "active": False, "activity_uM": 0, "potentiation": 0},
    {"id": "28", "smiles": "C(O)[C@H]1OC(=O)[C@H](O)[C@@H]1O", "active": False, "activity_uM": 0, "potentiation": 0},
    {"id": "29", "smiles": "C([C@@H](O)C1OC(=O)[C@@H](O)[C@H]1O)=O", "active": False, "activity_uM": 0, "potentiation": 0},
    {"id": "30", "smiles": "C(O)[C@@H](O)[C@H]1OC2O[C@@H](O[C@@H]2[C@H]1O)C(Cl)(Cl)Cl", "active": False, "activity_uM": 0, "potentiation": 0},
]

# Receptor PDBQT files (pre-converted from CIF)
RECEPTORS = {
    '2to3': OUTPUT_DIR / "receptor_2to3.pdbqt",
    '3to2': OUTPUT_DIR / "receptor_3to2.pdbqt",
}

# Candidate sites with centers and box sizes
SITES = {
    23: {'center': [0.0, 0.0, 20.0], 'box_size': 20, 'desc': 'ECD alpha9/alpha10'},
    5:  {'center': [10.0, 5.0, 15.0], 'box_size': 20, 'desc': 'ECD alpha9/alpha10'},
    21: {'center': [-5.0, 10.0, 25.0], 'box_size': 20, 'desc': 'ECD alpha9/alpha10'},
    34: {'center': [0.0, -10.0, -5.0], 'box_size': 22, 'desc': 'TMD M2 pore'},
    35: {'center': [5.0, -5.0, -10.0], 'box_size': 22, 'desc': 'TMD M2'},
    18: {'center': [-10.0, 0.0, 20.0], 'box_size': 20, 'desc': 'ECD alpha10/alpha10'},
    25: {'center': [5.0, 5.0, 30.0], 'box_size': 20, 'desc': 'ECD alpha9/alpha10'},
    11: {'center': [25.0, -20.0, 20.0], 'box_size': 25, 'desc': 'ECD alpha10/alpha9'},
}

def smiles_to_mol(smiles):
    """Convert SMILES to RDKit mol with 3D coords."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    result = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    if result != 0:
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3(), useRandomCoords=True)
    AllChem.MMFFOptimizeMolecule(mol)
    return mol

AD_LIGAND = {'C':'C','N':'NA','O':'OA','S':'SA','H':'HD','P':'P','F':'F','Cl':'Cl','Br':'Br','I':'I','FE':'FE','ZN':'ZN','CA':'CA','MG':'MG'}

def element_from_pdb_name(aname):
    s = aname.strip()
    if not s:
        return 'C'
    while s and s[0].isdigit():
        s = s[1:]
    if not s:
        return 'C'
    two = s[:2]
    if two in AD_LIGAND:
        return two
    one = s[0]
    if one in AD_LIGAND:
        return one
    return 'C'

def mol_to_pdbqt(mol, output_path):
    """Write mol to PDBQT format with strict columns."""
    sdf_path = str(output_path).replace('.pdbqt', '.sdf')
    pdb_path = str(output_path).replace('.pdbqt', '.pdb')

    writer = Chem.SDWriter(sdf_path)
    writer.write(mol)
    writer.close()

    supplier = Chem.SDMolSupplier(sdf_path)
    for m in supplier:
        if m is not None:
            Chem.MolToPDBFile(m, pdb_path)
            break

    atom_lines = []
    with open(pdb_path, 'r') as fin:
        for line in fin:
            if not line.startswith(('ATOM', 'HETATM')):
                continue
            record   = line[0:6]
            serial_s = line[6:11].strip()
            name     = line[12:16]
            altloc   = line[16:17]
            resname  = line[17:20]
            chain    = line[21:22]
            resseq_s = line[22:26].strip()
            icode    = line[26:27] if len(line) > 26 else ' '
            x_str = line[30:38].strip() if len(line) > 37 else '0.0'
            y_str = line[38:46].strip() if len(line) > 45 else '0.0'
            z_str = line[46:54].strip() if len(line) > 53 else '0.0'
            occ_s = line[54:60].strip() if len(line) > 59 else '1.00'
            bfac_s = line[60:66].strip() if len(line) > 65 else '0.00'

            serial = int(serial_s) if serial_s.isdigit() else 0
            resseq = int(resseq_s) if resseq_s.isdigit() else 0
            x = float(x_str)
            y = float(y_str)
            z = float(z_str)
            occ = float(occ_s) if occ_s else 1.0
            bfac = float(bfac_s) if bfac_s else 0.0

            elem = element_from_pdb_name(name)
            ad_type = AD_LIGAND.get(elem, elem)

            out = f"{record:<6.6}{serial:>5d} {name:<4.4}{altloc}{resname:<3.3} {chain}{resseq:>4d}{icode}   {x:8.3f}{y:8.3f}{z:8.3f}{occ:6.2f}{bfac:6.2f}          {ad_type:>2s}\n"
            atom_lines.append(out)

    with open(output_path, 'w') as fout:
        fout.write("ROOT\n")
        for al in atom_lines:
            fout.write(al)
        fout.write("ENDROOT\n")
        fout.write("TORSDOF 0\n")
    
    return True

def prepare_receptor_pdbqt(cif_path, output_pdbqt):
    """Convert receptor CIF to PDBQT."""
    pdb_path = str(output_pdbqt).replace('.pdbqt', '.pdb')
    
    # Try obabel
    try:
        import subprocess
        result = subprocess.run(['obabel', str(cif_path), '-opdb', '-O', pdb_path],
                              capture_output=True, text=True, timeout=120)
        if result.returncode == 0 and os.path.exists(pdb_path):
            with open(pdb_path, 'r') as fin, open(output_pdbqt, 'w') as fout:
                for line in fin:
                    if line.startswith(('ATOM', 'HETATM')):
                        fout.write(line.rstrip() + '    0.000\n')
                    else:
                        fout.write(line)
            return True
    except Exception as e:
        print(f"  obabel failed: {e}")
    
    # Try with rdkit
    try:
        from rdkit import Chem
        mol = Chem.MolFromPDBFile(str(cif_path))
        if mol is not None:
            Chem.MolToPDBFile(mol, pdb_path)
            with open(pdb_path, 'r') as fin, open(output_pdbqt, 'w') as fout:
                for line in fin:
                    if line.startswith(('ATOM', 'HETATM')):
                        fout.write(line.rstrip() + '    0.000\n')
                    else:
                        fout.write(line)
            return True
    except Exception as e:
        print(f"  rdkit failed: {e}")
    
    return False

def dock_single(v, ligand_pdbqt_path, center, box_size):
    """Run a single docking."""
    try:
        v.compute_vina_maps(center=center, box_size=[box_size]*3)
        energy = v.dock(exhaustiveness=8, n_poses=9)
        
        # Get best pose energy
        energies = v.energies()
        best_energy = float(energies[0][0]) if len(energies) > 0 else None
        
        # Get best pose
        out_path = str(ligand_pdbqt_path).replace('.pdbqt', '_docked.pdbqt')
        v.write_poses(out_path, n_poses=1, overwrite=True)
        
        return best_energy, out_path
    except Exception as e:
        print(f"  Docking failed: {e}")
        return None, None

def write_csv(path, rows):
    if not rows:
        return
    fieldnames = list(rows[0].keys())
    with open(path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

def main():
    print("="*70)
    print("FULL DOCKING: ALL 30 COMPOUNDS × 8 SITES × 2 STOICHIOMETRIES")
    print("="*70)
    
    if not HAS_VINA:
        print("FATAL: vina not importable")
        return
    
    if not HAS_RDKIT:
        print("FATAL: rdkit not importable")
        return
    
    results = []
    summary = {'total': 0, 'success': 0, 'failed': 0}
    
    for stoich, receptor_path in RECEPTORS.items():
        print(f"\n=== Stoichiometry: {stoich} ===")
        
        # Check receptor
        if not receptor_path.exists():
            print(f"  Receptor not found: {receptor_path}")
            print(f"  Attempting conversion...")
            
            # Find CIF source
            cif_candidates = [
                SCRATCH / 'af3_outputs' / f'a9a10_{stoich}_APO' / f'a9a10_{stoich}_APO_model.cif',
                SCRATCH / 'af3_outputs' / f'a9a10_{stoich}_LASC_ANION' / f'a9a10_{stoich}_LASC_ANION_model.cif',
            ]
            
            converted = False
            for cif in cif_candidates:
                if cif.exists():
                    if prepare_receptor_pdbqt(cif, receptor_path):
                        print(f"  Converted {cif.name} to PDBQT")
                        converted = True
                        break
            
            if not converted:
                print(f"  SKIPPING {stoich} - no receptor available")
                continue
        
        # Initialize Vina
        v = Vina(sf_name='vina')
        try:
            v.set_receptor(str(receptor_path))
        except Exception as e:
            print(f"  Failed to set receptor: {e}")
            continue
        
        for compound in COMPOUNDS:
            # Prepare ligand
            ligand_pdbqt = OUTPUT_DIR / f"ligand_{compound['id']}.pdbqt"
            
            if not ligand_pdbqt.exists():
                mol = smiles_to_mol(compound['smiles'])
                if mol is None:
                    print(f"  Compound {compound['id']}: SMILES failed")
                    summary['failed'] += 1
                    continue
                mol_to_pdbqt(mol, ligand_pdbqt)
            
            if not ligand_pdbqt.exists():
                print(f"  Compound {compound['id']}: PDBQT not created")
                summary['failed'] += 1
                continue
            
            try:
                v.set_ligand_from_file(str(ligand_pdbqt))
            except Exception as e:
                print(f"  Compound {compound['id']}: set_ligand failed: {e}")
                summary['failed'] += 1
                continue
            
            for site_id, site_info in SITES.items():
                summary['total'] += 1
                
                energy, pose_path = dock_single(
                    v, ligand_pdbqt, 
                    site_info['center'], site_info['box_size']
                )
                
                if energy is not None:
                    summary['success'] += 1
                    results.append({
                        'compound_id': compound['id'],
                        'smiles': compound['smiles'],
                        'is_active': compound['active'],
                        'activity_uM': compound['activity_uM'],
                        'potentiation_pct': compound['potentiation'],
                        'stoichiometry': stoich,
                        'site_id': site_id,
                        'site_description': site_info['desc'],
                        'binding_energy': round(energy, 4),
                        'has_pose': os.path.exists(pose_path) if pose_path else False,
                    })
                else:
                    summary['failed'] += 1
                
                if summary['total'] % 20 == 0:
                    print(f"  Progress: {summary['total']} total, "
                          f"{summary['success']} success, {summary['failed']} failed")
    
    # Save results
    if results:
        write_csv(OUTPUT_DIR / "full_docking_results.csv", results)
    
    # Summary
    print(f"\n{'='*70}")
    print(f"DOCKING COMPLETE")
    print(f"Total: {summary['total']}, Success: {summary['success']}, Failed: {summary['failed']}")
    if results:
        print(f"Results saved: {OUTPUT_DIR / 'full_docking_results.csv'}")
    
    # Save summary
    with open(OUTPUT_DIR / "docking_summary.json", 'w') as f:
        json.dump(summary, f, indent=2)
    
    # Quick analysis
    if results:
        active_energies = [r['binding_energy'] for r in results if r['is_active']]
        inactive_energies = [r['binding_energy'] for r in results if not r['is_active']]
        
        if active_energies:
            print(f"\nActive compounds: mean={sum(active_energies)/len(active_energies):.3f}")
        if inactive_energies:
            print(f"Inactive compounds: mean={sum(inactive_energies)/len(inactive_energies):.3f}")
        
        # Best binding per compound
        by_compound = defaultdict(list)
        for r in results:
            by_compound[r['compound_id']].append(r)
        
        print(f"\nBest binding energy per compound (Site 23):")
        for cid in sorted(by_compound.keys(), key=lambda x: int(x)):
            site23 = [r for r in by_compound[cid] if r['site_id'] == 23]
            if site23:
                best = min(site23, key=lambda x: x['binding_energy'])
                active = "ACTIVE" if best['is_active'] else "inactive"
                print(f"  Compound {cid}: {best['binding_energy']:.3f} kcal/mol ({active})")

if __name__ == "__main__":
    main()
