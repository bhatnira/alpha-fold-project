#!/usr/bin/env python3
"""
Full Interaction Fingerprint Analysis for all compounds and sites.
Analyzes protein-ligand contacts from docking and co-folding results.
"""
import csv
import json
import os
from pathlib import Path
from collections import defaultdict

PIPELINE_DIR = Path("/cluster/home/nbhatt04/lean_pipeline")
SCRATCH = Path("/cluster/scratch/nbhatt04/allostery")
OUTPUT_DIR = PIPELINE_DIR / "06_interaction_fingerprints"
OUTPUT_DIR.mkdir(exist_ok=True)

def parse_pdb_contacts(pdb_path, ligand_residue="LIG", cutoff=4.0):
    """Parse PDB file for protein-ligand contacts within cutoff."""
    contacts = []
    protein_atoms = []
    ligand_atoms = []
    
    try:
        with open(pdb_path, 'r') as f:
            for line in f:
                if line.startswith(('ATOM', 'HETATM')):
                    resname = line[17:20].strip()
                    resseq = int(line[22:26].strip())
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    atom_name = line[12:16].strip()
                    chain = line[21:22].strip()
                    
                    atom = {
                        'chain': chain,
                        'resname': resname,
                        'resseq': resseq,
                        'atom_name': atom_name,
                        'x': x, 'y': y, 'z': z,
                    }
                    
                    if resname == ligand_residue:
                        ligand_atoms.append(atom)
                    else:
                        protein_atoms.append(atom)
        
        # Calculate distances
        for pa in protein_atoms:
            for la in ligand_atoms:
                dx = pa['x'] - la['x']
                dy = pa['y'] - la['y']
                dz = pa['z'] - la['z']
                dist = (dx*dx + dy*dy + dz*dz) ** 0.5
                
                if dist <= cutoff:
                    contacts.append({
                        'chain': pa['chain'],
                        'resname': pa['resname'],
                        'resseq': pa['resseq'],
                        'atom': pa['atom_name'],
                        'distance': round(dist, 2),
                        'contact_type': classify_contact(pa, la, dist),
                    })
    except Exception as e:
        print(f"  Error parsing {pdb_path}: {e}")
    
    return contacts

def classify_contact(protein_atom, ligand_atom, distance):
    """Classify contact type based on atom types and distance."""
    pa_name = protein_atom['atom_name']
    la_name = ligand_atom['atom_name']
    
    # H-bond: donor-acceptor within 3.5A
    if distance < 3.5 and ('N' in pa_name or 'O' in pa_name) and ('N' in la_name or 'O' in la_name):
        return 'H-bond'
    
    # Salt bridge: charged residues within 4.0A
    if distance < 4.0 and protein_atom['resname'] in ('ARG', 'LYS', 'ASP', 'GLU'):
        return 'salt_bridge'
    
    # Aromatic: aromatic residues within 4.5A
    if distance < 4.5 and protein_atom['resname'] in ('PHE', 'TRP', 'TYR', 'HIS'):
        return 'aromatic'
    
    # Hydrophobic
    if distance < 4.5 and protein_atom['resname'] in ('ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PRO'):
        return 'hydrophobic'
    
    return 'van_der_waals'

def build_ifp_matrix(contacts_list):
    """Build interaction fingerprint matrix from contact lists."""
    # Collect all unique residues
    all_residues = set()
    for contacts in contacts_list:
        for c in contacts:
            key = f"{c['chain']}:{c['resname']}:{c['resseq']}"
            all_residues.add(key)
    
    all_residues = sorted(all_residues)
    
    # Build binary IFP matrix
    matrix = []
    for contacts in contacts_list:
        fp = {}
        for res in all_residues:
            fp[res] = 0
        for c in contacts:
            key = f"{c['chain']}:{c['resname']}:{c['resseq']}"
            fp[key] = 1
        matrix.append(fp)
    
    return all_residues, matrix

def main():
    print("="*70)
    print("FULL INTERACTION FINGERPRINT ANALYSIS")
    print("="*70)
    
    # Load docking results
    docking_file = PIPELINE_DIR / "09_docking" / "full_docking_results.csv"
    if not docking_file.exists():
        print("No docking results found. Run full_docking.py first.")
        return
    
    docking = []
    with open(docking_file, 'r') as f:
        reader = csv.DictReader(f)
        docking = list(reader)
    
    print(f"Loaded {len(docking)} docking results")
    
    # Analyze each docked pose
    all_ifps = []
    for d in docking:
        pdb_path = PIPELINE_DIR / "09_docking" / f"docked_{d['compound_id']}_{d['stoichiometry']}_site{d['site_id']}.pdb"
        
        contacts = []
        if pdb_path.exists():
            contacts = parse_pdb_contacts(pdb_path)
        
        ifps = {
            'compound_id': d['compound_id'],
            'stoichiometry': d['stoichiometry'],
            'site_id': d['site_id'],
            'binding_energy': d['binding_energy'],
            'n_contacts': len(contacts),
            'h_bonds': len([c for c in contacts if c['contact_type'] == 'H-bond']),
            'salt_bridges': len([c for c in contacts if c['contact_type'] == 'salt_bridge']),
            'aromatic': len([c for c in contacts if c['contact_type'] == 'aromatic']),
            'hydrophobic': len([c for c in contacts if c['contact_type'] == 'hydrophobic']),
            'vdw': len([c for c in contacts if c['contact_type'] == 'van_der_waals']),
            'contact_residues': ';'.join(sorted(set(
                f"{c['chain']}:{c['resname']}:{c['resseq']}" for c in contacts
            ))),
        }
        all_ifps.append(ifps)
    
    # Save IFP summary
    if all_ifps:
        fieldnames = list(all_ifps[0].keys())
        with open(OUTPUT_DIR / "full_ifp_summary.csv", 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(all_ifps)
    
    # Analyze contact frequency
    contact_freq = defaultdict(int)
    for ifp in all_ifps:
        if ifp['contact_residues']:
            for res in ifp['contact_residues'].split(';'):
                contact_freq[res] += 1
    
    top_contacts = sorted(contact_freq.items(), key=lambda x: x[1], reverse=True)[:30]
    
    print(f"\nTop 30 contact residues:")
    for res, freq in top_contacts:
        print(f"  {res}: {freq} contacts")
    
    # Save contact frequency
    with open(OUTPUT_DIR / "contact_frequency.csv", 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['residue', 'frequency'])
        for res, freq in top_contacts:
            writer.writerow([res, freq])
    
    # Active vs inactive comparison
    active_ifps = [i for i in all_ifps if any(d['is_active'] == 'True' for d in docking if d['compound_id'] == i['compound_id'])]
    inactive_ifps = [i for i in all_ifps if not any(d['is_active'] == 'True' for d in docking if d['compound_id'] == i['compound_id'])]
    
    print(f"\nActive IFPs: {len(active_ifps)}")
    print(f"Inactive IFPs: {len(inactive_ifps)}")
    
    if active_ifps:
        avg_hbonds_active = sum(i['h_bonds'] for i in active_ifps) / len(active_ifps)
        print(f"  Avg H-bonds (active): {avg_hbonds_active:.1f}")
    if inactive_ifps:
        avg_hbonds_inactive = sum(i['h_bonds'] for i in inactive_ifps) / len(inactive_ifps)
        print(f"  Avg H-bonds (inactive): {avg_hbonds_inactive:.1f}")
    
    print(f"\nIFP analysis complete: {OUTPUT_DIR / 'full_ifp_summary.csv'}")

if __name__ == "__main__":
    main()
