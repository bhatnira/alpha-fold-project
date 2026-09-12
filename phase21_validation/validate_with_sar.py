#!/usr/bin/env python3
"""Validate computational model with experimental SAR data.

Docks active compounds from modulator-dataset-a9a10.csv to Site 23
(ascorbate binding site) and correlates with experimental IC50.
"""

import csv, sys, os, subprocess, math
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
OUT = Path(__file__).resolve().parent
RESULTS = OUT / "results"
RESULTS.mkdir(exist_ok=True)

# AutoDock atom type mapping
AUTOATOM_MAP = {
    'C': 'C', 'CA': 'C', 'CB': 'C', 'CG': 'C', 'CD': 'C', 'CE': 'C', 'CZ': 'C',
    'CG1': 'C', 'CG2': 'C', 'CD1': 'C', 'CD2': 'C', 'CE1': 'C', 'CE2': 'C',
    'N': 'NA', 'ND1': 'NA', 'ND2': 'NA', 'NE': 'NA', 'NE1': 'NA', 'NE2': 'NA',
    'NH1': 'NA', 'NH2': 'NA', 'NZ': 'NA', 'NZ2': 'NA', 'NE2': 'NA',
    'O': 'OA', 'OD1': 'OA', 'OD2': 'OA', 'OE1': 'OA', 'OE2': 'OA', 'OG': 'OA',
    'OG1': 'OA', 'OH': 'OA', 'OXT': 'OA',
    'S': 'SA', 'SD': 'SA', 'SG': 'SA',
    'H': 'HD', 'HA': 'HD', 'HB': 'HD', 'HG': 'HD', 'HE': 'HD', 'HZ': 'HD',
    'HB1': 'HD', 'HB2': 'HD', 'HB3': 'HD', 'HG1': 'HD', 'HG2': 'HD',
    'F': 'F', 'CL': 'Cl', 'BR': 'Br', 'I': 'I',
}

RESIDUE_CHARGES = {
    'ALA': {'N': -0.4157, 'CA': 0.0337, 'C': 0.5973, 'O': -0.5679, 'CB': -0.0101},
    'ARG': {'N': -0.4157, 'CA': 0.0337, 'C': 0.5973, 'O': -0.5679,
            'CB': 0.0337, 'CG': 0.0337, 'CD': 0.0337, 'NE': -0.3274,
            'CZ': 0.8076, 'NH1': -0.5877, 'NH2': -0.5877},
    'ASN': {'N': -0.4157, 'CA': 0.0337, 'C': 0.5973, 'O': -0.5679,
            'CB': 0.0337, 'CG': 0.6755, 'OD1': -0.5877, 'ND2': -0.6755},
    'ASP': {'N': -0.4157, 'CA': 0.0337, 'C': 0.5973, 'O': -0.5679,
            'CB': 0.0337, 'CG': 0.6755, 'OD1': -0.6755, 'OD2': -0.6755},
    'CYS': {'N': -0.4157, 'CA': 0.0337, 'C': 0.5973, 'O': -0.5679, 'CB': 0.0337, 'SG': -0.2029},
    'GLU': {'N': -0.4157, 'CA': 0.0337, 'C': 0.5973, 'O': -0.5679,
            'CB': 0.0337, 'CG': 0.0337, 'CD': 0.6755, 'OE1': -0.6755, 'OE2': -0.6755},
    'GLN': {'N': -0.4157, 'CA': 0.0337, 'C': 0.5973, 'O': -0.5679,
            'CB': 0.0337, 'CG': 0.0337, 'CD': 0.6755, 'OE1': -0.5877, 'NE2': -0.6755},
    'GLY': {'N': -0.4157, 'CA': 0.0337, 'C': 0.5973, 'O': -0.5679},
    'HIS': {'N': -0.4157, 'CA': 0.0337, 'C': 0.5973, 'O': -0.5679,
            'CB': 0.0337, 'CG': 0.1879, 'ND1': -0.5877, 'CE1': 0.1366, 'NE2': -0.5877, 'CD2': 0.1879},
    'ILE': {'N': -0.4157, 'CA': 0.0337, 'C': 0.5973, 'O': -0.5679, 'CB': 0.0337, 'CG1': 0.0337, 'CG2': -0.0101, 'CD1': -0.0101},
    'LEU': {'N': -0.4157, 'CA': 0.0337, 'C': 0.5973, 'O': -0.5679, 'CB': 0.0337, 'CG': 0.0337, 'CD1': -0.0101, 'CD2': -0.0101},
    'LYS': {'N': -0.4157, 'CA': 0.0337, 'C': 0.5973, 'O': -0.5679,
            'CB': 0.0337, 'CG': 0.0337, 'CD': 0.0337, 'CE': 0.0337, 'NZ': -0.2917},
    'MET': {'N': -0.4157, 'CA': 0.0337, 'C': 0.5973, 'O': -0.5679, 'CB': 0.0337, 'CG': 0.0337, 'SD': -0.1022, 'CE': -0.0101},
    'PHE': {'N': -0.4157, 'CA': 0.0337, 'C': 0.5973, 'O': -0.5679, 'CB': 0.0337,
            'CG': 0.0189, 'CD1': -0.0189, 'CD2': -0.0189, 'CE1': -0.0189, 'CE2': -0.0189, 'CZ': 0.0189},
    'PRO': {'N': -0.2111, 'CA': 0.0458, 'C': 0.5718, 'O': -0.5284, 'CB': 0.0337, 'CG': 0.0337, 'CD': 0.0337},
    'SER': {'N': -0.4157, 'CA': 0.0337, 'C': 0.5973, 'O': -0.5679, 'CB': 0.0337, 'OG': -0.6755},
    'THR': {'N': -0.4157, 'CA': 0.0337, 'C': 0.5973, 'O': -0.5679, 'CB': 0.0337, 'OG1': -0.6755, 'CG2': -0.0101},
    'TRP': {'N': -0.4157, 'CA': 0.0337, 'C': 0.5973, 'O': -0.5679, 'CB': 0.0337,
            'CG': 0.0189, 'CD1': -0.0308, 'CD2': 0.0337, 'NE1': -0.3616, 'CE2': 0.0617, 'CE3': -0.0337, 'CZ2': -0.0337, 'CZ3': -0.0337, 'CH2': -0.0337},
    'TYR': {'N': -0.4157, 'CA': 0.0337, 'C': 0.5973, 'O': -0.5679, 'CB': 0.0337,
            'CG': 0.0189, 'CD1': -0.0189, 'CD2': -0.0189, 'CE1': -0.0189, 'CE2': -0.0189, 'CZ': 0.0617, 'OH': -0.6755},
    'VAL': {'N': -0.4157, 'CA': 0.0337, 'C': 0.5973, 'O': -0.5679, 'CB': 0.0337, 'CG1': -0.0101, 'CG2': -0.0101},
}


def get_autodock_type(atom_name):
    name = atom_name.strip()
    if name in AUTOATOM_MAP:
        return AUTOATOM_MAP[name]
    if name:
        first = name[0].upper()
        if first == 'C': return 'C'
        elif first == 'N': return 'NA'
        elif first == 'O': return 'OA'
        elif first == 'S': return 'SA'
        elif first == 'H': return 'HD'
        elif first == 'F': return 'F'
        elif first == 'B': return 'Br'
        elif first == 'I': return 'I'
    return 'C'


def get_charge(resname, atom_name):
    resname = resname.strip().upper()
    atom_name = atom_name.strip()
    if resname in RESIDUE_CHARGES and atom_name in RESIDUE_CHARGES[resname]:
        return RESIDUE_CHARGES[resname][atom_name]
    atype = get_autodock_type(atom_name)
    if atype == 'NA': return -0.5
    elif atype == 'OA': return -0.5
    elif atype == 'HD': return 0.5
    return 0.0


def fmt_pdbqt_line(serial, atom_name, altloc, resname, chain, resseq, icode,
                    x, y, z, occ, tf, charge, atype):
    return (
        f"ATOM  "
        f"{serial:>5d}"
        f" "
        f"{atom_name:<4s}"
        f"{altloc}"
        f"{resname:>3s}"
        f" "
        f"{chain}"
        f"{resseq:>4d}"
        f"{icode}"
        f"   "
        f"{x:>8.3f}"
        f"{y:>8.3f}"
        f"{z:>8.3f}"
        f"{occ:>6.2f}"
        f"{tf:>6.2f}"
        f"     "
        f"{charge:>5.2f}"
        f" "
        f"{atype:<2s}"
    )


def prepare_receptor_pdbqt(pdb_path, output_path):
    serial = 0
    altloc_prev = {}
    with open(pdb_path) as fin, open(output_path, 'w') as fout:
        for line in fin:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            try:
                atom_name = line[12:16]
                altloc = line[16] if len(line) > 16 else ' '
                resname = line[17:20]
                chain = line[21] if len(line) > 21 and line[21].strip() else 'A'
                resseq = int(line[22:26].strip())
                icode = line[26] if len(line) > 26 else ' '
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except (ValueError, IndexError):
                continue

            key = (chain, resseq, atom_name.strip())
            if altloc not in (' ', 'A', '1'):
                if key in altloc_prev:
                    continue
            altloc_prev[key] = altloc

            serial += 1
            atype = get_autodock_type(atom_name)
            charge = get_charge(resname, atom_name)

            fout.write(fmt_pdbqt_line(serial, atom_name, altloc, resname,
                                      chain, resseq, icode, x, y, z,
                                      1.0, 0.0, charge, atype) + "\n")
    return serial


def prepare_ligand_pdbqt(smiles, output_path):
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)
        AllChem.ComputeGasteigerCharges(mol)

        conf = mol.GetConformer()
        n_atoms = conf.GetNumAtoms()

        with open(output_path, 'w') as f:
            f.write(f'REMARK  Generated from SMILES: {smiles}\n')
            f.write('ROOT\n')
            for i in range(n_atoms):
                atom = mol.GetAtomWithIdx(i)
                pos = conf.GetAtomPosition(i)
                elem = atom.GetSymbol()
                gc = float(atom.GetDoubleProp('_GasteigerCharge'))
                if math.isnan(gc):
                    gc = 0.0

                atype_map = {'C': 'C', 'N': 'NA', 'O': 'OA', 'S': 'SA', 'H': 'HD',
                             'F': 'F', 'Cl': 'Cl', 'Br': 'Br', 'I': 'I'}
                atype = atype_map.get(elem, 'C')

                aname = f'{elem}{i+1:2d}' if len(elem) <= 1 else f'{elem}{i+1}'
                if len(aname) < 4:
                    aname = aname.ljust(4)
                elif len(aname) > 4:
                    aname = aname[:4]

                serial = i + 1
                f.write(f'ATOM  {serial:>5d} {aname:<4s} LIG A   1    '
                        f'{pos.x:>8.3f}{pos.y:>8.3f}{pos.z:>8.3f}  1.00  0.00'
                        f'     {gc:>5.2f} {atype:<2s}\n')
            f.write('ENDROOT\n')
            f.write('TORSDOF 1\n')
        return True
    except Exception as e:
        print(f"    Ligand prep error: {e}")
        return False


def run_vina_docking(receptor_pdbqt, ligand_pdbqt, center, box_size, output_path, exhaustiveness=16, n_poses=10):
    try:
        from vina import Vina
        v = Vina(sf_name='vina')
        v.set_receptor(receptor_pdbqt)
        v.set_ligand_from_file(ligand_pdbqt)
        v.compute_vina_maps(center=center, box_size=box_size)
        v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)
        energies = v.energies()
        v.write_poses(str(output_path), n_poses=n_poses, overwrite=True)
        best_energy = float(energies[0][0]) if len(energies) > 0 else None
        return True, best_energy
    except Exception as e:
        print(f"Vina error: {e}")
        return False, None


def load_active_compounds(csv_path):
    compounds = []
    with open(csv_path) as f:
        for row in csv.DictReader(f):
            ic50 = float(row['Activity (uM)'])
            if ic50 > 0:
                compounds.append({
                    'id': row['Identifier'],
                    'smiles': row['Smiles'],
                    'ic50': ic50,
                    'potentiation': float(row['%Potentiation']),
                })
    return compounds


def main():
    print("Model Validation with Experimental SAR Data")
    print("=" * 60)
    
    # Load active compounds
    csv_path = ROOT / "modulator-dataset-a9a10.csv"
    compounds = load_active_compounds(csv_path)
    print(f"Loaded {len(compounds)} active compounds")
    
    # Receptor (Site 23 - ascorbate binding site)
    receptor = ROOT / "data/deliverable/05_structures/2to3/a9a10_2to3_apo_reference.pdb"
    if not receptor.exists():
        receptor = ROOT / "data/deliverable/05_structures/3to2/a9a10_3to2_apo_reference.pdb"
    
    rec_pdbqt = RESULTS / "site23_receptor.pdbqt"
    print(f"\nPreparing receptor...")
    prepare_receptor_pdbqt(str(receptor), str(rec_pdbqt))
    
    # Site 23 center (from SAR dataset)
    center = (-11.65, 38.0, 22.41)
    box_size = (25, 25, 25)
    
    results = []
    for comp in compounds:
        lig_pdbqt = RESULTS / f"compound_{comp['id']}.pdbqt"
        print(f"\nCompound {comp['id']}: IC50={comp['ic50']:.2f} μM")
        
        if not lig_pdbqt.exists() or lig_pdbqt.stat().st_size < 10:
            print(f"  Preparing ligand...")
            if not prepare_ligand_pdbqt(comp['smiles'], str(lig_pdbqt)):
                print(f"  Failed to prepare ligand")
                continue
        
        out_pdb = RESULTS / f"compound_{comp['id']}_docked.pdb"
        print(f"  Docking to Site 23...", end=" ", flush=True)
        success, energy = run_vina_docking(
            str(rec_pdbqt), str(lig_pdbqt),
            center, box_size, str(out_pdb)
        )
        
        if success:
            results.append({
                'id': comp['id'],
                'smiles': comp['smiles'],
                'ic50': comp['ic50'],
                'potentiation': comp['potentiation'],
                'docking_energy': energy,
                'pIC50': -math.log10(comp['ic50'] * 1e-6) if comp['ic50'] > 0 else 0,
            })
            print(f"OK  E={energy:.2f} kcal/mol")
        else:
            print("FAILED")
    
    # Save results
    if results:
        csv_out = OUT / "validation_docking_results.csv"
        with open(csv_out, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=results[0].keys())
            writer.writeheader()
            writer.writerows(results)
        print(f"\nResults saved: {csv_out}")
        
        # Correlation analysis
        print("\n" + "=" * 60)
        print("CORRELATION ANALYSIS")
        print("=" * 60)
        
        # Calculate correlation
        n = len(results)
        if n >= 3:
            x = [r['docking_energy'] for r in results]
            y = [r['pIC50'] for r in results]
            
            mean_x = sum(x) / n
            mean_y = sum(y) / n
            
            cov_xy = sum((xi - mean_x) * (yi - mean_y) for xi, yi in zip(x, y)) / n
            std_x = (sum((xi - mean_x)**2 for xi in x) / n) ** 0.5
            std_y = (sum((yi - mean_y)**2 for yi in y) / n) ** 0.5
            
            if std_x > 0 and std_y > 0:
                r = cov_xy / (std_x * std_y)
                r2 = r ** 2
                print(f"Pearson r: {r:.4f}")
                print(f"R²: {r2:.4f}")
                print(f"n = {n} compounds")
                
                # Rank correlation
                x_rank = [sorted(x).index(xi) + 1 for xi in x]
                y_rank = [sorted(y).index(yi) + 1 for yi in y]
                d_sq = sum((xr - yr)**2 for xr, yr in zip(x_rank, y_rank))
                spearman_rho = 1 - (6 * d_sq) / (n * (n**2 - 1))
                print(f"Spearman ρ: {spearman_rho:.4f}")
        
        # Summary table
        print("\n" + "=" * 60)
        print("RESULTS SUMMARY")
        print("=" * 60)
        print(f"{'ID':<6} {'IC50 (μM)':<12} {'pIC50':<10} {'Docking':<12} {'Potentiation'}")
        print("-" * 60)
        for r in sorted(results, key=lambda x: x['ic50']):
            print(f"{r['id']:<6} {r['ic50']:<12.2f} {r['pIC50']:<10.2f} {r['docking_energy']:<12.2f} {r['potentiation']:.0f}%")
    
    print("\nDone.")


if __name__ == "__main__":
    main()
