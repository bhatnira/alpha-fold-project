#!/usr/bin/env python3
"""Phase V — Ensemble docking with AutoDock Vina.

Docks ligands into candidate pockets using the receptor ensemble.
Uses rigid docking (Vina 1.2.7 Python binding does not support flex).
"""

import csv, sys, os, subprocess, math
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
OUT  = Path(__file__).resolve().parent
RESULTS = OUT / "results"
RESULTS.mkdir(exist_ok=True)


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
    """Format exactly 78-column PDBQT ATOM line."""
    return (
        f"ATOM  "                     # 1-6
        f"{serial:>5d}"               # 7-11
        f" "                          # 12
        f"{atom_name:<4s}"            # 13-16
        f"{altloc}"                   # 17
        f"{resname:>3s}"              # 18-20
        f" "                          # 21
        f"{chain}"                    # 22
        f"{resseq:>4d}"               # 23-26
        f"{icode}"                    # 27
        f"   "                        # 28-30
        f"{x:>8.3f}"                  # 31-38
        f"{y:>8.3f}"                  # 39-46
        f"{z:>8.3f}"                  # 47-54
        f"{occ:>6.2f}"                # 55-60
        f"{tf:>6.2f}"                 # 61-66
        f"     "                      # 67-71
        f"{charge:>5.2f}"             # 72-76
        f" "                          # 77
        f"{atype:<2s}"                # 78
    )


def prepare_receptor_pdbqt(pdb_path, output_path):
    """Convert PDB to valid 78-column PDBQT with proper AutoDock atom types."""
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
    fout_name = Path(output_path).name
    print(f"    {fout_name}: {serial} atoms written")
    return serial


def prepare_ligand_pdbqt(smiles, output_path):
    """Convert SMILES to PDBQT via RDKit (3D embed + Gasteiger charges)."""
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
                sym = atom.GetAtomicNum()
                elem = atom.GetSymbol()
                gc = float(atom.GetDoubleProp('_GasteigerCharge'))
                if math.isnan(gc):
                    gc = 0.0

                atype_map = {'C': 'C', 'N': 'NA', 'O': 'OA', 'S': 'SA', 'H': 'HD',
                             'F': 'F', 'I': 'I'}
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


def run_vina_docking(receptor_pdbqt, ligand_pdbqt, center, box_size, output_path, exhaustiveness=8, n_poses=5):
    """Run AutoDock Vina rigid docking."""
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


def main():
    print("Phase V — Ensemble docking (rigid, Vina 1.2.7)")
    print("=" * 60)

    receptors = [
        ROOT / "data/deliverable/05_structures/2to3/a9a10_2to3_apo_reference.pdb",
        ROOT / "data/deliverable/05_structures/3to2/a9a10_3to2_apo_reference.pdb",
    ]

    ligands = [
        {"name": "ascorbate", "smiles": "OC[C@H](O)[C@H]1OC(=O)C(O)=C1O"},
        {"name": "acetate", "smiles": "CC(=O)[O-]"},
        {"name": "O-ethyl_ascorbate", "smiles": "CCOC(=O)C(=O)OC[C@@H](O)[C@@H]1OC(=O)C(O)=C1O"},
    ]

    sites_csv = ROOT / "data/deliverable/05_structures/2to3/sites_2to3.csv"
    sites = []
    if sites_csv.exists():
        with open(sites_csv) as f:
            for row in csv.DictReader(f):
                sites.append(row)
        print(f"Loaded {len(sites)} sites from sites_2to3.csv")
    else:
        print("No sites CSV, using default center")
        sites = [{"site": "default", "copy": "1", "centroid_x": "-24.2",
                   "centroid_y": "30.97", "centroid_z": "21.77",
                   "centroid_spread_A": "8"}]

    results = []
    for rec_path in receptors:
        if not rec_path.exists():
            print(f"  Skipping {rec_path.name} (not found)")
            continue

        rec_name = rec_path.stem
        rec_pdbqt = RESULTS / f"{rec_name}.pdbqt"

        print(f"\nPreparing receptor {rec_name}...")
        prepare_receptor_pdbqt(str(rec_path), str(rec_pdbqt))

        for lig in ligands:
            lig_pdbqt = RESULTS / f"{lig['name']}.pdbqt"
            if not lig_pdbqt.exists() or lig_pdbqt.stat().st_size < 10:
                print(f"  Preparing ligand {lig['name']}...")
                prepare_ligand_pdbqt(lig['smiles'], str(lig_pdbqt))

            if not lig_pdbqt.exists() or lig_pdbqt.stat().st_size < 10:
                print(f"  Ligand {lig['name']} preparation failed, skipping")
                continue

            for site_row in sites:
                site_id = site_row.get("site", "default")
                copy = site_row.get("copy", "1")
                try:
                    cx = float(site_row["centroid_x"])
                    cy = float(site_row["centroid_y"])
                    cz = float(site_row["centroid_z"])
                except (ValueError, KeyError):
                    continue
                try:
                    spread = float(site_row.get("centroid_spread_A", 8))
                except ValueError:
                    spread = 8
                box_dim = max(15, spread * 2 + 10)
                box_size = (box_dim, box_dim, box_dim)

                out_pdb = RESULTS / f"{rec_name}_site{site_id}_copy{copy}_{lig['name']}_docked.pdb"

                print(f"  Docking {lig['name']} -> site {site_id} copy {copy}...", end=" ", flush=True)
                success, energy_val = run_vina_docking(
                    str(rec_pdbqt), str(lig_pdbqt),
                    (cx, cy, cz), box_size, str(out_pdb)
                )
                results.append({
                    "receptor": rec_name,
                    "site_id": site_id,
                    "copy": copy,
                    "ligand": lig['name'],
                    "smiles": lig['smiles'],
                    "center_x": cx,
                    "center_y": cy,
                    "center_z": cz,
                    "box_size": box_dim,
                    "success": success,
                    "binding_energy": energy_val,
                    "output": str(out_pdb) if success else "",
                })
                print(f"{'OK' if success else 'FAILED'}" + (f"  E={energy_val:.2f}" if energy_val else ""))

    if results:
        csv_path = OUT / "docking_results.csv"
        with open(csv_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=results[0].keys())
            writer.writeheader()
            writer.writerows(results)
        print(f"\nResults: {csv_path}")

    n_ok = sum(1 for r in results if r["success"])
    print(f"\n{n_ok}/{len(results)} dockings successful")
    print("Done.")


if __name__ == "__main__":
    main()
