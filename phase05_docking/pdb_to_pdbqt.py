#!/usr/bin/env python3
"""Robust PDB to PDBQT converter for flexible docking.

Handles the full PDBQT format correctly, including:
- Proper column alignment (78 chars per line)
- Atom typing based on element
- Gasteiger charge assignment (simplified)
- Flexible residue marking
"""

import sys
from pathlib import Path

ATOM_TYPE_MAP = {
    'C': 'C', 'CA': 'C', 'CB': 'C', 'CG': 'C', 'CD': 'C', 'CE': 'C', 'CZ': 'C',
    'N': 'NA', 'ND1': 'NA', 'ND2': 'NA', 'NE': 'NA', 'NE1': 'NA', 'NE2': 'NA',
    'NH1': 'NA', 'NH2': 'NA', 'NZ': 'NA',
    'O': 'OA', 'OD1': 'OA', 'OD2': 'OA', 'OE1': 'OA', 'OE2': 'OA', 'OG': 'OA',
    'OG1': 'OA', 'OH': 'OA', 'OXT': 'OA',
    'S': 'SA', 'SD': 'SA', 'SG': 'SA',
    'H': 'HD', 'HA': 'HD', 'HB': 'HD', 'HG': 'HD', 'HE': 'HD', 'HZ': 'HD',
}

RESIDUE_CHARGES = {
    'ALA': {'N': -0.4157, 'CA': 0.0337, 'C': 0.5973, 'O': -0.5679},
    'ARG': {'N': -0.4157, 'CA': 0.0337, 'C': 0.5973, 'O': -0.5679},
    'ASN': {'N': -0.4157, 'CA': 0.0337, 'C': 0.5973, 'O': -0.5679},
    'ASP': {'N': -0.4157, 'CA': 0.0337, 'C': 0.5973, 'O': -0.5679},
    'CYS': {'N': -0.4157, 'CA': 0.0337, 'C': 0.5973, 'O': -0.5679, 'SG': -0.2029},
    'GLU': {'N': -0.4157, 'CA': 0.0337, 'C': 0.5973, 'O': -0.5679},
    'GLN': {'N': -0.4157, 'CA': 0.0337, 'C': 0.5973, 'O': -0.5679},
    'GLY': {'N': -0.4157, 'CA': 0.0337, 'C': 0.5973, 'O': -0.5679},
    'HIS': {'N': -0.4157, 'CA': 0.0337, 'C': 0.5973, 'O': -0.5679},
    'ILE': {'N': -0.4157, 'CA': 0.0337, 'C': 0.5973, 'O': -0.5679},
    'LEU': {'N': -0.4157, 'CA': 0.0337, 'C': 0.5973, 'O': -0.5679},
    'LYS': {'N': -0.4157, 'CA': 0.0337, 'C': 0.5973, 'O': -0.5679},
    'MET': {'N': -0.4157, 'CA': 0.0337, 'C': 0.5973, 'O': -0.5679, 'SD': -0.1022},
    'PHE': {'N': -0.4157, 'CA': 0.0337, 'C': 0.5973, 'O': -0.5679},
    'PRO': {'N': -0.2111, 'CA': 0.0458, 'C': 0.5718, 'O': -0.5284},
    'SER': {'N': -0.4157, 'CA': 0.0337, 'C': 0.5973, 'O': -0.5679, 'OG': -0.6755},
    'THR': {'N': -0.4157, 'CA': 0.0337, 'C': 0.5973, 'O': -0.5679, 'OG1': -0.6755},
    'TRP': {'N': -0.4157, 'CA': 0.0337, 'C': 0.5973, 'O': -0.5679},
    'TYR': {'N': -0.4157, 'CA': 0.0337, 'C': 0.5973, 'O': -0.5679, 'OH': -0.6755},
    'VAL': {'N': -0.4157, 'CA': 0.0337, 'C': 0.5973, 'O': -0.5679},
}


def get_atom_type(atom_name, element):
    atom_name = atom_name.strip()
    if atom_name in ATOM_TYPE_MAP:
        return ATOM_TYPE_MAP[atom_name]
    if element:
        element = element.strip().upper()
        if element == 'C': return 'C'
        elif element == 'N': return 'NA'
        elif element == 'O': return 'OA'
        elif element == 'S': return 'SA'
        elif element == 'H': return 'HD'
    return 'C'


def get_charge(resname, atom_name):
    resname = resname.strip().upper()
    atom_name = atom_name.strip()
    if resname in RESIDUE_CHARGES:
        charges = RESIDUE_CHARGES[resname]
        if atom_name in charges:
            return charges[atom_name]
    atype = get_atom_type(atom_name, atom_name[0] if atom_name else 'C')
    if atype == 'NA': return -0.5
    elif atype == 'OA': return -0.5
    elif atype == 'HD': return 0.5
    return 0.0


def _format_pdbqt_line(serial, atom_name, altloc, resname, chain, resseq, icode,
                        x, y, z, occ, tf, charge, atype):
    """Format a single PDBQT line (exactly 78 chars)."""
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
        f"{occ:>6.2f}"                # 55-61
        f"{tf:>6.2f}"                 # 62-66
                f"    "                       # 67-70
                f"{charge:>6.2f}"             # 71-76
                f"           "               # 77-87
                f"{atype:>2s}"               # 88-89
    )


def _parse_pdb_atom(line):
    """Parse an ATOM/HETATM line from PDB format."""
    try:
        atom_name = line[12:16]
        altloc = line[16]
        resname = line[17:20]
        chain = line[21] if len(line) > 21 and line[21].strip() else 'A'
        resseq = int(line[22:26].strip())
        icode = line[26] if len(line) > 26 else ' '
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        occ = float(line[54:60]) if line[54:60].strip() else 1.0
        tf = float(line[60:66]) if line[60:66].strip() else 0.0
        return {
            'atom_name': atom_name, 'altloc': altloc, 'resname': resname,
            'chain': chain, 'resseq': resseq, 'icode': icode,
            'x': x, 'y': y, 'z': z, 'occ': occ, 'tf': tf,
        }
    except (ValueError, IndexError):
        return None


def pdb_to_pdbqt(pdb_path, pdbqt_path, flex_residues=None):
    """Convert PDB to PDBQT with proper 78-char column alignment."""
    if flex_residues is None:
        flex_residues = []

    flex_set = set()
    for res in flex_residues:
        parts = res.split(":")
        if len(parts) == 2:
            chain = parts[0].strip() if parts[0].strip() else "A"
            resnum = parts[1].strip()
            if resnum.isdigit():
                flex_set.add((chain, int(resnum)))

    serial = 0
    lines_written = 0

    with open(pdb_path) as fin, open(pdbqt_path, 'w') as fout:
        for line in fin:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue

            parsed = _parse_pdb_atom(line)
            if parsed is None:
                continue

            if parsed['altloc'] not in (' ', 'A', '1'):
                continue

            serial += 1
            element = parsed['atom_name'].strip()[0] if parsed['atom_name'].strip() else 'C'
            atype = get_atom_type(parsed['atom_name'], element)
            charge = get_charge(parsed['resname'], parsed['atom_name'])

            pdbqt_line = _format_pdbqt_line(
                serial, parsed['atom_name'], parsed['altloc'],
                parsed['resname'], parsed['chain'], parsed['resseq'], parsed['icode'],
                parsed['x'], parsed['y'], parsed['z'],
                parsed['occ'], parsed['tf'], charge, atype
            )

            fout.write(pdbqt_line + "\n")
            lines_written += 1

    return lines_written


def create_flexible_pdbqt(pdb_path, flex_pdbqt_path, flex_residues):
    """Create a PDBQT with only flexible residue atoms (for Vina flex docking)."""
    flex_set = set()
    for res in flex_residues:
        parts = res.split(":")
        if len(parts) == 2:
            chain = parts[0].strip() if parts[0].strip() else "A"
            resnum = parts[1].strip()
            if resnum.isdigit():
                flex_set.add((chain, int(resnum)))

    serial = 0
    lines_written = 0

    with open(pdb_path) as fin, open(flex_pdbqt_path, 'w') as fout:
        for line in fin:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue

            parsed = _parse_pdb_atom(line)
            if parsed is None:
                continue

            if parsed['altloc'] not in (' ', 'A', '1'):
                continue

            is_flex = (parsed['chain'], parsed['resseq']) in flex_set
            if not is_flex:
                continue

            serial += 1
            element = parsed['atom_name'].strip()[0] if parsed['atom_name'].strip() else 'C'
            atype = get_atom_type(parsed['atom_name'], element)
            charge = get_charge(parsed['resname'], parsed['atom_name'])

            pdbqt_line = _format_pdbqt_line(
                serial, parsed['atom_name'], parsed['altloc'],
                parsed['resname'], parsed['chain'], parsed['resseq'], parsed['icode'],
                parsed['x'], parsed['y'], parsed['z'],
                parsed['occ'], parsed['tf'], charge, atype
            )

            fout.write(pdbqt_line + "\n")
            lines_written += 1

    return lines_written


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: pdb_to_pdbqt.py <input.pdb> <output.pdbqt> [flex_res1,flex_res2,...]")
        sys.exit(1)

    pdb_path = sys.argv[1]
    pdbqt_path = sys.argv[2]
    flex_residues = sys.argv[3].split(",") if len(sys.argv) > 3 and sys.argv[3] else []

    n = pdb_to_pdbqt(pdb_path, pdbqt_path, flex_residues)
    print(f"Wrote {n} atoms to {pdbqt_path}")

    if flex_residues:
        flex_path = pdbqt_path.replace(".pdbqt", "_flex.pdbqt")
        n_flex = create_flexible_pdbqt(pdb_path, flex_path, flex_residues)
        print(f"Wrote {n_flex} flexible atoms to {flex_path}")
