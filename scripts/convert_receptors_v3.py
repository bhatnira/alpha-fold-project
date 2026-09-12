#!/usr/bin/env python3
"""Convert CIF receptors to PDBQT using strict PDBQT column format."""
import sys, os
sys.path.insert(0, os.path.expanduser("~/.local/lib/python3.9/site-packages"))

from openbabel import openbabel as ob
from pathlib import Path

SCRATCH = Path("/cluster/scratch/nbhatt04/allostery")
OUTPUT = Path("/cluster/home/nbhatt04/lean_pipeline/09_docking")
OUTPUT.mkdir(exist_ok=True)

# AutoDock atom type map
AD_TYPE = {
    'C': 'C', 'N': 'NA', 'O': 'OA', 'S': 'SA', 'H': 'HD',
    'FE': 'FE', 'CA': 'CA', 'MG': 'MG', 'ZN': 'ZN',
    'BR': 'Br', 'CL': 'Cl', 'F': 'F', 'I': 'I',
}

def pdb_atom_line_to_pdbqt(line):
    """Convert a single PDB ATOM/HETATM line to PDBQT format with strict columns."""
    if not line.startswith(('ATOM', 'HETATM')):
        return None
    
    # Extract PDB fields (1-indexed, Python 0-indexed)
    record   = line[0:6].rstrip()       # 1-6
    serial   = line[6:11].rstrip()      # 7-11
    name     = line[12:16]              # 13-16 (keep spacing)
    altloc   = line[16:17]              # 17
    resname  = line[17:20]              # 18-20
    chain    = line[21:22]              # 22
    resseq   = line[22:26].rstrip()     # 23-26
    icode    = line[26:27]              # 27
    x        = line[30:38]              # 31-38
    y        = line[38:46]              # 39-46
    z        = line[46:54]              # 47-54
    occ      = line[54:60]              # 55-60
    bfactor  = line[60:66]              # 61-66
    
    # Determine element from atom name
    atom_name_stripped = name.strip()
    element = ''
    for c in atom_name_stripped:
        if c.isalpha():
            element += c
        else:
            break
    if not element:
        element = 'C'
    
    # Map to AutoDock type
    ad_type = AD_TYPE.get(element.upper(), element)
    
    # Build PDBQT line with strict column format
    # Columns 1-6:   Record (ATOM  )
    # Columns 7-11:  Serial (right-justified)
    # Column  12:    Space
    # Columns 13-16: Atom name
    # Column  17:    AltLoc
    # Columns 18-20: Residue name
    # Column  21:    Space
    # Column  22:    Chain
    # Columns 23-26: ResSeq (right-justified)
    # Column  27:    ICode
    # Columns 28-30: Spaces
    # Columns 31-38: X (8.3f)
    # Columns 39-46: Y (8.3f)
    # Columns 47-54: Z (8.3f)
    # Columns 55-60: Occupancy (6.2f)
    # Columns 61-66: B-factor (6.2f)
    # Columns 67-76: Spaces
    # Columns 77-78: Atom type (right-justified)
    
    line_out = (
        f"{record:<6}"
        f"{serial:>5}"
        f" "
        f"{name}"
        f"{altloc}"
        f"{resname}"
        f" "
        f"{chain}"
        f"{resseq:>4}"
        f"{icode}"
        f"   "
        f"{x}"
        f"{y}"
        f"{z}"
        f"{occ}"
        f"{bfactor}"
        f"          "
        f"{ad_type:>2}"
        f"\n"
    )
    
    return line_out

def convert_cif_to_pdbqt(cif_path, pdbqt_path):
    tmp_pdb = str(pdbqt_path) + ".tmp.pdb"
    
    conv = ob.OBConversion()
    conv.SetInFormat("cif")
    conv.SetOutFormat("pdb")
    
    mol = ob.OBMol()
    success = conv.ReadFile(mol, str(cif_path))
    
    if not success or mol.NumAtoms() == 0:
        print(f"  ERROR: Could not read {cif_path}")
        return False
    
    print(f"  Read {mol.NumAtoms()} atoms")
    conv.WriteFile(mol, tmp_pdb)
    
    with open(tmp_pdb, 'r') as fin, open(pdbqt_path, 'w') as fout:
        for line in fin:
            if line.startswith(('ATOM', 'HETATM')):
                pdbqt_line = pdb_atom_line_to_pdbqt(line)
                if pdbqt_line:
                    fout.write(pdbqt_line)
            elif line.startswith(('TER', 'END')):
                fout.write(line)
    
    os.unlink(tmp_pdb)
    
    atom_count = sum(1 for l in open(pdbqt_path) if l.startswith(('ATOM', 'HETATM')))
    print(f"  Wrote {atom_count} atoms to PDBQT")
    return True

for stoich in ["2to3", "3to2"]:
    cif = SCRATCH / f"af3_outputs/a9a10_{stoich}_APO/a9a10_{stoich}_APO_model.cif"
    pdbqt = OUTPUT / f"receptor_{stoich}.pdbqt"
    
    if not cif.exists():
        print(f"No CIF for {stoich}")
        continue
    
    print(f"Converting: {cif.name}")
    convert_cif_to_pdbqt(cif, pdbqt)

print("\nDone!")
