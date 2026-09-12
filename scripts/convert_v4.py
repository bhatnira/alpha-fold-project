#!/usr/bin/env python3
"""Convert CIF to PDBQT with exact PDBQT column alignment."""
import sys, os, re
sys.path.insert(0, os.path.expanduser("~/.local/lib/python3.9/site-packages"))

from openbabel import openbabel as ob
from pathlib import Path

SCRATCH = Path("/cluster/scratch/nbhatt04/allostery")
OUTPUT = Path("/cluster/home/nbhatt04/lean_pipeline/09_docking")
OUTPUT.mkdir(exist_ok=True)

VALID_ELEMENTS = {'H','C','N','O','S','P','F','Cl','Br','I','FE','ZN','CA','MG','MN','CU','CO','NA','SE','MO'}

AD = {'C':'C','N':'NA','O':'OA','S':'SA','H':'HD','FE':'FE','CA':'CA','MG':'MG','ZN':'ZN','BR':'Br','CL':'Cl','F':'F','I':'I',
      'P':'P','MN':'MN','CU':'CU','CO':'CO','NA':'NA','SE':'SE'}

def element_from_name(aname):
    s = aname.strip()
    if not s:
        return 'C'
    # In PDB, atom names like " CA " mean carbon alpha, " CB " means carbon beta
    # The element is determined by the first letter(s) that form a valid element
    # For standard amino acids: N, CA, C, O, CB, CG, CD, CE, CZ, OG, SG, etc.
    # The element is always the first character (or first two for two-letter elements)
    
    # Strip leading digits (like "1HB" -> "HB")
    while s and s[0].isdigit():
        s = s[1:]
    
    if not s:
        return 'C'
    
    # Try two-letter element first
    two = s[:2]
    if two in VALID_ELEMENTS:
        return two
    # Try single letter
    one = s[0]
    if one in VALID_ELEMENTS:
        return one
    
    return 'C'

def write_pdbqt_from_pdb(pdb_path, pdbqt_path):
    """Read PDB ATOM/HETATM lines and write strict PDBQT format."""
    with open(pdb_path, 'r') as fin, open(pdbqt_path, 'w') as fout:
        for line in fin:
            if not line.startswith(('ATOM', 'HETATM')):
                if line.startswith(('TER', 'END')):
                    fout.write(line.rstrip() + '\n')
                continue

            # Parse PDB fixed-width columns (0-indexed)
            record   = line[0:6]        # "ATOM  " or "HETATM"
            serial_s = line[6:11].strip()
            name     = line[12:16]      # 4 chars, keep as-is
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

            elem = element_from_name(name)
            ad_type = AD.get(elem.upper(), elem)

            # Strict PDBQT format using fixed-width formatting
            # PDBQT columns (1-indexed):
            # 1-6: record, 7-11: serial, 12: space
            # 13-16: name, 17: altloc, 18-20: resname
            # 21: space, 22: chain, 23-26: resseq, 27: icode
            # 28-30: spaces, 31-38: x, 39-46: y, 47-54: z
            # 55-60: occ, 61-66: bfac, 67-76: spaces, 77-78: adtype
            out = f"{record:<6.6}{serial:>5d} {name:<4.4}{altloc}{resname:<3.3} {chain}{resseq:>4d}{icode}   {x:8.3f}{y:8.3f}{z:8.3f}{occ:6.2f}{bfac:6.2f}          {ad_type:>2s}\n"
            fout.write(out)

for stoich in ["2to3", "3to2"]:
    cif = SCRATCH / f"af3_outputs/a9a10_{stoich}_APO/a9a10_{stoich}_APO_model.cif"
    pdbqt = OUTPUT / f"receptor_{stoich}.pdbqt"
    
    if not cif.exists():
        print(f"No CIF for {stoich}")
        continue
    
    # Step 1: CIF -> PDB using OpenBabel
    tmp_pdb = str(pdbqt) + ".tmp.pdb"
    conv = ob.OBConversion()
    conv.SetInFormat("cif")
    conv.SetOutFormat("pdb")
    mol = ob.OBMol()
    ok = conv.ReadFile(mol, str(cif))
    if not ok or mol.NumAtoms() == 0:
        print(f"  ERROR reading {cif}")
        continue
    print(f"{stoich}: Read {mol.NumAtoms()} atoms from CIF")
    conv.WriteFile(mol, tmp_pdb)
    
    # Step 2: PDB -> strict PDBQT
    write_pdbqt_from_pdb(tmp_pdb, pdbqt)
    
    # Verify
    with open(pdbqt) as f:
        lines = [l for l in f if l.startswith(('ATOM','HETATM'))]
    print(f"  Wrote {len(lines)} atoms")
    print(f"  Sample: {lines[0].rstrip()}")
    
    # Check column positions
    sample = lines[0]
    print(f"  X at col 31-38: '{sample[30:38]}'")
    print(f"  Y at col 39-46: '{sample[38:46]}'")
    print(f"  Z at col 47-54: '{sample[46:54]}'")
    print(f"  Type at col 77-78: '{sample[76:78]}'")
    print(f"  Line length: {len(sample.rstrip())}")
    
    os.unlink(tmp_pdb)

print("\nDone!")
