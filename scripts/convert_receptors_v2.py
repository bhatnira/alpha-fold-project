#!/usr/bin/env python3
"""Convert CIF receptors to proper PDBQT format for AutoDock Vina."""
import sys, os
sys.path.insert(0, os.path.expanduser("~/.local/lib/python3.9/site-packages"))

from openbabel import openbabel as ob
from pathlib import Path

SCRATCH = Path("/cluster/scratch/nbhatt04/allostery")
OUTPUT = Path("/cluster/home/nbhatt04/lean_pipeline/09_docking")
OUTPUT.mkdir(exist_ok=True)

def convert_cif_to_pdbqt(cif_path, pdbqt_path):
    """Convert CIF to PDBQT with proper AutoDock atom types."""
    tmp_pdb = str(pdbqt_path) + ".tmp.pdb"
    
    conv = ob.OBConversion()
    conv.SetInFormat("cif")
    conv.SetOutFormat("pdb")
    
    mol = ob.OBMol()
    success = conv.ReadFile(mol, str(cif_path))
    
    if not success or mol.NumAtoms() == 0:
        print(f"  ERROR: Could not read {cif_path}")
        return False
    
    print(f"  Read {mol.NumAtoms()} atoms from CIF")
    conv.WriteFile(mol, tmp_pdb)
    
    # AutoDock atom type mapping
    ad_types = {
        'C': 'C', 'N': 'NA', 'O': 'OA', 'S': 'SA', 'H': 'HD',
        'FE': 'FE', 'CA': 'CA', 'MG': 'MG', 'ZN': 'ZN',
        'BR': 'Br', 'CL': 'Cl', 'F': 'F', 'I': 'I',
    }
    
    # Read PDB and create proper PDBQT
    with open(tmp_pdb, 'r') as fin, open(pdbqt_path, 'w') as fout:
        for line in fin:
            if line.startswith(('ATOM', 'HETATM')):
                # PDB format: columns 1-6 record, 7-11 serial, 13-16 name, 
                # 17 altLoc, 18-20 resName, 22 chain, 23-26 resSeq, 27 iCode
                # 31-38 x, 39-46 y, 47-54 z, 55-60 occupancy, 61-66 tempFactor
                record = line[0:6]
                serial = line[6:11]
                name = line[12:16]
                altLoc = line[16:17]
                resName = line[17:20]
                chain = line[21:22]
                resSeq = line[22:26]
                iCode = line[26:27]
                x = line[30:38]
                y = line[38:46]
                z = line[46:54]
                occupancy = line[54:60]
                tempFactor = line[60:66]
                
                # Determine element from atom name
                element = name.strip()[0]
                if element.isdigit():
                    element = name.strip()[1] if len(name.strip()) > 1 else 'C'
                
                # Map to AutoDock type
                ad_type = ad_types.get(element, element)
                
                # Write proper PDBQT line
                # Format: RECORD SER NAME ALT RES CHAIN RESSEQ ICODE X Y Z OCC ADTYPE CHARGE
                pdbqt_line = (
                    f"{record}{serial} {name}{altLoc}{resName} {chain}{resSeq}{iCode}"
                    f"    {x}{y}{z}{occupancy}{tempFactor}"
                    f"                    {ad_type}\n"
                )
                fout.write(pdbqt_line)
            elif line.startswith('TER'):
                fout.write(line)
            elif line.startswith('END'):
                fout.write(line)
    
    # Cleanup
    os.unlink(tmp_pdb)
    
    # Count atoms
    atom_count = sum(1 for line in open(pdbqt_path) if line.startswith(('ATOM', 'HETATM')))
    print(f"  Wrote PDBQT: {pdbqt_path} ({atom_count} atoms)")
    return True

for stoich in ["2to3", "3to2"]:
    cif = SCRATCH / f"af3_outputs/a9a10_{stoich}_APO/a9a10_{stoich}_APO_model.cif"
    pdbqt = OUTPUT / f"receptor_{stoich}.pdbqt"
    
    if not cif.exists():
        alt = list((SCRATCH / "af3_outputs").glob(f"a9a10_{stoich}*/*model.cif"))
        if alt:
            cif = alt[0]
        else:
            print(f"No CIF for {stoich}")
            continue
    
    print(f"Converting: {cif.name}")
    convert_cif_to_pdbqt(cif, pdbqt)

print("\nDone!")
