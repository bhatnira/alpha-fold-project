#!/usr/bin/env python3
"""Convert CIF receptor files to PDBQT for AutoDock Vina."""
import sys, os
sys.path.insert(0, os.path.expanduser("~/.local/lib/python3.9/site-packages"))

from openbabel import openbabel as ob
from pathlib import Path

SCRATCH = Path("/cluster/scratch/nbhatt04/allostery")
OUTPUT = Path("/cluster/home/nbhatt04/lean_pipeline/09_docking")
OUTPUT.mkdir(exist_ok=True)

def convert_cif_to_pdbqt(cif_path, pdbqt_path):
    pdb_path = str(pdbqt_path).replace(".pdbqt", ".pdb")
    
    conv = ob.OBConversion()
    conv.SetInFormat("cif")
    conv.SetOutFormat("pdb")
    
    mol = ob.OBMol()
    success = conv.ReadFile(mol, str(cif_path))
    
    if not success or mol.NumAtoms() == 0:
        print(f"  ERROR: Could not read {cif_path}")
        return False
    
    print(f"  Read {mol.NumAtoms()} atoms from CIF")
    conv.WriteFile(mol, pdb_path)
    print(f"  Wrote PDB: {pdb_path}")
    
    # PDB to PDBQT
    with open(pdb_path, "r") as fin, open(pdbqt_path, "w") as fout:
        for line in fin:
            if line.startswith(("ATOM", "HETATM")):
                fout.write(line.rstrip() + "    0.000\n")
            else:
                fout.write(line)
    
    print(f"  Wrote PDBQT: {pdbqt_path}")
    return True

for stoich in ["2to3", "3to2"]:
    cif = SCRATCH / f"af3_outputs/a9a10_{stoich}_APO/a9a10_{stoich}_APO_model.cif"
    pdbqt = OUTPUT / f"receptor_{stoich}.pdbqt"
    
    if pdbqt.exists():
        print(f"Already exists: {pdbqt}")
        continue
    
    if not cif.exists():
        alt = list((SCRATCH / "af3_outputs").glob(f"a9a10_{stoich}*/*model.cif"))
        if alt:
            cif = alt[0]
            print(f"Using: {cif}")
        else:
            print(f"No CIF for {stoich}")
            continue
    
    print(f"Converting: {cif.name} -> {pdbqt.name}")
    convert_cif_to_pdbqt(cif, pdbqt)

print("\nReceptor conversion complete!")
