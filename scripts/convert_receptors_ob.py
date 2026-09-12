#!/usr/bin/env python3
"""Convert CIF to PDBQT using OpenBabel's native PDBQT writer."""
import sys, os
sys.path.insert(0, os.path.expanduser("~/.local/lib/python3.9/site-packages"))

from openbabel import openbabel as ob
from pathlib import Path

SCRATCH = Path("/cluster/scratch/nbhatt04/allostery")
OUTPUT = Path("/cluster/home/nbhatt04/lean_pipeline/09_docking")
OUTPUT.mkdir(exist_ok=True)

for stoich in ["2to3", "3to2"]:
    cif = SCRATCH / f"af3_outputs/a9a10_{stoich}_APO/a9a10_{stoich}_APO_model.cif"
    pdbqt = OUTPUT / f"receptor_{stoich}.pdbqt"
    
    if not cif.exists():
        print(f"No CIF for {stoich}")
        continue
    
    print(f"Converting: {cif.name}")
    
    # Read CIF
    conv_in = ob.OBConversion()
    conv_in.SetInFormat("cif")
    mol = ob.OBMol()
    success = conv_in.ReadFile(mol, str(cif))
    
    if not success or mol.NumAtoms() == 0:
        print(f"  ERROR reading CIF")
        continue
    
    print(f"  Read {mol.NumAtoms()} atoms")
    
    # Add hydrogens
    mol.AddHydrogens()
    print(f"  After adding H: {mol.NumAtoms()} atoms")
    
    # Write PDBQT directly with OpenBabel's PDBQT writer
    conv_out = ob.OBConversion()
    conv_out.SetOutFormat("pdbqt")
    conv_out.WriteFile(mol, str(pdbqt))
    
    # Verify
    with open(pdbqt, 'r') as f:
        lines = f.readlines()
    atom_lines = [l for l in lines if l.startswith(('ATOM', 'HETATM'))]
    print(f"  PDBQT: {len(atom_lines)} atoms")
    if atom_lines:
        print(f"  First: {atom_lines[0].rstrip()}")
        print(f"  Last:  {atom_lines[-1].rstrip()}")

print("\nDone!")
