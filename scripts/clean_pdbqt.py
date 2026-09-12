#!/usr/bin/env python3
"""Clean receptor PDBQT files for Vina compatibility."""
from pathlib import Path

OUTPUT = Path("/cluster/home/nbhatt04/lean_pipeline/09_docking")

for stoich in ["2to3", "3to2"]:
    pdbqt = OUTPUT / f"receptor_{stoich}.pdbqt"
    if not pdbqt.exists():
        continue
    
    # Read and filter
    lines = []
    with open(pdbqt, 'r') as f:
        for line in f:
            # Vina only accepts ATOM, HETATM, END, ROOT, ENDROOT, BRANCH, ENDBRANCH, TORSDOF
            if line.startswith(('ATOM', 'HETATM', 'END', 'ROOT', 'ENDROOT', 'BRANCH', 'ENDBRANCH', 'TORSDOF', 'REMARK')):
                lines.append(line)
            elif line.strip() == '':
                continue
            else:
                # Skip non-standard lines (COMPND, HEADER, etc.)
                pass
    
    # Write cleaned file
    with open(pdbqt, 'w') as f:
        f.writelines(lines)
    
    print(f"Cleaned: {pdbqt} ({len(lines)} lines)")

print("Done!")
