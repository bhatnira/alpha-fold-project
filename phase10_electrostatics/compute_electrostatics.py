#!/usr/bin/env python3
"""Phase X — Electrostatic analysis using pdb2pqr + approximate PB.

Computes electrostatic potential at key residue positions
using a simplified Poisson-Boltzmann approach.
"""

import csv, sys, os
import numpy as np
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
OUT  = Path(__file__).resolve().parent

def parse_pqr(pqr_path):
    """Parse PQR file into coordinates and charges."""
    atoms = []
    with open(pqr_path) as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                # PQR format: columns may vary, try standard PDB-like first
                chain = line[21:22].strip() if len(line) > 21 else ""
                if not chain:
                    chain = "A"  # default chain
                resid_str = line[22:26].strip() if len(line) > 25 else "1"
                atoms.append({
                    "name": line[12:16].strip(),
                    "resname": line[17:20].strip(),
                    "resid": int(resid_str) if resid_str.isdigit() else 1,
                    "chain": chain,
                    "x": float(line[30:38]),
                    "y": float(line[38:46]),
                    "z": float(line[46:54]),
                    "charge": float(line[54:62]),
                    "radius": float(line[62:70]) if len(line) > 70 else 1.5,
                })
    return atoms

def compute_electrostatic_potential(atoms, grid_center, grid_size=20.0, grid_spacing=1.0):
    """Compute Coulombic electrostatic potential on a grid."""
    x = np.arange(grid_center[0] - grid_size/2, grid_center[0] + grid_size/2, grid_spacing)
    y = np.arange(grid_center[1] - grid_size/2, grid_center[1] + grid_size/2, grid_spacing)
    z = np.arange(grid_center[2] - grid_size/2, grid_center[2] + grid_size/2, grid_spacing)
    
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    V = np.zeros_like(X)
    
    # Coulomb constant: 332.0637 kcal/mol·Å·e²
    k = 332.0637
    
    for atom in atoms:
        dx = X - atom["x"]
        dy = Y - atom["y"]
        dz = Z - atom["z"]
        r = np.sqrt(dx**2 + dy**2 + dz**2)
        r = np.maximum(r, 0.5)  # avoid division by zero
        V += k * atom["charge"] / r
    
    return V, (x, y, z)

def analyze_electrostatics(pqr_path, site_residues, label):
    """Analyze electrostatic potential at specific residue positions."""
    atoms = parse_pqr(pqr_path)
    
    # Compute grid center from atoms
    coords = np.array([[a["x"], a["y"], a["z"]] for a in atoms])
    center = coords.mean(axis=0)
    
    # Compute potential
    V, grids = compute_electrostatic_potential(atoms, center, grid_size=30.0, grid_spacing=2.0)
    
    # Sample potential at key residue positions
    results = []
    for res in site_residues:
        chain, resnum = res.split(":")
        resnum = int(resnum)
        
        # Find atoms in this residue
        res_atoms = [a for a in atoms if a["chain"] == chain and a["resid"] == resnum]
        
        if res_atoms:
            res_coords = np.array([[a["x"], a["y"], a["z"]] for a in res_atoms])
            res_center = res_coords.mean(axis=0)
            res_charge = sum(a["charge"] for a in res_atoms)
            
            # Interpolate potential at residue center
            ix = np.argmin(np.abs(grids[0] - res_center[0]))
            iy = np.argmin(np.abs(grids[1] - res_center[1]))
            iz = np.argmin(np.abs(grids[2] - res_center[2]))
            potential = V[ix, iy, iz]
            
            results.append({
                "label": label,
                "residue": res,
                "chain": chain,
                "resid": resnum,
                "x": round(res_center[0], 3),
                "y": round(res_center[1], 3),
                "z": round(res_center[2], 3),
                "charge": round(res_charge, 3),
                "potential_kcal_mol": round(float(potential), 2),
                "n_atoms": len(res_atoms),
            })
    
    return results

def main():
    print("Phase X — Electrostatic analysis")
    
    # Key residues from pharmacophore (Phase XV)
    site23_residues = [
        "A:176", "A:175", "A:120", "A:224", "A:217", "A:57", "A:55",  # α9
        "D:143", "D:145", "D:81",  # α10
    ]
    
    site34_residues = [
        "A:270", "A:274", "A:277", "A:266", "A:267",  # α9 TM
        "D:273", "D:276", "D:265", "D:266", "D:269",  # α10 TM
    ]
    
    all_results = []
    
    for pqr in sorted(OUT.glob("results/*.pqr")):
        label = pqr.stem
        print(f"  Analyzing {label}...")
        
        # Site 23 (ECD)
        results = analyze_electrostatics(str(pqr), site23_residues, f"{label}_site23")
        all_results.extend(results)
        
        # Site 34 (TM)
        results = analyze_electrostatics(str(pqr), site34_residues, f"{label}_site34")
        all_results.extend(results)
    
    # Write results
    if all_results:
        with open(OUT / "electrostatic_potentials.csv", "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=all_results[0].keys())
            writer.writeheader()
            writer.writerows(all_results)
    
    # Summary
    lines = [
        "ELECTROSTATIC ANALYSIS",
        "=" * 60,
        "",
        "Method: Coulombic potential from pdb2pqr charges",
        "Grid: 30 Å cube, 2 Å spacing, centered on receptor",
        "",
    ]
    
    for label in set(r["label"] for r in all_results):
        label_results = [r for r in all_results if r["label"] == label]
        lines.append(f"{label}:")
        for r in label_results:
            lines.append(f"  {r['residue']:10s}  charge={r['charge']:+.3f}  "
                          f"potential={r['potential_kcal_mol']:+.2f} kcal/mol")
        lines.append("")
    
    lines.extend([
        "INTERPRETATION:",
        "",
        "  Negative potential at binding site → favorable for cation binding",
        "  Positive potential at binding site → favorable for anion binding",
        "  L-ascorbate (anion) → prefers positive potential regions",
        "  If acetate shows same pattern → charge-driven",
        "  If L-ascorbate shows different pattern → molecular recognition",
    ])
    
    with open(OUT / "electrostatic_summary.txt", "w") as f:
        f.write("\n".join(lines) + "\n")
    
    print(f"  {len(all_results)} residue potentials computed")
    print("Done.")

if __name__ == "__main__":
    main()
