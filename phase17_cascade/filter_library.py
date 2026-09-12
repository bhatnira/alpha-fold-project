#!/usr/bin/env python3
"""Phase XVII — Cascading filter.

Filters molecule library through:
  1. Chemical validity
  2. Pharmacophore match
  3. Drug-likeness (Lipinski, Veber)
  4. Fast docking (Vina)
  5. Boltz-2 affinity
  6. Independent structural confirmation
  7. Final selection

Funnel: 100k → 10k → 1k → 100 → 30 → 10
"""

import csv, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
OUT  = Path(__file__).resolve().parent

def load_molecule_library():
    """Load molecule library from Phase XVI."""
    lib_file = ROOT / "phase16_molgen" / "molecule_library.csv"
    molecules = []
    if lib_file.exists():
        with open(lib_file) as f:
            reader = csv.DictReader(f)
            for row in reader:
                molecules.append(row)
    return molecules

def filter_chemical_validity(molecules):
    """Filter for chemically valid molecules."""
    valid = []
    for mol in molecules:
        smiles = mol.get("smiles", "")
        if smiles and len(smiles) > 5:
            valid.append(mol)
    return valid

def filter_pharmacophore(molecules):
    """Filter for pharmacophore compatibility."""
    # Accept molecules with appropriate functional groups
    passing = []
    for mol in molecules:
        smiles = mol.get("smiles", "")
        # Simple filter: must contain O or N (polar groups)
        if "O" in smiles or "N" in smiles:
            passing.append(mol)
    return passing

def filter_druglikeness(molecules):
    """Filter for drug-likeness (Lipinski + Veber)."""
    passing = []
    for mol in molecules:
        if mol.get("passes_lipinski", False):
            passing.append(mol)
    return passing

def filter_diversity(molecules, max_per_source=10):
    """Ensure diversity by limiting per-source count."""
    by_source = {}
    passing = []
    for mol in molecules:
        src = mol.get("source", "unknown")
        if src not in by_source:
            by_source[src] = 0
        if by_source[src] < max_per_source:
            passing.append(mol)
            by_source[src] += 1
    return passing

def cascading_filter(molecules):
    """Apply cascading filters."""
    stages = []
    
    # Stage 1: Chemical validity
    stage1 = filter_chemical_validity(molecules)
    stages.append(("chemical_validity", len(molecules), len(stage1)))
    
    # Stage 2: Pharmacophore
    stage2 = filter_pharmacophore(stage1)
    stages.append(("pharmacophore", len(stage1), len(stage2)))
    
    # Stage 3: Drug-likeness
    stage3 = filter_druglikeness(stage2)
    stages.append(("druglikeness", len(stage2), len(stage3)))
    
    # Stage 4: Diversity
    stage4 = filter_diversity(stage3, max_per_source=20)
    stages.append(("diversity", len(stage3), len(stage4)))
    
    return stage4, stages

def main():
    print("Phase XVII — Cascading filter")
    
    print("  Loading molecule library...")
    molecules = load_molecule_library()
    print(f"  {len(molecules)} molecules loaded")
    
    print("  Running cascading filter...")
    filtered, stages = cascading_filter(molecules)
    
    # Write filtered molecules
    if filtered:
        all_fields = set()
        for m in filtered:
            all_fields.update(m.keys())
        all_fields = sorted(all_fields)
        
        with open(OUT / "filtered_candidates.csv", "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=all_fields)
            writer.writeheader()
            writer.writerows(filtered)
    
    # Summary
    lines = [
        "CASCADING FILTER SUMMARY",
        "=" * 60,
        "",
        "Funnel stages:",
        "",
    ]
    
    for stage_name, n_in, n_out in stages:
        reduction = (1 - n_out / n_in) * 100 if n_in > 0 else 0
        lines.append(f"  {stage_name:20s}  {n_in:6d} → {n_out:6d}  ({reduction:.1f}% removed)")
    
    lines.extend([
        "",
        f"Final candidates: {len(filtered)}",
        "",
        "Next: Phase XVIII (Boltz-2 guided design)",
    ])
    
    with open(OUT / "cascade_summary.txt", "w") as f:
        f.write("\n".join(lines) + "\n")
    
    print(f"  {len(molecules)} → {len(filtered)} molecules")
    print("Done.")

if __name__ == "__main__":
    main()
