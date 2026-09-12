#!/usr/bin/env python3
"""Quick status check for the lean pipeline."""

import os
from pathlib import Path

ROOT = Path(__file__).resolve().parent

def check_phase(phase_dir, name):
    """Check if a phase has outputs."""
    if not phase_dir.exists():
        return f"  {name:35s}  MISSING"
    
    files = list(phase_dir.glob("*"))
    csvs = [f for f in files if f.suffix == ".csv"]
    txts = [f for f in files if f.suffix == ".txt"]
    pmls = [f for f in files if f.suffix == ".pml"]
    other = [f for f in files if f.suffix not in (".csv", ".txt", ".pml", ".py")]
    
    total_size = sum(f.stat().st_size for f in files if f.is_file())
    
    parts = []
    if csvs:
        parts.append(f"{len(csvs)} CSV")
    if txts:
        parts.append(f"{len(txts)} TXT")
    if pmls:
        parts.append(f"{len(pmls)} PML")
    if other:
        parts.append(f"{len(other)} other")
    
    desc = ", ".join(parts) if parts else "empty"
    return f"  {name:35s}  {total_size/1024:7.1f} KB  {desc}"

def main():
    print("LEAN PIPELINE STATUS")
    print("=" * 70)
    print()
    
    phases = [
        ("phase01_sar", "Phase I — SAR dataset"),
        ("phase02_receptor", "Phase II — Receptor ensemble"),
        ("phase03_pockets", "Phase III — Blind pocket discovery"),
        ("phase04_convergence", "Phase IV — AF3 + Boltz-2 convergence"),
        ("phase05_docking", "Phase V — Ensemble docking"),
        ("phase06_fingerprints", "Phase VI — Interaction fingerprints"),
        ("phase07_sar_model", "Phase VII — Structural SAR"),
        ("phase08_boltz_affinity", "Phase VIII — Boltz-2 affinity"),
        ("phase09_specificity", "Phase IX — Specificity/falsification"),
        ("phase10_electrostatics", "Phase X — Electrostatic analysis"),
        ("phase11_robustness", "Phase XI — Robustness without MD"),
        ("phase12_subtype", "Phase XII — α9α10 vs α7 vs α4β2"),
        ("phase13_ryanodine", "Phase XIII — Ryanodine"),
        ("phase14_evidence", "Phase XIV — Integrated evidence"),
        ("phase15_pharmacophore", "Phase XV — Pharmacophore"),
        ("phase16_molgen", "Phase XVI — Molecule generation"),
        ("phase17_cascade", "Phase XVII — Cascading filter"),
        ("phase18_boltz_design", "Phase XVIII — Boltz-2 guided design"),
        ("phase19_selectivity", "Phase XIX — α9α10 selectivity"),
        ("phase20_ranking", "Phase XX — Multi-objective ranking"),
        ("phase21_validation", "Phase XXI — Prospective validation"),
        ("phase22_closed_loop", "Phase XXII — Closed-loop"),
    ]
    
    completed = 0
    for dirname, name in phases:
        phase_dir = ROOT / dirname
        print(check_phase(phase_dir, name))
        if any(phase_dir.glob("*.csv")):
            completed += 1
    
    print()
    print(f"Completed: {completed} / {len(phases)} phases")
    print()
    
    # Data sources
    print("DATA SOURCES:")
    data_dir = ROOT / "data"
    for item in sorted(data_dir.iterdir()):
        if item.is_symlink():
            target = os.readlink(item)
            exists = item.exists()
            print(f"  {item.name:35s}  {'OK' if exists else 'BROKEN'}  → {target}")
        elif item.is_file():
            print(f"  {item.name:35s}  {item.stat().st_size/1024:.1f} KB")
    
    print()
    
    # Key outputs
    print("KEY OUTPUTS:")
    for phase_dir in ["phase01_sar", "phase04_convergence", "phase06_fingerprints",
                      "phase09_specificity", "phase11_robustness", "phase14_evidence",
                      "phase15_pharmacophore"]:
        p = ROOT / phase_dir
        if p.exists():
            csvs = sorted(p.glob("*.csv"))
            for f in csvs[:2]:
                print(f"  {phase_dir}/{f.name}")

if __name__ == "__main__":
    main()
