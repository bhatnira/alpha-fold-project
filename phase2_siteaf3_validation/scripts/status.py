#!/usr/bin/env python3
"""Status check for Phase 2 SiteAF3 validation."""

from pathlib import Path

PIPELINE = Path("/cluster/home/nbhatt04/lean_pipeline")
PHASE2 = PIPELINE / "phase2_siteaf3_validation"

print("PHASE 2 STATUS: SiteAF3 Validation")
print("=" * 60)

# Check current jobs
print("\nCURRENTLY RUNNING:")
import subprocess
try:
    result = subprocess.run(["squeue", "-u", "nbhatt04", "--Format=jobid,name,timeleft,state"],
                          capture_output=True, text=True, timeout=5)
    lines = result.stdout.strip().split("\n")
    for line in lines[1:]:  # skip header
        parts = line.split()
        if len(parts) >= 4:
            print(f"  {parts[0]:<10s} {parts[1]:<20s} {parts[2]:<12s} {parts[3]}")
except Exception:
    print("  (cannot check queue)")

# Check inputs
print("\nGENERATED INPUTS:")
n_hotspot = len(list((PHASE2 / "hotspot_pocket_pdbs").glob("*.pdb")))
n_json = len(list((PHASE2 / "siteaf3_inputs").glob("*.json")))
print(f"  Hotspot/pocket PDBs: {n_hotspot}")
print(f"  SiteAF3 JSON configs: {n_json}")

# Check scripts
print("\nSCRIPTS:")
scripts = list((PHASE2 / "scripts").glob("*"))
for s in sorted(scripts):
    print(f"  {s.name}")

# Check results
print("\nRESULTS:")
for subdir in ["cofolding_analysis", "siteaf3"]:
    d = PHASE2 / "results" / subdir
    if d.exists():
        files = list(d.iterdir())
        print(f"  {subdir}: {len(files)} files")
    else:
        print(f"  {subdir}: NOT YET")

print("\n" + "=" * 60)
print("TO START: sbatch phase2_siteaf3_validation/scripts/master_phase2.sh")
