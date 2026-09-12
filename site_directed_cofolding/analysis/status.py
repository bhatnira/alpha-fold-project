#!/usr/bin/env python3
"""Quick status check for site-directed cofolding study."""

from pathlib import Path

BASE = Path("/cluster/home/nbhatt04/lean_pipeline/site_directed_cofolding")
RESULTS = BASE / "results"

groups = {
    "state_models": ["2to3_open", "2to3_desensitized", "3to2_open", "3to2_desensitized"],
    "ach_states": ["2to3_ach_bound", "3to2_ach_bound"],
}

print("SITE-DIRECTED COFOLDING STATUS")
print("=" * 60)

total_done = 0
total_pending = 0

for group, expected in groups.items():
    group_dir = RESULTS / group
    done = 0
    for exp in expected:
        # Check for any confidence JSON in the result tree
        result_dir = group_dir / f"a9a10_{exp}"
        if list(result_dir.rglob("confidence_*.json")):
            done += 1
            status = "DONE"
        elif list(result_dir.rglob("*.pdb")):
            status = "RUNNING"
        elif result_dir.exists():
            status = "PENDING"
        else:
            status = "NOT_STARTED"
        print(f"  {group}/{exp:30s}  {status}")
    total_done += done
    total_pending += len(expected) - done

# Check binary/ternary
for complex_type in ["binary", "ternary", "full_panel"]:
    group_dir = RESULTS / complex_type
    if group_dir.exists():
        dirs = [d for d in group_dir.iterdir() if d.is_dir()]
        done = sum(1 for d in dirs if list(d.rglob("confidence_*.json")))
        running = sum(1 for d in dirs if not list(d.rglob("confidence_*.json")) and list(d.rglob("*.pdb")))
        pending = len(dirs) - done - running
        print(f"  {complex_type:20s}  {done:3d} DONE  {running:3d} RUNNING  {pending:3d} PENDING")
        total_done += done
        total_pending += pending

print(f"\nTotal: {total_done} done, {total_pending} pending")
