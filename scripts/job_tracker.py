#!/usr/bin/env python3
"""Job tracker for the lean pipeline.

Every computational job records:
  - molecule ID / site ID
  - receptor ID
  - stoichiometry
  - model (AF3/Boltz2/dock)
  - seed
  - software version
  - parameters
  - GPU/CPU
  - wall time
  - output path
  - success/failure
  - confidence metrics
"""

import csv, os, time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TRACKER_FILE = ROOT / "job_tracker.csv"

FIELDS = [
    "timestamp", "phase", "job_id", "site_id", "molecule_id", "receptor_id",
    "stoichiometry", "method", "seed", "sample", "software", "version",
    "parameters", "gpu_cpu", "wall_time_s", "output_path", "success",
    "confidence", "notes"
]

def log_job(phase, job_id, **kwargs):
    """Append a job record to the tracker."""
    row = {"timestamp": time.strftime("%Y-%m-%d %H:%M:%S"), "phase": phase, "job_id": job_id}
    for f in FIELDS:
        if f not in row:
            row[f] = kwargs.get(f, "")
    
    file_exists = TRACKER_FILE.exists()
    with open(TRACKER_FILE, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS)
        if not file_exists:
            writer.writeheader()
        writer.writerow(row)

def list_jobs(phase=None):
    """List all jobs, optionally filtered by phase."""
    if not TRACKER_FILE.exists():
        return []
    
    jobs = []
    with open(TRACKER_FILE) as f:
        reader = csv.DictReader(f)
        for row in reader:
            if phase is None or row["phase"] == phase:
                jobs.append(row)
    return jobs

def summary():
    """Print a summary of all jobs."""
    jobs = list_jobs()
    if not jobs:
        print("No jobs recorded yet.")
        return
    
    by_phase = {}
    for j in jobs:
        p = j["phase"]
        if p not in by_phase:
            by_phase[p] = {"total": 0, "success": 0, "failed": 0}
        by_phase[p]["total"] += 1
        if j["success"] == "True":
            by_phase[p]["success"] += 1
        else:
            by_phase[p]["failed"] += 1
    
    print(f"Total jobs: {len(jobs)}")
    print(f"{'Phase':15s} {'Total':>6s} {'OK':>6s} {'Fail':>6s}")
    print("-" * 40)
    for p, counts in sorted(by_phase.items()):
        print(f"{p:15s} {counts['total']:6d} {counts['success']:6d} {counts['failed']:6d}")

if __name__ == "__main__":
    summary()
