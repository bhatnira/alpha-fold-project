#!/usr/bin/env python3
"""
run_all.py - Master script that runs all steps of the mechanistic SAR
pipeline in order. Handles errors gracefully and prints progress.
"""
import subprocess
import sys
import time
import json
from pathlib import Path

PROJECT = Path("/cluster/home/nbhatt04/lean_pipeline")
MECHANISTIC = PROJECT / "mechanistic_sar"
SCRIPTS = MECHANISTIC / "scripts"
RESULTS = MECHANISTIC / "results"
STATUS_FILE = RESULTS / "pipeline_status.json"

# Ordered list of scripts to run
STEPS = [
    ("01_audit", "Data Audit"),
    ("02_prepare_data", "Prepare Data"),
    ("03_binding_features", "Binding Features"),
    ("04_structural_features", "Structural Features"),
    ("05_state_analysis", "State Analysis"),
    ("06_nma", "Normal Mode Analysis"),
    ("07_network", "Network Analysis"),
    ("08_consensus", "Consensus Analysis"),
    ("09_build_fingerprints", "Chemical Fingerprints"),
    ("10_train_models", "Model Training"),
    ("11_xai", "Explainability"),
    ("12_stability", "Stability Analysis"),
    ("13_ablation", "Ablation Study"),
    ("14_statistics", "Statistical Tests"),
    ("15_generate_figures", "Figure Generation"),
]


def run_step(script_name, step_name, results_dir):
    """Run a single pipeline step."""
    script_path = SCRIPTS / f"{script_name}.py"
    if not script_path.exists():
        print(f"  SKIPPED: {script_path} not found")
        return False, 0.0

    print(f"\n{'='*70}")
    print(f"Running: {step_name} ({script_name}.py)")
    print(f"{'='*70}")

    start_time = time.time()
    try:
        result = subprocess.run(
            [sys.executable, str(script_path)],
            cwd=str(MECHANISTIC),
            timeout=3600,  # 1 hour timeout per step
            capture_output=True,
            text=True,
        )
        elapsed = time.time() - start_time

        if result.returncode == 0:
            # Print last 20 lines of output
            lines = result.stdout.strip().split("\n")
            for line in lines[-20:]:
                print(f"  {line}")
            print(f"\n  COMPLETED in {elapsed:.1f}s")
            return True, elapsed
        else:
            print(f"  FAILED with return code {result.returncode}")
            if result.stdout:
                lines = result.stdout.strip().split("\n")
                for line in lines[-10:]:
                    print(f"  {line}")
            if result.stderr:
                lines = result.stderr.strip().split("\n")
                for line in lines[-10:]:
                    print(f"  ERROR: {line}")
            return False, elapsed

    except subprocess.TimeoutExpired:
        elapsed = time.time() - start_time
        print(f"  TIMEOUT after {elapsed:.1f}s")
        return False, elapsed
    except Exception as e:
        elapsed = time.time() - start_time
        print(f"  EXCEPTION: {e}")
        return False, elapsed


def main():
    print("=" * 70)
    print("MECHANISTIC SAR PIPELINE - Full Run")
    print("=" * 70)
    print(f"Project: {PROJECT}")
    print(f"Scripts: {SCRIPTS}")
    print(f"Results: {RESULTS}")
    print(f"Steps: {len(STEPS)}")

    RESULTS.mkdir(parents=True, exist_ok=True)

    # Parse command line arguments
    start_step = 0
    end_step = len(STEPS)
    skip_steps = set()

    if len(sys.argv) > 1:
        # Can specify: run_all.py [start] [end] or run_all.py --skip 05,06
        for arg in sys.argv[1:]:
            if arg.startswith("--skip"):
                idx = sys.argv.index(arg)
                if "=" in arg:
                    skip_str = arg.split("=")[1]
                elif idx + 1 < len(sys.argv):
                    skip_str = sys.argv[idx + 1]
                else:
                    skip_str = ""
                for s in skip_str.split(","):
                    s = s.strip()
                    if s:
                        skip_steps.add(s)
            elif arg.isdigit():
                if start_step == 0:
                    start_step = int(arg) - 1  # 1-indexed
                else:
                    end_step = int(arg)

    print(f"\nSteps to run: {start_step + 1} to {end_step}")
    if skip_steps:
        print(f"Skipping: {skip_steps}")

    # Run all steps
    total_start = time.time()
    status = {}
    completed = 0
    failed = 0

    for i, (script_name, step_name) in enumerate(STEPS):
        if i < start_step or i >= end_step:
            continue
        if script_name[:2] in skip_steps or str(i + 1) in skip_steps:
            print(f"\nSKIPPING: {step_name}")
            status[script_name] = {"status": "skipped", "time": 0}
            continue

        success, elapsed = run_step(script_name, step_name, RESULTS)
        status[script_name] = {
            "status": "completed" if success else "failed",
            "time": round(elapsed, 1),
        }

        if success:
            completed += 1
        else:
            failed += 1

    total_elapsed = time.time() - total_start

    # Save status
    pipeline_status = {
        "total_steps": len(STEPS),
        "completed": completed,
        "failed": failed,
        "skipped": len(STEPS) - completed - failed,
        "total_time_seconds": round(total_elapsed, 1),
        "steps": status,
    }
    with open(STATUS_FILE, "w") as f:
        json.dump(pipeline_status, f, indent=2)

    # Print summary
    print("\n" + "=" * 70)
    print("PIPELINE SUMMARY")
    print("=" * 70)
    print(f"  Completed: {completed}/{len(STEPS)}")
    print(f"  Failed: {failed}/{len(STEPS)}")
    print(f"  Total time: {total_elapsed:.1f}s ({total_elapsed / 60:.1f} min)")

    if failed > 0:
        print("\n  Failed steps:")
        for script_name, s in status.items():
            if s["status"] == "failed":
                print(f"    - {script_name}")

    print(f"\nStatus saved to: {STATUS_FILE}")
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
