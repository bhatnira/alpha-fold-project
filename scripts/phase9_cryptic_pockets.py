#!/usr/bin/env python3
"""
Phase 9: Cryptic Pocket Analysis

Compare pocket volumes and geometry across receptor states:
- Apo resting (existing fpocket results)
- ACh-bound (from Phase 6 results when available)
- Open-like (from Phase 3-4 results when available)
- Desensitized-like (from Phase 3-4 results when available)

Classifies each candidate site as:
- constitutive (present in all states)
- state-dependent (present in specific states)
- ACh-dependent (created by ACh binding)
- cryptic (only in specific conformations)
"""
import json, csv, os
from pathlib import Path
from collections import defaultdict

ROOT = Path(__file__).resolve().parent.parent
OUTDIR = ROOT / "09_cryptic_pocket"
OUTDIR.mkdir(parents=True, exist_ok=True)

FPOCKET_DIR = ROOT / "phase03_pockets"
CANDIDATE_SITES = {
    "site23": {"interface": "alpha9(+)/alpha10(-) ECD", "residues": ["alpha9:176", "alpha9:175", "alpha10:143", "alpha10:145", "alpha9:120"]},
    "site25": {"interface": "alpha9(+)/alpha9(-) ECD", "residues": ["alpha9:176", "alpha9:175", "alpha9:120"]},
    "site34": {"interface": "M1-M2/M2 TM", "residues": ["alpha9:TMD", "alpha10:TMD"]},
    "site21": {"interface": "alpha10 ECD", "residues": ["alpha10:81", "alpha10:109"]},
}


def parse_fpocket_results(pocket_dir):
    """Parse fpocket output to get pocket volumes and residue composition."""
    pockets = []
    pocket_file = pocket_dir / "info.txt"
    if not pocket_file.exists():
        return pockets

    current_pocket = None
    with open(pocket_file) as f:
        for line in f:
            line = line.strip()
            if "Pocket" in line and ":" in line:
                if current_pocket:
                    pockets.append(current_pocket)
                parts = line.split(":")
                current_pocket = {"id": parts[0].strip(), "score": 0, "volume": 0, "residues": []}
            elif current_pocket and "Score" in line:
                try:
                    current_pocket["score"] = float(line.split(":")[-1].strip())
                except ValueError:
                    pass
            elif current_pocket and "Volume" in line:
                try:
                    current_pocket["volume"] = float(line.split(":")[-1].strip())
                except ValueError:
                    pass
            elif current_pocket and "RESD" in line:
                try:
                    resname = line.split()[2] if len(line.split()) > 2 else ""
                    current_pocket["residues"].append(resname)
                except (IndexError, ValueError):
                    pass
    if current_pocket:
        pockets.append(current_pocket)

    return pockets


def find_fpocket_dirs():
    """Find all fpocket result directories."""
    dirs = {}
    for d in FPOCKET_DIR.iterdir():
        if d.is_dir() and "out" in d.name:
            state = "apo"
            if "open" in d.name.lower():
                state = "open"
            elif "desens" in d.name.lower():
                state = "desensitized"
            elif "ach" in d.name.lower():
                state = "ach_bound"
            dirs[state] = d
    return dirs


def classify_pocket_constitutivity(state_pockets):
    """Classify whether a pocket is constitutive, state-dependent, or cryptic."""
    present_states = set(state_pockets.keys())
    all_states = {"apo", "open", "desensitized", "ach_bound"}

    if present_states == all_states:
        return "CONSTITUTIVE"
    elif present_states.issubset({"apo"}):
        return "APO_ONLY"
    elif "ach_bound" in present_states and "apo" not in present_states:
        return "ACH_DEPENDENT"
    elif present_states.issubset({"open", "desensitized"}):
        return "STATE_DEPENDENT"
    elif len(present_states) > 1:
        return "PARTIALLY_CONSTITUTIVE"
    else:
        return "UNCERTAIN"


def analyze_candidate_sites(fpocket_dirs, state_data):
    """For each candidate site, check if it appears in fpocket results."""
    analysis = {}
    for site_name, site_info in CANDIDATE_SITES.items():
        site_analysis = {
            "site": site_name,
            "interface": site_info["interface"],
            "residues": site_info["residues"],
            "state_pockets": {},
            "classification": "UNCERTAIN",
        }

        for state, fpocket_dir in fpocket_dirs.items():
            pockets = parse_fpocket_results(fpocket_dir)
            site_residues = set(r.split(":")[0] for r in site_info["residues"] if ":" in r)

            for pocket in pockets:
                pocket_residues = set(r[:3] for r in pocket["residues"])
                overlap = site_residues & pocket_residues
                if overlap or len(pocket["residues"]) > 0:
                    site_analysis["state_pockets"][state] = {
                        "fpocket_id": pocket["id"],
                        "volume": pocket["volume"],
                        "score": pocket["score"],
                        "n_residues": len(pocket["residues"]),
                    }
                    break

        site_analysis["classification"] = classify_pocket_constitutivity(
            {k: v for k, v in site_analysis["state_pockets"].items()}
        )
        analysis[site_name] = site_analysis

    return analysis


def main():
    print("Phase 9: Cryptic Pocket Analysis")
    print("=" * 60)

    fpocket_dirs = find_fpocket_dirs()
    print(f"Available fpocket states: {list(fpocket_dirs.keys())}")

    state_data = {}
    for state, fpocket_dir in fpocket_dirs.items():
        pockets = parse_fpocket_results(fpocket_dir)
        state_data[state] = {"n_pockets": len(pockets), "pockets": pockets}
        print(f"  {state}: {len(pockets)} pockets detected")

    analysis = analyze_candidate_sites(fpocket_dirs, state_data)

    print("\nCandidate Site Classification:")
    for site_name, info in analysis.items():
        print(f"  {site_name}: {info['classification']} (states: {list(info['state_pockets'].keys())})")

    report = {
        "analysis_type": "cryptic_pocket",
        "states_analyzed": list(fpocket_dirs.keys()),
        "n_candidate_sites": len(analysis),
        "results": analysis,
        "classification_summary": {
            k: sum(1 for v in analysis.values() if v["classification"] == k)
            for k in set(v["classification"] for v in analysis.values())
        },
        "note": "Currently only APO state has fpocket results. Open/desensitized/ACh-bound states need Phase 3-4 and Phase 6 results.",
    }

    report_path = OUTDIR / "cryptic_pocket_report.json"
    report_path.write_text(json.dumps(report, indent=2))
    print(f"\nReport: {report_path}")

    csv_path = OUTDIR / "cryptic_pocket_analysis.csv"
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "site", "interface", "classification", "states_present", "n_states"
        ])
        writer.writeheader()
        for site_name, info in analysis.items():
            writer.writerow({
                "site": site_name,
                "interface": info["interface"],
                "classification": info["classification"],
                "states_present": ";".join(info["state_pockets"].keys()),
                "n_states": len(info["state_pockets"]),
            })
    print(f"CSV: {csv_path}")


if __name__ == "__main__":
    main()
