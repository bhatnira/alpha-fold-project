#!/usr/bin/env python3
"""
Phase 12: Flexible Docking via Ensemble Approach

Since Vina 1.2.7 SWIG does not support flex receptor, use ensemble docking:
dock each compound against multiple receptor conformers and score by consensus.

Approach:
1. Use existing receptor models (2to3, 3to2 APO + open state if available)
2. For each compound, dock against all receptor conformers
3. Score by consensus across conformers
4. Compare active vs inactive discrimination
"""
import json, csv, os, subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUTDIR = ROOT / "09_docking" / "flexible_ensemble"
OUTDIR.mkdir(parents=True, exist_ok=True)

RECEPTORS = {
    "2to3_apo": ROOT / "09_docking" / "receptor_2to3.pdbqt",
    "3to2_apo": ROOT / "09_docking" / "receptor_3to2.pdbqt",
}

# Candidate sites with box centers (from existing docking results)
SITES = {
    "site23": {"center": [0, 0, 0], "description": "alpha9/alpha10 ECD interface"},
    "site25": {"center": [0, 0, 0], "description": "alpha9/alpha9 ECD interface"},
    "site34": {"center": [0, 0, 0], "description": "M1-M2/M2 TM region"},
}

LIGANDS = {}
ligands_dir = ROOT / "09_docking"
for i in range(1, 31):
    pdbqt = ligands_dir / f"ligand_{i}.pdbqt"
    if pdbqt.exists():
        LIGANDS[f"ligand_{i}"] = pdbqt


def get_site_boxes():
    """Extract docking box coordinates from existing results."""
    csv_path = ROOT / "data" / "FINAL_SITE_RANKING.csv"
    if csv_path.exists():
        with open(csv_path) as f:
            for row in csv.DictReader(f):
                site_id = row.get("site_id", "")
                if site_id in SITES:
                    try:
                        x = float(row.get("center_x", 0))
                        y = float(row.get("center_y", 0))
                        z = float(row.get("center_z", 0))
                        SITES[site_id]["center"] = [x, y, z]
                    except (ValueError, KeyError):
                        pass


def run_ensemble_docking():
    """Run docking across all receptor conformers for each site."""
    results = []

    for site_name, site_info in SITES.items():
        for lig_name, lig_pdbqt in LIGANDS.items():
            best_score = None
            best_receptor = None
            all_scores = {}

            for rec_name, rec_pdbqt in RECEPTORS.items():
                if not rec_pdbqt.exists():
                    continue

                center = site_info["center"]
                out_pdbqt = OUTDIR / f"{site_name}_{lig_name}_{rec_name}_docked.pdbqt"

                cmd = [
                    "vina",
                    "--receptor", str(rec_pdbqt),
                    "--ligand", str(lig_pdbqt),
                    "--center_x", str(center[0]),
                    "--center_y", str(center[1]),
                    "--center_z", str(center[2]),
                    "--size_x", "20",
                    "--size_y", "20",
                    "--size_z", "20",
                    "--exhaustiveness", "16",
                    "--num_poses", "5",
                    "--out", str(out_pdbqt),
                ]

                try:
                    result = subprocess.run(
                        cmd, capture_output=True, text=True, timeout=300
                    )
                    score = parse_vina_output(result.stdout)
                    if score is not None:
                        all_scores[rec_name] = score
                        if best_score is None or score < best_score:
                            best_score = score
                            best_receptor = rec_name
                except Exception as e:
                    pass

            consensus_score = (
                sum(all_scores.values()) / len(all_scores) if all_scores else None
            )
            results.append({
                "site": site_name,
                "ligand": lig_name,
                "best_score": best_score,
                "best_receptor": best_receptor,
                "consensus_score": consensus_score,
                "n_conformers": len(all_scores),
                "all_scores": all_scores,
            })

    return results


def parse_vina_output(stdout):
    """Extract best binding energy from Vina output."""
    for line in stdout.split("\n"):
        line = line.strip()
        if line and not line.startswith("-----") and not line.startswith("Writing"):
            parts = line.split()
            if len(parts) >= 4:
                try:
                    return float(parts[1])
                except ValueError:
                    continue
    return None


def main():
    print("Phase 12: Ensemble Flexible Docking")
    print("=" * 60)

    get_site_boxes()

    available_receptors = [k for k, v in RECEPTORS.items() if v.exists()]
    print(f"Receptor conformers: {available_receptors}")
    print(f"Ligands: {len(LIGANDS)}")
    print(f"Sites: {list(SITES.keys())}")

    results = run_ensemble_docking()

    out_csv = OUTDIR / "ensemble_docking_results.csv"
    with open(out_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "site", "ligand", "best_score", "best_receptor",
            "consensus_score", "n_conformers"
        ])
        writer.writeheader()
        for r in results:
            writer.writerow({k: r[k] for k in writer.fieldnames})

    print(f"\nResults: {out_csv}")
    print(f"Total dockings: {len(results)}")

    report = {
        "method": "ensemble_docking",
        "description": "Docking across multiple receptor conformers for robustness",
        "receptor_conformers": available_receptors,
        "n_ligands": len(LIGANDS),
        "n_sites": len(SITES),
        "n_total": len(results),
        "results_file": str(out_csv),
    }
    report_path = OUTDIR / "ensemble_docking_report.json"
    report_path.write_text(json.dumps(report, indent=2))
    print(f"Report: {report_path}")


if __name__ == "__main__":
    main()
