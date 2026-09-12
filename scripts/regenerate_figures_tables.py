#!/usr/bin/env python3
"""
Regenerate affected figures and tables with new gap-fill data.
Updates:
- Figure 3: State ensemble (add open/desensitized when available)
- Figure 13: Stoichiometry/state/ACh dependence
- Table 2: Receptor ensemble
- Table 7: Stoichiometry analysis
- Table 9: Falsification (add ensemble docking result)
- Table 10: Final ranking (add ensemble scores)
"""
import json, csv, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
FIGDIR = ROOT / "figures"
TABLEDIR = ROOT / "tables"
FIGDIR.mkdir(parents=True, exist_ok=True)
TABLEDIR.mkdir(parents=True, exist_ok=True)


def update_table10_final_ranking():
    """Update Table 10 with ensemble docking scores."""
    ensemble_csv = ROOT / "09_docking" / "flexible_ensemble" / "ensemble_docking_results.csv"
    if not ensemble_csv.exists():
        return

    ensemble = {}
    with open(ensemble_csv) as f:
        for row in csv.DictReader(f):
            key = (row["compound_id"], row["site_id"])
            ensemble[key] = float(row["mean_energy"])

    convergence_csv = ROOT / "14_stoichiometry_state" / "full_convergence_matrix.csv"
    if convergence_csv.exists():
        updated_csv = TABLEDIR / "Table10_final_ranking_v2.csv"
        with open(convergence_csv) as f:
            reader = csv.DictReader(f)
            fieldnames = list(reader.fieldnames) + ["ensemble_docking_score"]
            rows = []
            for row in reader:
                site = row.get("site_id", "")
                best_energy = None
                for (cpd, s), energy in ensemble.items():
                    if str(s) == str(site):
                        if best_energy is None or energy < best_energy:
                            best_energy = energy
                row["ensemble_docking_score"] = best_energy if best_energy else ""
                rows.append(row)

        with open(updated_csv, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for row in rows:
                writer.writerow(row)
        print(f"  Updated Table 10: {updated_csv}")


def update_table9_falsification():
    """Add ensemble docking falsification test to Table 9."""
    report = {
        "test": "F9_ensemble_docking_discrimination",
        "description": "Active vs inactive discrimination across ensemble docking",
        "result": "FAIL",
        "detail": "No site achieves active/inactive separation in ensemble docking",
        "active_mean_energy": -7.86,
        "inactive_mean_energy": -8.30,
        "note": "Inactive compounds consistently score better than active compounds",
    }

    table9_path = TABLEDIR / "Table9_falsification_report.txt"
    with open(table9_path, "a") as f:
        f.write(f"\n{report['test']}: {report['result']}\n")
        f.write(f"  {report['description']}\n")
        f.write(f"  {report['detail']}\n")
    print(f"  Updated Table 9: {table9_path}")


def generate_figure3_state_ensemble():
    """Generate updated Figure 3: State Ensemble status."""
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))

    states = ["Resting (APO)", "Open", "Desensitized", "ACh-bound"]
    status = ["AVAILABLE", "PENDING", "PENDING", "PENDING"]
    colors = ["#2ecc71", "#e74c3c", "#e74c3c", "#e74c3c"]

    bars = ax.barh(states, [1]*len(states), color=colors, alpha=0.7)
    for i, (bar, s) in enumerate(zip(bars, status)):
        ax.text(0.5, i, s, ha='center', va='center', fontweight='bold', fontsize=12,
                color='white' if s == 'AVAILABLE' else 'black')

    ax.set_xlim(0, 1)
    ax.set_xticks([])
    ax.set_title("Figure 3: Receptor State Ensemble Status", fontsize=14, fontweight='bold')
    ax.set_xlabel("Status")
    ax.invert_yaxis()

    plt.tight_layout()
    out_path = FIGDIR / "fig3_state_ensemble_status.png"
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Generated: {out_path}")


def generate_figure13_stoich_state_ach():
    """Generate updated Figure 13: Stoichiometry/State/ACh analysis."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    axes[0].bar(["2to3", "3to2"], [23, 21], color=["#3498db", "#e67e22"])
    axes[0].set_title("Sites per Stoichiometry")
    axes[0].set_ylabel("Number of Sites")

    axes[1].bar(["Resting", "Open", "Desensitized"], [46, 0, 0], color=["#2ecc71", "#e74c3c", "#e74c3c"])
    axes[1].set_title("State Coverage")
    axes[1].set_ylabel("Models Generated")

    axes[2].bar(["R", "R+PAM", "R+ACh", "R+ACh+PAM"], [480, 480, 0, 0],
                color=["#2ecc71", "#3498db", "#e74c3c", "#e74c3c"])
    axes[2].set_title("ACh Occupancy Coverage")
    axes[2].set_ylabel("Docking Results")

    plt.suptitle("Figure 13: Stoichiometry × State × ACh Analysis Status",
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    out_path = FIGDIR / "fig13_stoich_state_ach_status.png"
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Generated: {out_path}")


def generate_ensemble_docking_figure():
    """Generate ensemble docking discrimination figure."""
    ensemble_csv = ROOT / "09_docking" / "flexible_ensemble" / "ensemble_docking_results.csv"
    if not ensemble_csv.exists():
        return

    sites = defaultdict(lambda: {"active": [], "inactive": []})
    with open(ensemble_csv) as f:
        for row in csv.DictReader(f):
            site = int(float(row["site_id"]))
            energy = float(row["mean_energy"])
            if row["is_active"] == "True":
                sites[site]["active"].append(energy)
            else:
                sites[site]["inactive"].append(energy)

    fig, ax = plt.subplots(figsize=(10, 6))
    site_ids = sorted(sites.keys())
    x = np.arange(len(site_ids))
    width = 0.35

    active_means = [np.mean(sites[s]["active"]) if sites[s]["active"] else 0 for s in site_ids]
    inactive_means = [np.mean(sites[s]["inactive"]) if sites[s]["inactive"] else 0 for s in site_ids]

    ax.bar(x - width/2, active_means, width, label='Active', color='#2ecc71', alpha=0.7)
    ax.bar(x + width/2, inactive_means, width, label='Inactive', color='#e74c3c', alpha=0.7)

    ax.set_xlabel('Site ID')
    ax.set_ylabel('Mean Binding Energy (kcal/mol)')
    ax.set_title('Figure 7b: Ensemble Docking - Active vs Inactive by Site')
    ax.set_xticks(x)
    ax.set_xticklabels([f"Site {s}" for s in site_ids])
    ax.legend()
    ax.invert_yaxis()

    plt.tight_layout()
    out_path = FIGDIR / "fig7b_ensemble_docking.png"
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Generated: {out_path}")


if __name__ == "__main__":
    print("Regenerating affected figures and tables")
    print("=" * 60)

    from collections import defaultdict
    update_table10_final_ranking()
    update_table9_falsification()
    generate_figure3_state_ensemble()
    generate_figure13_stoich_state_ach()
    generate_ensemble_docking_figure()

    print("\nDone!")
