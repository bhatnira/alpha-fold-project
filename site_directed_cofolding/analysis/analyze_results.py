#!/usr/bin/env python3
"""
Analyze site-directed cofolding results.

Reads Boltz-2 confidence, affinity, and structural outputs.
Builds convergence matrix across states, stoichiometries, and sites.
Cross-validates with experimental SAR.
"""

import json
import csv
import os
from pathlib import Path
from collections import defaultdict

PIPELINE = Path("/cluster/home/nbhatt04/lean_pipeline")
BASE = PIPELINE / "site_directed_cofolding"
RESULTS = BASE / "results"
OUTPUT = BASE / "analysis"

# Compound activity mapping
COMPOUND_ACTIVITY = {
    "L_ASCORBATE": {"activity": "active", "potency_uM": 1797},
    "CPD12_POTENT": {"activity": "active", "potency_uM": 0.198},
    "CPD25_STRONG": {"activity": "active", "potency_uM": 2.63},
    "CPD1_MODERATE": {"activity": "active", "potency_uM": 1797},
    "CPD18_MODERATE": {"activity": "active", "potency_uM": 1288},
    "CPD3_WEAK": {"activity": "active", "potency_uM": 6077},
    "CPD8_INACTIVE": {"activity": "inactive", "potency_uM": 0},
    "CPD9_INACTIVE": {"activity": "inactive", "potency_uM": 0},
}

CANDIDATE_SITES = [23, 5, 21, 34, 33]
STOICHIOMETRIES = ["2to3", "3to2"]


def parse_confidence(result_dir):
    """Extract confidence metrics from a Boltz-2 result directory."""
    metrics = {}
    for p in Path(result_dir).rglob("confidence_*.json"):
        try:
            with open(p) as f:
                data = json.load(f)
            metrics = {
                "confidence_score": data.get("confidence_score"),
                "ptm": data.get("ptm"),
                "iptm": data.get("iptm"),
                "ligand_iptm": data.get("ligand_iptm"),
                "complex_plddt": data.get("complex_plddt"),
            }
            break
        except Exception:
            continue
    return metrics


def parse_affinity(result_dir):
    """Extract affinity predictions from a Boltz-2 result directory."""
    for p in Path(result_dir).rglob("affinity_*.json"):
        try:
            with open(p) as f:
                data = json.load(f)
            return {
                "affinity_pred": data.get("affinity_pred_value"),
                "affinity_prob": data.get("affinity_probability_binary"),
            }
        except Exception:
            continue
    return {}


def scan_results():
    """Scan all result directories and extract metrics."""
    records = []
    for group in ["state_models", "ach_states", "binary", "ternary", "full_panel"]:
        group_dir = RESULTS / group
        if not group_dir.exists():
            continue
        for result_dir in group_dir.iterdir():
            if not result_dir.is_dir():
                continue
            name = result_dir.name
            # Try to find boltz_results subdirectory
            boltz_dirs = list(result_dir.glob("boltz_results_*"))
            if boltz_dirs:
                actual = boltz_dirs[0] / "predictions"
            else:
                actual = result_dir

            conf = parse_confidence(actual)
            aff = parse_affinity(actual)

            if conf or aff:
                record = {"name": name, "group": group}
                record.update(conf)
                record.update(aff)
                records.append(record)
    return records


def build_convergence_matrix(records):
    """Build convergence matrix: site x compound x stoich x complex_type."""
    matrix = defaultdict(list)
    for r in records:
        name = r["name"]
        # Parse name components
        parts = name.split("_")
        matrix[name].append(r)
    return matrix


def sar_cross_validation(records):
    """Cross-validate predictions against experimental SAR."""
    results_by_compound = defaultdict(list)
    for r in records:
        name = r["name"]
        for cpd, info in COMPOUND_ACTIVITY.items():
            if cpd in name:
                results_by_compound[cpd].append(r)
                break

    print("\n" + "=" * 70)
    print("SAR CROSS-VALIDATION")
    print("=" * 70)
    print(f"\n{'Compound':<20s} {'Activity':<10s} {'Potency':<10s} {'n_models':<10s} {'avg_conf':<10s} {'avg_iptm':<10s}")
    print("-" * 70)

    for cpd in ["CPD12_POTENT", "CPD25_STRONG", "CPD18_MODERATE",
                 "L_ASCORBATE", "CPD1_MODERATE", "CPD3_WEAK",
                 "CPD8_INACTIVE", "CPD9_INACTIVE"]:
        info = COMPOUND_ACTIVITY[cpd]
        recs = results_by_compound.get(cpd, [])
        if recs:
            confs = [r["confidence_score"] for r in recs if r.get("confidence_score")]
            iptms = [r["ligand_iptm"] for r in recs if r.get("ligand_iptm")]
            avg_conf = sum(confs) / len(confs) if confs else float("nan")
            avg_iptm = sum(iptms) / len(iptms) if iptms else float("nan")
            print(f"{cpd:<20s} {info['activity']:<10s} {info['potency_uM']:<10.1f} {len(recs):<10d} {avg_conf:<10.3f} {avg_iptm:<10.3f}")
        else:
            print(f"{cpd:<20s} {info['activity']:<10s} {info['potency_uM']:<10.1f} {'NO_DATA':<10s}")


def site_comparison(records):
    """Compare sites across all conditions."""
    site_metrics = defaultdict(lambda: {"confs": [], "iptms": [], "affs": []})
    for r in records:
        name = r["name"]
        for site in CANDIDATE_SITES:
            if f"site{site}" in name:
                if r.get("confidence_score"):
                    site_metrics[site]["confs"].append(r["confidence_score"])
                if r.get("ligand_iptm"):
                    site_metrics[site]["iptms"].append(r["ligand_iptm"])
                if r.get("affinity_pred"):
                    site_metrics[site]["affs"].append(r["affinity_pred"])
                break

    print("\n" + "=" * 70)
    print("SITE COMPARISON")
    print("=" * 70)
    print(f"\n{'Site':<8s} {'n':<6s} {'avg_conf':<10s} {'avg_iptm':<10s} {'avg_aff':<10s}")
    print("-" * 50)

    for site in CANDIDATE_SITES:
        m = site_metrics[site]
        n = max(len(m["confs"]), len(m["iptms"]), len(m["affs"]))
        avg_c = sum(m["confs"]) / len(m["confs"]) if m["confs"] else float("nan")
        avg_i = sum(m["iptms"]) / len(m["iptms"]) if m["iptms"] else float("nan")
        avg_a = sum(m["affs"]) / len(m["affs"]) if m["affs"] else float("nan")
        print(f"{site:<8d} {n:<6d} {avg_c:<10.3f} {avg_i:<10.3f} {avg_a:<10.3f}")


def stoichiometry_comparison(records):
    """Compare 2to3 vs 3to2."""
    stoich_metrics = defaultdict(lambda: {"confs": [], "iptms": []})
    for r in records:
        name = r["name"]
        for s in STOICHIOMETRIES:
            if s in name:
                if r.get("confidence_score"):
                    stoich_metrics[s]["confs"].append(r["confidence_score"])
                if r.get("ligand_iptm"):
                    stoich_metrics[s]["iptms"].append(r["ligand_iptm"])
                break

    print("\n" + "=" * 70)
    print("STOICHIOMETRY COMPARISON")
    print("=" * 70)
    print(f"\n{'Stoich':<10s} {'n':<6s} {'avg_conf':<10s} {'avg_iptm':<10s}")
    print("-" * 40)

    for s in STOICHIOMETRIES:
        m = stoich_metrics[s]
        avg_c = sum(m["confs"]) / len(m["confs"]) if m["confs"] else float("nan")
        avg_i = sum(m["iptms"]) / len(m["iptms"]) if m["iptms"] else float("nan")
        print(f"{s:<10s} {len(m['confs']):<6d} {avg_c:<10.3f} {avg_i:<10.3f}")


def write_results_csv(records):
    """Write all results to CSV."""
    out_path = OUTPUT / "cofolding_results.csv"
    if not records:
        print("No results to write")
        return
    fieldnames = list(records[0].keys())
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(records)
    print(f"\nResults CSV: {out_path}")


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("Site-Directed Cofolding Analysis")
    print("=" * 70)

    records = scan_results()
    print(f"\nFound {len(records)} prediction results")

    if not records:
        print("No results found. Jobs may still be running.")
        print("Re-run this script after predictions complete.")
        return

    write_results_csv(records)
    sar_cross_validation(records)
    site_comparison(records)
    stoichiometry_comparison(records)

    print("\n" + "=" * 70)
    print("Analysis complete")
    print("=" * 70)


if __name__ == "__main__":
    main()
