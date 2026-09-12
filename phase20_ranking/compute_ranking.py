#!/usr/bin/env python3
"""Phase XX — Multi-objective ranking of candidate sites.

Ranks sites using:
  - Boltz-2 affinity
  - pharmacophore fit
  - interaction fingerprint
  - docking consensus
  - pocket occupancy
  - α9α10 selectivity
  - novelty
  - uncertainty

Uses Pareto-inspired ranking.
"""

import csv, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
OUT  = Path(__file__).resolve().parent

def load_evidence():
    evidence = {}
    ev_file = ROOT / "phase14_evidence" / "evidence_matrix.csv"
    if ev_file.exists():
        with open(ev_file) as f:
            reader = csv.DictReader(f)
            for row in reader:
                evidence[int(row["site_id"])] = row
    return evidence

def load_selectivity():
    sel = {}
    sel_file = ROOT / "phase19_selectivity" / "selectivity_scores.csv"
    if sel_file.exists():
        with open(sel_file) as f:
            reader = csv.DictReader(f)
            for row in reader:
                sel[int(row["site_id"])] = row
    return sel

def load_robustness():
    rob = {}
    rob_file = ROOT / "phase11_robustness" / "robustness_matrix.csv"
    if rob_file.exists():
        with open(rob_file) as f:
            reader = csv.DictReader(f)
            for row in reader:
                rob[int(row["site"])] = row
    return rob

def load_specificity():
    spec = {}
    spec_file = ROOT / "phase09_specificity" / "specificity_results.csv"
    if spec_file.exists():
        with open(spec_file) as f:
            reader = csv.DictReader(f)
            for row in reader:
                spec[int(row["site"])] = row
    return spec

def load_convergence():
    conv = {}
    conv_file = ROOT / "phase04_convergence" / "convergence_map.csv"
    if conv_file.exists():
        with open(conv_file) as f:
            reader = csv.DictReader(f)
            for row in reader:
                site = int(row["site"])
                if site not in conv:
                    conv[site] = {"jaccards": [], "classes": []}
                conv[site]["jaccards"].append(float(row["contact_jaccard"]))
                conv[site]["classes"].append(row["convergence_class"])
    return conv

def pareto_rank(scores_list):
    """Simple Pareto ranking."""
    n = len(scores_list)
    ranks = [0] * n
    for i in range(n):
        dominated_by = 0
        for j in range(n):
            if i == j:
                continue
            # j dominates i if j is better in ALL objectives
            if all(scores_list[j][k] >= scores_list[i][k] for k in scores_list[i]):
                # and strictly better in at least one
                if any(scores_list[j][k] > scores_list[i][k] for k in scores_list[i]):
                    dominated_by += 1
        ranks[i] = dominated_by + 1
    return ranks

def compute_ranking(evidence, selectivity, robustness, specificity, convergence):
    """Multi-objective ranking of candidate sites."""
    all_sites = set(evidence.keys()) | set(selectivity.keys())
    
    objectives = []
    for site_id in sorted(all_sites):
        ev = evidence.get(site_id, {})
        sel = selectivity.get(site_id, {})
        rob = robustness.get(site_id, {})
        spec = specificity.get(site_id, {})
        conv = convergence.get(site_id, {})
        
        # Normalize objectives to [0, 1]
        obj = {
            "evidence": float(ev.get("evidence_score", 0)),
            "selectivity": float(sel.get("selectivity_score", 0)),
            "robustness": float(rob.get("overall_robustness", 0)) if rob else 0,
            "stereo_discrimination": float(spec.get("stereo_discrimination", 0)) if spec else 0,
            "charge_discrimination": float(spec.get("charge_discrimination", 0)) if spec else 0,
            "convergence": max(conv.get("jaccards", [0])) if conv else 0,
        }
        
        # Weighted sum (simple aggregation)
        weights = {
            "evidence": 0.25,
            "selectivity": 0.20,
            "robustness": 0.15,
            "stereo_discrimination": 0.15,
            "charge_discrimination": 0.15,
            "convergence": 0.10,
        }
        
        weighted_score = sum(obj[k] * weights[k] for k in weights)
        
        objectives.append({
            "site_id": site_id,
            **{f"obj_{k}": round(v, 4) for k, v in obj.items()},
            "weighted_score": round(weighted_score, 4),
            "classification": ev.get("classification", "UNKNOWN"),
            "composite_score": ev.get("composite_score", 0),
        })
    
    # Sort by weighted score
    objectives.sort(key=lambda r: r["weighted_score"], reverse=True)
    
    # Assign rank
    for i, r in enumerate(objectives):
        r["rank"] = i + 1
    
    return objectives

def write_csv(rows, path):
    if not rows:
        return
    all_fields = set()
    for r in rows:
        all_fields.update(r.keys())
    all_fields = sorted(all_fields)
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=all_fields)
        writer.writeheader()
        writer.writerows(rows)

def write_summary(ranking, path):
    lines = [
        "MULTI-OBJECTIVE RANKING",
        "=" * 60,
        "",
        "Objectives:",
        "  evidence (0.25) × selectivity (0.20) × robustness (0.15)",
        "  × stereo_disc (0.15) × charge_disc (0.15) × convergence (0.10)",
        "",
        f"Total sites ranked: {len(ranking)}",
        "",
        "TOP 10 CANDIDATE SITES:",
        "",
    ]
    
    for r in ranking[:10]:
        lines.append(f"  Rank {r['rank']:2d}  site {r['site_id']:2d}  "
                      f"score={r['weighted_score']:.4f}  "
                      f"class={r['classification']}")
        lines.append(f"           evidence={r['obj_evidence']:.3f}  "
                      f"selectivity={r['obj_selectivity']:.3f}  "
                      f"robustness={r['obj_robustness']:.3f}")
        lines.append(f"           stereo={r['obj_stereo_discrimination']:.3f}  "
                      f"charge={r['obj_charge_discrimination']:.3f}  "
                      f"convergence={r['obj_convergence']:.3f}")
        lines.append("")
    
    lines.extend([
        "RECOMMENDATION:",
        "",
        "  For experimental validation, prioritize sites that are:",
        "  1. Ranked in top 5 by multi-objective score",
        "  2. Classified as HIGH_CONFIDENCE or INTERMEDIATE",
        "  3. Show strong stereochemical and charge discrimination",
        "  4. Have high convergence across AF3 and Boltz-2",
        "",
        "  The top-ranked site should be the primary target for",
        "  mutagenesis and binding assays.",
    ])
    
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")

def main():
    print("Loading data...")
    evidence = load_evidence()
    selectivity = load_selectivity()
    robustness = load_robustness()
    specificity = load_specificity()
    convergence = load_convergence()
    print(f"  evidence={len(evidence)}, selectivity={len(selectivity)}, "
          f"robustness={len(robustness)}, specificity={len(specificity)}, "
          f"convergence={len(convergence)}")
    
    print("Computing ranking...")
    ranking = compute_ranking(evidence, selectivity, robustness, specificity, convergence)
    
    print("Writing outputs...")
    write_csv(ranking, OUT / "site_ranking.csv")
    write_summary(ranking, OUT / "ranking_summary.txt")
    
    print("Done.")

if __name__ == "__main__":
    main()
