#!/usr/bin/env python3
"""Phase VIII — Boltz-2 affinity analysis.

Constructs ligand × site × receptor-state affinity matrix.
Compares predicted ranking with experimental activity.
"""

import csv, sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
OUT  = Path(__file__).resolve().parent

def load_site_ranking():
    sites = {}
    with open(DATA / "FINAL_SITE_RANKING.csv") as f:
        reader = csv.DictReader(f)
        for row in reader:
            sites[int(row["site_id"])] = {
                "rank": int(row["rank_primary"]),
                "composite_score": float(row["SiteScore"]),
                "n_boltz2": int(row["n_boltz2"]),
                "max_enrichment": float(row["max_ligand_enrichment"]),
            }
    return sites

def load_ligand_comparison():
    rows = []
    with open(DATA / "LIGAND_COMPARISON.csv") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append({
                "site": int(row["site"]),
                "ligand": row["ligand_class"],
                "n_models_boltz2": int(row["n_models_boltz2"]),
                "enrichment": float(row["enrichment"]),
                "mean_confidence": float(row["mean_confidence"]),
                "stoichiometry": row["stoichiometry"],
            })
    return rows

def build_affinity_matrix(sites, lig_comp):
    """Build ligand × site affinity matrix from Boltz-2 data."""
    matrix = []
    
    for lc in lig_comp:
        site_id = lc["site"]
        if site_id not in sites:
            continue
        
        # Use enrichment as proxy for relative affinity
        # (actual Boltz-2 affinity scores would be computed separately)
        affinity_score = lc["enrichment"] * lc["mean_confidence"]
        
        matrix.append({
            "site_id": site_id,
            "site_rank": sites[site_id]["rank"],
            "ligand": lc["ligand"],
            "stoichiometry": lc["stoichiometry"],
            "n_boltz2_models": lc["n_models_boltz2"],
            "enrichment": lc["enrichment"],
            "confidence": lc["mean_confidence"],
            "affinity_score": round(affinity_score, 4),
        })
    
    return matrix

def rank_by_affinity(matrix):
    """Rank ligand × site pairs by predicted affinity."""
    return sorted(matrix, key=lambda r: r["affinity_score"], reverse=True)

def compare_with_activity(matrix, sar_rows):
    """Compare predicted ranking with experimental activity tiers."""
    # Activity tiers from SAR
    activity_tiers = {
        "L-ASC": 1,
        "O-ETHYL": 2,
        "D-ASC": 3,
        "RYANODINE": 1,
        "ACETATE": 0,
    }
    
    # For each site, check if predicted affinity preserves activity order
    by_site = defaultdict(list)
    for r in matrix:
        by_site[r["site_id"]].append(r)
    
    comparisons = []
    for site_id, records in sorted(by_site.items()):
        # Sort by affinity
        records.sort(key=lambda r: r["affinity_score"], reverse=True)
        
        # Check if L-ASC > D-ASC
        l_asc = next((r for r in records if r["ligand"] == "L-ASC"), None)
        d_asc = next((r for r in records if r["ligand"] == "D-ASC"), None)
        
        preserves_stereo = True
        if l_asc and d_asc:
            if l_asc["affinity_score"] < d_asc["affinity_score"]:
                preserves_stereo = False
        
        # Check if L-ASC > ACETATE
        acetate = next((r for r in records if r["ligand"] == "ACETATE"), None)
        preserves_specificity = True
        if l_asc and acetate:
            if l_asc["affinity_score"] < acetate["affinity_score"]:
                preserves_specificity = False
        
        comparisons.append({
            "site_id": site_id,
            "n_ligands": len(records),
            "top_ligand": records[0]["ligand"] if records else "",
            "top_affinity": records[0]["affinity_score"] if records else 0,
            "preserves_stereo": preserves_stereo,
            "preserves_specificity": preserves_specificity,
            "ranking": ";".join(f"{r['ligand']}={r['affinity_score']:.3f}" for r in records[:5]),
        })
    
    return comparisons

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

def write_summary(matrix, comparisons, path):
    lines = [
        "BOLTZ-2 AFFINITY ANALYSIS",
        "=" * 60,
        "",
        "Method: enrichment × confidence as proxy for relative affinity",
        "Note: Actual Boltz-2 affinity scores require separate computation",
        "",
        f"Total ligand × site pairs: {len(matrix)}",
        "",
    ]
    
    # Top affinity pairs
    ranked = rank_by_affinity(matrix)
    lines.append("Top 10 predicted affinity pairs:")
    for r in ranked[:10]:
        lines.append(f"  site {r['site_id']:2d}  {r['ligand']:12s}  "
                      f"affinity={r['affinity_score']:.4f}  "
                      f"enrichment={r['enrichment']:.2f}")
    lines.append("")
    
    # SAR preservation
    stereo_ok = sum(1 for c in comparisons if c["preserves_stereo"])
    spec_ok = sum(1 for c in comparisons if c["preserves_specificity"])
    lines.append(f"SAR preservation:")
    lines.append(f"  Stereochemical order (L>D):  {stereo_ok}/{len(comparisons)} sites")
    lines.append(f"  Specificity (L>acetate):     {spec_ok}/{len(comparisons)} sites")
    
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")

def main():
    print("Loading data...")
    sites = load_site_ranking()
    lig_comp = load_ligand_comparison()
    print(f"  {len(sites)} sites, {len(lig_comp)} ligand×site records")
    
    print("Building affinity matrix...")
    matrix = build_affinity_matrix(sites, lig_comp)
    print(f"  {len(matrix)} records")
    
    print("Comparing with experimental activity...")
    comparisons = compare_with_activity(matrix, [])
    
    print("Writing outputs...")
    write_csv(matrix, OUT / "affinity_matrix.csv")
    write_csv(comparisons, OUT / "affinity_comparisons.csv")
    write_summary(matrix, comparisons, OUT / "affinity_summary.txt")
    
    print("Done.")

if __name__ == "__main__":
    main()
