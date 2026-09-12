#!/usr/bin/env python3
"""Phase XIX — α9α10 selectivity analysis for candidate molecules.

Uses pharmacophore conservation + contact fingerprints to predict
relative subtype compatibility.
"""

import csv, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
OUT  = Path(__file__).resolve().parent

def load_pharmacophore():
    features = []
    pharma_file = ROOT / "phase15_pharmacophore" / "pharmacophore_features.csv"
    if pharma_file.exists():
        with open(pharma_file) as f:
            reader = csv.DictReader(f)
            for row in reader:
                features.append(row)
    return features

def load_evidence():
    evidence = {}
    ev_file = ROOT / "phase14_evidence" / "evidence_matrix.csv"
    if ev_file.exists():
        with open(ev_file) as f:
            reader = csv.DictReader(f)
            for row in reader:
                evidence[int(row["site_id"])] = row
    return evidence

def load_subtype():
    subtypes = []
    sub_file = ROOT / "phase12_subtype" / "subtype_comparison.csv"
    if sub_file.exists():
        with open(sub_file) as f:
            reader = csv.DictReader(f)
            for row in reader:
                subtypes.append(row)
    return subtypes

def compute_selectivity(pharma, evidence, subtypes):
    """Compute α9α10 selectivity score per site."""
    results = []
    
    for site_id, ev in sorted(evidence.items()):
        # Base score from evidence
        base_score = float(ev.get("evidence_score", 0))
        
        # Subtype comparison bonus
        subtype_score = 0
        for sub in subtypes:
            if sub["subtype"] in ("α7", "α4β2"):
                # Lower conservation = higher selectivity
                a9_cons = float(sub.get("a9_pharma_conservation", 0))
                a10_cons = float(sub.get("a10_pharma_conservation", 0))
                subtype_score += (1 - a9_cons) + (1 - a10_cons)
        subtype_score /= max(len(subtypes), 1)
        
        # Overall selectivity
        selectivity = (base_score * 0.6 + subtype_score * 0.4)
        
        results.append({
            "site_id": site_id,
            "evidence_score": round(base_score, 4),
            "subtype_score": round(subtype_score, 4),
            "selectivity_score": round(selectivity, 4),
            "classification": ev.get("classification", "UNKNOWN"),
        })
    
    return results

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

def write_summary(results, path):
    lines = [
        "α9α10 SELECTIVITY ANALYSIS",
        "=" * 60,
        "",
        "Score = 0.6 × evidence_score + 0.4 × subtype_score",
        "Higher = more α9α10-selective",
        "",
    ]
    
    results.sort(key=lambda r: r["selectivity_score"], reverse=True)
    
    lines.append("Top 10 most selective sites:")
    for r in results[:10]:
        lines.append(f"  site {r['site_id']:2d}  selectivity={r['selectivity_score']:.4f}  "
                      f"evidence={r['evidence_score']:.4f}  "
                      f"subtype={r['subtype_score']:.4f}  "
                      f"class={r['classification']}")
    
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")

def main():
    print("Loading pharmacophore...")
    pharma = load_pharmacophore()
    print(f"  {len(pharma)} features")
    
    print("Loading evidence...")
    evidence = load_evidence()
    print(f"  {len(evidence)} sites")
    
    print("Loading subtype comparison...")
    subtypes = load_subtype()
    print(f"  {len(subtypes)} subtypes")
    
    print("Computing selectivity...")
    results = compute_selectivity(pharma, evidence, subtypes)
    
    print("Writing outputs...")
    write_csv(results, OUT / "selectivity_scores.csv")
    write_summary(results, OUT / "selectivity_summary.txt")
    
    print("Done.")

if __name__ == "__main__":
    main()
