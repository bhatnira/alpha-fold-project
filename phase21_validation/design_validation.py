#!/usr/bin/env python3
"""Phase XXI — Prospective validation design.

Freezes computational ranking and designs experimental validation plan.
"""

import csv, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
OUT  = Path(__file__).resolve().parent

def load_ranking():
    ranking = []
    rank_file = ROOT / "phase20_ranking" / "site_ranking.csv"
    if rank_file.exists():
        with open(rank_file) as f:
            reader = csv.DictReader(f)
            for row in reader:
                ranking.append(row)
    return ranking

def load_mutagenesis():
    preds = []
    mut_file = DATA / "MUTAGENESIS_PREDICTIONS.csv"
    if mut_file.exists():
        with open(mut_file) as f:
            reader = csv.DictReader(f)
            for row in reader:
                preds.append(row)
    return preds

def design_validation(ranking, mutagenesis):
    """Design experimental validation plan."""
    # Top 3 sites
    top_sites = ranking[:3] if len(ranking) >= 3 else ranking
    
    # High-priority mutations from mutagenesis
    high_priority = [m for m in mutagenesis if m.get("priority") == "high"]
    
    validation_plan = {
        "computational_ranking_frozen": True,
        "top_sites": [],
        "mutation_plan": [],
        "assay_plan": [],
    }
    
    for site in top_sites:
        site_id = int(site["site_id"])
        site_mutations = [m for m in high_priority if int(m["site"]) == site_id]
        
        validation_plan["top_sites"].append({
            "site_id": site_id,
            "rank": site["rank"],
            "evidence_score": site.get("weighted_score", 0),
            "classification": site.get("classification", "UNKNOWN"),
            "n_mutations_to_test": len(site_mutations),
            "key_residues": [f"{m['subunit']}:{m['residue_identity']}{m['residue_number']}" 
                            for m in site_mutations[:5]],
        })
    
    # Mutation plan
    for m in high_priority[:15]:
        validation_plan["mutation_plan"].append({
            "site": m["site"],
            "mutation": f"{m['subunit']}:{m['residue_identity']}{m['residue_number']}→Ala",
            "anchor_type": m["anchor_type"],
            "hypothesis": m["discriminates_hypothesis"],
            "predicted_effect": m["prediction_if_A_specific_recognition"],
            "priority": m["priority"],
        })
    
    # Assay plan
    validation_plan["assay_plan"] = [
        {
            "assay": "Radioligand binding",
            "target": "α9α10 nAChR",
            "ligand": "L-ascorbate",
            "measure": "Kd, Bmax",
            "mutations": "Top 5 high-priority mutations",
        },
        {
            "assay": "Calcium flux",
            "target": "α9α10 nAChR",
            "ligand": "L-ascorbate, D-ascorbate, acetate",
            "measure": "IC50, EC50, efficacy",
            "mutations": "Full mutation panel",
        },
        {
            "assay": "Electrophysiology (patch-clamp)",
            "target": "α9α10 nAChR",
            "ligand": "L-ascorbate + ACh",
            "measure": "PAM/NCAM activity, concentration-response",
            "mutations": "Key contact residues",
        },
        {
            "assay": "Competition binding",
            "target": "α9α10 vs α7 vs α4β2",
            "ligand": "Top 3 candidate compounds",
            "measure": "Selectivity ratio",
            "mutations": "None (wild-type comparison)",
        },
    ]
    
    return validation_plan

def write_json(data, path):
    import json
    with open(path, "w") as f:
        json.dump(data, f, indent=2)

def write_summary(plan, path):
    lines = [
        "PROSPECTIVE VALIDATION PLAN",
        "=" * 60,
        "",
        "FROZEN COMPUTATIONAL RANKING:",
        f"  Top 3 sites: {[s['site_id'] for s in plan['top_sites']]}",
        "",
        "MUTAGENESIS PLAN:",
        f"  {len(plan['mutation_plan'])} mutations to test",
        "",
    ]
    
    for m in plan["mutation_plan"][:10]:
        lines.append(f"  {m['mutation']:30s}  {m['anchor_type']:20s}  {m['hypothesis']}")
    lines.append("")
    
    lines.append("ASSAY PLAN:")
    for a in plan["assay_plan"]:
        lines.append(f"  {a['assay']}")
        lines.append(f"    Target: {a['target']}")
        lines.append(f"    Ligands: {a['ligand']}")
        lines.append(f"    Measure: {a['measure']}")
        lines.append(f"    Mutations: {a['mutations']}")
        lines.append("")
    
    lines.extend([
        "SUCCESS CRITERIA:",
        "",
        "  Hit rate: ≥20% of top-ranked compounds show α9α10 activity",
        "  Enrichment: top-ranked compounds outperform random by ≥5x",
        "  Ranking: Spearman ρ ≥ 0.6 between predicted and experimental",
        "  Selectivity: top compounds show ≥10x α9α10 vs α7",
        "",
        "FAILURE MODES:",
        "",
        "  If top site shows no binding → re-evaluate pocket definition",
        "  If acetate reproduces pattern → downgrade molecular recognition",
        "  If stereochemistry doesn't matter → revise pharmacophore",
    ])
    
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")

def main():
    print("Loading ranking...")
    ranking = load_ranking()
    print(f"  {len(ranking)} sites")
    
    print("Loading mutagenesis predictions...")
    mutagenesis = load_mutagenesis()
    print(f"  {len(mutagenesis)} mutations")
    
    print("Designing validation plan...")
    plan = design_validation(ranking, mutagenesis)
    
    print("Writing outputs...")
    write_json(plan, OUT / "validation_plan.json")
    write_summary(plan, OUT / "validation_summary.txt")
    
    print("Done.")

if __name__ == "__main__":
    main()
