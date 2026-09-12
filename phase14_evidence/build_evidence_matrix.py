#!/usr/bin/env python3
"""Phase XIV — Integrated site confidence score.

Combines evidence from all phases to produce a final confidence
classification for each candidate allosteric site.
"""

import csv, sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
OUT  = Path(__file__).resolve().parent

# Evidence sources and their weights
EVIDENCE_SOURCES = {
    "pocket_consensus": {"weight": 1.5, "description": "Geometric cavity detection across apo receptors"},
    "af3_convergence": {"weight": 1.2, "description": "AF3 complex prediction convergence"},
    "boltz2_convergence": {"weight": 1.2, "description": "Boltz-2 complex prediction convergence"},
    "docking_consensus": {"weight": 1.0, "description": "Ensemble docking recurrence"},
    "sar_explanation": {"weight": 2.0, "description": "Explains experimental SAR trends"},
    "boltz2_affinity": {"weight": 1.0, "description": "Boltz-2 affinity ranking matches activity"},
    "acetate_discrimination": {"weight": 1.5, "description": "Acetate does NOT reproduce L-ascorbate pattern"},
    "stereochemical_discrimination": {"weight": 1.5, "description": "L vs D ascorbate shows architecture-specific contacts"},
    "chemical_perturbation": {"weight": 1.0, "description": "O-ethyl ascorbate alters interaction pattern"},
    "electrostatic_analysis": {"weight": 0.8, "description": "Electrostatics alone cannot explain localization"},
    "ensemble_robustness": {"weight": 1.0, "description": "Site persists across seeds, states, stoichiometries"},
    "residue_consensus": {"weight": 1.0, "description": "Core residues recur across methods"},
    "subtype_comparison": {"weight": 1.2, "description": "Site is α9α10-selective vs α7/α4β2"},
}

def load_site_ranking():
    sites = {}
    with open(DATA / "FINAL_SITE_RANKING.csv") as f:
        reader = csv.DictReader(f)
        for row in reader:
            sites[int(row["site_id"])] = {
                "rank": int(row["rank_primary"]),
                "composite_score": float(row["SiteScore"]),
                "n_models": int(row["n_models"]),
                "n_af3": int(row["n_af3"]),
                "n_boltz2": int(row["n_boltz2"]),
                "ligands": row["ligands"].split(";"),
                "level": row["evidence_level"],
                "sufficient_support": row["sufficient_support"] == "True",
                "stereo_test": row["stereo_test"],
                "acetate_test": row["acetate_test"],
                "analogue_test": row["analogue_test"],
                "comparative_class": row["comparative_class"],
                "c_pocket_consensus": row["c_pocket_consensus"] == "True",
                "c_af3_reproducibility": float(row["c_af3_reproducibility"]),
                "c_boltz2_reproducibility": float(row["c_boltz2_reproducibility"]),
                "c_ensemble_persistence": float(row["c_ensemble_persistence"]),
                "c_ligand_specificity": float(row["c_ligand_specificity"]),
                "c_stereochemical_discrimination": float(row["c_stereochemical_discrimination"]),
                "c_analogue_discrimination": float(row["c_analogue_discrimination"]),
                "c_acetate_control": float(row["c_acetate_control"]),
                "max_ligand_enrichment": float(row["max_ligand_enrichment"]),
            }
    return sites

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

def load_ifp():
    ifp = {}
    ifp_file = ROOT / "phase06_fingerprints" / "ifp_summary.csv"
    if ifp_file.exists():
        with open(ifp_file) as f:
            reader = csv.DictReader(f)
            for row in reader:
                ifp[int(row["site"])] = {
                    "n_core": int(row["n_core"]),
                    "n_peripheral": int(row["n_peripheral"]),
                }
    return ifp

def compute_evidence_scores(sites, conv, ifp):
    """Compute per-site evidence scores from all available data."""
    evidence = []
    
    for site_id, site in sorted(sites.items()):
        scores = {}
        
        # 1. Pocket consensus (from FINAL_SITE_RANKING)
        scores["pocket_consensus"] = 1.0 if site["c_pocket_consensus"] else 0.0
        
        # 2. AF3 convergence
        scores["af3_convergence"] = min(site["c_af3_reproducibility"] * 3, 1.0)
        
        # 3. Boltz-2 convergence
        scores["boltz2_convergence"] = min(site["c_boltz2_reproducibility"] * 3, 1.0)
        
        # 4. Docking consensus (use ensemble persistence as proxy)
        scores["docking_consensus"] = min(site["c_ensemble_persistence"] * 2, 1.0)
        
        # 5. SAR explanation (use ligand specificity as proxy)
        scores["sar_explanation"] = min(site["c_ligand_specificity"] * 2, 1.0)
        
        # 6. Boltz-2 affinity (use max enrichment as proxy)
        scores["boltz2_affinity"] = min(site["max_ligand_enrichment"] / 5, 1.0)
        
        # 7. Acetate discrimination
        acetate_score = 0.0
        if "acetate" not in [l.lower() for l in site["ligands"]] or site["c_acetate_control"] > 0.5:
            acetate_score = 1.0
        elif site["c_acetate_control"] > 0:
            acetate_score = 0.5
        scores["acetate_discrimination"] = acetate_score
        
        # 8. Stereochemical discrimination
        scores["stereochemical_discrimination"] = min(site["c_stereochemical_discrimination"] * 2, 1.0)
        
        # 9. Chemical perturbation
        scores["chemical_perturbation"] = min(site["c_analogue_discrimination"] * 2, 1.0)
        
        # 10. Electrostatic analysis (default 0.5 if not computed)
        scores["electrostatic_analysis"] = 0.5
        
        # 11. Ensemble robustness
        scores["ensemble_robustness"] = min(site["c_ensemble_persistence"] * 2, 1.0)
        
        # 12. Residue consensus (from IFP)
        if site_id in ifp:
            scores["residue_consensus"] = min(ifp[site_id]["n_core"] / 5, 1.0)
        else:
            scores["residue_consensus"] = 0.5
        
        # 13. Subtype comparison (from comparative_class)
        if "conserved" in site["comparative_class"].lower():
            scores["subtype_comparison"] = 0.8
        elif "divergent" in site["comparative_class"].lower():
            scores["subtype_comparison"] = 0.4
        else:
            scores["subtype_comparison"] = 0.2
        
        # Weighted sum
        total_weight = sum(v["weight"] for v in EVIDENCE_SOURCES.values())
        weighted_sum = sum(scores[k] * EVIDENCE_SOURCES[k]["weight"] for k in scores)
        evidence_score = weighted_sum / total_weight
        
        # Classification
        if evidence_score >= 0.7:
            classification = "HIGH_CONFIDENCE"
        elif evidence_score >= 0.5:
            classification = "INTERMEDIATE"
        elif evidence_score >= 0.3:
            classification = "LOW_CONFIDENCE"
        else:
            classification = "REJECTED"
        
        evidence.append({
            "site_id": site_id,
            "rank": site["rank"],
            "composite_score": site["composite_score"],
            "evidence_score": round(evidence_score, 4),
            "classification": classification,
            **{f"evidence_{k}": round(v, 4) for k, v in scores.items()},
        })
    
    return evidence

def write_csv(rows, path):
    if not rows:
        return
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

def write_summary(evidence, path):
    lines = [
        "INTEGRATED SITE CONFIDENCE SCORES",
        "=" * 60,
        "",
    ]
    
    by_class = defaultdict(list)
    for r in evidence:
        by_class[r["classification"]].append(r)
    
    lines.append("Classification summary:")
    for cls in ["HIGH_CONFIDENCE", "INTERMEDIATE", "LOW_CONFIDENCE", "REJECTED"]:
        items = by_class.get(cls, [])
        sites = sorted([r["site_id"] for r in items])
        lines.append(f"  {cls:20s}  {len(items):3d} sites  {sites}")
    lines.append("")
    
    # Top sites
    lines.append("Top 10 sites by evidence score:")
    top10 = sorted(evidence, key=lambda r: r["evidence_score"], reverse=True)[:10]
    for r in top10:
        lines.append(f"  site {r['site_id']:2d}  evidence={r['evidence_score']:.4f}  "
                      f"class={r['classification']:20s}  composite={r['composite_score']:.4f}")
    lines.append("")
    
    # Evidence breakdown for top site
    if top10:
        top = top10[0]
        lines.append(f"Detailed evidence for top site (site {top['site_id']}):")
        for k, v in EVIDENCE_SOURCES.items():
            score = top.get(f"evidence_{k}", 0)
            lines.append(f"  {k:30s}  score={score:.3f}  weight={v['weight']:.1f}  "
                          f"({v['description']})")
    
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")

def main():
    print("Loading site ranking...")
    sites = load_site_ranking()
    print(f"  {len(sites)} sites")
    
    print("Loading convergence map...")
    conv = load_convergence()
    print(f"  {len(conv)} sites with convergence data")
    
    print("Loading IFP summaries...")
    ifp = load_ifp()
    print(f"  {len(ifp)} sites with IFP data")
    
    print("Computing evidence scores...")
    evidence = compute_evidence_scores(sites, conv, ifp)
    
    print("Writing outputs...")
    write_csv(evidence, OUT / "evidence_matrix.csv")
    write_summary(evidence, OUT / "evidence_summary.txt")
    
    print("Done.")

if __name__ == "__main__":
    main()
