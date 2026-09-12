#!/usr/bin/env python
"""
PHASE 26-29: Candidate Site Competition + Falsification + Model Freeze
Master Prompt Compliance: Sections 30-33

Integrates all evidence to select the best candidate PAM site.
Tests falsification criteria.
Freezes the model.
"""

import json
import csv
import numpy as np
from pathlib import Path
from collections import defaultdict

PIPELINE_DIR = Path("/cluster/home/nbhatt04/lean_pipeline")
PROJECT_DIR = PIPELINE_DIR / "pam_project"
OUTPUT_DIR = PROJECT_DIR / "08_frozen_sites"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Load AF3 discovery results
AF3_RESULTS = PROJECT_DIR / "06_af3_discovery/af3_site_discovery.json"

# Key residues from literature (using mature protein numbering)
# AF3 models use full-length numbering which is +1 for alpha9 due to signal peptide
# So W176 in literature = W175 in AF3, S175 = S174, Y120 = Y119, Y224 = Y223, Y217 = Y216
# Alpha10 numbering is consistent between literature and AF3
KEY_RESIDUES = {
    "alpha9": {"W175": "anchor", "S174": "H-bond", "Y119": "H-bond", "Y223": "hydrophobic", "Y216": "H-bond", "W54": "aromatic"},
    "alpha10": {"D145": "electrostatic", "R83": "electrostatic", "W81": "H-bond", "R143": "electrostatic"},
}


def load_af3_results():
    """Load AF3 site discovery results."""
    with open(AF3_RESULTS, "r") as f:
        return json.load(f)


def build_candidate_site_matrix(af3_data):
    """Build comprehensive candidate site competition matrix."""
    
    clusters = af3_data["clusters"]
    
    # Take top 10 clusters
    top_clusters = clusters[:10]
    
    matrix = []
    
    for i, cluster in enumerate(top_clusters):
        # Check for Site 23 residues (both numbering schemes, both code formats)
        top_residues = cluster["top_residues"]
        
        # Site 23 key residues in AF3 numbering with three-letter codes
        site23_keys = [
            # One-letter codes
            "W175", "S174", "Y119", "Y223", "Y216", "W54",
            "D145", "R83", "W81", "R143",
            # Three-letter codes (AF3 format)
            "TRP175", "SER174", "TYR119", "TYR223", "TYR216", "TRP54",
            "ASP145", "ARG83", "TRP81", "ARG143",
        ]
        
        site23_residues_found = []
        for res_key in site23_keys:
            for r in top_residues:
                if res_key in r:
                    site23_residues_found.append(res_key)
                    break
        
        # Determine interface type
        interface = cluster["interface"]
        if "A" in interface and "C" in interface:
            interface_type = "alpha9(+)/alpha10(-)"
        elif "C" in interface and "A" in interface:
            interface_type = "alpha10(+)/alpha9(-)"
        elif interface.startswith("intra"):
            interface_type = "intra-subunit"
        else:
            interface_type = "multi-subunit"
        
        # Score components
        n_models = cluster["n_models"]
        n_stoich = cluster["recurrence_across_stoichiometries"]
        n_ligands = cluster["recurrence_across_ligands"]
        
        # Site 23 residue overlap score
        site23_score = len(site23_residues_found) / 9.0  # 9 key residues
        
        # Recurrence score (normalized)
        recurrence_score = min(n_models / 36.0, 1.0)  # 36 is max observed
        
        # Multi-stoichiometry bonus
        stoich_bonus = 1.0 if n_stoich > 1 else 0.5
        
        # Multi-ligand bonus
        ligand_bonus = min(n_ligands / 5.0, 1.0)
        
        # Composite score
        composite = (
            0.35 * site23_score +
            0.25 * recurrence_score +
            0.15 * stoich_bonus +
            0.15 * ligand_bonus +
            0.10 * (1.0 if "DERIVATIVE" in str(cluster["ligand_types"]) else 0.0)
        )
        
        # Classification
        if composite > 0.6 and site23_score > 0.3:
            classification = "HIGH_CANDIDATE"
        elif composite > 0.4:
            classification = "MEDIUM_CANDIDATE"
        else:
            classification = "LOW_CANDIDATE"
        
        entry = {
            "rank": i + 1,
            "cluster_id": cluster["cluster_id"],
            "n_models": n_models,
            "interface": interface,
            "interface_type": interface_type,
            "stoichiometries": cluster["stoichiometries"],
            "ligand_types": cluster["ligand_types"],
            "top_residues": top_residues[:10],
            "site23_residues_found": site23_residues_found,
            "site23_score": float(site23_score),
            "recurrence_score": float(recurrence_score),
            "stoich_bonus": float(stoich_bonus),
            "ligand_bonus": float(ligand_bonus),
            "composite_score": float(composite),
            "classification": classification,
            "centroid": cluster["centroid"],
        }
        
        matrix.append(entry)
    
    # Sort by composite score
    matrix.sort(key=lambda x: -x["composite_score"])
    
    # Re-rank
    for i, entry in enumerate(matrix):
        entry["final_rank"] = i + 1
    
    return matrix


def falsification_tests(matrix, af3_data):
    """Apply falsification criteria to leading candidate."""
    
    tests = []
    
    # Test 1: Site 23 residue coverage
    top = matrix[0] if matrix else None
    if top:
        tests.append({
            "test": "Site 23 residue coverage",
            "criterion": "Leading site must contain >= 3 of 9 key Site 23 residues",
            "result": len(top["site23_residues_found"]) >= 3,
            "score": f"{len(top['site23_residues_found'])}/9",
            "verdict": "PASS" if len(top["site23_residues_found"]) >= 3 else "FAIL",
        })
    
    # Test 2: Cross-stoichiometry recurrence
    if top:
        n_stoich = len(top["stoichiometries"])
        tests.append({
            "test": "Cross-stoichiometry recurrence",
            "criterion": "Leading site must appear in both 2to3 and 3to2",
            "result": n_stoich > 1,
            "score": f"{n_stoich} stoichiometries",
            "verdict": "PASS" if n_stoich > 1 else "FAIL",
        })
    
    # Test 3: Cross-ligand recurrence
    if top:
        n_ligands = len(top["ligand_types"])
        tests.append({
            "test": "Cross-ligand recurrence",
            "criterion": "Leading site must accommodate >= 3 different ligand types",
            "result": n_ligands >= 3,
            "score": f"{n_ligands} ligand types",
            "verdict": "PASS" if n_ligands >= 3 else "FAIL",
        })
    
    # Test 4: Active compound recurrence
    if top:
        has_derivative = "DERIVATIVE_EXPERIMENTAL" in top["ligand_types"]
        tests.append({
            "test": "Active compound recurrence",
            "criterion": "Leading site must accommodate the DERIVATIVE_EXPERIMENTAL ligand",
            "result": has_derivative,
            "score": "Yes" if has_derivative else "No",
            "verdict": "PASS" if has_derivative else "FAIL",
        })
    
    # Test 5: Interface plausibility
    if top:
        interface_ok = top["interface_type"] in ["alpha9(+)/alpha10(-)", "alpha10(+)/alpha9(-)", "multi-subunit"]
        tests.append({
            "test": "Interface plausibility",
            "criterion": "Leading site must be at a biologically plausible interface",
            "result": interface_ok,
            "score": top["interface_type"],
            "verdict": "PASS" if interface_ok else "FAIL",
        })
    
    # Test 6: Alternative site competition
    if len(matrix) > 1:
        second = matrix[1]
        margin = top["composite_score"] - second["composite_score"] if top else 0
        tests.append({
            "test": "Alternative site competition",
            "criterion": "Leading site must outscore second-best by > 0.1",
            "result": margin > 0.1,
            "score": f"Margin: {margin:.3f}",
            "verdict": "PASS" if margin > 0.1 else "FAIL",
        })
    
    # Test 7: Minimum model recurrence
    if top:
        tests.append({
            "test": "Minimum model recurrence",
            "criterion": "Leading site must have >= 5 AF3 models",
            "result": top["n_models"] >= 5,
            "score": f"{top['n_models']} models",
            "verdict": "PASS" if top["n_models"] >= 5 else "FAIL",
        })
    
    # Overall
    n_pass = sum(1 for t in tests if t["verdict"] == "PASS")
    n_total = len(tests)
    
    return {
        "tests": tests,
        "n_pass": n_pass,
        "n_total": n_total,
        "overall": "PASS" if n_pass >= n_total * 0.7 else "CONDITIONAL" if n_pass >= n_total * 0.5 else "FAIL",
    }


def freeze_model(matrix, falsification):
    """Freeze the model specification."""
    
    top = matrix[0] if matrix else None
    
    frozen_model = {
        "model_version": "1.0",
        "freeze_date": "2026-09-05",
        "leading_candidate": {
            "cluster_id": top["cluster_id"] if top else None,
            "interface": top["interface"] if top else None,
            "interface_type": top["interface_type"] if top else None,
            "stoichiometries": top["stoichiometries"] if top else None,
            "top_residues": top["top_residues"] if top else None,
            "site23_residues": top["site23_residues_found"] if top else None,
            "composite_score": top["composite_score"] if top else None,
        },
        "falsification": falsification,
        "frozen_parameters": {
            "receptor_models": "AF3 predicted models (2to3 and 3to2)",
            "states": "Apo/ligand-bound (no open/desensitized states available)",
            "candidate_sites": [m["cluster_id"] for m in matrix[:3]],
            "docking_grids": "To be defined based on frozen site centroids",
            "scoring_rules": "Multi-method convergence + SAR consistency",
        },
        "caveats": [
            "Only apo/ligand-bound states modeled (no open/desensitized states)",
            "ACh ternary complexes not yet modeled",
            "Boltz-2 affinity predictions pending",
            "Flexible docking not yet performed",
            "Dataset has 30 compounds (not 28 as in master prompt)",
        ],
    }
    
    return frozen_model


def main():
    print("=" * 70)
    print("PHASE 26-29: CANDIDATE SITE COMPETITION + FALSIFICATION")
    print("=" * 70)
    
    # Load AF3 results
    af3_data = load_af3_results()
    print(f"Loaded {af3_data['n_total_models']} AF3 models in {af3_data['n_clusters']} clusters")
    
    # Build competition matrix
    print("\n--- Building Candidate Site Competition Matrix ---")
    matrix = build_candidate_site_matrix(af3_data)
    
    print(f"\n{'='*80}")
    print(f"{'Rank':<5} {'Cluster':<10} {'Models':<8} {'Interface':<20} {'Site23':<10} {'Recurrence':<12} {'Composite':<10} {'Class'}")
    print(f"{'='*80}")
    
    for entry in matrix[:10]:
        print(f"{entry['final_rank']:<5} {entry['cluster_id']:<10} {entry['n_models']:<8} {entry['interface_type']:<20} {entry['site23_score']:.2f}     {entry['recurrence_score']:.2f}        {entry['composite_score']:.3f}    {entry['classification']}")
    
    # Falsification tests
    print("\n--- Falsification Tests ---")
    falsification = falsification_tests(matrix, af3_data)
    
    for test in falsification["tests"]:
        status = "✓" if test["verdict"] == "PASS" else "✗"
        print(f"  {status} {test['test']}: {test['score']} -> {test['verdict']}")
    
    print(f"\n  Overall: {falsification['n_pass']}/{falsification['n_total']} PASS -> {falsification['overall']}")
    
    # Freeze model
    print("\n--- Freezing Model ---")
    frozen = freeze_model(matrix, falsification)
    
    print(f"  Model version: {frozen['model_version']}")
    print(f"  Leading candidate: Cluster {frozen['leading_candidate']['cluster_id']}")
    print(f"  Interface: {frozen['leading_candidate']['interface_type']}")
    print(f"  Composite score: {frozen['leading_candidate']['composite_score']:.3f}")
    print(f"  Falsification: {frozen['falsification']['overall']}")
    
    print(f"\n  Caveats:")
    for c in frozen["caveats"]:
        print(f"    - {c}")
    
    # Save outputs
    with open(OUTPUT_DIR / "candidate_site_matrix.json", "w") as f:
        json.dump(matrix, f, indent=2, default=str)
    
    with open(OUTPUT_DIR / "falsification_report.json", "w") as f:
        json.dump(falsification, f, indent=2, default=str)
    
    with open(OUTPUT_DIR / "frozen_model.json", "w") as f:
        json.dump(frozen, f, indent=2, default=str)
    
    # Save CSV
    with open(OUTPUT_DIR / "candidate_site_matrix.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "final_rank", "cluster_id", "n_models", "interface_type",
            "site23_score", "recurrence_score", "composite_score", "classification"
        ])
        writer.writeheader()
        for entry in matrix:
            writer.writerow({k: entry[k] for k in writer.fieldnames})
    
    print(f"\nResults saved to: {OUTPUT_DIR}")
    
    return matrix, falsification, frozen


if __name__ == "__main__":
    main()
