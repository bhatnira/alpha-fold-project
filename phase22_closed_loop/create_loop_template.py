#!/usr/bin/env python3
"""Phase XXII — Closed-loop design template.

After experimental validation:
  1. Update SAR with new data
  2. Update pharmacophore
  3. Update interaction rules
  4. Generate new molecules
  5. Run Boltz-2
  6. Run docking
  7. Check selectivity
  8. Test experimentally
  9. Repeat
"""

import csv, sys, json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
OUT  = Path(__file__).resolve().parent

def create_loop_template():
    """Create the closed-loop iteration template."""
    template = {
        "iteration": 0,
        "status": "awaiting_experimental_data",
        "input": {
            "experimental_results": "path/to/experimental/binding/data.csv",
            "updated_sar": "path/to/updated/sar/dataset.csv",
            "new_compounds_tested": [],
        },
        "steps": [
            {
                "step": 1,
                "name": "Update SAR",
                "action": "Add experimental results to SAR dataset",
                "input": "experimental_results.csv",
                "output": "updated_sar_dataset.csv",
                "command": "python phase01_sar/build_sar.py --update experimental_results.csv",
            },
            {
                "step": 2,
                "name": "Update Pharmacophore",
                "action": "Refine pharmacophore based on new SAR",
                "input": "updated_sar_dataset.csv",
                "output": "refined_pharmacophore.pml",
                "command": "python phase15_pharmacophore/build_pharmacophore.py --refine",
            },
            {
                "step": 3,
                "name": "Update Interaction Rules",
                "action": "Recalculate IFPs and SAR rules",
                "input": "updated_sar_dataset.csv + new complex structures",
                "output": "updated_sar_rules.csv",
                "command": "python phase06_fingerprints/build_ifp.py && python phase07_sar_model/derive_sar_rules.py",
            },
            {
                "step": 4,
                "name": "Generate New Molecules",
                "action": "Generate next round of candidates",
                "input": "refined_pharmacophore.pml",
                "output": "new_molecule_library.csv",
                "command": "python phase16_molgen/generate_molecules.py --pharmacophore refined_pharmacophore.pml",
                "note": "Requires RDKit",
            },
            {
                "step": 5,
                "name": "Run Boltz-2",
                "action": "Predict structures and affinities for new candidates",
                "input": "new_molecule_library.csv",
                "output": "boltz2_predictions.csv",
                "command": "boltz2 predict --input new_molecule_library.csv",
                "note": "Requires Boltz-2 API access",
            },
            {
                "step": 6,
                "name": "Run Docking",
                "action": "Dock new candidates into target sites",
                "input": "new_molecule_library.csv + receptor ensemble",
                "output": "docking_results.csv",
                "command": "python phase05_docking/run_docking.py",
                "note": "Requires AutoDock Vina",
            },
            {
                "step": 7,
                "name": "Check Selectivity",
                "action": "Evaluate α9α10 selectivity of new candidates",
                "input": "boltz2_predictions.csv + docking_results.csv",
                "output": "selectivity_scores.csv",
                "command": "python phase19_selectivity/compute_selectivity.py",
            },
            {
                "step": 8,
                "name": "Rank Candidates",
                "action": "Multi-objective ranking of new candidates",
                "input": "all new data",
                "output": "updated_ranking.csv",
                "command": "python phase20_ranking/compute_ranking.py",
            },
            {
                "step": 9,
                "name": "Experimental Testing",
                "action": "Test top 10-30 candidates experimentally",
                "input": "updated_ranking.csv",
                "output": "experimental_results.csv",
                "note": "Requires wet lab access",
            },
            {
                "step": 10,
                "name": "Evaluate and Iterate",
                "action": "Assess hit rate, update model, decide next iteration",
                "input": "experimental_results.csv",
                "output": "iteration_summary.json",
                "command": "python phase22_closed_loop/evaluate_iteration.py",
            },
        ],
        "success_criteria": {
            "hit_rate": "≥20% of top-ranked compounds show α9α10 activity",
            "enrichment": "top compounds outperform random by ≥5x",
            "ranking_correlation": "Spearman ρ ≥ 0.6",
            "selectivity": "top compounds show ≥10x α9α10 vs α7",
        },
        "stop_criteria": {
            "max_iterations": 3,
            "convergence": "Hit rate ≥50% or ranking ρ ≥ 0.8",
            "resource_limit": "Budget or time exhausted",
        },
    }
    return template

def write_json(data, path):
    with open(path, "w") as f:
        json.dump(data, f, indent=2)

def write_summary(template, path):
    lines = [
        "CLOSED-LOOP DESIGN TEMPLATE",
        "=" * 60,
        "",
        "This template defines the iterative cycle for improving",
        "the computational model based on experimental results.",
        "",
        f"Iteration: {template['iteration']}",
        f"Status: {template['status']}",
        "",
        "STEPS:",
        "",
    ]
    
    for step in template["steps"]:
        lines.append(f"  {step['step']:2d}. {step['name']}")
        lines.append(f"      {step['action']}")
        lines.append(f"      Input:  {step['input']}")
        lines.append(f"      Output: {step['output']}")
        if "note" in step:
            lines.append(f"      Note:   {step['note']}")
        lines.append("")
    
    lines.extend([
        "SUCCESS CRITERIA:",
        "",
    ])
    for k, v in template["success_criteria"].items():
        lines.append(f"  {k}: {v}")
    lines.append("")
    
    lines.extend([
        "STOP CRITERIA:",
        "",
    ])
    for k, v in template["stop_criteria"].items():
        lines.append(f"  {k}: {v}")
    
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")

def main():
    print("Creating closed-loop template...")
    template = create_loop_template()
    
    print("Writing outputs...")
    write_json(template, OUT / "closed_loop_template.json")
    write_summary(template, OUT / "closed_loop_summary.txt")
    
    print("Done.")

if __name__ == "__main__":
    main()
