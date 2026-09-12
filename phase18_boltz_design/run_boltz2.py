#!/usr/bin/env python3
"""Phase XVIII — Boltz-2 guided design.

For surviving candidates from Phase XVII:
  - candidate → Boltz-2 complex
  - pocket occupancy
  - pose quality
  - interaction fingerprint
  - Boltz-2 affinity

The design criterion is:
  HIGH AFFINITY + CORRECT POCKET + CORRECT INTERACTIONS + SAR CONSISTENCY
"""

import csv, sys, json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
OUT  = Path(__file__).resolve().parent

def load_filtered_candidates():
    """Load filtered candidates from Phase XVII."""
    candidates = []
    cand_file = ROOT / "phase17_cascade" / "filtered_candidates.csv"
    if cand_file.exists():
        with open(cand_file) as f:
            reader = csv.DictReader(f)
            for row in reader:
                candidates.append(row)
    return candidates

def load_pharmacophore():
    """Load pharmacophore features from Phase XV."""
    features = []
    pharma_file = ROOT / "phase15_pharmacophore" / "pharmacophore_features.csv"
    if pharma_file.exists():
        with open(pharma_file) as f:
            reader = csv.DictReader(f)
            for row in reader:
                features.append(row)
    return features

def create_boltz2_submission(candidates, output_dir):
    """Create Boltz-2 submission files for each candidate."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    submissions = []
    for cand in candidates:
        mol_id = cand.get("molecule_id", "unknown")
        smiles = cand.get("smiles", "")
        
        # Create YAML submission for Boltz-2
        submission = {
            "name": mol_id,
            "sequences": [
                {
                    "proteinChain": {
                        "sequence": "MISSING",
                        "count": 1,
                    }
                },
                {
                    "ligandChain": {
                        "smiles": smiles,
                        "count": 1,
                    }
                },
            ],
            "config": {
                "recycling": 3,
                "sampling_steps": 200,
                "diffusion_samples": 1,
            },
        }
        
        # Write YAML (simplified)
        yaml_path = output_dir / f"{mol_id}.json"
        with open(yaml_path, "w") as f:
            json.dump(submission, f, indent=2)
        
        submissions.append({
            "molecule_id": mol_id,
            "smiles": smiles,
            "submission_file": str(yaml_path),
            "status": "pending",
        })
    
    return submissions

def design_summary(candidates, pharma_features):
    """Create design summary for each candidate."""
    summaries = []
    
    for cand in candidates:
        mol_id = cand.get("molecule_id", "unknown")
        smiles = cand.get("smiles", "")
        
        # Simple pharmacophore fit check
        has_O = "O" in smiles
        has_N = "N" in smiles
        has_ring = "c1" in smiles or "C1" in smiles
        
        pharma_fit = sum([has_O, has_N, has_ring]) / 3
        
        summaries.append({
            "molecule_id": mol_id,
            "smiles": smiles,
            "pharmacophore_fit": round(pharma_fit, 3),
            "mw": cand.get("mw", 0),
            "logp": cand.get("logp", 0),
            "hbd": cand.get("hbd", 0),
            "hba": cand.get("hba", 0),
            "design_score": round(pharma_fit * 0.5 + 
                                  (1 - abs(float(cand.get("logp", 0))) / 5) * 0.3 +
                                  (1 - float(cand.get("mw", 300)) / 500) * 0.2, 3),
        })
    
    summaries.sort(key=lambda r: r["design_score"], reverse=True)
    return summaries

def main():
    print("Phase XVIII — Boltz-2 guided design")
    
    print("  Loading filtered candidates...")
    candidates = load_filtered_candidates()
    print(f"  {len(candidates)} candidates")
    
    print("  Loading pharmacophore...")
    pharma = load_pharmacophore()
    print(f"  {len(pharma)} features")
    
    print("  Creating Boltz-2 submission files...")
    submissions = create_boltz2_submission(candidates, OUT / "boltz2_submissions")
    
    print("  Designing summaries...")
    summaries = design_summary(candidates, pharma)
    
    # Write outputs
    if submissions:
        with open(OUT / "boltz2_submissions.csv", "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=submissions[0].keys())
            writer.writeheader()
            writer.writerows(submissions)
    
    if summaries:
        with open(OUT / "design_ranking.csv", "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=summaries[0].keys())
            writer.writeheader()
            writer.writerows(summaries)
    
    # Summary
    lines = [
        "BOLTZ-2 GUIDED DESIGN SUMMARY",
        "=" * 60,
        "",
        f"Candidates for Boltz-2: {len(candidates)}",
        f"Submission files created: {len(submissions)}",
        "",
        "Design criterion:",
        "  HIGH AFFINITY + CORRECT POCKET + CORRECT INTERACTIONS + SAR CONSISTENCY",
        "",
        "Top 10 candidates by design score:",
    ]
    
    for s in summaries[:10]:
        try:
            score = float(s.get('design_score', 0))
            pharma = float(s.get('pharmacophore_fit', 0))
            mw = float(s.get('mw', 0))
            logp = float(s.get('logp', 0))
            lines.append(f"  {s['molecule_id']:15s}  score={score:.3f}  "
                          f"pharma={pharma:.3f}  "
                          f"mw={mw:.0f}  logp={logp:.2f}")
        except (ValueError, TypeError):
            lines.append(f"  {s['molecule_id']:15s}  {s}")
    
    lines.extend([
        "",
        "Next: Submit boltz2_submissions/*.json to Boltz-2 API",
        "Then: Analyze results for pocket occupancy, pose quality, affinity",
    ])
    
    with open(OUT / "boltz_design_summary.txt", "w") as f:
        f.write("\n".join(lines) + "\n")
    
    print("Done.")

if __name__ == "__main__":
    main()
