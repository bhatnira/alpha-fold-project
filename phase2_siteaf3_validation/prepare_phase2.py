#!/usr/bin/env python3
"""
Phase 2 Prepare: Generate SiteAF3 inputs + analysis scripts
Runs AFTER current Boltz-2 cofolding job (3332533) completes.

Produces:
  1. Hotspot/pocket PDBs for 5 candidate sites
  2. SiteAF3 JSON configs for all site × compound × stoichiometry
  3. Analysis script for Boltz-2 cofolding results
  4. Master orchestration SLURM script
"""

import csv
import json
import os
from pathlib import Path

# ============================================================
# PATHS
# ============================================================

PIPELINE = Path("/cluster/home/nbhatt04/lean_pipeline")
PHASE2 = PIPELINE / "phase2_siteaf3_validation"
SITEAF3_DIR = PHASE2 / "siteaf3_inputs"
RESULTS_DIR = PHASE2 / "results"
SCRIPTS_DIR = PHASE2 / "scripts"
HOTSPOT_DIR = PHASE2 / "hotspot_pocket_pdbs"
RECEPTOR_2TO3 = PIPELINE / "09_docking" / "receptor_2to3.pdb"
RECEPTOR_3TO2 = PIPELINE / "09_docking" / "receptor_3to2.pdb"
COFOLDING_RESULTS = PIPELINE / "site_directed_cofolding" / "results"
COFOLDING_YAMLS = PIPELINE / "site_directed_cofolding" / "yaml_inputs"

# ============================================================
# SITE DEFINITIONS — from convergence + evidence matrix
# ============================================================

SITES = {
    23: {
        "interface": "alpha9_plus_alpha10_minus",
        "evidence_score": 0.627,
        "hotspot_residues": [
            ("alpha9", 176), ("alpha10", 145), ("alpha10", 83),
            ("alpha9", 120), ("alpha10", 81), ("alpha9", 224),
        ],
        "pocket_residues": [
            ("alpha9", 176), ("alpha9", 175), ("alpha9", 177),
            ("alpha9", 120), ("alpha9", 121), ("alpha9", 119),
            ("alpha9", 224), ("alpha9", 223), ("alpha9", 225),
            ("alpha10", 145), ("alpha10", 143), ("alpha10", 146),
            ("alpha10", 83), ("alpha10", 81), ("alpha10", 84),
            ("alpha10", 195), ("alpha10", 62),
        ],
        "af3_convergence": 0.14,
        "boltz2_convergence": 1.0,
        "docking": 1.0,
        "sar": 0.9,
    },
    21: {
        "interface": "alpha10_plus_alpha9_minus",
        "evidence_score": 0.535,
        "hotspot_residues": [
            ("alpha10", 103), ("alpha10", 105), ("alpha10", 30),
            ("alpha9", 178), ("alpha9", 47),
        ],
        "pocket_residues": [
            ("alpha10", 103), ("alpha10", 105), ("alpha10", 104),
            ("alpha10", 30), ("alpha10", 31), ("alpha10", 29),
            ("alpha9", 178), ("alpha9", 177), ("alpha9", 179),
            ("alpha9", 47), ("alpha9", 46), ("alpha9", 48),
            ("alpha9", 82),
        ],
        "af3_convergence": 0.33,
        "boltz2_convergence": 0.04,
        "docking": 1.0,
        "sar": 1.0,
    },
    5: {
        "interface": "alpha9_plus_alpha10_minus",
        "evidence_score": 0.511,
        "hotspot_residues": [
            ("alpha9", 176), ("alpha10", 145), ("alpha10", 83),
            ("alpha9", 120), ("alpha10", 81),
        ],
        "pocket_residues": [
            ("alpha9", 176), ("alpha9", 175), ("alpha9", 177),
            ("alpha9", 120), ("alpha9", 121),
            ("alpha10", 145), ("alpha10", 143), ("alpha10", 146),
            ("alpha10", 83), ("alpha10", 81), ("alpha10", 84),
        ],
        "af3_convergence": 0.89,
        "boltz2_convergence": 0.0,
        "docking": 1.0,
        "sar": 1.0,
    },
    34: {
        "interface": "alpha9_plus_alpha9_minus",
        "evidence_score": 0.477,
        "hotspot_residues": [
            ("alpha10", 273), ("alpha10", 276),
            ("alpha9", 270), ("alpha9", 274), ("alpha9", 277),
        ],
        "pocket_residues": [
            ("alpha10", 273), ("alpha10", 276), ("alpha10", 274),
            ("alpha9", 270), ("alpha9", 274), ("alpha9", 277),
            ("alpha9", 271), ("alpha9", 273),
        ],
        "af3_convergence": 0.74,
        "boltz2_convergence": 0.04,
        "docking": 1.0,
        "sar": 1.0,
    },
    7: {
        "interface": "alpha9_plus_alpha10_minus",
        "evidence_score": 0.476,
        "hotspot_residues": [
            ("alpha9", 176), ("alpha10", 145), ("alpha10", 83),
            ("alpha9", 120),
        ],
        "pocket_residues": [
            ("alpha9", 176), ("alpha9", 175), ("alpha9", 177),
            ("alpha9", 120), ("alpha9", 121),
            ("alpha10", 145), ("alpha10", 143), ("alpha10", 146),
            ("alpha10", 83), ("alpha10", 81),
        ],
        "af3_convergence": 0.25,
        "boltz2_convergence": 0.0,
        "docking": 1.0,
        "sar": 1.0,
    },
}

# ============================================================
# COMPOUNDS
# ============================================================

COMPOUNDS = {
    "L_ASCORBATE": "OC[C@H](O)[C@H]1OC(=O)C(O)=C1O",
    "CPD12_POTENT": "C(#C)COC1[C@@H]([C@@H]2OC(C)(C)OC2)OC(=O)C=1O",
    "CPD25_STRONG": "C(Br)[C@H]([C@H]1OC(=O)C(O)=C1O)O",
    "CPD1_MODERATE": "C1([C@@H]([C@H](O)CO)OC(=O)C=1O)O",
    "CPD18_MODERATE": "C1C=CC(COC2[C@@H]([C@H](O)CO)OC(=O)C=2O)=CC=1",
    "CPD3_WEAK": "CC1(O[C@H]([C@H]2OC(=O)C(O)=C2O)CO1)C",
    "CPD8_INACTIVE": "CCCOC1C(C2OC(C)(C)OC2)OC(=O)C=1O",
    "CPD9_INACTIVE": "CCCCOC1[C@@H]([C@@H]2OC(C)(C)OC2)OC(=O)C=1O",
}

# Chain mapping: which chains are alpha9 vs alpha10 in each stoichiometry
# 2to3: A,B=alpha10; C,D,E=alpha9
# 3to2: A,B,C=alpha10; D,E=alpha9
CHAIN_MAP = {
    "2to3": {"alpha10": ["A", "B"], "alpha9": ["C", "D", "E"]},
    "3to2": {"alpha10": ["A", "B", "C"], "alpha9": ["D", "E"]},
}


def extract_residue_from_pdb(pdb_path, subunit, residue_number, chain_for_subunit):
    """Extract a single residue from PDB as a string line."""
    lines = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                chain = line[21].strip()
                res_seq = int(line[22:26].strip())
                if chain == chain_for_subunit and res_seq == residue_number:
                    lines.append(line.rstrip())
    return lines


def generate_site_pdb(site_id, site_info, receptor_pdb, stoich, output_dir):
    """Generate hotspot and pocket PDB files for a site."""
    chain_map = CHAIN_MAP[stoich]

    hotspot_lines = []
    pocket_lines = []

    for subunit, resnum in site_info["hotspot_residues"]:
        chains = chain_map.get(subunit, [])
        for chain in chains:
            res_lines = extract_residue_from_pdb(receptor_pdb, subunit, resnum, chain)
            hotspot_lines.extend(res_lines)

    for subunit, resnum in site_info["pocket_residues"]:
        chains = chain_map.get(subunit, [])
        for chain in chains:
            res_lines = extract_residue_from_pdb(receptor_pdb, subunit, resnum, chain)
            pocket_lines.extend(res_lines)

    # Write hotspot PDB
    hotspot_path = output_dir / f"site{site_id}_{stoich}_hotspot.pdb"
    with open(hotspot_path, "w") as f:
        f.write("HEADER HOTSPOT RESIDUES FOR SITEAF3\n")
        for line in hotspot_lines:
            f.write(line + "\n")
        f.write("END\n")

    # Write pocket PDB
    pocket_path = output_dir / f"site{site_id}_{stoich}_pocket.pdb"
    with open(pocket_path, "w") as f:
        f.write("HEADER POCKET RESIDUES FOR SITEAF3\n")
        for line in pocket_lines:
            f.write(line + "\n")
        f.write("END\n")

    return hotspot_path, pocket_path


def generate_siteaf3_json(site_id, stoich, compound_name, smiles, hotspot_path, pocket_path, receptor_pdb, seeds, output_dir):
    """Generate SiteAF3 JSON config."""
    config = {
        "name": f"siteaf3_site{site_id}_{stoich}_{compound_name}",
        "receptor": [{
            "rec_struct_path": str(receptor_pdb),
            "fixed_chain_id": CHAIN_MAP[stoich]["alpha9"] + CHAIN_MAP[stoich]["alpha10"],
            "hotspot_path": str(hotspot_path),
            "pocket_path": str(pocket_path),
        }],
        "ligand": [{
            "small_molecule": {
                "id": "LIG",
                "smiles": smiles,
            }
        }],
        "modelSeeds": seeds,
    }

    json_path = output_dir / f"site{site_id}_{stoich}_{compound_name}.json"
    with open(json_path, "w") as f:
        json.dump(config, f, indent=2)
    return json_path


def create_analysis_script():
    """Create the analysis script for Boltz-2 cofolding results."""
    script = '''#!/usr/bin/env python3
"""
Analyze Boltz-2 cofolding results from job 3332533 (site21, 30 compounds).
Compare active vs inactive binding.
"""

import json
import csv
from pathlib import Path
from collections import defaultdict

PIPELINE = Path("/cluster/home/nbhatt04/lean_pipeline")
COFOLDING = PIPELINE / "site_directed_cofolding" / "results"
OUTPUT = PIPELINE / "phase2_siteaf3_validation" / "analysis"

COMPOUNDS = {
    1: {"smiles": "C1([C@@H]([C@H](O)CO)OC(=O)C=1O)O", "activity": "active", "uM": 1797},
    2: {"smiles": "C(O)C(C1OC(=O)C(O)=C1O)O", "activity": "active", "uM": 1316},
    3: {"smiles": "CC1(O[C@H]([C@H]2OC(=O)C(O)=C2O)CO1)C", "activity": "active", "uM": 6077},
    12: {"smiles": "C(#C)COC1[C@@H]([C@@H]2OC(C)(C)OC2)OC(=O)C=1O", "activity": "active", "uM": 0.198},
    18: {"smiles": "C1C=CC(COC2[C@@H]([C@H](O)CO)OC(=O)C=2O)=CC=1", "activity": "active", "uM": 1288},
    24: {"smiles": "CCCOC1C(=O)O[C@H]([C@H](O)CO)C=1O", "activity": "active", "uM": 1202},
    25: {"smiles": "C(Br)[C@H]([C@H]1OC(=O)C(O)=C1O)O", "activity": "active", "uM": 2.63},
}


def parse_confidence(result_dir):
    for p in Path(result_dir).rglob("confidence_*.json"):
        try:
            with open(p) as f:
                return json.load(f)
        except Exception:
            continue
    return {}


def parse_affinity(result_dir):
    for p in Path(result_dir).rglob("affinity_*.json"):
        try:
            with open(p) as f:
                return json.load(f)
        except Exception:
            continue
    return {}


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    
    print("=" * 70)
    print("BOLTZ-2 COFOLDING ANALYSIS — Site21, 30 Compounds")
    print("=" * 70)
    
    # Scan binary and ternary results
    results = {"binary": [], "ternary": []}
    
    for complex_type in ["binary", "ternary"]:
        site_dir = COFOLDING / complex_type
        if not site_dir.exists():
            print(f"  {complex_type}: NOT FOUND")
            continue
        
        for result_dir in sorted(site_dir.iterdir()):
            if not result_dir.is_dir():
                continue
            name = result_dir.name
            if "site21" not in name:
                continue
            
            # Find boltz_results subdirectory
            boltz_dirs = list(result_dir.glob("boltz_results_*"))
            if boltz_dirs:
                actual = boltz_dirs[0] / "predictions"
            else:
                actual = result_dir
            
            conf = parse_confidence(actual)
            aff = parse_affinity(actual)
            
            if conf or aff:
                record = {
                    "name": name,
                    "complex_type": complex_type,
                    "confidence": conf.get("confidence_score"),
                    "iptm": conf.get("iptm"),
                    "ligand_iptm": conf.get("ligand_iptm"),
                    "complex_plddt": conf.get("complex_plddt"),
                    "affinity_pred": aff.get("affinity_pred_value"),
                    "affinity_prob": aff.get("affinity_probability_binary"),
                }
                results[complex_type].append(record)
    
    # Print summary
    for ct in ["binary", "ternary"]:
        recs = results[ct]
        print(f"\\n{ct.upper()}: {len(recs)} predictions")
        if not recs:
            continue
        
        actives = [r for r in recs if any(f"cpd{c}_" in r["name"] or f"cpd{c}." in r["name"] for c in COMPOUNDS if COMPOUNDS[c]["activity"] == "active")]
        inactives = [r for r in recs if any(f"cpd{c}_" in r["name"] or f"cpd{c}." in r["name"] for c in COMPOUNDS if COMPOUNDS[c]["activity"] == "inactive")]
        
        active_confs = [r["confidence"] for r in actives if r.get("confidence")]
        inactive_confs = [r["confidence"] for r in inactives if r.get("confidence")]
        
        active_iptms = [r["ligand_iptm"] for r in actives if r.get("ligand_iptm")]
        inactive_iptms = [r["ligand_iptm"] for r in inactives if r.get("ligand_iptm")]
        
        if active_confs:
            print(f"  Actives (n={len(actives)}): avg_conf={sum(active_confs)/len(active_confs):.3f}")
        if inactive_confs:
            print(f"  Inactives (n={len(inactives)}): avg_conf={sum(inactive_confs)/len(inactive_confs):.3f}")
        if active_iptms:
            print(f"  Actives iptm: {sum(active_iptms)/len(active_iptms):.3f}")
        if inactive_iptms:
            print(f"  Inactives iptm: {sum(inactive_iptms)/len(inactive_iptms):.3f}")
    
    # Write CSV
    all_records = results["binary"] + results["ternary"]
    if all_records:
        csv_path = OUTPUT / "cofolding_site21_results.csv"
        with open(csv_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=all_records[0].keys())
            writer.writeheader()
            writer.writerows(all_records)
        print(f"\\nCSV: {csv_path}")
    
    print("\\n" + "=" * 70)


if __name__ == "__main__":
    main()
'''
    script_path = SCRIPTS_DIR / "analyze_cofolding_results.py"
    with open(script_path, "w") as f:
        f.write(script)
    os.chmod(script_path, 0o755)
    return script_path


def create_master_slurm():
    """Create master SLURM orchestration script."""
    script = '''#!/bin/bash
# =============================================================================
# PHASE 2: SiteAF3 Validation — Master Orchestration
# =============================================================================
# Runs AFTER current jobs complete (3332533, 3327732, 3327733)
# Submit this script when ready.
# =============================================================================

set -e

PIPELINE=/cluster/home/nbhatt04/lean_pipeline
PHASE2=$PIPELINE/phase2_siteaf3_validation
SCRIPTS=$PHASE2/scripts
LOGS=$PHASE2/logs
mkdir -p $LOGS

echo "=========================================="
echo "Phase 2: SiteAF3 Validation"
echo "Started: $(date)"
echo "=========================================="

# Step 1: Analyze current Boltz-2 cofolding results
echo ""
echo "[Step 1] Analyzing Boltz-2 cofolding results..."
python3 $SCRIPTS/analyze_cofolding_results.py
echo "  Analysis complete"

# Step 2: Generate SiteAF3 hotspot/pocket PDBs
echo ""
echo "[Step 2] Generating hotspot/pocket PDBs..."
python3 $SCRIPTS/generate_hotspot_pocket.py
echo "  Hotspot/pocket PDBs generated"

# Step 3: Generate SiteAF3 JSON configs
echo ""
echo "[Step 3] Generating SiteAF3 JSON configs..."
python3 $SCRIPTS/generate_siteaf3_configs.py
echo "  JSON configs generated"

# Step 4: Submit SiteAF3 runs (if SiteAF3 is installed)
if command -v run_SiteAF3.py &> /dev/null || [ -f "$PHASE2/siteaf3_setup_done.flag" ]; then
    echo ""
    echo "[Step 4] Submitting SiteAF3 runs..."
    sbatch $SCRIPTS/submit_siteaf3.sh
    echo "  SiteAF3 jobs submitted"
else
    echo ""
    echo "[Step 4] SiteAF3 not yet configured — skipping"
    echo "  Run: python3 $SCRIPTS/setup_siteaf3.py first"
fi

echo ""
echo "=========================================="
echo "Phase 2 setup complete: $(date)"
echo "=========================================="
'''
    script_path = SCRIPTS_DIR / "master_phase2.sh"
    with open(script_path, "w") as f:
        f.write(script)
    os.chmod(script_path, 0o755)
    return script_path


def main():
    print("=" * 70)
    print("Phase 2 Prepare: SiteAF3 Inputs + Scripts")
    print("=" * 70)

    # Create directories
    for d in [PHASE2, SITEAF3_DIR, RESULTS_DIR, SCRIPTS_DIR, HOTSPOT_DIR,
              RESULTS_DIR / "siteaf3", RESULTS_DIR / "cofolding_analysis"]:
        d.mkdir(parents=True, exist_ok=True)

    # Generate hotspot/pocket PDBs
    print("\n[1] Generating hotspot/pocket PDBs...")
    for stoich in ["2to3", "3to2"]:
        receptor = RECEPTOR_2TO3 if stoich == "2to3" else RECEPTOR_3TO2
        for site_id, site_info in SITES.items():
            hotspot, pocket = generate_site_pdb(
                site_id, site_info, receptor, stoich, HOTSPOT_DIR
            )
            print(f"  Site {site_id} ({stoich}): hotspot={hotspot.name}, pocket={pocket.name}")

    # Generate SiteAF3 JSON configs
    print("\n[2] Generating SiteAF3 JSON configs...")
    seeds = [42, 123, 456]
    json_count = 0
    for stoich in ["2to3", "3to2"]:
        receptor = RECEPTOR_2TO3 if stoich == "2to3" else RECEPTOR_3TO2
        for site_id in SITES:
            for cpd_name, smiles in COMPOUNDS.items():
                hotspot = HOTSPOT_DIR / f"site{site_id}_{stoich}_hotspot.pdb"
                pocket = HOTSPOT_DIR / f"site{site_id}_{stoich}_pocket.pdb"
                json_path = generate_siteaf3_json(
                    site_id, stoich, cpd_name, smiles,
                    hotspot, pocket, receptor, seeds, SITEAF3_DIR
                )
                json_count += 1
    print(f"  Generated {json_count} JSON configs (5 sites × 8 compounds × 2 stoich × 1)")

    # Create analysis script
    print("\n[3] Creating analysis script...")
    create_analysis_script()
    print(f"  Script: {SCRIPTS_DIR / 'analyze_cofolding_results.py'}")

    # Create master SLURM script
    print("\n[4] Creating master SLURM script...")
    create_master_slurm()
    print(f"  Script: {SCRIPTS_DIR / 'master_phase2.sh'}")

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    n_hotspot = len(list(HOTSPOT_DIR.glob("*.pdb")))
    print(f"  Hotspot/pocket PDBs: {n_hotspot} files")
    print(f"  SiteAF3 JSON configs: {json_count} files")
    print(f"  Scripts: {len(list(SCRIPTS_DIR.glob('*.py')))} Python + {len(list(SCRIPTS_DIR.glob('*.sh')))} shell")
    print()
    print("NEXT STEPS:")
    print("  1. Wait for job 3332533 to complete")
    print("  2. Run: python3 phase2_siteaf3_validation/scripts/analyze_cofolding_results.py")
    print("  3. Install SiteAF3 if not done: cd phase2_siteaf3_validation/siteaf3_inputs && pip install -e ../SiteAF3/")
    print("  4. Submit: sbatch phase2_siteaf3_validation/scripts/master_phase2.sh")


if __name__ == "__main__":
    main()
