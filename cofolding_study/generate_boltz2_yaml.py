#!/usr/bin/env python
"""
Generate Boltz-2 affinity prediction YAML inputs for ALL 30 compounds.

Tests Hypothesis 1: Do all analogs share the same/near binding site?
Tests Hypothesis 2: Does Boltz-2 affinity correlate with experimental potency?

For each compound at top candidate sites (23, 5, 21):
  - Binary complex (receptor + PAM)
  - Ternary complex (receptor + ACh + PAM)

Output: YAML files ready for `boltz predict` with affinity scoring.
"""

import csv
import yaml
import os
from pathlib import Path

# Configuration
PIPELINE_DIR = Path("/cluster/home/nbhatt04/lean_pipeline")
OUTPUT_DIR = PIPELINE_DIR / "cofolding_study" / "boltz2_inputs"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Receptor sequences (verified against UniProt Q9UGM1 / Q9GZZ6 on 2026-09-10)
# α9 subunit sequence (human, Q9UGM1, 479 aa)
ALPHA9_SEQUENCE = "MNWSHSCISFCWIYFAASRLRAAETADGKYAQKLFNDLFEDYSNALRPVEDTDKVLNVTLQITLSQIKDMDERNQILTAYLWIRQIWHDAYLTWDRDQYDGLDSIRIPSDLVWRPDIVLYNKADDESSEPVNTNVVLRYDGLITWDAPAITKSSCVVDVTYFPFDNQQCNLTFGSWTYNGNQVDIFNALDSGDLSDFIEDVEWEVHGMPAVKNVISYGCCSEPYPDVTFTLLLKRRSSFYIVNLLIPCVLISFLAPLSFYLPAASGEKVSLGVTILLAMTVFQLMVAEIMPASENVPLIGKYYIATMALITASTALTIMVMNIHFCGAEARPVPHWARVVILKYMSRVLFVYDVGESCLSPHHSRERDHLTKVYSKLPESNLKAARNKDLSRKKDMNKRLKNDLGCQGKNPQEAESYCAQYKVLTRNIEYIAKCLKDHKATNSKGSEWKKVAKVIDRFFMWIFFIMVFVMTILIIARAD"
# α10 subunit sequence (human, Q9GZZ6, 450 aa)
ALPHA10_SEQUENCE = "MGLRSHHLSLGLLLLFLLPAECLGAEGRLALKLFRDLFANYTSALRPVADTDQTLNVTLEVTLSQIIDMDERNQVLTLYLWIRQEWTDAYLRWDPNAYGGLDAIRIPSSLVWRPDIVLYNKADAQPPGSASTNVVLRHDGAVRWDAPAITRSSCRVDVAAFPFDAQHCGLTFGSWTHGGHQLDVRPRGAAASLADFVENVEWRVLGMPARRRVLTYGCCSEPYPDVTFTLLLRRRAAAYVCNLLLPCVLISLLAPLAFHLPADSGEKVSLGVTVLLALTVFQLLLAESMPPAESVPLIGKYYMATMTMVTFSTALTILIMNLHYCGPSVRPVPAWARALLLGHLARGLCVRERGEPCGQSRPPELSPSPQSPEGGAGPPAGPCHEPRCLCRQEALLHHVATIANTFRSHRAAQRCHEDWKRLARVMDRFFLAIFFSMALVMSLLVLVQAL"

# Receptor stoichiometry of the α9α10 pentamer: 2 α9 + 3 α10
RECEPTOR_STOICHIOMETRY = [
    ("A", ALPHA9_SEQUENCE),
    ("B", ALPHA9_SEQUENCE),
    ("C", ALPHA10_SEQUENCE),
    ("D", ALPHA10_SEQUENCE),
    ("E", ALPHA10_SEQUENCE),
]

# ACh SMILES
ACH_SMILES = "CC[N+](C)(C)CCO"

# Top candidate sites from phase23 ranking
TOP_SITES = [23, 5, 21]

# Load compound dataset
CSV_PATH = PIPELINE_DIR / "modulator-dataset-a9a10.csv"


def load_compounds():
    """Load all 30 compounds from CSV."""
    compounds = []
    with open(CSV_PATH, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            compounds.append({
                "identifier": int(row["Identifier"]),
                "smiles": row["Smiles"],
                "activity_uM": float(row["Activity (uM)"]),
                "potentiation_pct": float(row["%Potentiation"]),
                "activity_class": classify_compound(
                    float(row["Activity (uM)"]),
                    float(row["%Potentiation"])
                ),
            })
    return compounds


def classify_compound(activity_uM, potentiation_pct):
    """Classify compound by activity."""
    if activity_uM == 0 and potentiation_pct == 0:
        return "inactive"
    elif activity_uM < 1:
        return "highly_potent"
    elif activity_uM < 10:
        return "potent"
    elif activity_uM < 2000:
        return "moderate"
    else:
        return "weak"


def create_boltz2_yaml(compound, site_id, include_ach=False, output_dir=None):
    """Create Boltz-2 YAML input for affinity prediction."""
    
    # Build sequences
    sequences = []
    
    # Protein chains (α9α10 heteropentamer: 2 α9 + 3 α10)
    for chain_id, chain_seq in RECEPTOR_STOICHIOMETRY:
        sequences.append({
            "protein": {
                "id": chain_id,
                "sequence": chain_seq,
                "msa": "empty",
            }
        })
    
    # PAM ligand
    sequences.append({
        "ligand": {
            "id": "PAM",
            "smiles": compound["smiles"],
        }
    })
    
    # Optionally add ACh for ternary complex
    if include_ach:
        sequences.append({
            "ligand": {
                "id": "ACH",
                "smiles": ACH_SMILES,
            }
        })
    
    # Build YAML content
    yaml_content = {
        "version": 1,
        "sequences": sequences,
        "properties": [
            {
                "affinity": {
                    "binder": "PAM",
                }
            }
        ],
    }
    
    # Determine output path
    complex_type = "ternary" if include_ach else "binary"
    filename = f"site{site_id}_cpd{compound['identifier']}_{complex_type}.yaml"
    
    if output_dir is None:
        output_dir = OUTPUT_DIR / f"site{site_id}"
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    filepath = output_dir / filename
    with open(filepath, "w") as f:
        yaml.dump(yaml_content, f, default_flow_style=False, sort_keys=False)
    
    return filepath


def main():
    print("=" * 60)
    print("Boltz-2 Affinity YAML Generator")
    print("=" * 60)
    
    # Load compounds
    compounds = load_compounds()
    print(f"Loaded {len(compounds)} compounds from {CSV_PATH}")
    
    # Summary by class
    classes = {}
    for c in compounds:
        cls = c["activity_class"]
        classes[cls] = classes.get(cls, 0) + 1
    print(f"Classification: {classes}")
    print(f"Top sites: {TOP_SITES}")
    print()
    
    all_yaml_files = []
    
    for site_id in TOP_SITES:
        print(f"\n--- Site {site_id} ---")
        site_dir = OUTPUT_DIR / f"site{site_id}"
        
        for compound in compounds:
            # Binary (receptor + PAM)
            path = create_boltz2_yaml(compound, site_id, include_ach=False, output_dir=site_dir)
            all_yaml_files.append(path)
            
            # Ternary (receptor + ACh + PAM)
            path = create_boltz2_yaml(compound, site_id, include_ach=True, output_dir=site_dir)
            all_yaml_files.append(path)
        
        print(f"  Generated {len(compounds) * 2} YAML files ({len(compounds)} binary + {len(compounds)} ternary)")
    
    # Create manifest
    manifest = {
        "total_yamls": len(all_yaml_files),
        "compounds": len(compounds),
        "sites": TOP_SITES,
        "complex_types": ["binary", "ternary"],
        "files": [str(f) for f in all_yaml_files],
    }
    
    manifest_path = OUTPUT_DIR / "manifest.json"
    with open(manifest_path, "w") as f:
        yaml.dump(manifest, f, default_flow_style=False)
    
    # Create compound summary table
    summary_path = OUTPUT_DIR / "compound_summary.csv"
    with open(summary_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "identifier", "smiles", "activity_uM", "potentiation_pct",
            "activity_class", "n_binary", "n_ternary"
        ])
        for c in compounds:
            writer.writerow([
                c["identifier"], c["smiles"], c["activity_uM"],
                c["potentiation_pct"], c["activity_class"],
                len(TOP_SITES), len(TOP_SITES)
            ])
    
    print(f"\n{'=' * 60}")
    print(f"Total YAML files: {len(all_yaml_files)}")
    print(f"  Binary: {len(all_yaml_files) // 2}")
    print(f"  Ternary: {len(all_yaml_files) // 2}")
    print(f"Manifest: {manifest_path}")
    print(f"Summary: {summary_path}")


if __name__ == "__main__":
    main()
