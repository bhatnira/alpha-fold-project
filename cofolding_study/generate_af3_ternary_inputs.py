#!/usr/bin/env python
"""
Generate AF3 cofolding inputs for ternary complexes:
  receptor (α9α10 pentamer) + ACh + PAM analog

Tests Hypothesis 1: Do all analogs share the same/near binding site?
Tests Hypothesis 2: Does interaction correlate with experimental potency?

Representative panel:
  - Compound 12: highly potent (0.198 uM)
  - Compound 25: potent (2.63 uM)
  - Compound 1: moderate (1797 uM)
  - Compound 18: moderate (1288 uM)
  - Compound 3: weak (6077 uM)
  - Compound 8: inactive (0 uM)
  - Compound 9: inactive (0 uM)
  - L-ascorbic acid (parent reference)
  - ACh alone (binary control)

Each compound gets: binary (receptor+PAM) and ternary (receptor+ACh+PAM)
"""

import json
import os
from pathlib import Path

# Configuration
PIPELINE_DIR = Path("/cluster/home/nbhatt04/lean_pipeline")
OUTPUT_DIR = PIPELINE_DIR / "cofolding_study" / "af3_inputs"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Receptor sequences (verified against UniProt Q9UGM1 / Q9GZZ6 on 2026-09-10)
# α9 subunit sequence (human, Q9UGM1, 479 aa)
ALPHA9_SEQ = "MNWSHSCISFCWIYFAASRLRAAETADGKYAQKLFNDLFEDYSNALRPVEDTDKVLNVTLQITLSQIKDMDERNQILTAYLWIRQIWHDAYLTWDRDQYDGLDSIRIPSDLVWRPDIVLYNKADDESSEPVNTNVVLRYDGLITWDAPAITKSSCVVDVTYFPFDNQQCNLTFGSWTYNGNQVDIFNALDSGDLSDFIEDVEWEVHGMPAVKNVISYGCCSEPYPDVTFTLLLKRRSSFYIVNLLIPCVLISFLAPLSFYLPAASGEKVSLGVTILLAMTVFQLMVAEIMPASENVPLIGKYYIATMALITASTALTIMVMNIHFCGAEARPVPHWARVVILKYMSRVLFVYDVGESCLSPHHSRERDHLTKVYSKLPESNLKAARNKDLSRKKDMNKRLKNDLGCQGKNPQEAESYCAQYKVLTRNIEYIAKCLKDHKATNSKGSEWKKVAKVIDRFFMWIFFIMVFVMTILIIARAD"

# α10 subunit sequence (human, Q9GZZ6, 450 aa)
ALPHA10_SEQ = "MGLRSHHLSLGLLLLFLLPAECLGAEGRLALKLFRDLFANYTSALRPVADTDQTLNVTLEVTLSQIIDMDERNQVLTLYLWIRQEWTDAYLRWDPNAYGGLDAIRIPSSLVWRPDIVLYNKADAQPPGSASTNVVLRHDGAVRWDAPAITRSSCRVDVAAFPFDAQHCGLTFGSWTHGGHQLDVRPRGAAASLADFVENVEWRVLGMPARRRVLTYGCCSEPYPDVTFTLLLRRRAAAYVCNLLLPCVLISLLAPLAFHLPADSGEKVSLGVTVLLALTVFQLLLAESMPPAESVPLIGKYYMATMTMVTFSTALTILIMNLHYCGPSVRPVPAWARALLLGHLARGLCVRERGEPCGQSRPPELSPSPQSPEGGAGPPAGPCHEPRCLCRQEALLHHVATIANTFRSHRAAQRCHEDWKRLARVMDRFFLAIFFSMALVMSLLVLVQAL"

# ACh SMILES
ACH_SMILES = "CC[N+](C)(C)CCO"

# Representative compounds
COMPOUNDS = {
    "L_ASCORBATE_PARENT": {
        "smiles": "OC[C@H](O)[C@H]1OC(=O)C(O)=C1O",
        "activity_uM": None,
        "class": "parent_reference",
        "description": "L-ascorbic acid parent"
    },
    "CPD12_HIGHLY_POTENT": {
        "smiles": "C(#C)COC1[C@@H]([C@@H]2OC(C)(C)OC2)OC(=O)C=1O",
        "activity_uM": 0.1981,
        "class": "highly_potent",
        "description": "Most potent analog (0.198 uM)"
    },
    "CPD25_POTENT": {
        "smiles": "C(Br)[C@H]([C@H]1OC(=O)C(O)=C1O)O",
        "activity_uM": 2.63,
        "class": "potent",
        "description": "Potent analog (2.63 uM)"
    },
    "CPD1_MODERATE": {
        "smiles": "C1([C@@H]([C@H](O)CO)OC(=O)C=1O)O",
        "activity_uM": 1797,
        "class": "moderate",
        "description": "Moderate activity (1797 uM)"
    },
    "CPD18_MODERATE": {
        "smiles": "C1C=CC(COC2[C@@H]([C@H](O)CO)OC(=O)C=2O)=CC=1",
        "activity_uM": 1288,
        "class": "moderate",
        "description": "Moderate activity (1288 uM)"
    },
    "CPD3_WEAK": {
        "smiles": "CC1(O[C@H]([C@H]2OC(=O)C(O)=C2O)CO1)C",
        "activity_uM": 6077,
        "class": "weak",
        "description": "Weak activity (6077 uM)"
    },
    "CPD8_INACTIVE": {
        "smiles": "CCCOC1C(C2OC(C)(C)OC2)OC(=O)C=1O",
        "activity_uM": 0,
        "class": "inactive",
        "description": "Inactive compound"
    },
    "CPD9_INACTIVE": {
        "smiles": "CCCCOC1[C@@H]([C@@H]2OC(C)(C)OC2)OC(=O)C=1O",
        "activity_uM": 0,
        "class": "inactive",
        "description": "Inactive compound"
    },
}

# Stoichiometries to test
STOICHIOMETRIES = {
    "2to3": {"alpha9": 2, "alpha10": 3},
    "3to2": {"alpha9": 3, "alpha10": 2},
}

# Model seeds for reproducibility
SEEDS = [42, 123, 456, 789, 1011]


def create_af3_input_json(name, sequences, seeds, output_path):
    """Create AF3 input data JSON."""
    af3_input = {
        "name": name,
        "modelSeeds": seeds,
        "sequences": sequences,
    }
    with open(output_path, "w") as f:
        json.dump(af3_input, f, indent=2)
    return af3_input


def create_binary_input(compound_id, compound_info, stoich_name, stoich_info, seed_list):
    """Create AF3 input for receptor + PAM (binary complex)."""
    sequences = []
    
    # Add protein chains
    for subunit, count in stoich_info.items():
        seq = ALPHA9_SEQ if subunit == "alpha9" else ALPHA10_SEQ
        sequences.append({
            "proteinChain": {
                "sequence": seq,
                "count": count,
            }
        })
    
    # Add PAM ligand
    sequences.append({
        "ligand": {
            "smiles": compound_info["smiles"],
            "count": 1,
        }
    })
    
    name = f"a9a10_{stoich_name}_{compound_id}_binary"
    output_path = OUTPUT_DIR / f"{name}_data.json"
    create_af3_input_json(name, sequences, seed_list, output_path)
    return output_path


def create_ternary_input(compound_id, compound_info, stoich_name, stoich_info, seed_list):
    """Create AF3 input for receptor + ACh + PAM (ternary complex)."""
    sequences = []
    
    # Add protein chains
    for subunit, count in stoich_info.items():
        seq = ALPHA9_SEQ if subunit == "alpha9" else ALPHA10_SEQ
        sequences.append({
            "proteinChain": {
                "sequence": seq,
                "count": count,
            }
        })
    
    # Add ACh ligand
    sequences.append({
        "ligand": {
            "smiles": ACH_SMILES,
            "count": 1,
        }
    })
    
    # Add PAM ligand
    sequences.append({
        "ligand": {
            "smiles": compound_info["smiles"],
            "count": 1,
        }
    })
    
    name = f"a9a10_{stoich_name}_{compound_id}_ternary"
    output_path = OUTPUT_DIR / f"{name}_data.json"
    create_af3_input_json(name, sequences, seed_list, output_path)
    return output_path


def create_ach_only_input(stoich_name, stoich_info, seed_list):
    """Create AF3 input for receptor + ACh only (binary control)."""
    sequences = []
    
    for subunit, count in stoich_info.items():
        seq = ALPHA9_SEQ if subunit == "alpha9" else ALPHA10_SEQ
        sequences.append({
            "proteinChain": {
                "sequence": seq,
                "count": count,
            }
        })
    
    sequences.append({
        "ligand": {
            "smiles": ACH_SMILES,
            "count": 1,
        }
    })
    
    name = f"a9a10_{stoich_name}_ACH_ONLY"
    output_path = OUTPUT_DIR / f"{name}_data.json"
    create_af3_input_json(name, sequences, seed_list, output_path)
    return output_path


def main():
    print("=" * 60)
    print("AF3 Ternary Complex Input Generator")
    print("=" * 60)
    print(f"Output directory: {OUTPUT_DIR}")
    print(f"Compounds: {len(COMPOUNDS)}")
    print(f"Stoichiometries: {list(STOICHIOMETRIES.keys())}")
    print(f"Seeds per job: {len(SEEDS)}")
    print()
    
    all_inputs = []
    
    for stoich_name, stoich_info in STOICHIOMETRIES.items():
        print(f"\n--- Stoichiometry: {stoich_name} ---")
        
        # ACh-only control
        path = create_ach_only_input(stoich_name, stoich_info, SEEDS)
        all_inputs.append(("ACH_ONLY", stoich_name, "ach_only", path))
        print(f"  Created: {path.name}")
        
        for compound_id, compound_info in COMPOUNDS.items():
            # Binary complex (receptor + PAM)
            path = create_binary_input(compound_id, compound_info, stoich_name, stoich_info, SEEDS)
            all_inputs.append((compound_id, stoich_name, "binary", path))
            print(f"  Created: {path.name}")
            
            # Ternary complex (receptor + ACh + PAM)
            path = create_ternary_input(compound_id, compound_info, stoich_name, stoich_info, SEEDS)
            all_inputs.append((compound_id, stoich_name, "ternary", path))
            print(f"  Created: {path.name}")
    
    # Create manifest
    manifest = {
        "total_jobs": len(all_inputs),
        "compounds": list(COMPOUNDS.keys()),
        "stoichiometries": list(STOICHIOMETRIES.keys()),
        "seeds": SEEDS,
        "jobs": [
            {
                "compound": c,
                "stoichiometry": s,
                "complex_type": t,
                "input_file": str(p),
            }
            for c, s, t, p in all_inputs
        ],
    }
    
    manifest_path = OUTPUT_DIR / "manifest.json"
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)
    
    print(f"\n{'=' * 60}")
    print(f"Total AF3 jobs generated: {len(all_inputs)}")
    print(f"Manifest: {manifest_path}")
    print(f"\nBreakdown:")
    print(f"  Binary (receptor+PAM): {sum(1 for _,_,t,_ in all_inputs if t=='binary')}")
    print(f"  Ternary (receptor+ACh+PAM): {sum(1 for _,_,t,_ in all_inputs if t=='ternary')}")
    print(f"  ACh-only control: {sum(1 for _,_,t,_ in all_inputs if t=='ach_only')}")


if __name__ == "__main__":
    main()
