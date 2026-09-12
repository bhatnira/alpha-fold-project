#!/usr/bin/env python3
"""
Tier 2.2: Generate Open State Templates via Boltz-2
Creates receptor structures in open (activated) state for PAM binding analysis.
Uses Boltz-2 YAML input format for structure prediction.
"""
import csv
import os
import yaml
from pathlib import Path

ROOT = Path("/cluster/home/nbhatt04/lean_pipeline")
OUTDIR = ROOT / "02_receptor_ensemble" / "open_states"
OUTDIR.mkdir(parents=True, exist_ok=True)

PDB_DIR = ROOT / "phase03_pockets" / "a9a10_2to3_apo_out"


def extract_sequences_from_pdb(pdb_path):
    """Extract protein sequences from PDB file."""
    aa3to1 = {
        'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C',
        'GLU':'E','GLN':'Q','GLY':'G','HIS':'H','ILE':'I',
        'LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P',
        'SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V',
    }
    chains = {}
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM"):
                ch = line[21]
                resseq = int(line[22:26])
                resname = line[17:20].strip()
                if ch not in chains:
                    chains[ch] = {}
                if resseq not in chains[ch]:
                    chains[ch][resseq] = resname
    sequences = {}
    for ch, resdict in chains.items():
        seq = ''.join(aa3to1.get(resdict[r], 'X') for r in sorted(resdict.keys()))
        sequences[ch] = seq
    return sequences


def create_boltz2_yaml_open_state(sequences, output_yaml, ligand_smiles=None):
    """Create Boltz-2 YAML input for open state structure prediction."""
    # Build chains list
    chains_list = []
    chain_ids = sorted(sequences.keys())
    
    # Map chains: A/B = alpha9 (2 copies), C/D/E = alpha10 (3 copies)
    alpha9_seq = None
    alpha10_seq = None
    
    for ch_id, seq in sequences.items():
        if alpha9_seq is None:
            alpha9_seq = seq
        elif alpha10_seq is None and seq != alpha9_seq:
            alpha10_seq = seq
    
    if alpha9_seq is None or alpha10_seq is None:
        print(f"  WARNING: Could not distinguish alpha9/alpha10 sequences")
        return None
    
    # 2to3 stoichiometry: 2x alpha9 + 3x alpha10
    chain_config = [
        {"protein": {"id": "A", "sequence": alpha9_seq}},
        {"protein": {"id": "B", "sequence": alpha9_seq}},
        {"protein": {"id": "C", "sequence": alpha10_seq}},
        {"protein": {"id": "D", "sequence": alpha10_seq}},
        {"protein": {"id": "E", "sequence": alpha10_seq}},
    ]
    
    content = {
        "version": 1,
        "name": "a9a10_2to3_open",
        "sequences": chain_config,
    }
    
    if ligand_smiles:
        content["sequences"].append({
            "ligand": {"id": "PAM", "smiles": ligand_smiles}
        })
    
    with open(output_yaml, 'w') as f:
        yaml.dump(content, f, default_flow_style=False, sort_keys=False)
    
    return output_yaml


def main():
    print("=" * 60)
    print("Tier 2.2: Open State Template Generation")
    print("=" * 60)
    
    # Find PDB files
    pdb_files = list(PDB_DIR.glob("*.pdb"))
    print(f"\nFound {len(pdb_files)} PDB files in {PDB_DIR}")
    
    for pdb_path in pdb_files:
        print(f"\nProcessing: {pdb_path.name}")
        sequences = extract_sequences_from_pdb(pdb_path)
        print(f"  Extracted {len(sequences)} chains")
        for ch_id, seq in sorted(sequences.items()):
            print(f"  Chain {ch_id}: {len(seq)} residues")
    
    # Create YAML for 2to3 open state
    print("\n--- Creating Boltz-2 Input: 2to3 Open State ---")
    
    # Use first PDB that has both chains
    for pdb_path in pdb_files:
        sequences = extract_sequences_from_pdb(pdb_path)
        if len(sequences) >= 2:
            yaml_path = OUTDIR / "boltz2_2to3_open.yaml"
            result = create_boltz2_yaml_open_state(sequences, yaml_path)
            if result:
                print(f"  Created: {yaml_path}")
                break
    
    # PAM compound 12 SMILES
    pam_smiles = "C#CCOC1=C(O)C(=O)O[C@@H]1[C@H]1COC(C)(C)O1"
    
    # Create YAML with PAM for ternary complex
    print("\n--- Creating Boltz-2 Input: Ternary Complex (2to3 + PAM) ---")
    for pdb_path in pdb_files:
        sequences = extract_sequences_from_pdb(pdb_path)
        if len(sequences) >= 2:
            yaml_path = OUTDIR / "boltz2_2to3_ternary.yaml"
            result = create_boltz2_yaml_open_state(sequences, yaml_path, pam_smiles)
            if result:
                print(f"  Created: {yaml_path}")
                break
    
    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"\nGenerated files in: {OUTDIR}/")
    print("\nTo run Boltz-2 predictions:")
    print(f"  boltz predict {OUTDIR}/boltz2_2to3_open.yaml --out_dir {OUTDIR}/predictions/")
    print(f"  boltz predict {OUTDIR}/boltz2_2to3_ternary.yaml --out_dir {OUTDIR}/predictions/")
    print("\nOr submit via ColabFold web interface.")
    print("=" * 60)


if __name__ == "__main__":
    main()
