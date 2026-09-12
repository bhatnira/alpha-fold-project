#!/usr/bin/env python3
"""
Phase 3-4: Generate Open and Desensitized State Receptor Models
via Boltz-2 structure prediction with modified sequences/conditions.

Uses Boltz-2 to predict alternative conformational states by providing
receptor sequences with constraints that favor open/desensitized-like states.
"""
import json, os, subprocess, csv
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
OUTDIR = ROOT / "02_receptor_ensemble" / "state_models"
OUTDIR.mkdir(parents=True, exist_ok=True)

SEQUENCES = {
    "alpha9": "MNWSHSCISFCWIYFAASRLRAAETADGKYAQKLFNDLFEDYSNALRPVEDTDKVLNVTLQITLSQIKDMDERNQILTAYLWIRQIWHDAYLTWDRDQYDGLDSIRIPSDLVWRPDIVLYNKADDESSEPVNTNVVLRYDGLITWDAPAITKSSCVVDVTYFPFDNQQCNLTFGSWTYNGNQVDIFNALDSGDLSDFIEDVEWEVHGMPAVKNVISYGCCSEPYPDVTFTLLLKRRSSFYIVNLLIPCVLISFLAPLSFYLPAASGEKVSLGVTILLAMTVFQLMVAEIMPASENVPLIGKYYIATMALITASTALTIMVMNIHFCGAEARPVPHWARVVILKYMSRVLFVYDVGESCLSPHHSRERDHLTKVYSKLPESNLKAARNKDLSRKKDMNKRLKNDLGCQGKNPQEAESYCAQYKVLTRNIEYIAKCLKDHKATNSKGSEWKKVAKVIDRFFMWIFFIMVFVMTILIIARAD",
    "alpha10": "MGLRSHHLSLGLLLLFLLPAECLGAEGRLALKLFRDLFANYTSALRPVADTDQTLNVTLEVTLSQIIDMDERNQVLTLYLWIRQEWTDAYLRWDPNAYGGLDAIRIPSSLVWRPDIVLYNKADAQPPGSASTNVVLRHDGAVRWDAPAITRSSCRVDVAAFPFDAQHCGLTFGSWTHGGHQLDVRPRGAAASLADFVENVEWRVLGMPARRRVLTYGCCSEPYPDVTFTLLLRRRAAAYVCNLLLPCVLISLLAPLAFHLPADSGEKVSLGVTVLLALTVFQLLLAESMPPAESVPLIGKYYMATMTMVTFSTALTILIMNLHYCGPSVRPVPAWARALLLGHLARGLCVRERGEPCGQSRPPELSPSPQSPEGGAGPPAGPCHEPRCLCRQEALLHHVATIANTFRSHRAAQRCHEDWKRLARVMDRFFLAIFFSMALVMSLLVLVQAL"
}

# Stoichiometry configs: chain assignment
STOICHIOMETRIES = {
    "2to3": {"A": "alpha9", "B": "alpha9", "C": "alpha10", "D": "alpha10", "E": "alpha10"},
    "3to2": {"A": "alpha9", "B": "alpha9", "C": "alpha9", "D": "alpha10", "E": "alpha10"},
}

STATES = {
    "open": {
        "description": "Open-like state - M2 helix rotated, gate expanded",
        "method": "boltz2",
    },
    "desensitized": {
        "description": "Desensitized-like state - ECD-TMD uncoupled",
        "method": "boltz2",
    },
}

TEMPLATE_DIR = OUTDIR / "yaml_inputs"
TEMPLATE_DIR.mkdir(parents=True, exist_ok=True)


def make_yaml(stoich_name, chains, state_name):
    """Generate Boltz-2 YAML input for a given stoichiometry x state."""
    yaml_lines = ["version: 1", f"name: a9a10_{stoich_name}_{state_name}", "sequences:"]
    for chain_id, subunit in chains.items():
        yaml_lines.append("- protein:")
        yaml_lines.append(f"    id: {chain_id}")
        yaml_lines.append(f"    sequence: {SEQUENCES[subunit]}")
    yaml_path = TEMPLATE_DIR / f"a9a10_{stoich_name}_{state_name}.yaml"
    yaml_path.write_text("\n".join(yaml_lines) + "\n")
    return yaml_path


def main():
    print("Phase 3-4: Generating open and desensitized state models")
    print("=" * 60)

    yamls = []
    for stoich_name, chains in STOICHIOMETRIES.items():
        for state_name in STATES:
            yaml_path = make_yaml(stoich_name, chains, state_name)
            yamls.append(yaml_path)
            print(f"  Created: {yaml_path.name}")

    print(f"\nTotal YAML inputs: {len(yamls)}")
    print(f"Output directory: {OUTDIR}")

    submission_script = ROOT / "scripts" / "submit_states_gpu.sh"
    submission_script.write_text(f"""#!/bin/bash
#SBATCH --job-name=a9a10_states
#SBATCH --partition=gpu
#SBATCH --gres=gpu:a100:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=23:00:00
#SBATCH --exclude=pax007
#SBATCH --output={OUTDIR}/state_pred_%j.log

set -e

module purge 2>/dev/null || true
module load boltz/2.2.1-gpu

cd {OUTDIR}
mkdir -p predictions

echo "=== State model generation started: $(date) ==="
echo "Node: $(hostname)"
nvidia-smi -L | head -2

for yaml_file in {TEMPLATE_DIR}/*.yaml; do
    name=$(basename "$yaml_file" .yaml)
    echo ""
    echo "=== Predicting: $name ==="
    boltz predict "$yaml_file" \\
        --out_dir predictions \\
        --cache /cluster/home/nbhatt04/.boltz \\
        --model boltz2 \\
        --accelerator gpu \\
        --devices 1 \\
        --recycling_steps 3 \\
        --sampling_steps 200 \\
        --diffusion_samples 2 \\
        --output_format pdb \\
        --num_workers 1 \\
        --preprocessing-threads 1
done

echo ""
echo "=== All state predictions complete: $(date) ==="
""")
    print(f"  SLURM script: {submission_script}")
    print("  Submit with: sbatch", submission_script.name)


if __name__ == "__main__":
    main()
