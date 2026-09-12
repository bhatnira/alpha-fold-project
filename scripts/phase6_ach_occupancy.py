#!/usr/bin/env python3
"""
Phase 6: Generate ACh Occupancy State Models

Model the four occupancy states:
1. R (receptor alone) - DONE (existing APO models)
2. R + PAM - DONE (existing AF3/Boltz-2 models)
3. R + ACh - TO DO
4. R + ACh + PAM - TO DO

Generates YAML inputs for Boltz-2 cofolding with ACh at the orthosteric site.
"""
import json, os
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUTDIR = ROOT / "05_ach_occupancy" / "yaml_inputs"
OUTDIR.mkdir(parents=True, exist_ok=True)

SEQUENCES = {
    "alpha9": "MNWSHSCISFCWIYFAASRLRAAETADGKYAQKLFNDLFEDYSNALRPVEDTDKVLNVTLQITLSQIKDMDERNQILTAYLWIRQIWHDAYLTWDRDQYDGLDSIRIPSDLVWRPDIVLYNKADDESSEPVNTNVVLRYDGLITWDAPAITKSSCVVDVTYFPFDNQQCNLTFGSWTYNGNQVDIFNALDSGDLSDFIEDVEWEVHGMPAVKNVISYGCCSEPYPDVTFTLLLKRRSSFYIVNLLIPCVLISFLAPLSFYLPAASGEKVSLGVTILLAMTVFQLMVAEIMPASENVPLIGKYYIATMALITASTALTIMVMNIHFCGAEARPVPHWARVVILKYMSRVLFVYDVGESCLSPHHSRERDHLTKVYSKLPESNLKAARNKDLSRKKDMNKRLKNDLGCQGKNPQEAESYCAQYKVLTRNIEYIAKCLKDHKATNSKGSEWKKVAKVIDRFFMWIFFIMVFVMTILIIARAD",
    "alpha10": "MGLRSHHLSLGLLLLFLLPAECLGAEGRLALKLFRDLFANYTSALRPVADTDQTLNVTLEVTLSQIIDMDERNQVLTLYLWIRQEWTDAYLRWDPNAYGGLDAIRIPSSLVWRPDIVLYNKADAQPPGSASTNVVLRHDGAVRWDAPAITRSSCRVDVAAFPFDAQHCGLTFGSWTHGGHQLDVRPRGAAASLADFVENVEWRVLGMPARRRVLTYGCCSEPYPDVTFTLLLRRRAAAYVCNLLLPCVLISLLAPLAFHLPADSGEKVSLGVTVLLALTVFQLLLAESMPPAESVPLIGKYYMATMTMVTFSTALTILIMNLHYCGPSVRPVPAWARALLLGHLARGLCVRERGEPCGQSRPPELSPSPQSPEGGAGPPAGPCHEPRCLCRQEALLHHVATIANTFRSHRAAQRCHEDWKRLARVMDRFFLAIFFSMALVMSLLVLVQAL"
}

ACH_SMILES = "CC(=O)OCC[N+](C)(C)C"

# Top PAM compounds for ternary modeling
PAM_COMPOUNDS = {
    "ascorbate": "OC[C@H](O)[C@H]1OC(=O)C(O)=C1O",
    "cpd12_best": "CC#CC(=O)OC[C@@H](O)[C@@H]1OC(=O)C(O)=C1O",
    "cpd25_strong": "OC[C@H](O)[C@H]1OC(=O)C(Br)=C1O",
}

STOICHIOMETRIES = {
    "2to3": {"A": "alpha9", "B": "alpha9", "C": "alpha10", "D": "alpha10", "E": "alpha10"},
    "3to2": {"A": "alpha9", "B": "alpha9", "C": "alpha9", "D": "alpha10", "E": "alpha10"},
}


def make_ternary_yaml(stoich_name, chains, pam_name, pam_smiles):
    """Generate YAML for R + ACh + PAM ternary complex."""
    yaml_lines = [
        "version: 1",
        f"name: a9a10_{stoich_name}_ach_{pam_name}",
        "sequences:",
    ]
    for chain_id, subunit in chains.items():
        yaml_lines.append("- protein:")
        yaml_lines.append(f"    id: {chain_id}")
        yaml_lines.append(f"    sequence: {SEQUENCES[subunit]}")
    yaml_lines.append("- ligand:")
    yaml_lines.append("    id: ACH")
    yaml_lines.append(f"    smiles: {ACH_SMILES}")
    yaml_lines.append("- ligand:")
    yaml_lines.append(f"    id: PAM")
    yaml_lines.append(f"    smiles: {pam_smiles}")

    yaml_path = OUTDIR / f"ternary_{stoich_name}_{pam_name}.yaml"
    yaml_path.write_text("\n".join(yaml_lines) + "\n")
    return yaml_path


def make_ach_only_yaml(stoich_name, chains):
    """Generate YAML for R + ACh binary complex."""
    yaml_lines = [
        "version: 1",
        f"name: a9a10_{stoich_name}_ach_only",
        "sequences:",
    ]
    for chain_id, subunit in chains.items():
        yaml_lines.append("- protein:")
        yaml_lines.append(f"    id: {chain_id}")
        yaml_lines.append(f"    sequence: {SEQUENCES[subunit]}")
    yaml_lines.append("- ligand:")
    yaml_lines.append("    id: ACH")
    yaml_lines.append(f"    smiles: {ACH_SMILES}")

    yaml_path = OUTDIR / f"ach_only_{stoich_name}.yaml"
    yaml_path.write_text("\n".join(yaml_lines) + "\n")
    return yaml_path


def main():
    print("Phase 6: Generating ACh occupancy YAML inputs")
    print("=" * 60)

    yamls = []
    for stoich_name, chains in STOICHIOMETRIES.items():
        yamls.append(make_ach_only_yaml(stoich_name, chains))
        print(f"  R+ACh: a9a10_{stoich_name}_ach_only.yaml")

        for pam_name, pam_smiles in PAM_COMPOUNDS.items():
            yamls.append(make_ternary_yaml(stoich_name, chains, pam_name, pam_smiles))
            print(f"  R+ACh+PAM: ternary_{stoich_name}_{pam_name}.yaml")

    print(f"\nTotal YAML inputs: {len(yamls)}")

    submission_script = ROOT / "scripts" / "submit_ach_occupancy_gpu.sh"
    submission_script.write_text(f"""#!/bin/bash
#SBATCH --job-name=a9a10_ach
#SBATCH --partition=gpu
#SBATCH --gres=gpu:a100:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=23:00:00
#SBATCH --exclude=pax007
#SBATCH --output={ROOT}/logs/ach_occupancy_%j.log

set -e

module purge 2>/dev/null || true
module load boltz/2.2.1-gpu

cd {OUTDIR}
mkdir -p predictions

echo "=== ACh occupancy modeling started: $(date) ==="
echo "Node: $(hostname)"
nvidia-smi -L | head -2

for yaml_file in {OUTDIR}/*.yaml; do
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
echo "=== ACh occupancy modeling complete: $(date) ==="
""")
    print(f"\n  SLURM script: {submission_script}")
    print("  Submit with: sbatch", submission_script.name)


if __name__ == "__main__":
    main()
