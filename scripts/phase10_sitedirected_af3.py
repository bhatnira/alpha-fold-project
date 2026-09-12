#!/usr/bin/env python3
"""
Phase 10: Site-Directed AF3 Modeling

After candidate sites have been discovered, use AF3/Boltz-2 to test whether
each candidate site can support stable ligand placement.

For each top candidate site, generate cofolding predictions with:
- Site 23 (top candidate, HIGH_CONFIDENCE)
- Site 25 (INTERMEDIATE)
- Site 34 (LOW_CONFIDENCE but STRONG convergence)
- Site 21 (LOW_CONFIDENCE but STRONG convergence)
"""
import json, os
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUTDIR = ROOT / "08_site_directed_af3" / "yaml_inputs"
OUTDIR.mkdir(parents=True, exist_ok=True)

SEQUENCES = {
    "alpha9": "MNWSHSCISFCWIYFAASRLRAAETADGKYAQKLFNDLFEDYSNALRPVEDTDKVLNVTLQITLSQIKDMDERNQILTAYLWIRQIWHDAYLTWDRDQYDGLDSIRIPSDLVWRPDIVLYNKADDESSEPVNTNVVLRYDGLITWDAPAITKSSCVVDVTYFPFDNQQCNLTFGSWTYNGNQVDIFNALDSGDLSDFIEDVEWEVHGMPAVKNVISYGCCSEPYPDVTFTLLLKRRSSFYIVNLLIPCVLISFLAPLSFYLPAASGEKVSLGVTILLAMTVFQLMVAEIMPASENVPLIGKYYIATMALITASTALTIMVMNIHFCGAEARPVPHWARVVILKYMSRVLFVYDVGESCLSPHHSRERDHLTKVYSKLPESNLKAARNKDLSRKKDMNKRLKNDLGCQGKNPQEAESYCAQYKVLTRNIEYIAKCLKDHKATNSKGSEWKKVAKVIDRFFMWIFFIMVFVMTILIIARAD",
    "alpha10": "MGLRSHHLSLGLLLLFLLPAECLGAEGRLALKLFRDLFANYTSALRPVADTDQTLNVTLEVTLSQIIDMDERNQVLTLYLWIRQEWTDAYLRWDPNAYGGLDAIRIPSSLVWRPDIVLYNKADAQPPGSASTNVVLRHDGAVRWDAPAITRSSCRVDVAAFPFDAQHCGLTFGSWTHGGHQLDVRPRGAAASLADFVENVEWRVLGMPARRRVLTYGCCSEPYPDVTFTLLLRRRAAAYVCNLLLPCVLISLLAPLAFHLPADSGEKVSLGVTVLLALTVFQLLLAESMPPAESVPLIGKYYMATMTMVTFSTALTILIMNLHYCGPSVRPVPAWARALLLGHLARGLCVRERGEPCGQSRPPELSPSPQSPEGGAGPPAGPCHEPRCLCRQEALLHHVATIANTFRSHRAAQRCHEDWKRLARVMDRFFLAIFFSMALVMSLLVLVQAL"
}

CANDIDATE_SITES = {
    "site23": {
        "description": "alpha9/alpha10 ECD interface (top candidate)",
        "stoichiometry": "2to3",
        "classification": "HIGH_CONFIDENCE",
    },
    "site25": {
        "description": "alpha9/alpha9 ECD interface",
        "stoichiometry": "2to3",
        "classification": "INTERMEDIATE",
    },
    "site34": {
        "description": "M1-M2/M2 TM region",
        "stoichiometry": "2to3",
        "classification": "LOW_CONFIDENCE",
    },
    "site21": {
        "description": "alpha10 ECD region",
        "stoichiometry": "2to3",
        "classification": "LOW_CONFIDENCE",
    },
}

LIGANDS = {
    "ascorbate": "OC[C@H](O)[C@H]1OC(=O)C(O)=C1O",
    "cpd12_potent": "CC#CC(=O)OC[C@@H](O)[C@@H]1OC(=O)C(O)=C1O",
    "acetate_control": "CC(=O)[O-]",
}

CHAIN_MAP = {"2to3": {"A": "alpha9", "B": "alpha9", "C": "alpha10", "D": "alpha10", "E": "alpha10"}}


def make_site_directed_yaml(site_name, site_info, ligand_name, ligand_smiles):
    stoich = site_info["stoichiometry"]
    chains = CHAIN_MAP[stoich]

    yaml_lines = [
        "version: 1",
        f"name: sitedirected_{site_name}_{ligand_name}",
        "sequences:",
    ]
    for chain_id, subunit in chains.items():
        yaml_lines.append("- protein:")
        yaml_lines.append(f"    id: {chain_id}")
        yaml_lines.append(f"    sequence: {SEQUENCES[subunit]}")
    yaml_lines.append("- ligand:")
    yaml_lines.append("    id: LIG")
    yaml_lines.append(f"    smiles: {ligand_smiles}")

    yaml_path = OUTDIR / f"sd_{site_name}_{ligand_name}.yaml"
    yaml_path.write_text("\n".join(yaml_lines) + "\n")
    return yaml_path


def main():
    print("Phase 10: Site-Directed AF3 Modeling")
    print("=" * 60)

    yamls = []
    for site_name, site_info in CANDIDATE_SITES.items():
        for ligand_name, ligand_smiles in LIGANDS.items():
            yaml_path = make_site_directed_yaml(site_name, site_info, ligand_name, ligand_smiles)
            yamls.append(yaml_path)
            print(f"  {site_name} + {ligand_name}")

    print(f"\nTotal YAML inputs: {len(yamls)}")

    submission_script = ROOT / "scripts" / "submit_sitedirected_gpu.sh"
    submission_script.write_text(f"""#!/bin/bash
#SBATCH --job-name=a9a10_sdaf3
#SBATCH --partition=gpu
#SBATCH --gres=gpu:a100:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=23:00:00
#SBATCH --exclude=pax007
#SBATCH --output={ROOT}/logs/sitedirected_%j.log

set -e

module purge 2>/dev/null || true
module load boltz/2.2.1-gpu

cd {OUTDIR}
mkdir -p predictions

echo "=== Site-directed AF3 started: $(date) ==="
echo "Node: $(hostname)"
nvidia-smi -L | head -2

for yaml_file in {OUTDIR}/sd_*.yaml; do
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
        --diffusion_samples 3 \\
        --output_format pdb \\
        --num_workers 1 \\
        --preprocessing-threads 1
done

echo ""
echo "=== Site-directed AF3 complete: $(date) ==="
""")
    print(f"\n  SLURM script: {submission_script}")
    print("  Submit with: sbatch", submission_script.name)


if __name__ == "__main__":
    main()
