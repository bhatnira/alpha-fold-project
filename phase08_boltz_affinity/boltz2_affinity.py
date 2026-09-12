#!/usr/bin/env python3
"""Boltz-2 affinity prediction for alpha9alpha10 nAChR.

Generates YAML inputs for Boltz-2, runs inference on GPU, and collects
affinity predictions (log10 IC50 in uM) for all ligand × site combinations.
"""

import csv, json, os, subprocess, sys, yaml
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PYTHON = str(Path.home() / ".venvs/pymol/bin/python")

LIGANDS = {
    "ascorbate": "OC[C@H](O)[C@H]1OC(=O)C(O)=C1O",
    "O-ethyl_ascorbate": "CCOC(=O)C(=O)OC[C@@H](O)[C@@H]1OC(=O)C(O)=C1O",
    "acetate": "CC(=O)[O-]",
    "ryanodine": "CC1(C)CC(C)(C)C(=O)OC2CC(C)CCC3C(C)CCC4(O)C5CC6OC(O)C(C)C6OC5C43C21",
}

# Receptor chains: A/B = alpha9, C/D/E = alpha10
PROTEIN_CHAINS = {
    "alpha9": None,   # filled at runtime from PDB
    "alpha10": None,
}


def extract_sequences(pdb_path):
    """Extract protein sequences from PDB."""
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


def make_yaml(protein_seq, ligand_name, ligand_smiles, output_yaml,
              protein_id="A", ligand_id="L", msa="empty"):
    """Create a Boltz-2 YAML input file with affinity prediction."""
    content = {
        "version": 1,
        "sequences": [
            {
                "protein": {
                    "id": protein_id,
                    "sequence": protein_seq,
                    "msa": msa,
                }
            },
            {
                "ligand": {
                    "id": ligand_id,
                    "smiles": ligand_smiles,
                }
            }
        ],
        "properties": [
            {"affinity": {"binder": ligand_id}}
        ]
    }
    with open(output_yaml, 'w') as f:
        yaml.dump(content, f, default_flow_style=False, sort_keys=False)
    return output_yaml


def run_boltz(yaml_path, out_dir, gpu_id=0):
    """Run boltz predict on a single YAML file."""
    env = os.environ.copy()
    env["CUDA_VISIBLE_DEVICES"] = str(gpu_id)
    cmd = [
        PYTHON, "-m", "boltz.main",
        "predict", str(yaml_path),
        "--out_dir", str(out_dir),
        "--model", "boltz2",
        "--accelerator", "gpu",
        "--devices", "1",
        "--recycling_steps", "3",
        "--sampling_steps", "200",
        "--diffusion_samples", "1",
        "--diffusion_samples_affinity", "5",
        "--sampling_steps_affinity", "200",
        "--output_format", "pdb",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=1800)
    return result


def parse_affinity(out_dir, yaml_stem):
    """Parse affinity output from Boltz-2."""
    pred_dir = Path(out_dir) / f"boltz_results_{yaml_stem}" / "predictions" / yaml_stem
    if not pred_dir.exists():
        # Try alternate naming
        for d in Path(out_dir).rglob("affinity*.json"):
            pred_dir = d.parent
            break
    if not pred_dir.exists():
        return None

    aff_file = pred_dir / "affinity_complex.json"
    if not aff_file.exists():
        # Search more broadly
        for f in Path(out_dir).rglob("affinity*.json"):
            aff_file = f
            break
    if not aff_file.exists():
        return None

    with open(aff_file) as f:
        data = json.load(f)
    return data


def main():
    print("Boltz-2 Affinity Prediction Pipeline")
    print("=" * 60)

    pdb_path = ROOT / "data/deliverable/05_structures/2to3/a9a10_2to3_apo_reference.pdb"
    dock_csv = ROOT / "phase05_docking/docking_results.csv"
    boltz_dir = ROOT / "phase08_boltz_affinity"
    boltz_dir.mkdir(exist_ok=True)

    # Extract sequences
    print("Extracting receptor sequences...")
    sequences = extract_sequences(str(pdb_path))
    print(f"  Chains: {', '.join(f'{ch}:{len(seq)}aa' for ch, seq in sorted(sequences.items()))}")

    # Use chain C (alpha10) as representative
    alpha10_seq = sequences.get("C", sequences.get("D", ""))
    alpha9_seq = sequences.get("A", sequences.get("B", ""))
    print(f"  alpha9  (chain A): {len(alpha9_seq)} residues")
    print(f"  alpha10 (chain C): {len(alpha10_seq)} residues")

    # Load docking results to identify top sites
    dock_results = []
    if dock_csv.exists():
        with open(dock_csv) as f:
            for row in csv.DictReader(f):
                dock_results.append(row)
        # Find best ligand per site by binding energy
        best_per_site = {}
        for r in dock_results:
            if r["success"] == "True" and r["binding_energy"]:
                site = r["site_id"]
                energy = float(r["binding_energy"])
                if site not in best_per_site or energy < best_per_site[site][1]:
                    best_per_site[site] = (r["ligand"], energy)
        # Sort sites by best energy
        ranked_sites = sorted(best_per_site.items(), key=lambda x: x[1][1])
        print(f"\nDocking results: {len(dock_results)} total, {len(ranked_sites)} sites")
        print("Top 5 sites by docking energy:")
        for site_id, (lig, energy) in ranked_sites[:5]:
            print(f"  Site {site_id}: {lig} E={energy:.2f} kcal/mol")

    # Generate YAML files for Boltz-2
    # Test with top 3 sites and all ligands
    yaml_dir = boltz_dir / "yaml_inputs"
    yaml_dir.mkdir(exist_ok=True)
    out_dir = boltz_dir / "boltz2_outputs"
    out_dir.mkdir(exist_ok=True)

    test_cases = []
    sites_to_test = [s[0] for s in ranked_sites[:3]] if ranked_sites else ["23"]

    for site_id in sites_to_test:
        for lig_name, lig_smiles in LIGANDS.items():
            yaml_name = f"site{site_id}_{lig_name}"
            yaml_path = yaml_dir / f"{yaml_name}.yaml"
            make_yaml(alpha10_seq, lig_name, lig_smiles, str(yaml_path))
            test_cases.append({
                "site_id": site_id,
                "ligand": lig_name,
                "smiles": lig_smiles,
                "yaml": str(yaml_path),
                "yaml_stem": yaml_name,
                "out_dir": str(out_dir),
            })

    print(f"\nGenerated {len(test_cases)} YAML files in {yaml_dir}")

    # Run Boltz-2 predictions
    results = []
    for i, tc in enumerate(test_cases):
        print(f"\n[{i+1}/{len(test_cases)}] Running Boltz-2: site {tc['site_id']} + {tc['ligand']}...")
        result = run_boltz(tc["yaml"], tc["out_dir"])

        aff_data = parse_affinity(tc["out_dir"], tc["yaml_stem"])
        if aff_data:
            pred_val = aff_data.get("affinity_pred_value", None)
            prob_bin = aff_data.get("affinity_probability_binary", None)
            print(f"  Affinity: pred={pred_val}, prob_binder={prob_bin}")
            results.append({
                "site_id": tc["site_id"],
                "ligand": tc["ligand"],
                "smiles": tc["smiles"],
                "affinity_pred_value": pred_val,
                "affinity_probability_binary": prob_bin,
                "boltz2_success": True,
            })
        else:
            print(f"  FAILED (rc={result.returncode})")
            if result.stderr:
                print(f"  stderr: {result.stderr[:200]}")
            results.append({
                "site_id": tc["site_id"],
                "ligand": tc["ligand"],
                "smiles": tc["smiles"],
                "affinity_pred_value": None,
                "affinity_probability_binary": None,
                "boltz2_success": False,
            })

    # Write results
    if results:
        csv_path = boltz_dir / "boltz2_affinity_results.csv"
        with open(csv_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=results[0].keys())
            writer.writeheader()
            writer.writerows(results)
        print(f"\nResults: {csv_path}")

        # Summary
        n_ok = sum(1 for r in results if r["boltz2_success"])
        print(f"\n{n_ok}/{len(results)} predictions successful")

        # Rank by affinity
        valid = [r for r in results if r["affinity_pred_value"] is not None]
        if valid:
            valid.sort(key=lambda r: r["affinity_pred_value"])
            print("\nRanked by affinity (lower = tighter binding):")
            for r in valid:
                print(f"  Site {r['site_id']:>2}  {r['ligand']:20s}  "
                      f"log10(IC50/uM)={r['affinity_pred_value']:.2f}  "
                      f"P(binder)={r['affinity_probability_binary']:.3f}")

    print("\nDone.")


if __name__ == "__main__":
    main()
