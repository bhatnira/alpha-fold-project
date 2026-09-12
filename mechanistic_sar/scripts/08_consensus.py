#!/usr/bin/env python3
"""
Step 8: Consensus - Compute pose consensus, structural consensus,
prediction disagreement across multiple models per compound.
Save to features/.
"""
import csv, json, os, math
from pathlib import Path
import numpy as np

PROJECT = Path("/cluster/home/nbhatt04/lean_pipeline")
MECHANISTIC = PROJECT / "mechanistic_sar"
RESULTS = MECHANISTIC / "results"
FEATURES = MECHANISTIC / "features"
CONSENSUS = FEATURES / "consensus"
CONSENSUS.mkdir(parents=True, exist_ok=True)

CACHE = RESULTS / ".cache"
CACHE.mkdir(parents=True, exist_ok=True)


def _load_cache(name):
    p = CACHE / f"{name}.json"
    if p.exists():
        with open(p) as f:
            return json.load(f)
    return None


def _save_cache(name, data):
    with open(CACHE / f"{name}.json", "w") as f:
        json.dump(data, f)


def parse_ca_atoms(pdb_path):
    """Extract C-alpha coordinates from PDB."""
    ca_atoms = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM") and "CA" in line[12:16]:
                try:
                    chain = line[21].strip()
                    res_seq = int(line[22:26].strip())
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    ca_atoms.append({
                        "chain": chain,
                        "res_seq": res_seq,
                        "x": x, "y": y, "z": z,
                    })
                except (ValueError, IndexError):
                    continue
    return ca_atoms


def parse_ligand_atoms(pdb_path):
    """Extract ligand (HETATM) atoms from PDB."""
    atoms = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("HETATM"):
                try:
                    name = line[12:16].strip()
                    res_name = line[17:20].strip()
                    chain = line[21].strip()
                    res_seq = int(line[22:26].strip())
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    if res_name not in ("HOH", "WAT", "SO4", "PO4", "GOL"):
                        atoms.append({
                            "name": name,
                            "res_name": res_name,
                            "chain": chain,
                            "res_seq": res_seq,
                            "x": x, "y": y, "z": z,
                        })
                except (ValueError, IndexError):
                    continue
    return atoms


def compute_structural_rmsd(atoms1, atoms2):
    """Compute RMSD between two sets of atoms (matched by index)."""
    if len(atoms1) != len(atoms2) or len(atoms1) == 0:
        return float("inf")
    coords1 = np.array([[a["x"], a["y"], a["z"]] for a in atoms1])
    coords2 = np.array([[a["x"], a["y"], a["z"]] for a in atoms2])
    diff = coords1 - coords2
    return float(np.sqrt(np.mean(np.sum(diff * diff, axis=1))))


def compute_ligand_rmsd(lig1, lig2):
    """Compute RMSD between ligand atoms (matched by atom name)."""
    if not lig1 or not lig2:
        return float("inf")

    names1 = [a["name"] for a in lig1]
    names2 = [a["name"] for a in lig2]
    common = sorted(set(names1) & set(names2))

    if len(common) < 3:
        return float("inf")

    coords1 = []
    coords2 = []
    for name in common:
        for a in lig1:
            if a["name"] == name:
                coords1.append([a["x"], a["y"], a["z"]])
                break
        for a in lig2:
            if a["name"] == name:
                coords2.append([a["x"], a["y"], a["z"]])
                break

    coords1 = np.array(coords1)
    coords2 = np.array(coords2)
    diff = coords1 - coords2
    return float(np.sqrt(np.mean(np.sum(diff * diff, axis=1))))


def group_models_by_compound():
    """
    Group PDB files by compound, so we can compare multiple models.
    Returns dict: compound_key -> list of PDB paths
    """
    groups = {}
    for root_dir in [
        PROJECT / "phase08_boltz_affinity" / "boltz2_sar_outputs",
        PROJECT / "phase08_boltz_affinity" / "boltz2_inactive_outputs",
    ]:
        if not root_dir.exists():
            continue
        for p in root_dir.glob("*_model_*.pdb"):
            name = p.stem
            # Extract compound key (remove _model_N suffix)
            parts = name.rsplit("_model_", 1)
            if len(parts) == 2:
                compound_key = parts[0]
            else:
                compound_key = name
            if compound_key not in groups:
                groups[compound_key] = []
            groups[compound_key].append(p)

    return groups


def load_binding_predictions():
    """Load binding energy predictions from docking CSV."""
    csv_path = PROJECT / "09_docking" / "full_docking_results.csv"
    if not csv_path.exists():
        return {}

    predictions = {}
    with open(csv_path) as f:
        for row in csv.DictReader(f):
            if row.get("stoichiometry") != "2to3":
                continue
            cid = row["compound_id"]
            site = row.get("site_id", "")
            energy = float(row.get("binding_energy", 0))
            if cid not in predictions:
                predictions[cid] = {}
            predictions[cid][site] = energy

    return predictions


def main():
    print("=" * 70)
    print("CONSENSUS ANALYSIS - Mechanistic SAR Analysis")
    print("=" * 70)

    print("\n[1/4] Grouping models by compound...")
    groups = group_models_by_compound()
    print(f"  Compounds with multiple models: {sum(1 for v in groups.values() if len(v) > 1)}")
    print(f"  Total compounds: {len(groups)}")

    print("\n[2/4] Loading binding predictions...")
    predictions = load_binding_predictions()
    print(f"  Compounds with docking: {len(predictions)}")

    print("\n[3/4] Computing consensus features...")
    consensus_features = {}

    for compound_key, pdb_paths in sorted(groups.items()):
        try:
            n_models = len(pdb_paths)

            # Parse all models
            all_ca = []
            all_lig = []
            for p in sorted(pdb_paths):
                ca = parse_ca_atoms(str(p))
                lig = parse_ligand_atoms(str(p))
                all_ca.append(ca)
                all_lig.append(lig)

            # Structural consensus: mean pairwise RMSD between models
            if n_models > 1:
                structural_rmsds = []
                for i in range(n_models):
                    for j in range(i + 1, n_models):
                        # Match by (chain, res_seq)
                        keys_i = {(a["chain"], a["res_seq"]): a for a in all_ca[i]}
                        keys_j = {(a["chain"], a["res_seq"]): a for a in all_ca[j]}
                        common_keys = sorted(set(keys_i.keys()) & set(keys_j.keys()))
                        if len(common_keys) > 10:
                            matched_i = [keys_i[k] for k in common_keys]
                            matched_j = [keys_j[k] for k in common_keys]
                            rmsd = compute_structural_rmsd(matched_i, matched_j)
                            structural_rmsds.append(rmsd)

                mean_structural_rmsd = float(np.mean(structural_rmsds)) if structural_rmsds else 0.0
                std_structural_rmsd = float(np.std(structural_rmsds)) if structural_rmsds else 0.0
                max_structural_rmsd = float(np.max(structural_rmsds)) if structural_rmsds else 0.0
            else:
                mean_structural_rmsd = 0.0
                std_structural_rmsd = 0.0
                max_structural_rmsd = 0.0

            # Pose consensus: ligand RMSD between models
            if n_models > 1 and any(len(l) > 0 for l in all_lig):
                ligand_rmsds = []
                for i in range(n_models):
                    for j in range(i + 1, n_models):
                        if all_lig[i] and all_lig[j]:
                            rmsd = compute_ligand_rmsd(all_lig[i], all_lig[j])
                            if rmsd < float("inf"):
                                ligand_rmsds.append(rmsd)

                mean_ligand_rmsd = float(np.mean(ligand_rmsds)) if ligand_rmsds else float("inf")
                pose_consensus = 1.0 / (1.0 + mean_ligand_rmsd) if mean_ligand_rmsd < float("inf") else 0.0
            else:
                mean_ligand_rmsd = 0.0
                pose_consensus = 1.0

            # Prediction disagreement across docking sites
            compound_id = compound_key.replace("sar_cpd", "").replace("inactive_cpd", "")
            # Try to extract numeric ID
            for prefix in ["sar_cpd", "inactive_cpd", "boltz_results_sar_cpd"]:
                compound_id = compound_key.replace(prefix, "")
                break

            site_energies = predictions.get(compound_id, {})
            if site_energies:
                energies = list(site_energies.values())
                prediction_disagreement = float(np.std(energies))
                mean_prediction = float(np.mean(energies))
            else:
                prediction_disagreement = 0.0
                mean_prediction = 0.0

            # Consensus score: combination of structural and pose agreement
            structural_consensus = 1.0 / (1.0 + mean_structural_rmsd) if mean_structural_rmsd > 0 else 1.0
            consensus_score = 0.5 * structural_consensus + 0.5 * pose_consensus

            consensus_features[compound_key] = {
                "n_models": n_models,
                "structural_rmsd_mean": mean_structural_rmsd,
                "structural_rmsd_std": std_structural_rmsd,
                "structural_rmsd_max": max_structural_rmsd,
                "ligand_rmsd_mean": mean_ligand_rmsd,
                "pose_consensus": pose_consensus,
                "structural_consensus": structural_consensus,
                "consensus_score": consensus_score,
                "prediction_disagreement": prediction_disagreement,
                "mean_prediction": mean_prediction,
            }

            print(f"  {compound_key}: {n_models} models, "
                  f"struct_rmsd={mean_structural_rmsd:.3f}, "
                  f"pose_cons={pose_consensus:.3f}")

        except Exception as e:
            print(f"  WARNING: Error processing {compound_key}: {e}")

    print(f"\n  Processed: {len(consensus_features)} compounds")

    print("\n[4/4] Saving consensus features...")
    # Save JSON
    with open(CONSENSUS / "consensus_features.json", "w") as f:
        json.dump(consensus_features, f, indent=2)

    # Save CSV
    if consensus_features:
        keys = sorted(consensus_features.keys())
        csv_keys = [k for k in consensus_features[keys[0]].keys()]
        with open(CONSENSUS / "consensus_features.csv", "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=csv_keys)
            writer.writeheader()
            for key in keys:
                writer.writerow(consensus_features[key])

    # Summary
    if consensus_features:
        scores = [v["consensus_score"] for v in consensus_features.values()]
        print(f"\n  Consensus score range: [{min(scores):.3f}, {max(scores):.3f}]")
        multi = sum(1 for v in consensus_features.values() if v["n_models"] > 1)
        print(f"  Multi-model compounds: {multi}")

    print(f"\nResults saved to: {CONSENSUS}")
    print("  consensus_features.json")
    print("  consensus_features.csv")


if __name__ == "__main__":
    main()
