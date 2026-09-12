#!/usr/bin/env python3
"""
Step 3: Binding Features - Extract binding features from docking CSV and
Boltz-2 affinity CSV. Save to features/binding/.
"""
import csv, json, os
from pathlib import Path
import numpy as np

PROJECT = Path("/cluster/home/nbhatt04/lean_pipeline")
MECHANISTIC = PROJECT / "mechanistic_sar"
RESULTS = MECHANISTIC / "results"
FEATURES = MECHANISTIC / "features"
BINDING = FEATURES / "binding"
BINDING.mkdir(parents=True, exist_ok=True)

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


def load_docking_scores():
    """Load docking scores from full_docking_results.csv."""
    cache = _load_cache("docking_scores")
    if cache is not None:
        return cache

    csv_path = PROJECT / "09_docking" / "full_docking_results.csv"
    results = {}
    if not csv_path.exists():
        print(f"  WARNING: {csv_path} not found")
        return results

    with open(csv_path) as f:
        for row in csv.DictReader(f):
            if row.get("stoichiometry") != "2to3":
                continue
            cid = row["compound_id"]
            if cid not in results:
                results[cid] = {
                    "compound_id": cid,
                    "smiles": row.get("smiles", ""),
                    "is_active": row.get("is_active") == "True",
                    "site_energies": [],
                    "best_energy": None,
                    "mean_energy": None,
                    "n_sites": 0,
                }
            try:
                e = float(row.get("binding_energy", 0))
                results[cid]["site_energies"].append(e)
                results[cid]["n_sites"] += 1
            except (ValueError, TypeError):
                pass

    # Compute summary statistics
    for cid, d in results.items():
        if d["site_energies"]:
            d["best_energy"] = float(np.min(d["site_energies"]))
            d["mean_energy"] = float(np.mean(d["site_energies"]))
            d["std_energy"] = float(np.std(d["site_energies"])) if len(d["site_energies"]) > 1 else 0.0
            d["median_energy"] = float(np.median(d["site_energies"]))
        else:
            d["best_energy"] = 0.0
            d["mean_energy"] = 0.0
            d["std_energy"] = 0.0
            d["median_energy"] = 0.0

    _save_cache("docking_scores", results)
    return results


def load_boltz2_features():
    """Load Boltz-2 affinity and confidence features."""
    cache = _load_cache("boltz2_features")
    if cache is not None:
        return cache

    # Load affinity predictions
    aff_path = PROJECT / "phase08_boltz_affinity" / "boltz2_full_sar_results.csv"
    aff_data = {}
    if aff_path.exists():
        with open(aff_path) as f:
            for row in csv.DictReader(f):
                cid = row["compound_id"]
                aff_data[cid] = {
                    "boltz2_affinity": float(row.get("affinity_value", 0) or 0),
                    "boltz2_probability": float(row.get("affinity_probability", 0) or 0),
                    "confidence_score": float(row.get("confidence_score", 0) or 0),
                    "iptm": float(row.get("iptm", 0) or 0),
                }

    # Load confidence scores (per-model)
    conf_path = PROJECT / "11_boltz2_analysis" / "boltz2_all_confidences.csv"
    conf_data = {}
    if conf_path.exists():
        with open(conf_path) as f:
            for row in csv.DictReader(f):
                name = row.get("name", "")
                # Extract compound info from name if possible
                conf_data[name] = {
                    "ptm": float(row.get("ptm", 0) or 0),
                    "iptm_conf": float(row.get("iptm", 0) or 0),
                    "ligand_iptm": float(row.get("ligand_iptm", 0) or 0),
                    "protein_iptm": float(row.get("protein_iptm", 0) or 0),
                    "complex_plddt": float(row.get("complex_plddt", 0) or 0),
                }

    # Merge into per-compound features
    results = {}
    for cid, aff in aff_data.items():
        results[cid] = aff.copy()

    # Compute mean confidence from per-model data if available
    # Group conf_data by compound-like prefix
    # For now, use the affinity confidence as the primary feature
    if not results:
        print("  WARNING: No Boltz-2 data found")

    _save_cache("boltz2_features", results)
    return results


def compute_binding_features(docking, boltz2):
    """Merge all binding features into a unified feature vector."""
    all_cids = sorted(set(list(docking.keys()) + list(boltz2.keys())))
    features = {}

    for cid in all_cids:
        d = docking.get(cid, {})
        b = boltz2.get(cid, {})
        feat = {
            "compound_id": cid,
            "docking_score": d.get("best_energy", 0.0) or 0.0,
            "docking_mean": d.get("mean_energy", 0.0) or 0.0,
            "docking_std": d.get("std_energy", 0.0) or 0.0,
            "boltz2_affinity": b.get("boltz2_affinity", 0.0) or 0.0,
            "boltz2_probability": b.get("boltz2_probability", 0.0) or 0.0,
            "confidence_score": b.get("confidence_score", 0.0) or 0.0,
            "iptm": b.get("iptm", 0.0) or 0.0,
        }
        # Derived features
        feat["docking_boltz_correlation"] = (
            feat["docking_score"] * feat["boltz2_affinity"]
        )
        features[cid] = feat

    return features


def main():
    print("=" * 70)
    print("BINDING FEATURES - Mechanistic SAR Analysis")
    print("=" * 70)

    print("\n[1/4] Loading docking scores...")
    docking = load_docking_scores()
    print(f"  Compounds with docking: {len(docking)}")
    if docking:
        energies = [d.get("best_energy", 0) for d in docking.values()]
        print(f"  Energy range: [{min(energies):.2f}, {max(energies):.2f}] kcal/mol")

    print("\n[2/4] Loading Boltz-2 features...")
    boltz2 = load_boltz2_features()
    print(f"  Compounds with Boltz-2: {len(boltz2)}")

    print("\n[3/4] Computing binding features...")
    features = compute_binding_features(docking, boltz2)
    print(f"  Total features computed: {len(features)}")

    print("\n[4/4] Saving binding features...")
    # Save as JSON
    with open(BINDING / "binding_features.json", "w") as f:
        json.dump(features, f, indent=2)

    # Save as CSV
    if features:
        cids = sorted(features.keys())
        fieldnames = ["compound_id"] + [k for k in features[cids[0]].keys() if k != "compound_id"]
        with open(BINDING / "binding_features.csv", "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for cid in cids:
                row = {k: v for k, v in features[cid].items()}
                writer.writerow(row)

    # Save feature matrix as numpy
    if features:
        cids = sorted(features.keys())
        feat_keys = [k for k in features[cids[0]].keys() if k != "compound_id"]
        matrix = np.array([[features[cid][k] for k in feat_keys] for cid in cids])
        np.save(BINDING / "binding_features.npy", matrix)
        with open(BINDING / "feature_names.json", "w") as f:
            json.dump(feat_keys, f)

    print(f"\nResults saved to: {BINDING}")
    print("  binding_features.json")
    print("  binding_features.csv")
    print("  binding_features.npy")
    print("  feature_names.json")


if __name__ == "__main__":
    main()
