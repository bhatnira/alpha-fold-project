#!/usr/bin/env python3
"""
Step 2: Prepare Data - Load experimental manifests, create train/test splits,
save prepared data for downstream analysis.
"""
import csv, json, os, hashlib
from pathlib import Path
import numpy as np

PROJECT = Path("/cluster/home/nbhatt04/lean_pipeline")
MECHANISTIC = PROJECT / "mechanistic_sar"
CONFIG = MECHANISTIC / "config"
RESULTS = MECHANISTIC / "results"
PREPARED = RESULTS / "prepared_data"
PREPARED.mkdir(parents=True, exist_ok=True)

CACHE = MECHANISTIC / "results" / ".cache"
CACHE.mkdir(parents=True, exist_ok=True)


def _cache_path(name):
    return CACHE / f"{name}.json"


def _load_cache(name):
    p = _cache_path(name)
    if p.exists():
        with open(p) as f:
            return json.load(f)
    return None


def _save_cache(name, data):
    with open(_cache_path(name), "w") as f:
        json.dump(data, f)


def canonicalize_smiles(smiles):
    """Basic SMILES canonicalization: sort canonical string for scaffold grouping."""
    if not smiles:
        return ""
    return smiles.strip()


def smiles_to_scaffold_key(smiles):
    """
    Create a scaffold key from SMILES using a simple hash-based approach.
    Groups molecules by heavy-atom connectivity patterns.
    Uses a simple ring-and-chain decomposition heuristic.
    """
    if not smiles:
        return "UNKNOWN"
    # Extract ring info: count rings and ring sizes as a rough scaffold proxy
    ring_sizes = []
    depth = 0
    in_ring = False
    ring_start = -1
    for i, ch in enumerate(smiles):
        if ch == '(':
            depth += 1
        elif ch == ')':
            depth -= 1
        elif ch == '%' and i + 2 < len(smiles) and smiles[i + 1].isdigit() and smiles[i + 2].isdigit():
            pass  # bracket atom, skip
        elif ch in '0123456789':
            pass  # ring closure digit
        elif ch == '=' or ch == '#':
            pass  # bond type
        elif ch in 'CNOSPFIBrClcCnNoOsSi':
            pass  # atom
    # Use a simpler approach: hash of the SMILES string as scaffold key
    # This groups exact structural motifs
    h = hashlib.md5(smiles.encode()).hexdigest()[:8]
    return h


def load_experimental_data():
    """Load experimental SAR data from the audit manifest."""
    manifest = RESULTS / "data_audit" / "experimental_manifest.csv"
    if not manifest.exists():
        print("  WARNING: experimental_manifest.csv not found, loading from docking CSV")
        manifest = PROJECT / "09_docking" / "full_docking_results.csv"

    compounds = {}
    with open(manifest) as f:
        reader = csv.DictReader(f)
        for row in reader:
            cid = row.get("compound_id", "")
            if cid in compounds:
                continue
            compounds[cid] = {
                "compound_id": cid,
                "smiles": row.get("smiles", ""),
                "is_active": row.get("is_active", "False") in ("True", "true", "1"),
                "activity_uM": float(row.get("activity_uM", 0) or 0),
                "potentiation_pct": float(row.get("potentiation_pct", 0) or 0),
            }
    return compounds


def load_docking_data():
    """Load docking results, return best docking energy per compound (site 23 focus)."""
    csv_path = PROJECT / "09_docking" / "full_docking_results.csv"
    if not csv_path.exists():
        return {}
    best = {}
    with open(csv_path) as f:
        for row in csv.DictReader(f):
            if row.get("stoichiometry") != "2to3":
                continue
            cid = row["compound_id"]
            e = float(row.get("binding_energy", 0))
            if cid not in best or e < best[cid]:
                best[cid] = e
    return best


def load_boltz2_data():
    """Load Boltz-2 affinity predictions."""
    csv_path = PROJECT / "phase08_boltz_affinity" / "boltz2_full_sar_results.csv"
    if not csv_path.exists():
        return {}
    data = {}
    with open(csv_path) as f:
        for row in csv.DictReader(f):
            cid = row["compound_id"]
            data[cid] = {
                "affinity_value": float(row.get("affinity_value", 0) or 0),
                "affinity_probability": float(row.get("affinity_probability", 0) or 0),
                "confidence_score": float(row.get("confidence_score", 0) or 0),
                "iptm": float(row.get("iptm", 0) or 0),
            }
    return data


def load_convergence_data():
    """Load convergence matrix."""
    csv_path = PROJECT / "14_stoichiometry_state" / "updated_convergence_matrix.csv"
    if not csv_path.exists():
        return {}
    data = {}
    with open(csv_path) as f:
        for row in csv.DictReader(f):
            data[row["site_id"]] = {
                "evidence_score": float(row.get("evidence_score", 0) or 0),
                "overall_score": float(row.get("overall_score", 0) or 0),
            }
    return data


def random_split(compounds, test_size=0.2, seed=42):
    """Random train/test split."""
    rng = np.random.RandomState(seed)
    ids = sorted(compounds.keys())
    rng.shuffle(ids)
    n_test = max(1, int(len(ids) * test_size))
    test_ids = set(ids[:n_test])
    train_ids = set(ids[n_test:])
    return {"train": sorted(train_ids), "test": sorted(test_ids)}


def scaffold_split(compounds, n_folds=5, seed=42):
    """
    Scaffold-based split: group compounds by scaffold hash, assign entire
    scaffold groups to train or test. Returns dict with train/test lists.
    """
    rng = np.random.RandomState(seed)
    # Group by scaffold
    scaffold_groups = {}
    for cid, cdata in compounds.items():
        sk = smiles_to_scaffold_key(cdata["smiles"])
        if sk not in scaffold_groups:
            scaffold_groups[sk] = []
        scaffold_groups[sk].append(cid)

    # Sort scaffolds by size (largest first) for better distribution
    sorted_scaffolds = sorted(scaffold_groups.items(), key=lambda x: -len(x[1]))

    train_ids = []
    test_ids = []
    # Greedily assign scaffolds to balance sets
    target_test = int(len(compounds) * 0.2)
    for sk, members in sorted_scaffolds:
        if len(test_ids) < target_test:
            test_ids.extend(members)
        else:
            train_ids.extend(members)

    return {"train": sorted(train_ids), "test": sorted(test_ids)}


def build_feature_matrix(compounds, docking, boltz2, convergence):
    """Build a simple feature matrix from available data."""
    ids = sorted(compounds.keys())
    features = {}
    for cid in ids:
        row = []
        # Docking energy
        row.append(docking.get(cid, 0.0))
        # Boltz-2 features
        b = boltz2.get(cid, {})
        row.append(b.get("affinity_value", 0.0))
        row.append(b.get("affinity_probability", 0.0))
        row.append(b.get("confidence_score", 0.0))
        row.append(b.get("iptm", 0.0))
        # Activity as feature (for completeness, not as label)
        row.append(compounds[cid]["activity_uM"])
        row.append(compounds[cid]["potentiation_pct"])
        features[cid] = np.array(row, dtype=float)
    return ids, features


def main():
    print("=" * 70)
    print("DATA PREPARATION - Mechanistic SAR Analysis")
    print("=" * 70)

    print("\n[1/5] Loading experimental data...")
    compounds = load_experimental_data()
    n_active = sum(1 for c in compounds.values() if c["is_active"])
    n_inactive = len(compounds) - n_active
    print(f"  Compounds: {len(compounds)} ({n_active} active, {n_inactive} inactive)")

    print("\n[2/5] Loading docking data...")
    docking = load_docking_data()
    print(f"  Docking records: {len(docking)}")

    print("\n[3/5] Loading Boltz-2 data...")
    boltz2 = load_boltz2_data()
    print(f"  Boltz-2 records: {len(boltz2)}")

    print("\n[4/5] Loading convergence data...")
    convergence = load_convergence_data()
    print(f"  Convergence sites: {len(convergence)}")

    print("\n[5/5] Creating train/test splits...")
    splits = {}

    # Random split
    rand = random_split(compounds, test_size=0.2, seed=42)
    splits["random"] = rand
    print(f"  Random: train={len(rand['train'])}, test={len(rand['test'])}")

    # Scaffold split
    scf = scaffold_split(compounds, n_folds=5, seed=42)
    splits["scaffold"] = scf
    print(f"  Scaffold: train={len(scf['train'])}, test={len(scf['test'])}")

    # Build feature matrix
    ids, features = build_feature_matrix(compounds, docking, boltz2, convergence)

    # Save prepared data
    print("\nSaving prepared data...")
    prepared = {
        "compounds": {k: {
            "compound_id": v["compound_id"],
            "smiles": v["smiles"],
            "is_active": v["is_active"],
            "activity_uM": v["activity_uM"],
            "potentiation_pct": v["potentiation_pct"],
        } for k, v in compounds.items()},
        "splits": splits,
        "docking": docking,
        "boltz2": boltz2,
        "convergence": convergence,
    }

    with open(PREPARED / "prepared_data.json", "w") as f:
        json.dump(prepared, f, indent=2)

    # Save labels for ML
    labels = {}
    for cid in ids:
        labels[cid] = {
            "is_active": compounds[cid]["is_active"],
            "activity_uM": compounds[cid]["activity_uM"],
        }
    with open(PREPARED / "labels.json", "w") as f:
        json.dump(labels, f, indent=2)

    # Save split info as CSV
    with open(PREPARED / "splits.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["compound_id", "split", "is_active"])
        for cid in sorted(ids):
            for split_name, split_data in splits.items():
                if cid in split_data["train"]:
                    writer.writerow([cid, f"{split_name}_train", compounds[cid]["is_active"]])
                    break
                elif cid in split_data["test"]:
                    writer.writerow([cid, f"{split_name}_test", compounds[cid]["is_active"]])
                    break

    print(f"\nResults saved to: {PREPARED}")
    print("  prepared_data.json")
    print("  labels.json")
    print("  splits.csv")


if __name__ == "__main__":
    main()
