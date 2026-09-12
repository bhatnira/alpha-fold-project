#!/usr/bin/env python3
"""
Step 13: Ablation Study - Remove feature families one at a time,
compare performance. Save to ablations/.
"""
import csv, json, os, math
from pathlib import Path
import numpy as np
from scipy import stats as sp_stats

PROJECT = Path("/cluster/home/nbhatt04/lean_pipeline")
MECHANISTIC = PROJECT / "mechanistic_sar"
RESULTS = MECHANISTIC / "results"
FEATURES = MECHANISTIC / "features"
ABLATIONS = RESULTS / "ablations"
ABLATIONS.mkdir(parents=True, exist_ok=True)


def load_features_and_labels():
    """Load features and labels organized by family."""
    feature_families = {
        "chemical": {},
        "binding": {},
    }
    labels = {}

    # Chemical
    desc_path = FEATURES / "chemical" / "chemical_descriptors.json"
    if desc_path.exists():
        with open(desc_path) as f:
            for cid, d in json.load(f).items():
                feature_families["chemical"][cid] = {
                    "mw": d.get("molecular_weight", 0),
                    "logp": d.get("logP", 0),
                    "hbd": d.get("hbd", 0),
                    "hba": d.get("hba", 0),
                    "n_rings": d.get("n_rings", 0),
                    "tpsa": d.get("tpsa", 0),
                }

    # Binding
    bind_path = FEATURES / "binding" / "binding_features.json"
    if bind_path.exists():
        with open(bind_path) as f:
            for cid, b in json.load(f).items():
                feature_families["binding"][cid] = {
                    "docking_score": b.get("docking_score", 0),
                    "boltz2_affinity": b.get("boltz2_affinity", 0),
                    "boltz2_probability": b.get("boltz2_probability", 0),
                    "confidence_score": b.get("confidence_score", 0),
                    "iptm": b.get("iptm", 0),
                }

    # Labels
    labels_path = RESULTS / "prepared_data" / "labels.json"
    if labels_path.exists():
        with open(labels_path) as f:
            labels = json.load(f)

    return feature_families, labels


def simple_cross_validate(X, y, n_folds=5, seed=42):
    """Simple cross-validation with correlation-based scoring."""
    rng = np.random.RandomState(seed)
    n_samples = len(y)

    pos_idx = np.where(y == 1)[0]
    neg_idx = np.where(y == 0)[0]
    rng.shuffle(pos_idx)
    rng.shuffle(neg_idx)

    n_folds_use = min(n_folds, len(pos_idx), len(neg_idx))
    if n_folds_use < 2:
        return {"auroc": 0.5, "spearman": 0.0, "rmse": 1.0}

    pos_folds = np.array_split(pos_idx, n_folds_use)
    neg_folds = np.array_split(neg_idx, n_folds_use)

    results = {"auroc": [], "spearman": [], "rmse": []}

    for fold in range(n_folds_use):
        test_idx = np.concatenate([pos_folds[fold], neg_folds[fold]])
        train_idx = np.concatenate([
            np.concatenate([pos_folds[i] for i in range(n_folds_use) if i != fold]),
            np.concatenate([neg_folds[i] for i in range(n_folds_use) if i != fold])
        ])

        X_train, X_test = X[train_idx], X[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]

        # Simple linear model
        mean_vals = np.mean(X_train, axis=0)
        std_vals = np.std(X_train, axis=0)
        std_vals = np.where(std_vals > 0, std_vals, 1.0)
        z_train = (X_train - mean_vals) / std_vals
        z_test = (X_test - mean_vals) / std_vals

        w = np.linalg.lstsq(z_train, y_train, rcond=None)[0]
        y_pred = z_test @ w

        # AUROC
        sorted_idx = np.argsort(-y_pred)
        y_sorted = y_test[sorted_idx]
        n_pos = np.sum(y_test == 1)
        n_neg = np.sum(y_test == 0)
        if n_pos > 0 and n_neg > 0:
            tpr = np.concatenate([[0], np.cumsum(y_sorted) / n_pos])
            fpr = np.concatenate([[0], np.cumsum(1 - y_sorted) / n_neg])
            auroc = abs(float(np.trapz(tpr, fpr)))
        else:
            auroc = 0.5

        # Spearman
        rho, _ = sp_stats.spearmanr(y_test, y_pred) if len(y_test) > 2 else (0.0, 1.0)
        rho = float(rho) if not np.isnan(rho) else 0.0

        # RMSE
        rmse = float(np.sqrt(np.mean((y_test - y_pred) ** 2)))

        results["auroc"].append(auroc)
        results["spearman"].append(rho)
        results["rmse"].append(rmse)

    return {k: float(np.mean(v)) for k, v in results.items()}


def main():
    print("=" * 70)
    print("ABLATION STUDY - Mechanistic SAR Analysis")
    print("=" * 70)

    print("\n[1/4] Loading features and labels...")
    feature_families, labels = load_features_and_labels()

    # Get common compounds
    all_cids = set()
    for family_data in feature_families.values():
        all_cids.update(family_data.keys())
    common = sorted(all_cids & set(labels.keys()))
    print(f"  Compounds: {len(common)}")

    if len(common) < 5:
        print("  ERROR: Not enough compounds. Exiting.")
        return

    y = np.array([1 if labels[cid].get("is_active", False) else 0 for cid in common])
    print(f"  Active: {np.sum(y == 1)}, Inactive: {np.sum(y == 0)}")

    # Build combined feature matrix
    family_names = sorted(feature_families.keys())
    family_matrices = {}

    for family in family_names:
        family_data = feature_families[family]
        feat_keys = sorted(next(iter(family_data.values())).keys()) if family_data else []
        X = np.zeros((len(common), len(feat_keys)))
        for i, cid in enumerate(common):
            for j, k in enumerate(feat_keys):
                X[i, j] = family_data.get(cid, {}).get(k, 0)
        family_matrices[family] = (feat_keys, X)

    # Full model
    X_full = np.hstack([fm[1] for fm in family_matrices.values()])
    print(f"\n  Full feature matrix: {X_full.shape}")

    print("\n[2/4] Baseline: full model performance...")
    baseline = simple_cross_validate(X_full, y, n_folds=5)
    print(f"  Full model AUROC: {baseline['auroc']:.3f}")
    print(f"  Full model Spearman: {baseline['spearman']:.3f}")

    print("\n[3/4] Running ablations...")
    ablation_results = {"full_model": baseline}

    for family in family_names:
        # Remove this family
        X_ablated = np.hstack([
            family_matrices[f][1] for f in family_names if f != family
        ])
        if X_ablated.shape[1] == 0:
            continue

        perf = simple_cross_validate(X_ablated, y, n_folds=5)
        delta_auroc = baseline["auroc"] - perf["auroc"]
        delta_spearman = baseline["spearman"] - perf["spearman"]

        ablation_results[f"no_{family}"] = {
            **perf,
            "delta_auroc": delta_auroc,
            "delta_spearman": delta_spearman,
            "features_removed": family,
            "n_features_removed": family_matrices[family][1].shape[1],
            "n_features_remaining": X_ablated.shape[1],
        }

        direction = "WORSE" if delta_auroc > 0 else "BETTER"
        print(f"  No {family}: AUROC={perf['auroc']:.3f} (delta={delta_auroc:+.3f}) "
              f"[{direction}]")

    # Single-family models
    for family in family_names:
        feat_keys, X_family = family_matrices[family]
        if X_family.shape[1] == 0:
            continue
        perf = simple_cross_validate(X_family, y, n_folds=5)
        ablation_results[f"only_{family}"] = perf
        print(f"  Only {family}: AUROC={perf['auroc']:.3f}")

    print("\n[4/4] Saving ablation results...")
    with open(ABLATIONS / "ablation_results.json", "w") as f:
        json.dump(ablation_results, f, indent=2)

    # Save CSV
    with open(ABLATIONS / "ablation_comparison.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "model", "auroc", "spearman", "rmse",
            "delta_auroc", "delta_spearman"
        ])
        writer.writeheader()
        for name, metrics in ablation_results.items():
            writer.writerow({
                "model": name,
                "auroc": metrics.get("auroc", 0),
                "spearman": metrics.get("spearman", 0),
                "rmse": metrics.get("rmse", 0),
                "delta_auroc": metrics.get("delta_auroc", 0),
                "delta_spearman": metrics.get("delta_spearman", 0),
            })

    print(f"\nResults saved to: {ABLATIONS}")
    print("  ablation_results.json")
    print("  ablation_comparison.csv")


if __name__ == "__main__":
    main()
