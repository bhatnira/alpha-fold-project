#!/usr/bin/env python3
"""
Step 12: Stability - Bootstrap analysis for explanation stability.
Save to stability/.
"""
import csv, json, os, math
from pathlib import Path
import numpy as np
from scipy import stats as sp_stats

PROJECT = Path("/cluster/home/nbhatt04/lean_pipeline")
MECHANISTIC = PROJECT / "mechanistic_sar"
RESULTS = MECHANISTIC / "results"
FEATURES = MECHANISTIC / "features"
STABILITY = RESULTS / "stability"
STABILITY.mkdir(parents=True, exist_ok=True)


def load_features_and_labels():
    """Load features and labels."""
    features = {}
    labels = {}

    # Chemical
    desc_path = FEATURES / "chemical" / "chemical_descriptors.json"
    if desc_path.exists():
        with open(desc_path) as f:
            for cid, d in json.load(f).items():
                features.setdefault(cid, {}).update({
                    "mw": d.get("molecular_weight", 0),
                    "logp": d.get("logP", 0),
                    "hbd": d.get("hbd", 0),
                    "hba": d.get("hba", 0),
                    "n_rings": d.get("n_rings", 0),
                    "tpsa": d.get("tpsa", 0),
                })

    # Binding
    bind_path = FEATURES / "binding" / "binding_features.json"
    if bind_path.exists():
        with open(bind_path) as f:
            for cid, b in json.load(f).items():
                features.setdefault(cid, {}).update({
                    "docking_score": b.get("docking_score", 0),
                    "boltz2_affinity": b.get("boltz2_affinity", 0),
                    "boltz2_probability": b.get("boltz2_probability", 0),
                    "confidence_score": b.get("confidence_score", 0),
                    "iptm": b.get("iptm", 0),
                })

    # Labels
    labels_path = RESULTS / "prepared_data" / "labels.json"
    if labels_path.exists():
        with open(labels_path) as f:
            labels = json.load(f)

    return features, labels


def bootstrap_permutation_importance(X, y, feature_names, n_bootstrap=100, n_repeats=5, seed=42):
    """Compute permutation importance with bootstrap resampling."""
    rng = np.random.RandomState(seed)
    n_samples = len(y)
    all_importances = np.zeros((n_bootstrap, len(feature_names)))

    for b in range(n_bootstrap):
        # Bootstrap sample
        boot_idx = rng.choice(n_samples, size=n_samples, replace=True)
        X_boot = X[boot_idx]
        y_boot = y[boot_idx]

        # Compute importance on this sample
        for feat_idx in range(len(feature_names)):
            # Simple score: correlation-based
            mean_vals = np.mean(X_boot, axis=0)
            std_vals = np.std(X_boot, axis=0)
            std_vals = np.where(std_vals > 0, std_vals, 1.0)
            z = (X_boot - mean_vals) / std_vals

            w = np.linalg.lstsq(z, y_boot, rcond=None)[0]
            base_pred = z @ w
            base_score = -np.mean((base_pred - y_boot) ** 2)

            scores = []
            for _ in range(n_repeats):
                X_perm = X_boot.copy()
                X_perm[:, feat_idx] = rng.permutation(X_perm[:, feat_idx])
                z_perm = (X_perm - mean_vals) / std_vals
                pred = z_perm @ w
                s = -np.mean((pred - y_boot) ** 2)
                scores.append(base_score - s)

            all_importances[b, feat_idx] = np.mean(scores)

    return all_importances


def compute_stability_metrics(all_importances, feature_names):
    """Compute stability metrics across bootstrap iterations."""
    n_features = all_importances.shape[1]
    stability = {}

    for i, name in enumerate(feature_names):
        values = all_importances[:, i]
        mean_imp = np.mean(values)
        std_imp = np.std(values)

        # Coefficient of variation
        cv = std_imp / mean_imp if mean_imp > 0 else float("inf")

        # Rank stability: how often is this feature in top 5?
        ranks = np.zeros_like(all_importances)
        for b in range(all_importances.shape[0]):
            ranks[b] = all_importances.shape[1] - sp_stats.rankdata(all_importances[b])
        mean_rank = np.mean(ranks[:, i])
        in_top5 = np.mean(ranks[:, i] < 5)

        stability[name] = {
            "mean_importance": float(mean_imp),
            "std_importance": float(std_imp),
            "cv": float(cv),
            "mean_rank": float(mean_rank),
            "fraction_in_top5": float(in_top5),
            "median_importance": float(np.median(values)),
            "ci_95_lower": float(np.percentile(values, 2.5)),
            "ci_95_upper": float(np.percentile(values, 97.5)),
        }

    return stability


def main():
    print("=" * 70)
    print("STABILITY ANALYSIS - Mechanistic SAR Analysis")
    print("=" * 70)

    print("\n[1/4] Loading features and labels...")
    features, labels = load_features_and_labels()

    common = sorted(set(features.keys()) & set(labels.keys()))
    print(f"  Compounds: {len(common)}")

    if len(common) < 5:
        print("  ERROR: Not enough compounds. Exiting.")
        return

    non_fp_keys = [k for k in features[common[0]].keys()
                   if not isinstance(features[common[0]][k], (list, np.ndarray))]

    X = np.array([[features[cid].get(k, 0) for k in non_fp_keys] for cid in common])
    y = np.array([1 if labels[cid].get("is_active", False) else 0 for cid in common])

    print(f"  Feature matrix: {X.shape}")

    print("\n[2/4] Running bootstrap permutation importance...")
    all_importances = bootstrap_permutation_importance(
        X, y, non_fp_keys, n_bootstrap=100, n_repeats=5, seed=42
    )
    print(f"  Bootstrap iterations: {all_importances.shape[0]}")

    print("\n[3/4] Computing stability metrics...")
    stability = compute_stability_metrics(all_importances, non_fp_keys)

    # Sort by mean importance
    ranked = sorted(stability.items(), key=lambda x: -x[1]["mean_importance"])
    print("  Top features (by mean importance):")
    for name, metrics in ranked[:5]:
        print(f"    {name}: {metrics['mean_importance']:.4f} "
              f"(CI: [{metrics['ci_95_lower']:.4f}, {metrics['ci_95_upper']:.4f}])")

    print("\n[4/4] Saving stability results...")
    results = {
        "stability_metrics": stability,
        "ranking": [
            {"feature": name, **metrics}
            for name, metrics in ranked
        ],
        "n_bootstrap": all_importances.shape[0],
        "n_features": len(non_fp_keys),
    }

    with open(STABILITY / "stability_results.json", "w") as f:
        json.dump(results, f, indent=2)

    # Save CSV
    with open(STABILITY / "stability_metrics.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "feature", "mean_importance", "std_importance", "cv",
            "mean_rank", "fraction_in_top5", "ci_95_lower", "ci_95_upper",
            "median_importance"
        ])
        writer.writeheader()
        for name, metrics in ranked:
            writer.writerow({"feature": name, **metrics})

    # Save raw bootstrap importances
    np.save(STABILITY / "bootstrap_importances.npy", all_importances)
    with open(STABILITY / "bootstrap_feature_names.json", "w") as f:
        json.dump(non_fp_keys, f)

    print(f"\nResults saved to: {STABILITY}")
    print("  stability_results.json")
    print("  stability_metrics.csv")
    print("  bootstrap_importances.npy")
    print("  bootstrap_feature_names.json")


if __name__ == "__main__":
    main()
