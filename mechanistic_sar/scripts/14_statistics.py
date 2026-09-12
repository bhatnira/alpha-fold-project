#!/usr/bin/env python3
"""
Step 14: Statistics - Statistical tests: paired bootstrap, permutation tests,
confidence intervals. Save to statistics/.
"""
import csv, json, os, math
from pathlib import Path
import numpy as np
from scipy import stats as sp_stats

PROJECT = Path("/cluster/home/nbhatt04/lean_pipeline")
MECHANISTIC = PROJECT / "mechanistic_sar"
RESULTS = MECHANISTIC / "results"
STATS_DIR = RESULTS / "statistics"
STATS_DIR.mkdir(parents=True, exist_ok=True)


def load_results():
    """Load model results and ablation results."""
    results = {}

    # Model results
    model_path = RESULTS / "models" / "model_results.json"
    if model_path.exists():
        with open(model_path) as f:
            results["models"] = json.load(f)

    # Ablation results
    abl_path = RESULTS / "ablations" / "ablation_results.json"
    if abl_path.exists():
        with open(abl_path) as f:
            results["ablations"] = json.load(f)

    # XAI results
    xai_path = RESULTS / "explainability" / "xai_results.json"
    if xai_path.exists():
        with open(xai_path) as f:
            results["xai"] = json.load(f)

    # Stability results
    stab_path = RESULTS / "stability" / "stability_results.json"
    if stab_path.exists():
        with open(stab_path) as f:
            results["stability"] = json.load(f)

    return results


def paired_bootstrap_test(scores1, scores2, n_bootstrap=10000, seed=42):
    """
    Paired bootstrap test for comparing two models.
    Tests H0: mean(scores1) == mean(scores2)
    """
    rng = np.random.RandomState(seed)
    scores1 = np.array(scores1)
    scores2 = np.array(scores2)

    if len(scores1) != len(scores2) or len(scores1) < 3:
        return {"p_value": 1.0, "mean_diff": 0.0, "ci_lower": 0.0, "ci_upper": 0.0}

    observed_diff = np.mean(scores1 - scores2)
    n = len(scores1)

    bootstrap_diffs = np.zeros(n_bootstrap)
    for i in range(n_bootstrap):
        idx = rng.choice(n, size=n, replace=True)
        bootstrap_diffs[i] = np.mean(scores1[idx] - scores2[idx])

    # Two-tailed p-value
    p_value = np.mean(np.abs(bootstrap_diffs) >= np.abs(observed_diff))

    ci_lower = np.percentile(bootstrap_diffs, 2.5)
    ci_upper = np.percentile(bootstrap_diffs, 97.5)

    return {
        "p_value": float(p_value),
        "mean_diff": float(observed_diff),
        "ci_lower": float(ci_lower),
        "ci_upper": float(ci_upper),
        "significant_005": p_value < 0.05,
        "significant_001": p_value < 0.01,
    }


def permutation_test(x, y, n_permutations=10000, seed=42):
    """
    Permutation test for association between feature x and label y.
    """
    rng = np.random.RandomState(seed)
    x = np.array(x, dtype=float)
    y = np.array(y, dtype=float)

    observed_corr = abs(float(np.corrcoef(x, y)[0, 1])) if len(x) > 2 else 0.0

    perm_corrs = np.zeros(n_permutations)
    for i in range(n_permutations):
        y_perm = rng.permutation(y)
        perm_corrs[i] = abs(float(np.corrcoef(x, y_perm)[0, 1])) if len(x) > 2 else 0.0

    p_value = np.mean(perm_corrs >= observed_corr)

    return {
        "observed_correlation": observed_corr,
        "p_value": float(p_value),
        "mean_perm_correlation": float(np.mean(perm_corrs)),
        "std_perm_correlation": float(np.std(perm_corrs)),
        "significant_005": p_value < 0.05,
        "significant_001": p_value < 0.01,
    }


def compute_confidence_intervals(data, confidence=0.95, method="bootstrap", n_bootstrap=1000, seed=42):
    """Compute confidence intervals for mean."""
    data = np.array(data, dtype=float)
    if len(data) < 3:
        return {"mean": float(np.mean(data)), "ci_lower": float(np.mean(data)),
                "ci_upper": float(np.mean(data))}

    mean_val = float(np.mean(data))

    if method == "bootstrap":
        rng = np.random.RandomState(seed)
        boot_means = np.zeros(n_bootstrap)
        for i in range(n_bootstrap):
            sample = rng.choice(data, size=len(data), replace=True)
            boot_means[i] = np.mean(sample)
        alpha = (1 - confidence) / 2
        ci_lower = float(np.percentile(boot_means, alpha * 100))
        ci_upper = float(np.percentile(boot_means, (1 - alpha) * 100))
    else:
        se = float(np.std(data) / np.sqrt(len(data)))
        t_val = sp_stats.t.ppf((1 + confidence) / 2, len(data) - 1)
        ci_lower = mean_val - t_val * se
        ci_upper = mean_val + t_val * se

    return {
        "mean": mean_val,
        "ci_lower": ci_lower,
        "ci_upper": ci_upper,
        "std": float(np.std(data)),
        "n": len(data),
    }


def effect_size_cohens_d(group1, group2):
    """Compute Cohen's d effect size."""
    g1, g2 = np.array(group1, dtype=float), np.array(group2, dtype=float)
    n1, n2 = len(g1), len(g2)
    if n1 < 2 or n2 < 2:
        return 0.0
    pooled_std = math.sqrt(((n1 - 1) * np.var(g1, ddof=1) + (n2 - 1) * np.var(g2, ddof=1)) / (n1 + n2 - 2))
    if pooled_std == 0:
        return 0.0
    return float((np.mean(g1) - np.mean(g2)) / pooled_std)


def main():
    print("=" * 70)
    print("STATISTICAL TESTS - Mechanistic SAR Analysis")
    print("=" * 70)

    print("\n[1/5] Loading results...")
    results = load_results()

    print("\n[2/5] Paired bootstrap tests...")
    bootstrap_results = {}

    # Compare basic vs augmented features
    if "models" in results:
        basic = results["models"].get("basic_features", {}).get("cv_results", {})
        augmented = results["models"].get("augmented_features", {}).get("cv_results", {})

        basic_folds = results["models"].get("basic_features", {}).get("fold_results", {})
        aug_folds = results["models"].get("augmented_features", {}).get("fold_results", {})

        if basic_folds.get("auroc") and aug_folds.get("auroc"):
            test = paired_bootstrap_test(
                basic_folds["auroc"], aug_folds["auroc"], n_bootstrap=10000
            )
            bootstrap_results["basic_vs_augmented_auroc"] = test
            print(f"  Basic vs Augmented AUROC: diff={test['mean_diff']:.3f}, "
                  f"p={test['p_value']:.4f}")

    print("\n[3/5] Permutation tests...")
    permutation_results = {}

    # Load features for permutation test
    features_path = FEATURES_DIR / "binding" / "binding_features.json"
    labels_path = RESULTS / "prepared_data" / "labels.json"
    if features_path.exists() and labels_path.exists():
        with open(features_path) as f:
            features = json.load(f)
        with open(labels_path) as f:
            labels = json.load(f)

        common = sorted(set(features.keys()) & set(labels.keys()))
        if len(common) >= 5:
            for feat_name in ["docking_score", "boltz2_affinity", "confidence_score"]:
                x = np.array([features[cid].get(feat_name, 0) for cid in common])
                y = np.array([1 if labels[cid].get("is_active", False) else 0 for cid in common])
                result = permutation_test(x, y, n_permutations=5000)
                permutation_results[feat_name] = result
                print(f"  {feat_name}: r={result['observed_correlation']:.3f}, "
                      f"p={result['p_value']:.4f}")

    print("\n[4/5] Confidence intervals...")
    ci_results = {}

    if "ablations" in results:
        ablations = results["ablations"]
        if "full_model" in ablations:
            full_auroc = ablations["full_model"].get("auroc", 0)
            ci_results["full_model_auroc"] = compute_confidence_intervals(
                [full_auroc], confidence=0.95, method="parametric"
            )
            print(f"  Full model AUROC: {ci_results['full_model_auroc']['mean']:.3f} "
                  f"CI: [{ci_results['full_model_auroc']['ci_lower']:.3f}, "
                  f"{ci_results['full_model_auroc']['ci_upper']:.3f}]")

    # Effect sizes
    effect_sizes = {}
    if "ablations" in results:
        ablations = results["ablations"]
        full = ablations.get("full_model", {}).get("auroc", 0)
        for key, val in ablations.items():
            if key.startswith("no_"):
                family = val.get("features_removed", key)
                ablated_auroc = val.get("auroc", 0)
                effect_sizes[f"remove_{family}"] = {
                    "cohens_d": effect_size_cohens_d([full], [ablated_auroc]),
                    "delta_auroc": full - ablated_auroc,
                }

    print("\n[5/5] Saving statistical results...")
    stat_results = {
        "paired_bootstrap": bootstrap_results,
        "permutation_tests": permutation_results,
        "confidence_intervals": ci_results,
        "effect_sizes": effect_sizes,
    }

    with open(STATS_DIR / "statistical_results.json", "w") as f:
        json.dump(stat_results, f, indent=2, default=lambda x: bool(x) if hasattr(x, 'item') else x)

    # Summary CSV
    with open(STATS_DIR / "test_summary.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["test", "metric", "value", "p_value", "significant"])
        for name, res in bootstrap_results.items():
            writer.writerow(["paired_bootstrap", name, res["mean_diff"], res["p_value"],
                           res["significant_005"]])
        for name, res in permutation_results.items():
            writer.writerow(["permutation", name, res["observed_correlation"], res["p_value"],
                           res["significant_005"]])

    print(f"\nResults saved to: {STATS_DIR}")
    print("  statistical_results.json")
    print("  test_summary.csv")


if __name__ == "__main__":
    FEATURES_DIR = MECHANISTIC / "features"
    main()
