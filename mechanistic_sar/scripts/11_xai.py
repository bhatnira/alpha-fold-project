#!/usr/bin/env python3
"""
Step 11: XAI - Compute feature importance via permutation importance
and correlation analysis. Save to explainability/.
"""
import csv, json, os, math
from pathlib import Path
import numpy as np
from scipy import stats as sp_stats

PROJECT = Path("/cluster/home/nbhatt04/lean_pipeline")
MECHANISTIC = PROJECT / "mechanistic_sar"
RESULTS = MECHANISTIC / "results"
FEATURES = MECHANISTIC / "features"
XAI = RESULTS / "explainability"
XAI.mkdir(parents=True, exist_ok=True)


def load_features_and_labels():
    """Load all features and labels."""
    features = {}
    labels = {}

    # Chemical descriptors
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

    # Binding features
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

    # Structural features
    struct_path = FEATURES / "structural" / "structural_features.json"
    if struct_path.exists():
        with open(struct_path) as f:
            for name, s in json.load(f).items():
                # Map name to compound ID
                cid = name.replace("sar_cpd", "").replace("inactive_cpd", "")
                for prefix in ["boltz_results_sar_cpd", "boltz_results_inactive_cpd",
                               "sar_cpd", "inactive_cpd"]:
                    cid = cid.replace(prefix, "")
                    break
                if cid:
                    features.setdefault(cid, {}).update({
                        "contact_count": s.get("contact_count", 0),
                        "hbond_count": s.get("hbond_count", 0),
                        "min_distance": s.get("min_distance", 999),
                    })

    # State analysis
    state_path = FEATURES / "structural" / "state_analysis.json"
    if state_path.exists():
        with open(state_path) as f:
            for name, s in json.load(f).items():
                cid = name.replace("sar_cpd", "").replace("inactive_cpd", "")
                for prefix in ["boltz_results_sar_cpd", "boltz_results_inactive_cpd",
                               "sar_cpd", "inactive_cpd"]:
                    cid = cid.replace(prefix, "")
                    break
                if cid:
                    features.setdefault(cid, {}).update({
                        "ca_rmsd": s.get("ca_rmsd", 0),
                        "pocket_rmsd": s.get("pocket_rmsd", 0),
                        "open_state_score": s.get("open_state_score", 0),
                    })

    # Labels
    labels_path = RESULTS / "prepared_data" / "labels.json"
    if labels_path.exists():
        with open(labels_path) as f:
            labels = json.load(f)

    return features, labels


def permutation_importance(X, y, feature_names, n_repeats=10, seed=42):
    """Compute permutation importance for each feature."""
    rng = np.random.RandomState(seed)

    # Simple model: logistic-like using correlation-based scoring
    def score(X_sub, y_sub):
        if X_sub.shape[1] == 0:
            return 0.0
        # Use correlation with first principal component as proxy
        mean_vals = np.mean(X_sub, axis=0)
        std_vals = np.std(X_sub, axis=0)
        std_vals = np.where(std_vals > 0, std_vals, 1.0)
        z = (X_sub - mean_vals) / std_vals
        # Simple linear model score
        w = np.linalg.lstsq(z, y_sub, rcond=None)[0]
        y_pred = z @ w
        # AUROC-like score
        sorted_idx = np.argsort(-y_pred)
        y_sorted = y_sub[sorted_idx]
        n_pos = np.sum(y_sub == 1)
        n_neg = np.sum(y_sub == 0)
        if n_pos == 0 or n_neg == 0:
            return 0.5
        tpr = np.cumsum(y_sorted) / n_pos
        fpr = np.cumsum(1 - y_sorted) / n_neg
        tpr = np.concatenate([[0], tpr])
        fpr = np.concatenate([[0], fpr])
        return abs(float(np.trapz(tpr, fpr)))

    base_score = score(X, y)
    importances = np.zeros(X.shape[1])

    for feat_idx in range(X.shape[1]):
        scores = []
        for _ in range(n_repeats):
            X_perm = X.copy()
            X_perm[:, feat_idx] = rng.permutation(X_perm[:, feat_idx])
            s = score(X_perm, y)
            scores.append(base_score - s)
        importances[feat_idx] = float(np.mean(scores))

    # Normalize
    importances = np.maximum(importances, 0)
    total = np.sum(importances)
    if total > 0:
        importances /= total

    return importances, base_score


def correlation_analysis(X, y, feature_names):
    """Compute correlation of each feature with the label."""
    correlations = {}
    for i, name in enumerate(feature_names):
        r, p = sp_stats.pearsonr(X[:, i], y)
        correlations[name] = {
            "pearson_r": float(r) if not np.isnan(r) else 0.0,
            "p_value": float(p) if not np.isnan(p) else 1.0,
            "abs_r": abs(float(r)) if not np.isnan(r) else 0.0,
        }
    return correlations


def main():
    print("=" * 70)
    print("EXPLAINABILITY (XAI) - Mechanistic SAR Analysis")
    print("=" * 70)

    print("\n[1/4] Loading features and labels...")
    features, labels = load_features_and_labels()

    common = sorted(set(features.keys()) & set(labels.keys()))
    print(f"  Compounds: {len(common)}")

    if len(common) < 5:
        print("  ERROR: Not enough compounds. Exiting.")
        return

    # Build feature matrix
    non_fp_keys = [k for k in features[common[0]].keys()
                   if not isinstance(features[common[0]][k], (list, np.ndarray))]

    X = np.array([[features[cid].get(k, 0) for k in non_fp_keys] for cid in common])
    y = np.array([1 if labels[cid].get("is_active", False) else 0 for cid in common])

    print(f"  Feature matrix: {X.shape}")
    print(f"  Active: {np.sum(y == 1)}, Inactive: {np.sum(y == 0)}")

    print("\n[2/4] Permutation importance...")
    importances, base_score = permutation_importance(X, y, non_fp_keys, n_repeats=10)
    print(f"  Base AUROC: {base_score:.3f}")

    importance_ranking = sorted(
        zip(non_fp_keys, importances),
        key=lambda x: -x[1]
    )
    print("  Top features:")
    for name, imp in importance_ranking[:5]:
        print(f"    {name}: {imp:.4f}")

    print("\n[3/4] Correlation analysis...")
    correlations = correlation_analysis(X, y, non_fp_keys)

    # Sort by absolute correlation
    corr_ranking = sorted(
        correlations.items(),
        key=lambda x: -x[1]["abs_r"]
    )
    print("  Top correlated features:")
    for name, corr in corr_ranking[:5]:
        print(f"    {name}: r={corr['pearson_r']:.3f}, p={corr['p_value']:.4f}")

    print("\n[4/4] Saving explainability results...")
    xai_results = {
        "permutation_importance": {
            name: float(imp) for name, imp in importance_ranking
        },
        "base_auroc": base_score,
        "correlations": correlations,
        "importance_ranking": [
            {"feature": name, "importance": float(imp)}
            for name, imp in importance_ranking
        ],
        "correlation_ranking": [
            {"feature": name, **corr}
            for name, corr in corr_ranking
        ],
    }

    with open(XAI / "xai_results.json", "w") as f:
        json.dump(xai_results, f, indent=2)

    # Save CSV
    with open(XAI / "feature_importance.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["feature", "importance", "pearson_r", "p_value"])
        writer.writeheader()
        for name, imp in importance_ranking:
            corr = correlations.get(name, {})
            writer.writerow({
                "feature": name,
                "importance": imp,
                "pearson_r": corr.get("pearson_r", 0),
                "p_value": corr.get("p_value", 1),
            })

    print(f"\nResults saved to: {XAI}")
    print("  xai_results.json")
    print("  feature_importance.csv")


if __name__ == "__main__":
    main()
