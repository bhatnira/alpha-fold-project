#!/usr/bin/env python3
"""
Step 10: Train Models - Train simple ML models using numpy/scipy:
implement Random Forest from scratch (decision tree ensemble), evaluate
with cross-validation, compute metrics (AUROC, Spearman, RMSE).
Save to models/.
"""
import csv, json, os, math, time
from pathlib import Path
import numpy as np
from scipy import stats as sp_stats

PROJECT = Path("/cluster/home/nbhatt04/lean_pipeline")
MECHANISTIC = PROJECT / "mechanistic_sar"
RESULTS = MECHANISTIC / "results"
MODELS = RESULTS / "models"
MODELS.mkdir(parents=True, exist_ok=True)

CACHE = RESULTS / ".cache"
CACHE.mkdir(parents=True, exist_ok=True)


# ============================================================
# Simple Decision Tree Implementation
# ============================================================

class SimpleDecisionTree:
    """Decision tree for classification and regression."""

    def __init__(self, max_depth=5, min_samples_split=5, min_samples_leaf=2):
        self.max_depth = max_depth
        self.min_samples_split = min_samples_split
        self.min_samples_leaf = min_samples_leaf
        self.tree = None

    def fit(self, X, y):
        self.tree = self._build_tree(X, y, depth=0)
        return self

    def _build_tree(self, X, y, depth):
        n_samples = len(y)
        n_classes = len(np.unique(y))

        # Stopping criteria
        if (depth >= self.max_depth or n_samples < self.min_samples_split
                or n_classes == 1):
            if len(y) > 0:
                # Leaf node: return most common class (classification) or mean (regression)
                if self._is_classification:
                    values, counts = np.unique(y, return_counts=True)
                    return {"leaf": True, "value": values[np.argmax(counts)]}
                else:
                    return {"leaf": True, "value": float(np.mean(y))}
            return {"leaf": True, "value": 0}

        # Find best split
        best_feature, best_threshold, best_score = None, None, float("inf")
        for feature_idx in range(X.shape[1]):
            thresholds = np.unique(X[:, feature_idx])
            for threshold in thresholds:
                left_mask = X[:, feature_idx] <= threshold
                right_mask = ~left_mask
                if (np.sum(left_mask) < self.min_samples_leaf or
                        np.sum(right_mask) < self.min_samples_leaf):
                    continue

                score = self._compute_split_score(y[left_mask], y[right_mask])
                if score < best_score:
                    best_score = score
                    best_feature = feature_idx
                    best_threshold = threshold

        if best_feature is None:
            if len(y) > 0:
                values, counts = np.unique(y, return_counts=True)
                return {"leaf": True, "value": values[np.argmax(counts)]}
            return {"leaf": True, "value": 0}

        left_mask = X[:, best_feature] <= best_threshold
        right_mask = ~left_mask

        return {
            "leaf": False,
            "feature": best_feature,
            "threshold": best_threshold,
            "left": self._build_tree(X[left_mask], y[left_mask], depth + 1),
            "right": self._build_tree(X[right_mask], y[right_mask], depth + 1),
        }

    def _compute_split_score(self, y_left, y_right):
        """Compute Gini impurity for classification split."""
        def gini(y):
            if len(y) == 0:
                return 0
            _, counts = np.unique(y, return_counts=True)
            probs = counts / len(y)
            return 1.0 - np.sum(probs ** 2)

        n = len(y_left) + len(y_right)
        return (len(y_left) * gini(y_left) + len(y_right) * gini(y_right)) / n

    def predict(self, X):
        return np.array([self._predict_one(x, self.tree) for x in X])

    def _predict_one(self, x, node):
        if node["leaf"]:
            return node["value"]
        if x[node["feature"]] <= node["threshold"]:
            return self._predict_one(x, node["left"])
        else:
            return self._predict_one(x, node["right"])

    def predict_proba(self, X):
        """For classification: return class probabilities."""
        preds = self.predict(X)
        return preds


# ============================================================
# Simple Random Forest Implementation
# ============================================================

class SimpleRandomForest:
    """Random Forest from scratch using numpy."""

    def __init__(self, n_estimators=100, max_depth=5, min_samples_split=5,
                 max_features="sqrt", random_state=42):
        self.n_estimators = n_estimators
        self.max_depth = max_depth
        self.min_samples_split = min_samples_split
        self.max_features = max_features
        self.random_state = random_state
        self.trees = []
        self.feature_indices = []

    def fit(self, X, y):
        rng = np.random.RandomState(self.random_state)
        n_samples, n_features = X.shape

        if self.max_features == "sqrt":
            n_feats = max(1, int(math.sqrt(n_features)))
        elif self.max_features == "log2":
            n_feats = max(1, int(math.log2(n_features)))
        else:
            n_features_use = int(self.max_features) if isinstance(self.max_features, (int, float)) else n_features
            n_feats = min(n_features_use, n_features)

        self.trees = []
        self.feature_indices = []

        for i in range(self.n_estimators):
            # Bootstrap sample
            boot_idx = rng.choice(n_samples, size=n_samples, replace=True)
            X_boot = X[boot_idx]
            y_boot = y[boot_idx]

            # Feature subspace
            feat_idx = rng.choice(n_features, size=n_feats, replace=False)
            self.feature_indices.append(feat_idx)

            # Build tree
            tree = SimpleDecisionTree(
                max_depth=self.max_depth,
                min_samples_split=self.min_samples_split,
            )
            tree._is_classification = True
            tree.fit(X_boot[:, feat_idx], y_boot)
            self.trees.append(tree)

        return self

    def predict(self, X):
        predictions = []
        for tree, feat_idx in zip(self.trees, self.feature_indices):
            pred = tree.predict(X[:, feat_idx])
            predictions.append(pred)
        predictions = np.array(predictions)
        return np.mean(predictions, axis=0)

    def predict_proba(self, X):
        return self.predict(X)

    def feature_importances(self, X, y):
        """Compute permutation importance."""
        n_features = X.shape[1]
        base_score = self._accuracy_score(y, self.predict(X))
        importances = np.zeros(n_features)

        for i in range(n_features):
            X_perm = X.copy()
            rng = np.random.RandomState(42)
            X_perm[:, i] = rng.permutation(X_perm[:, i])
            perm_score = self._accuracy_score(y, self.predict(X_perm))
            importances[i] = base_score - perm_score

        importances = np.maximum(importances, 0)
        total = np.sum(importances)
        if total > 0:
            importances /= total
        return importances

    def _accuracy_score(self, y_true, y_pred):
        y_pred_binary = (y_pred >= 0.5).astype(int)
        return float(np.mean(y_true == y_pred_binary))


# ============================================================
# Metrics
# ============================================================

def compute_auroc(y_true, y_scores):
    """Compute AUROC using numpy."""
    y_true = np.array(y_true)
    y_scores = np.array(y_scores)

    # Sort by score descending
    sorted_idx = np.argsort(-y_scores)
    y_true_sorted = y_true[sorted_idx]

    n_pos = np.sum(y_true == 1)
    n_neg = np.sum(y_true == 0)

    if n_pos == 0 or n_neg == 0:
        return 0.5

    # Compute TPR and FPR
    tpr = np.cumsum(y_true_sorted) / n_pos
    fpr = np.cumsum(1 - y_true_sorted) / n_neg

    # Add origin
    tpr = np.concatenate([[0], tpr])
    fpr = np.concatenate([[0], fpr])

    # AUROC = area under curve using trapezoidal rule
    auroc = np.trapz(tpr, fpr)
    return float(abs(auroc))


def compute_spearman(x, y):
    """Compute Spearman rank correlation."""
    x = np.array(x)
    y = np.array(y)
    if len(x) < 3:
        return 0.0
    rho, _ = sp_stats.spearmanr(x, y)
    return float(rho) if not np.isnan(rho) else 0.0


def compute_rmse(y_true, y_pred):
    """Compute RMSE."""
    return float(np.sqrt(np.mean((np.array(y_true) - np.array(y_pred)) ** 2)))


def compute_pearson(x, y):
    """Compute Pearson correlation."""
    if len(x) < 3:
        return 0.0
    r, _ = sp_stats.pearsonr(x, y)
    return float(r) if not np.isnan(r) else 0.0


# ============================================================
# Cross-validation
# ============================================================

def cross_validate(X, y, n_folds=5, seed=42):
    """Stratified k-fold cross-validation."""
    rng = np.random.RandomState(seed)
    n_samples = len(y)

    # Stratify
    pos_idx = np.where(y == 1)[0]
    neg_idx = np.where(y == 0)[0]
    rng.shuffle(pos_idx)
    rng.shuffle(neg_idx)

    pos_folds = np.array_split(pos_idx, n_folds)
    neg_folds = np.array_split(neg_idx, n_folds)

    results = {"auroc": [], "spearman": [], "rmse": [], "pearson": []}

    for fold in range(n_folds):
        test_idx = np.concatenate([pos_folds[fold], neg_folds[fold]])
        train_idx = np.concatenate([np.concatenate([pos_folds[i] for i in range(n_folds) if i != fold]),
                                     np.concatenate([neg_folds[i] for i in range(n_folds) if i != fold])])

        X_train, X_test = X[train_idx], X[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]

        # Train model
        model = SimpleRandomForest(n_estimators=50, max_depth=5, random_state=seed)
        model.fit(X_train, y_train)

        # Predict
        y_pred = model.predict(X_test)

        # Metrics
        results["auroc"].append(compute_auroc(y_test, y_pred))
        results["spearman"].append(compute_spearman(y_test, y_pred))
        results["rmse"].append(compute_rmse(y_test, y_pred))
        results["pearson"].append(compute_pearson(y_test, y_pred))

    # Average results
    avg_results = {k: float(np.mean(v)) for k, v in results.items()}
    std_results = {k: float(np.std(v)) for k, v in results.items()}

    return avg_results, std_results, results


def load_all_features():
    """Load all feature families and merge into unified matrix."""
    features = {}
    feature_names = {}

    # Chemical features
    fp_path = FEATURES_DIR / "chemical" / "fingerprints.npy"
    if fp_path.exists():
        fps = np.load(str(fp_path))
        compounds_path = FEATURES_DIR / "chemical" / "fingerprint_compounds.json"
        with open(compounds_path) as f:
            fp_cids = json.load(f)
        for i, cid in enumerate(fp_cids):
            if cid not in features:
                features[cid] = {}
            features[cid]["chemical_fp"] = fps[i]

    # Chemical descriptors
    desc_path = FEATURES_DIR / "chemical" / "chemical_descriptors.json"
    if desc_path.exists():
        with open(desc_path) as f:
            descs = json.load(f)
        for cid, d in descs.items():
            if cid not in features:
                features[cid] = {}
            features[cid]["mw"] = d.get("molecular_weight", 0)
            features[cid]["logp"] = d.get("logP", 0)
            features[cid]["hbd"] = d.get("hbd", 0)
            features[cid]["hba"] = d.get("hba", 0)
            features[cid]["n_rings"] = d.get("n_rings", 0)
            features[cid]["tpsa"] = d.get("tpsa", 0)

    # Binding features
    bind_path = FEATURES_DIR / "binding" / "binding_features.json"
    if bind_path.exists():
        with open(bind_path) as f:
            binds = json.load(f)
        for cid, b in binds.items():
            if cid not in features:
                features[cid] = {}
            for k in ["docking_score", "boltz2_affinity", "boltz2_probability",
                       "confidence_score", "iptm"]:
                features[cid][k] = b.get(k, 0)

    # Labels
    labels_path = RESULTS / "prepared_data" / "labels.json"
    labels = {}
    if labels_path.exists():
        with open(labels_path) as f:
            labels = json.load(f)

    return features, labels


def main():
    print("=" * 70)
    print("MODEL TRAINING - Mechanistic SAR Analysis")
    print("=" * 70)

    print("\n[1/5] Loading features...")
    features, labels = load_all_features()
    print(f"  Compounds with features: {len(features)}")
    print(f"  Compounds with labels: {len(labels)}")

    # Get common compounds
    common = sorted(set(features.keys()) & set(labels.keys()))
    print(f"  Common compounds: {len(common)}")

    if len(common) < 5:
        print("  ERROR: Not enough compounds for modeling. Exiting.")
        return

    print("\n[2/5] Building feature matrix...")
    # Build matrix: use non-fingerprint features first
    non_fp_keys = [k for k in features[common[0]].keys()
                   if k != "chemical_fp" and not isinstance(features[common[0]][k], np.ndarray)]

    X_basic = np.array([[features[cid].get(k, 0) for k in non_fp_keys] for cid in common])
    y = np.array([1 if labels[cid].get("is_active", False) else 0 for cid in common])

    print(f"  Basic feature matrix: {X_basic.shape}")
    print(f"  Active: {np.sum(y == 1)}, Inactive: {np.sum(y == 0)}")

    # Also build fingerprint-augmented matrix if available
    if "chemical_fp" in features[common[0]]:
        fp_matrix = np.array([features[cid]["chemical_fp"] for cid in common])
        X_augmented = np.hstack([X_basic, fp_matrix])
        print(f"  Augmented feature matrix: {X_augmented.shape}")
    else:
        X_augmented = X_basic

    print("\n[3/5] Cross-validation (basic features)...")
    avg_basic, std_basic, fold_basic = cross_validate(X_basic, y, n_folds=5, seed=42)
    print(f"  AUROC: {avg_basic['auroc']:.3f} +/- {std_basic['auroc']:.3f}")
    print(f"  Spearman: {avg_basic['spearman']:.3f} +/- {std_basic['spearman']:.3f}")
    print(f"  RMSE: {avg_basic['rmse']:.3f} +/- {std_basic['rmse']:.3f}")

    print("\n[4/5] Cross-validation (augmented features)...")
    avg_aug, std_aug, fold_aug = cross_validate(X_augmented, y, n_folds=5, seed=42)
    print(f"  AUROC: {avg_aug['auroc']:.3f} +/- {std_aug['auroc']:.3f}")
    print(f"  Spearman: {avg_aug['spearman']:.3f} +/- {std_aug['spearman']:.3f}")
    print(f"  RMSE: {avg_aug['rmse']:.3f} +/- {std_aug['rmse']:.3f}")

    print("\n[5/5] Training final model and saving...")
    # Train final model on all data
    final_model = SimpleRandomForest(n_estimators=100, max_depth=5, random_state=42)
    final_model.fit(X_augmented, y)

    # Feature importance
    importances = final_model.feature_importances(X_augmented, y)

    # Save model results
    model_results = {
        "basic_features": {
            "feature_names": non_fp_keys,
            "cv_results": avg_basic,
            "cv_std": std_basic,
            "fold_results": fold_basic,
        },
        "augmented_features": {
            "n_features": X_augmented.shape[1],
            "cv_results": avg_aug,
            "cv_std": std_aug,
            "fold_results": fold_aug,
        },
        "feature_importance": {
            non_fp_keys[i]: float(importances[i]) for i in range(min(len(non_fp_keys), len(importances)))
        },
        "n_compounds": len(common),
        "n_active": int(np.sum(y == 1)),
        "n_inactive": int(np.sum(y == 0)),
    }

    with open(MODELS / "model_results.json", "w") as f:
        json.dump(model_results, f, indent=2)

    # Save feature importance as CSV
    importance_data = []
    for i, name in enumerate(non_fp_keys):
        if i < len(importances):
            importance_data.append({"feature": name, "importance": float(importances[i])})
    importance_data.sort(key=lambda x: -x["importance"])

    with open(MODELS / "feature_importance.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["feature", "importance"])
        writer.writeheader()
        for row in importance_data:
            writer.writerow(row)

    print(f"\n  Final AUROC: {avg_aug['auroc']:.3f}")
    print(f"  Top features:")
    for row in importance_data[:5]:
        print(f"    {row['feature']}: {row['importance']:.4f}")

    print(f"\nResults saved to: {MODELS}")
    print("  model_results.json")
    print("  feature_importance.csv")


if __name__ == "__main__":
    FEATURES_DIR = MECHANISTIC / "features"
    main()
