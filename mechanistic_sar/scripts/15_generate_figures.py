#!/usr/bin/env python3
"""
Step 15: Generate Figures - Generate publication figures using matplotlib.
Save to figures/.
"""
import csv, json, os, math
from pathlib import Path
import numpy as np

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    HAS_MPL = True
except ImportError:
    HAS_MPL = False
    print("WARNING: matplotlib not available, will skip figure generation")

PROJECT = Path("/cluster/home/nbhatt04/lean_pipeline")
MECHANISTIC = PROJECT / "mechanistic_sar"
RESULTS = MECHANISTIC / "results"
FEATURES = MECHANISTIC / "features"
FIGURES = RESULTS / "figures"
FIGURES.mkdir(parents=True, exist_ok=True)


def load_all_results():
    """Load all analysis results."""
    data = {}
    for name, path in [
        ("binding_features", FEATURES / "binding" / "binding_features.json"),
        ("chemical_descriptors", FEATURES / "chemical" / "chemical_descriptors.json"),
        ("structural_features", FEATURES / "structural" / "structural_features.json"),
        ("state_analysis", FEATURES / "structural" / "state_analysis.json"),
        ("nma_features", FEATURES / "nma" / "nma_features.json"),
        ("network_features", FEATURES / "network" / "network_features.json"),
        ("consensus_features", FEATURES / "consensus" / "consensus_features.json"),
        ("model_results", RESULTS / "models" / "model_results.json"),
        ("feature_importance", RESULTS / "models" / "feature_importance.csv"),
        ("xai_results", RESULTS / "explainability" / "xai_results.json"),
        ("stability_results", RESULTS / "stability" / "stability_results.json"),
        ("ablation_results", RESULTS / "ablations" / "ablation_results.json"),
        ("statistical_results", RESULTS / "statistics" / "statistical_results.json"),
        ("labels", RESULTS / "prepared_data" / "labels.json"),
    ]:
        if path.exists():
            try:
                if path.suffix == ".json":
                    with open(path) as f:
                        data[name] = json.load(f)
                elif path.suffix == ".csv":
                    with open(path) as f:
                        reader = csv.DictReader(f)
                        data[name] = list(reader)
            except Exception:
                pass
    return data


def fig1_binding_features(data):
    """Figure 1: Binding feature distributions by activity."""
    if not HAS_MPL:
        return
    binding = data.get("binding_features", {})
    labels = data.get("labels", {})
    if not binding or not labels:
        return

    features_to_plot = ["docking_score", "boltz2_affinity", "confidence_score", "iptm"]
    available = [f for f in features_to_plot if any(f in v for v in binding.values())]
    if not available:
        return

    fig, axes = plt.subplots(1, len(available), figsize=(4 * len(available), 4))
    if len(available) == 1:
        axes = [axes]

    for ax, feat in zip(axes, available):
        active_vals = []
        inactive_vals = []
        for cid, b in binding.items():
            val = b.get(feat, None)
            if val is None:
                continue
            if labels.get(cid, {}).get("is_active", False):
                active_vals.append(float(val))
            else:
                inactive_vals.append(float(val))

        if active_vals:
            ax.hist(active_vals, bins=10, alpha=0.6, label="Active", color="steelblue")
        if inactive_vals:
            ax.hist(inactive_vals, bins=10, alpha=0.6, label="Inactive", color="coral")
        ax.set_xlabel(feat)
        ax.set_ylabel("Count")
        ax.legend(fontsize=8)

    plt.tight_layout()
    plt.savefig(FIGURES / "fig1_binding_features.png", dpi=150, bbox_inches="tight")
    plt.savefig(FIGURES / "fig1_binding_features.pdf", bbox_inches="tight")
    plt.close()
    print("  Saved fig1_binding_features.png/pdf")


def fig2_feature_importance(data):
    """Figure 2: Feature importance ranking."""
    if not HAS_MPL:
        return
    importance_data = data.get("feature_importance", [])
    if not importance_data:
        return

    names = [d.get("feature", "") for d in importance_data[:10]]
    values = [float(d.get("importance", 0)) for d in importance_data[:10]]

    fig, ax = plt.subplots(figsize=(8, 5))
    y_pos = np.arange(len(names))
    ax.barh(y_pos, values, color="steelblue")
    ax.set_yticks(y_pos)
    ax.set_yticklabels(names, fontsize=9)
    ax.set_xlabel("Importance")
    ax.set_title("Feature Importance (Permutation)")
    ax.invert_yaxis()
    plt.tight_layout()
    plt.savefig(FIGURES / "fig2_feature_importance.png", dpi=150, bbox_inches="tight")
    plt.savefig(FIGURES / "fig2_feature_importance.pdf", bbox_inches="tight")
    plt.close()
    print("  Saved fig2_feature_importance.png/pdf")


def fig3_ablation(data):
    """Figure 3: Ablation study results."""
    if not HAS_MPL:
        return
    ablations = data.get("ablation_results", {})
    if not ablations:
        return

    models = []
    aurocs = []
    for name, metrics in ablations.items():
        models.append(name.replace("_", " ").title())
        aurocs.append(metrics.get("auroc", 0))

    if not models:
        return

    fig, ax = plt.subplots(figsize=(10, 5))
    x_pos = np.arange(len(models))
    bars = ax.bar(x_pos, aurocs, color="steelblue", alpha=0.8)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(models, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("AUROC")
    ax.set_title("Ablation Study: Feature Family Removal")
    ax.set_ylim(0, 1)
    plt.tight_layout()
    plt.savefig(FIGURES / "fig3_ablation.png", dpi=150, bbox_inches="tight")
    plt.savefig(FIGURES / "fig3_ablation.pdf", bbox_inches="tight")
    plt.close()
    print("  Saved fig3_ablation.png/pdf")


def fig4_stability(data):
    """Figure 4: Feature importance stability."""
    if not HAS_MPL:
        return
    stability = data.get("stability_results", {})
    if not stability:
        return

    ranking = stability.get("ranking", [])
    if not ranking:
        return

    names = [r["feature"] for r in ranking[:10]]
    means = [r["mean_importance"] for r in ranking[:10]]
    ci_lower = [r["ci_95_lower"] for r in ranking[:10]]
    ci_upper = [r["ci_95_upper"] for r in ranking[:10]]

    fig, ax = plt.subplots(figsize=(8, 5))
    y_pos = np.arange(len(names))
    errors_lower = [m - l for m, l in zip(means, ci_lower)]
    errors_upper = [u - m for m, u in zip(means, ci_upper)]
    ax.barh(y_pos, means, xerr=[errors_lower, errors_upper],
            color="steelblue", alpha=0.8, capsize=3)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(names, fontsize=9)
    ax.set_xlabel("Importance")
    ax.set_title("Feature Importance Stability (Bootstrap 95% CI)")
    ax.invert_yaxis()
    plt.tight_layout()
    plt.savefig(FIGURES / "fig4_stability.png", dpi=150, bbox_inches="tight")
    plt.savefig(FIGURES / "fig4_stability.pdf", bbox_inches="tight")
    plt.close()
    print("  Saved fig4_stability.png/pdf")


def fig5_structural_features(data):
    """Figure 5: Structural feature distributions."""
    if not HAS_MPL:
        return
    structural = data.get("structural_features", {})
    if not structural:
        return

    features_to_plot = ["contact_count", "hbond_count", "min_distance"]
    available = []
    for f in features_to_plot:
        for v in structural.values():
            if f in v:
                available.append(f)
                break
    available = available[:3]

    if not available:
        return

    fig, axes = plt.subplots(1, len(available), figsize=(4 * len(available), 4))
    if len(available) == 1:
        axes = [axes]

    for ax, feat in zip(axes, available):
        vals = [float(v.get(feat, 0)) for v in structural.values() if feat in v]
        if vals:
            ax.hist(vals, bins=15, color="steelblue", alpha=0.7)
        ax.set_xlabel(feat)
        ax.set_ylabel("Count")

    plt.tight_layout()
    plt.savefig(FIGURES / "fig5_structural_features.png", dpi=150, bbox_inches="tight")
    plt.savefig(FIGURES / "fig5_structural_features.pdf", bbox_inches="tight")
    plt.close()
    print("  Saved fig5_structural_features.png/pdf")


def fig6_nma_fluctuations(data):
    """Figure 6: NMA eigenvalue spectra."""
    if not HAS_MPL:
        return
    nma = data.get("nma_features", {})
    if not nma:
        return

    fig, ax = plt.subplots(figsize=(8, 5))
    for name, features in list(nma.items())[:5]:
        evals = features.get("eigenvalues", [])
        if evals:
            ax.plot(range(1, len(evals) + 1), evals, "o-", markersize=4,
                    label=name[:20], alpha=0.7)

    ax.set_xlabel("Mode Index")
    ax.set_ylabel("Eigenvalue")
    ax.set_title("NMA Eigenvalue Spectra")
    ax.legend(fontsize=7)
    plt.tight_layout()
    plt.savefig(FIGURES / "fig6_nma_eigenvalues.png", dpi=150, bbox_inches="tight")
    plt.savefig(FIGURES / "fig6_nma_eigenvalues.pdf", bbox_inches="tight")
    plt.close()
    print("  Saved fig6_nma_eigenvalues.png/pdf")


def fig7_network_features(data):
    """Figure 7: Network properties."""
    if not HAS_MPL:
        return
    network = data.get("network_features", {})
    if not network:
        return

    fig, axes = plt.subplots(1, 3, figsize=(12, 4))

    features = ["network_efficiency", "mean_betweenness", "n_communities"]
    for ax, feat in zip(axes, features):
        vals = [float(v.get(feat, 0)) for v in network.values() if feat in v]
        if vals:
            ax.hist(vals, bins=10, color="steelblue", alpha=0.7)
        ax.set_xlabel(feat)
        ax.set_ylabel("Count")

    plt.tight_layout()
    plt.savefig(FIGURES / "fig7_network_features.png", dpi=150, bbox_inches="tight")
    plt.savefig(FIGURES / "fig7_network_features.pdf", bbox_inches="tight")
    plt.close()
    print("  Saved fig7_network_features.png/pdf")


def fig8_state_analysis(data):
    """Figure 8: State analysis RMSD distribution."""
    if not HAS_MPL:
        return
    state = data.get("state_analysis", {})
    if not state:
        return

    fig, ax = plt.subplots(figsize=(6, 4))
    rmsds = [float(v.get("ca_rmsd", 0)) for v in state.values()]
    if rmsds:
        ax.hist(rmsds, bins=15, color="steelblue", alpha=0.7)
        ax.axvline(np.mean(rmsds), color="red", linestyle="--",
                   label=f"Mean: {np.mean(rmsds):.2f}")
        ax.legend()
    ax.set_xlabel("CA-RMSD (Angstrom)")
    ax.set_ylabel("Count")
    ax.set_title("Structural RMSD Distribution vs Reference")
    plt.tight_layout()
    plt.savefig(FIGURES / "fig8_state_analysis.png", dpi=150, bbox_inches="tight")
    plt.savefig(FIGURES / "fig8_state_analysis.pdf", bbox_inches="tight")
    plt.close()
    print("  Saved fig8_state_analysis.png/pdf")


def fig9_correlation_heatmap(data):
    """Figure 9: Correlation heatmap of features."""
    if not HAS_MPL:
        return

    # Collect all numeric features
    binding = data.get("binding_features", {})
    if not binding:
        return

    feat_names = ["docking_score", "boltz2_affinity", "boltz2_probability",
                  "confidence_score", "iptm"]
    available = []
    for f in feat_names:
        for v in binding.values():
            if f in v:
                available.append(f)
                break

    if len(available) < 2:
        return

    cids = sorted(binding.keys())
    matrix = np.array([[binding[cid].get(f, 0) for f in available] for cid in cids])

    corr = np.corrcoef(matrix.T)

    fig, ax = plt.subplots(figsize=(6, 5))
    im = ax.imshow(corr, cmap="RdBu_r", vmin=-1, vmax=1)
    ax.set_xticks(range(len(available)))
    ax.set_yticks(range(len(available)))
    ax.set_xticklabels(available, rotation=45, ha="right", fontsize=8)
    ax.set_yticklabels(available, fontsize=8)
    plt.colorbar(im, ax=ax, shrink=0.8)
    ax.set_title("Feature Correlation Heatmap")
    plt.tight_layout()
    plt.savefig(FIGURES / "fig9_correlation_heatmap.png", dpi=150, bbox_inches="tight")
    plt.savefig(FIGURES / "fig9_correlation_heatmap.pdf", bbox_inches="tight")
    plt.close()
    print("  Saved fig9_correlation_heatmap.png/pdf")


def fig10_summary(data):
    """Figure 10: Summary dashboard."""
    if not HAS_MPL:
        return

    fig = plt.figure(figsize=(14, 10))
    gs = gridspec.GridSpec(3, 3, figure=fig, hspace=0.4, wspace=0.4)

    # Panel A: Model performance
    ax_a = fig.add_subplot(gs[0, 0])
    models = data.get("model_results", {})
    if models:
        basic_auroc = models.get("basic_features", {}).get("cv_results", {}).get("auroc", 0)
        aug_auroc = models.get("augmented_features", {}).get("cv_results", {}).get("auroc", 0)
        ax_a.bar(["Basic", "Augmented"], [basic_auroc, aug_auroc],
                color=["steelblue", "coral"])
        ax_a.set_ylabel("AUROC")
        ax_a.set_title("A. Model Performance")
        ax_a.set_ylim(0, 1)

    # Panel B: Feature importance
    ax_b = fig.add_subplot(gs[0, 1])
    importance = data.get("feature_importance", [])
    if importance:
        top5 = importance[:5]
        names = [d["feature"] for d in top5]
        vals = [float(d["importance"]) for d in top5]
        ax_b.barh(range(len(names)), vals, color="steelblue")
        ax_b.set_yticks(range(len(names)))
        ax_b.set_yticklabels(names, fontsize=7)
        ax_b.set_title("B. Top Features")
        ax_b.invert_yaxis()

    # Panel C: Ablation
    ax_c = fig.add_subplot(gs[0, 2])
    ablations = data.get("ablation_results", {})
    if ablations:
        names = list(ablations.keys())[:5]
        vals = [ablations[n].get("auroc", 0) for n in names]
        ax_c.bar(range(len(names)), vals, color="steelblue", alpha=0.7)
        ax_c.set_xticks(range(len(names)))
        ax_c.set_xticklabels([n[:10] for n in names], rotation=45, fontsize=6)
        ax_c.set_ylabel("AUROC")
        ax_c.set_title("C. Ablation")

    plt.savefig(FIGURES / "fig10_summary_dashboard.png", dpi=150, bbox_inches="tight")
    plt.savefig(FIGURES / "fig10_summary_dashboard.pdf", bbox_inches="tight")
    plt.close()
    print("  Saved fig10_summary_dashboard.png/pdf")


def main():
    print("=" * 70)
    print("FIGURE GENERATION - Mechanistic SAR Analysis")
    print("=" * 70)

    print("\n[1/2] Loading results...")
    data = load_all_results()
    print(f"  Loaded {len(data)} result files")

    print("\n[2/2] Generating figures...")
    fig1_binding_features(data)
    fig2_feature_importance(data)
    fig3_ablation(data)
    fig4_stability(data)
    fig5_structural_features(data)
    fig6_nma_fluctuations(data)
    fig7_network_features(data)
    fig8_state_analysis(data)
    fig9_correlation_heatmap(data)
    fig10_summary(data)

    print(f"\nAll figures saved to: {FIGURES}")
    print(f"  Total PNG files: {len(list(FIGURES.glob('*.png')))}")
    print(f"  Total PDF files: {len(list(FIGURES.glob('*.pdf')))}")


if __name__ == "__main__":
    main()
