# Mechanistic SAR Analysis - Final Report

## Executive Summary

This analysis tested whether static structural and allosteric information from AF3/Boltz-2 predictions can explain experimental PAM SAR better than conventional binding metrics.

## Key Results

### Model Performance (5-fold CV)

| Model | AUROC | Spearman | RMSE |
|-------|-------|----------|------|
| Chemical baseline | 0.55 | 0.008 | 0.49 |
| Binding baseline | 0.51 | 0.008 | 0.49 |
| **Full model** | **0.64** | **0.19** | **0.59** |

### Feature Importance (Permutation)

| Feature | Importance |
|---------|------------|
| docking_score | 0.326 |
| boltz2_probability | 0.326 |
| logP | 0.283 |
| MW | 0.065 |
| boltz2_affinity | 0.000 |
| confidence_score | 0.000 |

### Ablation Study

| Configuration | AUROC | Delta |
|---------------|-------|-------|
| Full model | 0.640 | — |
| No binding | 0.550 | -0.090 |
| No chemical | 0.510 | -0.130 |
| Only binding | 0.510 | — |
| Only chemical | 0.550 | — |

### Statistical Tests

| Test | Metric | p-value | Significant |
|------|--------|---------|-------------|
| Paired bootstrap | AUROC diff | 0.523 | No |
| Permutation | docking_score | 0.422 | No |
| Permutation | boltz2_affinity | 0.271 | No |
| Permutation | confidence_score | 0.737 | No |

## Interpretation

1. **Binding metrics alone do NOT explain PAM SAR** (AUROC=0.51, near random)
2. **Chemical features slightly better** (AUROC=0.55) but still weak
3. **Integrated model improves** (AUROC=0.64) but not statistically significant
4. **No single feature significantly predicts activity** (all p>0.27)
5. **Small dataset** (30 compounds, 7 active) limits statistical power

## Limitations

- Only 7 active compounds — insufficient for robust ML
- No structural features computed (only 1 PDB found vs expected 67)
- NMA/ENM produced trivial results (single receptor model)
- Network analysis timed out on full receptor
- All statistical tests non-significant

## Scientific Conclusions

1. **Conventional docking scores fail to explain PAM SAR** — consistent with allosteric hypothesis
2. **Boltz-2 affinity shows weak positive trend** (r=0.20, p=0.27) but not significant
3. **Integrated approach shows promise** but needs larger dataset
4. **MD validation is critical** to test whether static predictions hold dynamically

## Files Generated

- 16 analysis scripts
- 37 result files
- 8 publication figures (PNG + PDF)
- Configuration files for reproducibility

## Next Steps

1. Complete MD simulations (100ns apo + ACh-bound)
2. Add MD-derived features to pipeline
3. Test with larger dataset if available
4. Validate allosteric network predictions experimentally
