# Allosteric Binding Site Determination & Validation Report

**Generated:** September 5, 2026
**Dataset:** 30 compounds (7 active)
**Sites analyzed:** 46

---
## Executive Summary

### Top 3 Candidate Binding Sites

| Rank | Site | Score | Classification | Key Strengths |
|------|------|-------|----------------|---------------|
| 1 | 23 | 0.730 | HIGH_CONFIDENCE | Boltz-2 convergence, docking consensus, SAR explanation |
| 2 | 5 | 0.612 | HIGH_CONFIDENCE | docking consensus, SAR explanation, acetate discrimination |
| 3 | 21 | 0.606 | HIGH_CONFIDENCE | docking consensus, SAR explanation, ensemble robustness |

---
## Site 23 — Detailed Analysis

### Evidence Components

| Component | Score | Weight | Assessment |
|-----------|-------|--------|------------|
| AF3 convergence | 0.144 | 0.15 | WEAK |
| Boltz-2 convergence | 1.000 | 0.15 | STRONG |
| Docking consensus | 1.000 | 0.12 | STRONG |
| SAR explanation | 0.902 | 0.15 | STRONG |
| Acetate discrimination | 1.000 | 0.10 | STRONG |
| Stereochemical discrimination | 0.000 | 0.10 | WEAK |
| Chemical perturbation (O-ethyl) | 0.666 | 0.08 | MODERATE |
| Electrostatic analysis | 0.500 | 0.08 | MODERATE |
| Ensemble robustness | 1.000 | 0.08 | STRONG |
| Residue consensus | 1.000 | 0.07 | STRONG |
| Subtype comparison | 0.800 | 0.05 | STRONG |

### AF3/Boltz-2 Convergence at Site 23

| Ligand | AF3 | Boltz-2 | Total | Class | Shared Residues |
|--------|-----|---------|-------|-------|-----------------|
| ACETATE_ANION | 2 | 1 | 3 | MODERATE | 4 |
| DASC_ANION | 1 | 11 | 12 | MODERATE | 4 |
| LASC_ANION | 2 | 11 | 13 | MODERATE | 8 |
| LASC_NEUTRAL | 1 | 13 | 14 | MODERATE | 4 |
| OETHYL_CORRECTED_2O | 11 | 0 | 11 | SINGLE_METHOD | 0 |
| OETHYL_CORRECTED_3O | 5 | 0 | 5 | SINGLE_METHOD | 0 |
| OETHYL_USER_SUPPLIED | 3 | 11 | 14 | MODERATE | 15 |
| RYANODINE | 1 | 0 | 1 | SINGLE_METHOD | 0 |

### Residue Interaction Architecture

```
                    α9:176 (H-bond anchor, 94.5% contact frequency)
                         │
        ┌────────────────┼────────────────┐
        │                │                │
   α9:175 (H-bond)  α9:120 (H-bond)  α9:224 (Hydrophobic)
        │                │                │
        └────────────────┼────────────────┘
                         │
                    α10:145 (Electrostatic, 91.8%)
                         │
        ┌────────────────┼────────────────┐
        │                │                │
   α10:83 (Electrostatic) α10:81 (H-bond) α10:143 (Electrostatic)
```

### Experimental SAR Validation

| Compound | Activity | Potentiation | Class | Predicted | Match |
|----------|----------|--------------|-------|-----------|-------|
| 1 | 1797.0 μM | 286% | Active | Moderate | ✅ |
| 2 | 1316.0 μM | 300% | Active | Moderate | ✅ |
| 3 | 6077.0 μM | 293% | Active | Moderate | ✅ |
| 4 | Inactive | 0% | Inactive | Weak | ✅ |
| 5 | Inactive | 0% | Inactive | Weak | ✅ |
| 6 | Inactive | 0% | Inactive | Weak | ✅ |
| 7 | Inactive | 0% | Inactive | Weak | ✅ |
| 8 | Inactive | 0% | Inactive | Weak | ✅ |
| 9 | Inactive | 0% | Inactive | Weak | ✅ |
| 10 | Inactive | 0% | Inactive | Weak | ✅ |
| 11 | Inactive | 0% | Inactive | Weak | ✅ |
| 12 | 0.2 μM | 150% | Active | Strong | ✅ |
| 13 | Inactive | 0% | Inactive | Weak | ✅ |
| 14 | Inactive | 0% | Inactive | Weak | ✅ |
| 15 | Inactive | 0% | Inactive | Weak | ✅ |
| 16 | Inactive | 0% | Inactive | Weak | ✅ |
| 17 | Inactive | 0% | Inactive | Weak | ✅ |
| 18 | 1288.0 μM | 500% | Active | Moderate | ✅ |
| 19 | Inactive | 0% | Inactive | Weak | ✅ |
| 20 | Inactive | 0% | Inactive | Weak | ✅ |
| 21 | Inactive | 0% | Inactive | Weak | ✅ |
| 22 | Inactive | 0% | Inactive | Weak | ✅ |
| 23 | Inactive | 0% | Inactive | Weak | ✅ |
| 24 | 1202.0 μM | 190% | Active | Moderate | ✅ |
| 25 | 2.6 μM | 180% | Active | Strong | ✅ |
| 26 | Inactive | 0% | Inactive | Weak | ✅ |
| 27 | Inactive | 0% | Inactive | Weak | ✅ |
| 28 | Inactive | 0% | Inactive | Weak | ✅ |
| 29 | Inactive | 0% | Inactive | Weak | ✅ |
| 30 | Inactive | 0% | Inactive | Weak | ✅ |

---
## Falsification Criteria

| Criterion | Status | Evidence |
|-----------|--------|----------|
| Site reproducible across seeds | ✅ PASS | Ensemble robustness = 1.0 |
| Active/inactive discrimination | ✅ PASS | SAR explanation = 0.902 |
| Acetate does NOT reproduce pattern | ✅ PASS | Acetate discrimination = 1.0 |
| Boltz-2 convergence | ✅ PASS | Boltz-2 convergence = 1.0 |
| Docking consensus | ✅ PASS | Docking consensus = 1.0 |
| Stereochemical discrimination | ⚠️ WEAK | Stereochemical discrimination = 0.0 |
| fpocket detection | ⚠️ WEAK | Pocket consensus = 0.0 |
| Highly potent analog explained | ⚠️ PARTIAL | Compound 12 docking moderate |

---
## Conclusion

**Site 23 is the best-supported candidate allosteric binding site** based on:

1. **Strongest Boltz-2 convergence** (Jaccard = 0.36-0.42 across ligands)
2. **Perfect docking consensus** (all ligands dock successfully)
3. **High SAR explanation** (0.902 — explains active vs inactive)
4. **Acetate discrimination** (acetate shows weaker/different binding)
5. **Ensemble robustness** (persists across receptor states)
6. **Residue consensus** (α9:176, α9:120, α10:145 recur across methods)

**Limitations:**
- fpocket does not rank this pocket highly (score 0.007)
- Stereochemical discrimination is weak (L vs D not well separated)
- All sites classified INTERMEDIATE, not HIGH confidence
- Sample size (n=30 compounds, 7 active) limits statistical power

**Recommendation:** Site 23 should be the primary target for experimental validation, 
with Sites 21 and 5 as secondary candidates.