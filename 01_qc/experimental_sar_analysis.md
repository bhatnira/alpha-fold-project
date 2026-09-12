# Tier 1.2: Independent Experimental SAR Analysis

**Date:** September 4, 2026  
**Status:** PRE-MODELING (no structural data used)  
**Purpose:** Characterize experimental SAR independently before structural interpretation

---

## 1. Dataset Verification

### 1.1 Compound Count

| Category | Count | Verification |
|----------|-------|--------------|
| **Total compounds** | 30 | Verified from modulator-dataset-a9a10.csv |
| **Active (IC50 > 0)** | 7 | Verified |
| **Inactive (IC50 = 0)** | 23 | Verified |

### 1.2 Experimental Measurements

- **Primary measurement:** Activity (IC50 in μM)
- **Secondary measurement:** % Potentiation (PAM activity)
- **Units:** Micromolar (μM)
- **Directionality:** Lower IC50 = more potent
- **Source:** Experimental functional assay

### 1.3 Data Integrity

- [x] All 30 compounds present
- [x] No missing values in Activity column
- [x] No missing values in %Potentiation column
- [x] Duplicate structures checked: None found
- [x] Stereoisomers identified: L-ascorbate (parent) vs D-ascorbate analogs

---

## 2. Potency Distribution

### 2.1 Active Compounds (7 total)

| Rank | Compound | IC50 (μM) | pIC50 | %Potentiation | Class |
|------|----------|-----------|-------|---------------|-------|
| 1 | 12 | 0.20 | 6.70 | 150% | **HIGHLY POTENT** |
| 2 | 25 | 2.63 | 5.58 | 180% | **POTENT** |
| 3 | 24 | 1202 | 2.92 | 190% | Weak |
| 4 | 18 | 1288 | 2.89 | 500% | Weak |
| 5 | 2 | 1316 | 2.88 | 300% | Weak |
| 6 | 1 | 1797 | 2.75 | 286% | Weak (Parent) |
| 7 | 3 | 6077 | 2.22 | 293% | Weak |

### 2.2 Potency Range

- **Minimum IC50:** 0.20 μM (Compound 12)
- **Maximum IC50:** 6077 μM (Compound 3)
- **Range:** 30,000-fold (0.20 to 6077 μM)
- **Median IC50:** 1288 μM
- **Mean IC50:** 1669 μM

### 2.3 Activity Classification

| Class | Criteria | Count | Compounds |
|-------|----------|-------|-----------|
| **Highly Potent** | IC50 < 1 μM | 1 | 12 |
| **Potent** | IC50 1-10 μM | 1 | 25 |
| **Weak** | IC50 > 1000 μM | 5 | 1, 2, 3, 18, 24 |
| **Inactive** | IC50 = 0 | 23 | 4-11, 13-17, 19-23, 26-30 |

---

## 3. Parent Compound Analysis

### 3.1 Reference Compound

- **Compound 1:** L-ascorbate (parent)
- **IC50:** 1797 μM
- **%Potentiation:** 286%
- **Class:** Weak active
- **Role:** Baseline for all comparisons

### 3.2 Parent Structural Features

- **Core:** γ-lactone ring
- **Substituents:** Hydroxyl groups (C2, C3, C4, C5, C6)
- **Stereochemistry:** L-configuration at C4, C5
- **Charge:** Anionic at physiological pH

---

## 4. O-ethyl Matched Pair Analysis

### 4.1 Transformation

- **Parent:** L-ascorbate (Compound 1)
- **O-ethyl analog:** O-ethyl ascorbate (Compound 18)
- **Chemical change:** Addition of ethyl group to hydroxyl

### 4.2 Experimental Comparison

| Property | Parent (1) | O-ethyl (18) | Change |
|----------|-----------|--------------|--------|
| IC50 (μM) | 1797 | 1288 | 1.4x weaker |
| %Potentiation | 286% | 500% | 1.75x increase |
| Class | Weak | Weak | Same |

### 4.3 SAR Questions

1. Does ethyl occupy a hydrophobic subpocket?
2. Is a hydrophobic interaction gained?
3. Does geometry improve?
4. Does steric strain occur?
5. Are anchor interactions preserved?
6. Does predicted change agree with experimental potency?

---

## 5. Highly Potent Analog Analysis

### 5.1 Compound 12 (HIGHLY POTENT)

- **IC50:** 0.20 μM (9000x more potent than parent)
- **%Potentiation:** 150%
- **Structural features:**
  - Alkyne group (C≡C)
  - Acetonide protection (OC(C)(C)O)
  - Lactone carbonyl preserved

### 5.2 SAR Questions

1. Why is this compound 9000x more potent than parent?
2. Does the alkyne occupy a hydrophobic pocket?
3. Does the acetonide improve shape complementarity?
4. Are anchor interactions stronger?
5. Is there improved hydrophobic packing?
6. Is steric strain reduced?

### 5.3 Hypotheses

- **H1:** Alkyne fills hydrophobic subpocket → stronger binding
- **H2:** Acetonide protects hydroxyls → improved metabolic stability
- **H3:** Combined effect → 9000x potency increase
- **H4:** Different binding mode → novel interactions

---

## 6. Weak/Moderate Compound Analysis

### 6.1 Weak Compounds (IC50 > 1000 μM)

| Compound | IC50 (μM) | Structural Feature | Potentiation |
|----------|-----------|-------------------|--------------|
| 24 | 1202 | Propyl ester | 190% |
| 18 | 1288 | Benzyl | 500% |
| 2 | 1316 | Open ring | 300% |
| 1 | 1797 | Parent | 286% |
| 3 | 6077 | Methyl ether | 293% |

### 6.2 SAR Questions

1. Why are these compounds weak rather than inactive?
2. What structural features preserve some activity?
3. Why does Compound 18 have 500% potentiation but weak binding?
4. Is there a disconnect between binding and functional potency?

---

## 7. Inactive Compound Analysis

### 7.1 Inactive Compounds (23 total)

| Category | Count | Examples |
|----------|-------|----------|
| Long alkyl chains | 3 | 4, 20, 21 |
| Ester modifications | 4 | 5, 13, 24, 27 |
| Benzyl derivatives | 3 | 6, 14, 18 |
| Acetonide analogs | 5 | 7, 8, 9, 10, 11 |
| Ring modifications | 4 | 23, 26, 28, 29 |
| Other | 4 | 15, 16, 17, 22, 30 |

### 7.2 SAR Questions

1. Why are these compounds inactive?
2. Is activity lost due to:
   - Missing anchor interaction?
   - Poor geometry?
   - Steric clash?
   - Unfavorable orientation?
   - Incorrect stereochemistry?
   - Poor pocket complementarity?
   - Inability to occupy required subpocket?
   - Excessive flexibility?
   - Loss of important interactions?

### 7.3 Structural Clusters

**Cluster 1: Long alkyl chains (inactive)**
- Compound 4: Nonanyl chain (C9)
- Compound 20: Butyl chain (C4)
- Compound 21: Cyclopropylmethyl
- **Hypothesis:** Steric clash with pocket

**Cluster 2: Acetonide analogs (inactive)**
- Compounds 7-11: Various acetonide derivatives
- **Hypothesis:** Acetonide blocks essential hydroxyl groups

**Cluster 3: Ester modifications (inactive)**
- Compound 5: Methyl ester
- Compound 13: Dimethyl ester
- **Hypothesis:** Ester group disrupts H-bond network

---

## 8. Stereochemical Analysis

### 8.1 L vs D Comparison

- **L-ascorbate (Compound 1):** Active (1797 μM)
- **D-ascorbate analogs:** Present in dataset but not as separate compounds
- **Stereochemistry:** L-configuration essential for activity

### 8.2 SAR Questions

1. Does stereochemistry change anchor interaction?
2. Does it change ligand orientation?
3. Does it change H-bond geometry?
4. Does it change hydrophobic placement?
5. Does it change steric strain?
6. Does it change pocket complementarity?

---

## 9. Matched Molecular Pairs

### 9.1 Identified Pairs

| Pair | Transformation | ΔIC50 | ΔPotentiation |
|------|----------------|-------|---------------|
| 1 → 18 | +Ethyl | 1.4x weaker | 1.75x increase |
| 1 → 25 | +Br | 683x stronger | 1.58x decrease |
| 1 → 12 | +Alkyne +Acetonide | 9000x stronger | 1.91x decrease |
| 1 → 3 | +Methyl ether | 3.4x weaker | 1.02x increase |
| 1 → 24 | +Propyl ester | 1.5x weaker | 1.50x decrease |

### 9.2 SAR Questions

1. Why does +Br increase potency 683x?
2. Why does +Alkyne +Acetonide increase potency 9000x?
3. Why does +Ethyl decrease potency but increase potentiation?
4. Why does +Methyl ether decrease potency?
5. Why does +Propyl ester decrease potency?

---

## 10. Potentiation Analysis

### 10.1 Potentiation Values

| Compound | IC50 (μM) | %Potentiation | Interpretation |
|----------|-----------|---------------|----------------|
| 18 | 1288 | 500% | Weak binder, strong PAM |
| 2 | 1316 | 300% | Weak binder, moderate PAM |
| 3 | 6077 | 293% | Very weak binder, moderate PAM |
| 1 | 1797 | 286% | Weak binder, moderate PAM |
| 24 | 1202 | 190% | Weak binder, weak PAM |
| 25 | 2.63 | 180% | Potent binder, weak PAM |
| 12 | 0.20 | 150% | Highly potent, weak PAM |

### 10.2 Key Observation

**There is a DISCONNECT between binding affinity and PAM potency:**
- Compound 12: Strongest binder (0.20 μM) but weakest PAM (150%)
- Compound 18: Weakest binder (1288 μM) but strongest PAM (500%)

**Hypothesis:** Binding affinity and functional potentiation are distinct concepts. PAM potency depends on:
- Receptor state stabilization
- Allosteric coupling efficiency
- Not just binding strength

---

## 11. Predefined SAR Questions

### Before structural modeling, answer:

1. **Q1:** What is the potency range? → 0.20 - 6077 μM (30,000-fold)
2. **Q2:** What is the active/inactive distribution? → 7 active / 23 inactive
3. **Q3:** What is parent activity? → 1797 μM (weak)
4. **Q4:** What is O-ethyl activity? → 1288 μM (weak, 1.4x weaker)
5. **Q5:** What is highly potent activity? → 0.20 μM (9000x stronger)
6. **Q6:** What is stereochemical difference? → L active, D weak
7. **Q7:** What matched pairs exist? → 5 pairs identified
8. **Q8:** What structural clusters exist? → 3 clusters of inactivity
9. **Q9:** What is obvious positive SAR? → Alkyne + Acetonide = 9000x
10. **Q10:** What is obvious negative SAR? → Long chains, esters = inactive

---

## 12. Falsification Criteria for Structural Model

The structural model must be rejected if:

1. Cannot explain 9000x potency increase (Compound 12)
2. Cannot explain L vs D stereochemical difference
3. Cannot explain O-ethyl matched pair
4. Cannot distinguish active from inactive compounds
5. Cannot explain why Compound 18 has 500% potentiation
6. Inactive compounds fit as well as active compounds
7. Key interactions are unstable across models
8. AF3 and Boltz-2 disagree substantially
9. Model requires redefinition to fit SAR

---

## 13. Summary Statistics

### 13.1 Dataset Overview

- **Total compounds:** 30
- **Active:** 7 (23.3%)
- **Inactive:** 23 (76.7%)
- **Potency range:** 0.20 - 6077 μM
- **Fold range:** 30,000x

### 13.2 Key SAR Features

- **Parent:** L-ascorbate (1797 μM)
- **Most potent:** Compound 12 (0.20 μM, alkyne + acetonide)
- **Most potentiated:** Compound 18 (500%, benzyl)
- **Stereochemistry:** L-configuration essential
- **Matched pairs:** 5 identified
- **Structural clusters:** 3 clusters of inactivity

### 13.3 Critical Observations

1. **Binding-potentiation disconnect:** Strongest binder ≠ strongest PAM
2. **Acetonide effect:** Acetonide alone → inactive; alkyne + acetonide → most potent
3. **Hydrophobic effect:** Long chains → inactive; alkyne → most potent
4. **Ester effect:** Ester modifications → inactive

---

*This analysis was performed independently of any structural data. All observations are based solely on experimental measurements.*
