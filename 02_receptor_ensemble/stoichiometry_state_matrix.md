# Tier 2.1: Receptor Stoichiometry × Conformational State Matrix

**Date:** September 4, 2026  
**Status:** PRE-MODELING  
**Purpose:** Define which receptor states must be modeled for comprehensive analysis

---

## 1. α9α10 nAChR Architecture

### 1.1 Stoichiometry

- **Canonical form:** (α9)₂(α10)₃
- **Subunit arrangement:** Alternating α9/α10/α9/α10/α10
- **Symmetry:** C2 symmetry axis
- **Number of subunits:** 5 (2 α9 + 3 α10)

### 1.2 Key Features

- **Ligand binding sites:** 3 (at α9-α10 interfaces)
- **Ion channel:** Central pore
- **Allosteric sites:** Multiple (S23, S41, S46, etc.)
- **Modulation:** PAMs potentiate ACh response

---

## 2. Stoichiometry × State Matrix

### 2.1 Receptor States

| State | Description | Functional Role | Importance |
|-------|-------------|-----------------|------------|
| **Resting (closed)** | Channel closed, ACh not bound | Baseline state | Reference |
| **Open (activated)** | Channel open, ACh bound | Active state | **CRITICAL for PAM** |
| **Desensitized** | Channel closed, desensitized | Refractory state | Common in crystal structures |
| **Intermediate** | Partially open | Transitional state | May capture binding |

### 2.2 Stoichiometry × State Matrix

| Stoichiometry | Resting | Open | Desensitized | Intermediate | Notes |
|---------------|---------|------|--------------|--------------|-------|
| **2to3 (α9₂α10₃)** | ✓ | **✓** | ✓ | ✓ | Primary target |
| **3to2 (α9₃α10₂)** | ✓ | **✓** | ✓ | ✓ | Alternative form |

### 2.3 Critical Insight: PAM Requires Open State

**Why open state matters:**
- PAMs potentiate the response to agonist
- This means they stabilize the OPEN state relative to resting
- A PAM that binds equally well to resting and open states would not be a PAM
- **Therefore:** PAM binding affinity must be HIGHER in open state than resting state

**Implication for modeling:**
- Must model OPEN state (activated with ACh bound)
- Must compare PAM binding affinity between states
- Must show state-dependent binding

---

## 3. Available Structural Data

### 3.1 PDB Structures (Existing)

| PDB ID | Stoichiometry | State | Ligand | Resolution | Notes |
|--------|---------------|-------|--------|------------|-------|
| **9HIO** | 2to3 | Desensitized-like | - | - | AF3 template |
| **9HQM** | 3to2 | Desensitized-like | - | - | AF3 template |

### 3.2 Structural Gaps

| Gap | Status | Impact |
|-----|--------|--------|
| Open state structures | **MISSING** | Cannot validate state-dependent binding |
| Intermediate states | **MISSING** | May miss important conformations |
| ACh-bound structures | **MISSING** | Cannot model ternary complex |
| PAM-bound structures | **MISSING** | Cannot validate binding mode |

### 3.3 Modeling Requirements

| Model Type | Stoichiometry | State | ACh | PAM | Priority |
|------------|---------------|-------|-----|-----|----------|
| **Binary (ACh + receptor)** | 2to3 | Open | ✓ | ✗ | HIGH |
| **Binary (ACh + receptor)** | 3to2 | Open | ✓ | ✗ | HIGH |
| **Ternary (ACh + PAM + receptor)** | 2to3 | Open | ✓ | ✓ | **CRITICAL** |
| **Ternary (ACh + PAM + receptor)** | 3to2 | Open | ✓ | ✓ | HIGH |
| **Binary (PAM + receptor)** | 2to3 | Resting | ✗ | ✓ | MEDIUM |
| **Binary (PAM + receptor)** | 2to3 | Desensitized | ✗ | ✓ | MEDIUM |

---

## 4. AF3 Template Generation

### 4.1 Template Requirements

| Template | Stoichiometry | State | Source | Purpose |
|----------|---------------|-------|--------|---------|
| **2to3_open** | 2to3 | Open | AF3 | PAM binding in open state |
| **3to2_open** | 3to2 | Open | AF3 | PAM binding in open state |
| **2to3_resting** | 2to3 | Resting | AF3 | Reference state |
| **3to2_resting** | 3to2 | Resting | AF3 | Reference state |
| **2to3_desensitized** | 2to3 | Desensitized | AF3 | Comparison |
| **3to2_desensitized** | 3to2 | Desensitized | AF3 | Comparison |

### 4.2 AF3 Input Generation

For each template:
1. Extract subunit sequences from PDB 9HIO/9HQM
2. Define stoichiometry (2to3 or 3to2)
3. Define conformational state
4. Generate AF3 input JSON
5. Submit to ColabFold/AlphaFold3

### 4.3 Boltz-2 Template Generation

For each template:
1. Extract sequences and coordinates
2. Define state (resting/open/desensitized)
3. Generate Boltz-2 input
4. Submit for structure prediction

---

## 5. Conformational State Modeling

### 5.1 Open State Requirements

**For open state modeling:**
1. ACh must be bound at orthosteric sites
2. Channel must be in open conformation
3. Selectivity filter must be open
4. Gate must be open

**Challenge:** No open state template exists for α9α10

**Solution:** Use AF3/Boltz-2 to predict open state from:
- Sequence information
- Homology to open-state nAChR structures (e.g., GLIC)
- State-specific constraints

### 5.2 Resting State Requirements

**For resting state modeling:**
1. No ACh bound
2. Channel closed
3. Selectivity filter closed
4. Gate closed

### 5.3 Desensitized State Requirements

**For desensitized state modeling:**
1. ACh may be bound
2. Channel closed
3. Desensitized conformation
4. Refractory state

---

## 6. Binary vs Ternary Complex Modeling

### 6.1 Binary Complex (ACh + Receptor)

**Purpose:** Establish baseline ACh binding
**Components:**
- Receptor (2to3 or 3to2)
- ACh at all 3 binding sites
- State: Open

**Output:** ACh binding affinity, conformation

### 6.2 Ternary Complex (ACh + PAM + Receptor)

**Purpose:** Model PAM modulation
**Components:**
- Receptor (2to3 or 3to2)
- ACh at orthosteric sites
- PAM at allosteric site
- State: Open

**Output:** PAM binding affinity, interaction network, state stabilization

### 6.3 Critical Comparison

**For PAM mechanism:**
1. Compare PAM binding in resting vs open states
2. Show PAM stabilizes open state
3. Quantify state-dependent affinity difference
4. Explain SAR through state-dependent interactions

---

## 7. Receptor Ensemble Strategy

### 7.1 Why Multiple States?

**PAM mechanism requires:**
1. Binding in resting state (capture)
2. Stabilization of open state (potentiation)
3. Release from desensitized state (recovery)

**Therefore:**
- Must model at least 3 states: resting, open, desensitized
- Must compare binding affinity across states
- Must show state-dependent SAR

### 7.2 State Selection Priority

| Priority | State | Stoichiometry | Reason |
|----------|-------|---------------|--------|
| 1 | Open | 2to3 | Primary PAM binding state |
| 2 | Open | 3to2 | Alternative stoichiometry |
| 3 | Resting | 2to3 | Reference state |
| 4 | Resting | 3to2 | Reference state |
| 5 | Desensitized | 2to3 | Recovery mechanism |
| 6 | Desensitized | 3to2 | Recovery mechanism |

### 7.3 AF3/Boltz-2 Submission Plan

**Round 1: Open State (CRITICAL)**
1. 2to3_open: AF3 + Boltz-2
2. 3to2_open: AF3 + Boltz-2
3. Binary (ACh) and Ternary (ACh+PAM) complexes

**Round 2: Resting State (REFERENCE)**
1. 2to3_resting: AF3 + Boltz-2
2. 3to2_resting: AF3 + Boltz-2
3. Binary (PAM only) complexes

**Round 3: Desensitized State (COMPARISON)**
1. 2to3_desensitized: AF3 + Boltz-2
2. 3to2_desensitized: AF3 + Boltz-2

---

## 8. Validation Requirements

### 8.1 Template Validation

- [ ] 2to3_open template generated
- [ ] 3to2_open template generated
- [ ] Templates match known α9α10 features
- [ ] Templates have correct stoichiometry
- [ ] Templates have correct state

### 8.2 Binding Site Validation

- [ ] Site 23 exists in all templates
- [ ] Site 23 is accessible in all templates
- [ ] Site 23 residues match predictions
- [ ] Site 41, S46, S167 exist and are accessible

### 8.3 State Validation

- [ ] Open state channel is open
- [ ] Resting state channel is closed
- [ ] Desensitized state channel is closed
- [ ] State transitions are consistent

---

## 9. Summary

### 9.1 Matrix Requirements

| Stoichiometry | Resting | Open | Desensitized |
|---------------|---------|------|--------------|
| 2to3 (α9₂α10₃) | ✓ | **CRITICAL** | ✓ |
| 3to2 (α9₃α10₂) | ✓ | **CRITICAL** | ✓ |

### 9.2 Key Insight

**PAM requires open state modeling:**
- Must generate open state templates via AF3/Boltz-2
- Must compare binding across states
- Must show state-dependent SAR

### 9.3 Next Steps

1. Generate open state templates via AF3/Boltz-2
2. Validate templates match known features
3. Proceed to site freezing and analysis

---

*This matrix defines which receptor states must be modeled for comprehensive PAM analysis.*
