# REMAINING JOBS — Master Prompt Completion Tracker

**Created:** September 6, 2026
**Master Prompt Reference:** PNAS 2521048122 site-directed cofolding + α9α10 PAM discovery
**Last Updated:** September 13, 2026 (validation module 22_binding_site_validation scripts 01-14 fixed + ready to run)

---

## STATUS OVERVIEW

| Category | Done | Running/Queued | Remaining | Total |
|----------|------|----------------|-----------|-------|
| Core pipeline phases | 22/22 (smoke-tested) | 0 | re-verify | 22 |
| Receptor state models | 1/6 | 0 | 5 | 6 |
| ACh occupancy models | 0/2 predictions | 0 | 2 | 2 (re-run needed; cache bug fixed) |
| Site-directed cofolding | 0/160 | 0 | 160 | 160 |
| Full compound panel | 0/180 | 0 | 180 | 180 |
| AF3 outputs (existing) | 468 models | 0 | 0 | 468 |
| Boltz-2 outputs (existing) | 90 CIFs; 83/180 affinity jsons | 0 | re-run | invalid (α10-monomer input, FALLACIES_AUDIT C7) |
| MD simulations | 0/4 | 0 | 4 | 4 (all timed out; 0 DCD) |
| Analysis scripts | 3/8 | 0 | 5 | 8 |
| Reports/deliverables | 10/27 | 0 | 17 | 27 |
| Figures | 0/15 | 0 | 15 | 15 |
| Tables | 0/12 | 0 | 12 | 12 |

---

## PHASE MAP: MASTER PROMPT vs EXISTING PIPELINE

### Section Mapping

| Master Prompt Section | Pipeline Phase | Status | Notes |
|----------------------|----------------|--------|-------|
| §1 Non-negotiable principles | — | DOCUMENTED | Enforced in code comments |
| §2 Dataset | phase01_sar | DONE | 30 compounds, 7 active, 21 inactive |
| §3 Hypothesis space | — | PARTIAL | Stoichiometry/state/ACh/interface defined |
| §4 Overall workflow | — | MAPPED | This file |
| §5 Phase 1 — Data QC | phase01_sar | DONE | Discovery/held-out split |
| §6 Phase 2 — Experimental SAR | phase01_sar, 02_sar_analysis | DONE | SAR anchors identified |
| §7 Phase 3 — Receptor ensemble | phase02_receptor | PARTIAL | Only resting-like; open/desensitized MISSING |
| §8 Phase 4 — State assignment | 04_state_validation | PARTIAL | No pore/M2 analysis; states not classified |
| §9 Phase 5 — Structural QC | 04_state_validation | PARTIAL | Basic checks only |
| §10 Phase 6 — ACh occupancy | 05_ach_occupancy | PARTIAL | 2 runs, 0 predictions (cache FileExistsError, FIXED); re-run needed |
| §11 Phase 7 — Blind AF3 discovery | phase03_pockets, 04_convergence | DONE | fpocket + AF3/Boltz-2 convergence |
| §12 Phase 8 — Candidate clustering | phase04_convergence | DONE | 46 sites clustered |
| §13 Phase 9 — Cryptic pocket | 09_cryptic_pocket | PARTIAL | APO only; no open/desensitized/ACh comparison |
| §14 Phase 10 — Site-directed AF3 | 08_site_directed_af3 | PARTIAL | 2-chain reduced model; needs full pentamer |
| §15 Phase 11 — Standard docking | phase05_docking | DONE | 46 sites, 3 ligands, both stoich |
| §16 Phase 12 — Flexible docking | scripts/phase12_flexible_docking.py | EXISTS | Script exists, not run at scale |
| §17 Phase 13 — All-28 SAR | phase06_fingerprints, 07_sar_model | DONE | IFP + SAR rules |
| §18 Phase 14 — Active/inactive | phase09_specificity | DONE | F1-F9 falsification battery |
| §19 Phase 15 — Potency continuum | phase07_sar_model | DONE | SAR rules derived |
| §20 Phase 16 — High-potency stress | — | MISSING | Need compound 12 focus analysis |
| §21 Phase 17 — MMP analysis | — | MISSING | Need matched molecular pair script |
| §22 Phase 18 — Stereochemistry | phase09_specificity (F2) | PARTIAL | F2 FAILS — L/D not distinguished |
| §23 Phase 19 — Boltz-2 | phase08_boltz_affinity, 11_boltz2_analysis | PARTIAL | runs used α10-monomer input (FALLACIES_AUDIT C7); must re-run on α9α10 pentamer |
| §24 Phase 20 — Boltz-2 SAR | phase08_boltz_affinity | PARTIAL | correlation non-significant (r=0.20, p=0.27) and based on wrong-architecture input |
| §25 Phase 21 — Multi-method convergence | phase14_evidence | DONE | Evidence matrix built |
| §26 Phase 22 — Ternary modeling | phase08_af3, cofolding_study | PARTIAL | AF3 ternary runs exist (18); Boltz-2 ternary inputs fixed (2×α9+3×α10) but not yet re-run |
| §27 Phase 23 — Allosteric coupling | — | MISSING | PAM→TMD pathway not traced |
| §28 Phase 24 — Stoichiometry comparison | phase12_subtype | PARTIAL | Compares subtypes, not stoichiometries |
| §29 Phase 25 — State dependence | scripts/phase25_state_dependence.py | PARTIAL | Only resting data available |
| §30 Phase 26 — Site competition | phase14_evidence | DONE | Evidence matrix with weights |
| §31 Phase 27 — Falsification | phase09_specificity | DONE | F1-F9 battery: F2, F9 FAIL (7/9 PASS); F5 partial |
| §32 Phase 28 — Negative controls | — | MISSING | Need random-pocket comparison |
| §33 Phase 29 — Model freeze | 16_model_freeze | EXISTS | Template only |
| §34 Phase 30 — Prospective design | phase16_molgen, 17_cascade | DONE | Library + filter |
| §35 Phase 31 — Prospective screening | phase18_boltz_design | DONE | Design ranking |
| §36 Phase 32 — Blind test | scripts/phase32_blind_test.py | EXISTS | Script exists |
| §37 Phase 33 — Experimental validation | phase21_validation | DONE | Plan exists |
| §38 Final model | scripts/phase_q_t_final.py | EXISTS | Compilation script |
| §39 Required outputs | — | 10/27 DONE | See output tracking below |
| §40 Required figures | figures/ | 0/15 | None generated |
| §41 Required tables | tables/ | Partial | Some tables exist |

---

## CURRENTLY RUNNING / QUEUED JOBS

### Active SLURM Jobs (as of Sep 6, 2026)

| Job ID | Name | Node | Time Left | Purpose |
|--------|------|------|-----------|---------|
| 3327733 | a9a10_md100 | (ended) | time out | MD simulation (apo 2to3) — TIMED OUT, no DCD |
| 3327732 | a9a10_md100 | (ended) | time out | MD simulation (ach 2to3) — TIMED OUT, no DCD |
| 3332533 | boltz2_cofold | pax105 | 22h | Boltz-2 cofolding study |
| 3327734 | a9a10_md100 | — | — | MD simulation (queued, never ran); 0 DCD exist |

## ACTIVE SLURM JOBS (as of Sep 13, 2026)

| Job ID | Name | Purpose |
|--------|------|---------|
| 3633143 | sdcofold (Boltz-2 step3) | full_panel array (180 YAMLs, %2) — still RUNNING; validation run candidates to `afterok:3633143` |

## VALIDATION MODULE (22_binding_site_validation) — READY

All 14 scripts fixed (2026-09-13) and smoke-tested against live data; the seven rewritten
scripts (05,06,09,10,11,13,14) were de-risked by a trial run — all pass, including a fixed
tuple-key JSON bug in 09. Current interim classification: **B. Structurally supported (SiteAF3 pending)**.

**Run AFTER 3633143 completes:**
```bash
cd /cluster/home/nbhatt04/lean_pipeline/22_binding_site_validation
sbatch --dependency=afterok:3633143 submit_validation.sh
```
- STEP 1 = CPU analysis (`run_all.sh`, scripts 01-14, ~10-15 min). Regenerates all
  `outputs/*.json` + `outputs/14_validation_report.md` on the FINAL 180-condition panel.
- STEP 2 = SiteAF3 GPU array (80 configs). Gated on `SITEAF3_AF3_ENV`/`SITEAF3_MODEL_DIR`.
  **NOT configured** — AF3 conda env + model weights + MSA DBs are absent on this account,
  so STEP 2 prints a WARN and exits; report stays at class B until SiteAF3 runs.

### Queued: Site-Directed Cofolding (starts Sep 9)

| Step | Job Name | YAMLs | GPUs | Est. Duration | Begins |
|------|----------|-------|------|---------------|--------|
| 1 | sdcofold_step1 | 6 | 1 | ~3h | Sep 9 00:00 |
| 2 | sdcofold_step2 | 160 | 4 array | ~20h | After step 1 |
| 3 | sdcofold_step3 | 180 | 4 array | ~20h | After step 2 |

**Submit command:**
```bash
cd /cluster/home/nbhatt04/lean_pipeline/site_directed_cofolding/slurm_scripts
sbatch master_submit.sh
```

---

## REMAINING TASKS BY PRIORITY

### PRIORITY 1 — CRITICAL (blocks multiple downstream phases)

#### T1.1: Run Open/Desensitized Receptor State Predictions
- **What:** Boltz-2 cofolding of pentameric receptor in open and desensitized conformations
- **Why:** Enables cryptic pocket analysis, state-dependent binding, ternary modeling
- **Input:** `site_directed_cofolding/yaml_inputs/state_models/` (4 YAMLs)
- **Script:** `site_directed_cofolding/slurm_scripts/step1_state_models.sh`
- **Output:** `site_directed_cofolding/results/state_models/`
- **GPU time:** ~3h
- **Status:** READY TO SUBMIT

#### T1.2: Run ACh-Bound State Predictions
- **What:** Boltz-2 cofolding of receptor + ACh for both stoichiometries
- **Why:** Enables ternary modeling (R+ACh+PAM), ACh dependence analysis
- **Input:** `site_directed_cofolding/yaml_inputs/ach_states/` (2 YAMLs)
- **Script:** `site_directed_cofolding/slurm_scripts/step1_state_models.sh` (same script)
- **Output:** `site_directed_cofolding/results/ach_states/`
- **GPU time:** ~2h
- **Status:** READY TO SUBMIT

#### T1.3: Run Site-Directed Binary Cofolding
- **What:** Boltz-2 cofolding of receptor + PAM at 5 candidate sites, 8 compounds, 2 stoichiometries
- **Why:** Tests whether independently discovered sites support stable ligand placement
- **Input:** `site_directed_cofolding/yaml_inputs/binary/` (80 YAMLs)
- **Script:** `site_directed_cofolding/slurm_scripts/step2_site_directed.sh`
- **Output:** `site_directed_cofolding/results/binary/`
- **GPU time:** ~20h (4 GPUs parallel)
- **Status:** READY TO SUBMIT

#### T1.4: Run Site-Directed Ternary Cofolding
- **What:** Boltz-2 cofolding of receptor + ACh + PAM at 5 candidate sites
- **Why:** Tests ACh effect on PAM binding (master prompt Phase 22)
- **Input:** `site_directed_cofolding/yaml_inputs/ternary/` (80 YAMLs)
- **Script:** `site_directed_cofolding/slurm_scripts/step2_site_directed.sh` (same script)
- **Output:** `site_directed_cofolding/results/ternary/`
- **GPU time:** ~20h (4 GPUs parallel)
- **Status:** READY TO SUBMIT

#### T1.5: Run Full 30-Compound Panel
- **What:** Boltz-2 cofolding of all 30 compounds at top 3 sites, both stoichiometries
- **Why:** Tests active/inactive discrimination, potency ranking, SAR explanation
- **Input:** `site_directed_cofolding/yaml_inputs/full_panel/` (180 YAMLs)
- **Script:** `site_directed_cofolding/slurm_scripts/step3_full_panel.sh`
- **Output:** `site_directed_cofolding/results/full_panel/`
- **GPU time:** ~20h (4 GPUs parallel)
- **Status:** READY TO SUBMIT

---

### PRIORITY 2 — HIGH (required for complete scientific narrative)

#### T2.1: State Classification (Pore/M2 Analysis)
- **What:** Classify each receptor model as resting/open/desensitized using pore radius, M2 orientation, gate diameter
- **Why:** Master prompt §8 requires state assignment from structure, not ligand occupancy
- **Input:** Receptor PDBs from `site_directed_cofolding/results/state_models/`
- **Script:** NEEDS TO BE WRITTEN — `scripts/analyze_state_classification.py`
- **Output:** `04_state_validation/state_classification_report.json`
- **Dependencies:** T1.1 (state model predictions)
- **Est. effort:** 2-3 hours coding

#### T2.2: Structural QC for All Models
- **What:** Evaluate stereochemistry, clashes, missing residues, backbone geometry for every receptor model
- **Why:** Master prompt §9 requires QC before any analysis
- **Input:** All receptor PDBs from cofolding results
- **Script:** NEEDS TO BE WRITTEN — `scripts/structural_qc_all.py`
- **Output:** `04_state_validation/structural_qc_report.json`
- **Dependencies:** T1.1, T1.2
- **Est. effort:** 2-3 hours coding

#### T2.3: Cryptic Pocket Comparison
- **What:** Compare fpocket results across resting, open, desensitized, ACh-bound states
- **Why:** Master prompt §13 requires cryptic/state-dependent pocket classification
- **Input:** fpocket results from `03_pocket_discovery/` + new state fpocket runs
- **Script:** NEEDS TO BE WRITTEN — `scripts/cryptic_pocket_comparison.py`
- **Output:** `09_cryptic_pocket/cryptic_analysis_report.json`
- **Dependencies:** T1.1 (need open/desensitized structures for fpocket)
- **Est. effort:** 3-4 hours coding

#### T2.4: Ternary Complex Analysis (R vs R+ACh vs R+ACh+PAM)
- **What:** Compare PAM binding in ACh-free vs ACh-bound systems
- **Why:** Master prompt §26 requires ACh effect on PAM pocket geometry
- **Input:** Results from T1.2, T1.3, T1.4
- **Script:** NEEDS TO BE WRITTEN — `scripts/ternary_analysis.py`
- **Output:** `13_ternary_modeling/ternary_report.json`
- **Dependencies:** T1.2, T1.3, T1.4
- **Est. effort:** 4-5 hours coding

#### T2.5: Allosteric Coupling Pathway Analysis
- **What:** Trace PAM site → local residues → ECD/TMD coupling → M2 helix → pore
- **Why:** Master prompt §27 requires mechanistic hypothesis for PAM potentiation
- **Input:** Ternary complex structures from T1.4
- **Script:** NEEDS TO BE WRITTEN — `scripts/allosteric_coupling.py`
- **Output:** `13_ternary_modeling/allosteric_pathway.json`
- **Dependencies:** T1.4
- **Est. effort:** 5-6 hours coding

---

### PRIORITY 3 — MEDIUM (completes the scientific analysis)

#### T3.1: Flexible Docking at All Candidate Sites
- **What:** Re-dock all 30 compounds with receptor flexibility at each site
- **Why:** Master prompt §16 requires rigid vs flexible comparison
- **Input:** Receptor PDBs + compound SMILES
- **Script:** `scripts/phase12_flexible_docking.py` (EXISTS, needs scaling)
- **Output:** Flexible docking results CSV
- **Dependencies:** T1.1 (need multiple state structures)
- **Est. effort:** 1 hour coding + 10h GPU

#### T3.2: High-Potency Stress Test
- **What:** Specifically test whether compound 12 (0.198 μM) binds well at leading site
- **Why:** Master prompt §20 — strongest active is critical stress test
- **Input:** Compound 12 docking/cofolding results from T1.3
- **Script:** NEEDS TO BE WRITTEN — `scripts/high_potency_stress_test.py`
- **Output:** `tables/Table14_stress_test.json`
- **Dependencies:** T1.3
- **Est. effort:** 2 hours coding

#### T3.3: Matched Molecular Pair (MMP) Analysis
- **What:** Compare activity changes for structurally similar compound pairs
- **Why:** Master prompt §21 requires MMP-level SAR explanation
- **Input:** 30 compounds + docking/cofolding results
- **Script:** NEEDS TO BE WRITTEN — `scripts/mmp_analysis.py`
- **Output:** `10_sar_validation/mmp_report.json`
- **Dependencies:** T1.3 (need all 30 compound results)
- **Est. effort:** 3-4 hours coding

#### T3.4: Stereochemistry Deep Analysis
- **What:** Model L-ascorbate vs D-ascorbate binding poses explicitly
- **Why:** Master prompt §22; F2 falsification currently fails
- **Input:** L-ascorbate and D-ascorbate cofolding results
- **Script:** NEEDS TO BE WRITTEN — `scripts/stereochemistry_analysis.py`
- **Output:** `09_specificity/stereochemistry_report.json`
- **Dependencies:** T1.3
- **Est. effort:** 2-3 hours coding

#### T3.5: Candidate-Site Competition Matrix
- **What:** Build the full evidence matrix for all 5 candidate sites
- **Why:** Master prompt §30 requires head-to-head site comparison
- **Input:** All cofolding + docking + SAR results
- **Script:** NEEDS TO BE WRITTEN — `scripts/site_competition_matrix.py`
- **Output:** `14_evidence/site_competition_matrix.csv`
- **Dependencies:** T1.3, T1.4, T2.4
- **Est. effort:** 3-4 hours coding

#### T3.6: Falsification Battery (Expanded)
- **What:** Run all falsification tests on leading site hypothesis
- **Why:** Master prompt §31 requires active attempt to disprove
- **Input:** All results from cofolding study
- **Script:** NEEDS TO BE WRITTEN — `scripts/expanded_falsification.py`
- **Output:** `15_falsification/expanded_falsification_report.json`
- **Dependencies:** T1.3, T1.4, T3.5
- **Est. effort:** 4-5 hours coding

#### T3.7: Negative-Site Controls
- **What:** Compare leading site against random plausible pockets
- **Why:** Master prompt §32 requires defense against post-hoc selection
- **Input:** Random pocket models + cofolding results
- **Script:** NEEDS TO BE WRITTEN — `scripts/negative_site_controls.py`
- **Output:** `15_falsification/negative_controls.json`
- **Dependencies:** T1.3
- **Est. effort:** 3-4 hours coding

---

### PRIORITY 4 — LOW (documentation and deliverables)

#### T4.1: Stoichiometry × State × ACh Analysis
- **What:** Full factorial comparison of all conditions
- **Why:** Master prompt §28
- **Script:** NEEDS TO BE WRITTEN
- **Dependencies:** All T1 tasks

#### T4.2: Model Freeze
- **What:** Lock all parameters before prospective design
- **Why:** Master prompt §33
- **Script:** `16_model_freeze/` (template exists, needs population)

#### T4.3: Prospective Compound Design
- **What:** Design new compounds based on frozen model
- **Why:** Master prompt §34
- **Dependencies:** Model freeze

#### T4.4: Prospective Screening
- **What:** Screen purchasable molecules using frozen model
- **Why:** Master prompt §35
- **Dependencies:** Prospective design

#### T4.5: Blind Prospective Test
- **What:** Make predictions before seeing experimental results
- **Why:** Master prompt §36

#### T4.6: All 27 Required Outputs
- **What:** Generate all deliverables listed in master prompt §39
- **Status:** 10/27 complete
- **Missing outputs:**
  - Output 3: Receptor ensemble report
  - Output 4: State classification report
  - Output 5: ACh occupancy report
  - Output 7: Candidate-site table (expanded)
  - Output 8: Site-directed AF3 report
  - Output 10: Flexible docking report
  - Output 12: Active/inactive discrimination (expanded)
  - Output 13: Potency analysis
  - Output 14: High-potency stress test
  - Output 15: MMP analysis
  - Output 16: Stereochemical analysis
  - Output 18: ACh/PAM ternary analysis
  - Output 19: Stoichiometry × state × ACh analysis
  - Output 20: Allosteric communication analysis
  - Output 21: Candidate-site competition matrix
  - Output 22: Falsification analysis (expanded)
  - Output 23: Frozen model specification
  - Output 24: Prospective compound ranking
  - Output 25: Blind prospective predictions
  - Output 26: Experimental validation plan (updated)
  - Output 27: Final mechanistic model

#### T4.7: All 15 Required Figures
- **What:** Publication-quality figures
- **Status:** 0/15 generated
- **Script:** NEEDS TO BE WRITTEN — `scripts/generate_figures.py`

#### T4.8: All 12 Required Tables
- **What:** Formatted tables for publication
- **Status:** Partial (some exist in `tables/`)
- **Script:** NEEDS TO BE WRITTEN — `scripts/generate_tables.py`

---

## JOB SUBMISSION SEQUENCE

```
NOW (Sep 6):
  ├── Current jobs running (MD + Boltz-2 cofold)
  └── All YAMLs generated ✓

Sep 9 (after current jobs finish):
  ├── sbatch master_submit.sh
  │   ├── Step 1: state_models + ach_states (6 YAMLs, ~3h)
  │   ├── Step 2: binary + ternary (160 YAMLs, ~20h, depends on step 1)
  │   └── Step 3: full_panel (180 YAMLs, ~20h, depends on step 2)
  └── Total: ~43h on 4 GPUs

Sep 11-12 (after cofolding completes):
  ├── Run analysis scripts
  ├── Generate state classification (T2.1)
  ├── Generate structural QC (T2.2)
  └── Run cryptic pocket comparison (T2.3)

Sep 13-14:
  ├── Ternary analysis (T2.4)
  ├── Allosteric coupling (T2.5)
  ├── Flexible docking (T3.1)
  └── High-potency stress test (T3.2)

Sep 15-16:
  ├── MMP analysis (T3.3)
  ├── Stereochemistry analysis (T3.4)
  ├── Site competition matrix (T3.5)
  └── Expanded falsification (T3.6)

Sep 17-18:
  ├── Negative controls (T3.7)
  ├── Model freeze (T4.2)
  ├── All 27 outputs (T4.6)
  ├── All 15 figures (T4.7)
  └── All 12 tables (T4.8)

Sep 19+:
  ├── Prospective design (T4.3)
  ├── Prospective screening (T4.4)
  └── Blind test preparation (T4.5)
```

---

## GPU BUDGET ESTIMATE

| Task | GPUs | Duration | Total GPU-hours |
|------|------|----------|-----------------|
| Site-directed cofolding (steps 1-3) | 4 | 43h | 172 |
| MD simulations (remaining) | 3 | 48h | 144 |
| Flexible docking | 1 | 10h | 10 |
| Additional Boltz-2 (validation) | 1 | 20h | 20 |
| **TOTAL** | | | **346 GPU-hours** |

At 4 GPUs continuous: **~3.6 days**

---

## QUICK REFERENCE COMMANDS

```bash
# Check running jobs
squeue -u nbhatt04

# Submit site-directed cofolding
cd /cluster/home/nbhatt04/lean_pipeline/site_directed_cofolding/slurm_scripts
sbatch master_submit.sh

# Check cofolding progress
python3 /cluster/home/nbhatt04/lean_pipeline/site_directed_cofolding/analysis/status.py

# Run analysis after completion
python3 /cluster/home/nbhatt04/lean_pipeline/site_directed_cofolding/analysis/analyze_results.py

# Check pipeline status
python3 /cluster/home/nbhatt04/lean_pipeline/status.py

# Cancel all site-directed jobs
scancel --name=sdcofold_step*

# Check disk space
du -sh /cluster/home/nbhatt04/lean_pipeline/site_directed_cofolding/results/
du -sh /cluster/scratch/nbhatt04/allostery/
```

---

## KEY FILES AND DIRECTORIES

```
lean_pipeline/
├── site_directed_cofolding/          # NEW — main working directory
│   ├── generate_all_inputs.py        # YAML input generator
│   ├── yaml_inputs/
│   │   ├── state_models/             # 4 YAMLs (open/desensitized × 2 stoich)
│   │   ├── ach_states/               # 2 YAMLs (ACh-bound × 2 stoich)
│   │   ├── binary/                   # 80 YAMLs (5 sites × 8 cpds × 2 stoich)
│   │   ├── ternary/                  # 80 YAMLs (5 sites × 8 cpds × 2 stoich)
│   │   ├── full_panel/               # 180 YAMLs (30 cpds × 3 sites × 2 stoich)
│   │   └── manifest.json
│   ├── results/                      # Output directory (populated by jobs)
│   ├── slurm_scripts/
│   │   ├── master_submit.sh          # Run this to submit all jobs
│   │   ├── step1_state_models.sh
│   │   ├── step2_site_directed.sh
│   │   └── step3_full_panel.sh
│   └── analysis/
│       ├── analyze_results.py        # Post-completion analysis
│       └── status.py                 # Progress check
│
├── REMAINING_JOBS.md                 # THIS FILE
├── modulator-dataset-a9a10.csv       # 30-compound dataset
├── data/allostery/af3_outputs/       # 450 AF3 models (existing)
├── data/allostery/boltz2/outputs/    # 36 Boltz-2 sets (existing)
└── cofolding_study/                  # Previous cofolding (partial)
```

---

## KNOWN ISSUES

1. **F2 stereochemistry falsification FAILS** — L/D-ascorbate not distinguished computationally. This is a known weakness acknowledged in `VALIDATION_REPORT.md`. The stereochemistry deep analysis (T3.4) may resolve this.

2. **Broken data symlinks** — `data/campaign/`, `data/deliverable/`, and several CSV symlinks point to old paths in `alpha9alpha10_campaign/`. These need re-linking or the data copied locally.

3. **SMILES parse errors** — Several ascorbate derivatives with `OOH`, `ONH2`, `OCH3`, `OCF3` groups fail RDKit parsing (logged in `lean_pipeline_3266066.err`). These compounds are excluded from the 30-compound panel.

4. **Cofolding study incomplete** — `cofolding_study/results/boltz2/site5/` is empty. Site 23 only has 4 compounds. This is superseded by the new `site_directed_cofolding/` study.

5. **Alpha10 sequence duplicate** — In `generate_af3_ternary_inputs.py`, ALPHA9_SEQ and ALPHA10_SEQ are identical. This appears to be a copy-paste error; the actual sequences differ (visible in state model YAMLs where chains A/B use one sequence and C/D/E use another).
