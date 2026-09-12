# FALLACIES / FLAWS AUDIT — α9α10 nAChR PAM binding-site validation

Status: OPEN. These are NEGATIVE findings, not gains. They remain unresolved until each is explicitly
fixed and re-verified. Audited: 2026-09-10. Second full cross-check: 2026-09-10 (all claims re-checked
numerically; two prior findings were WRONG and are struck below).

> Bottom line: the only defensible conclusions are (a) F1/F3/F4/F8 pass, F5 is a PARTIAL pass (C9),
> F2 and F9 FAIL; (b) 34/46 sites preserve L>D and 38/46 L>acetate at the *data* level (not positivity);
> (c) BINDING_SITE_VALIDATION_REPORT.txt downgrades Site 23 to a "Class C hypothesis".
> Nothing below may be cited as support until fixed.

---

## CRITICAL

### C1. "Boltz-2 affinity" in the report is a mislabeled proxy + circular correlation  [PARTLY FIXED]
- `VALIDATION_REPORT.md:135-139` (§5.1 "Boltz-2 Affinity Predictions", column "Affinity (log10 IC50)") uses
  `affinity_matrix.csv`, where `affinity_score = enrichment × confidence` (a proxy, stated in `affinity_summary.txt`).
  67/94 rows have `n_boltz2_models=0` (e.g. RYANODINE affinity 7.1484 with 0 Boltz-2 models).
- The claimed "Ryanodine enrichment–affinity r²=0.87" is circular — affinity is built FROM enrichment.
- The report itself contradicts this: lines 225/234 say "Boltz-2 affinity predictions not yet computed for
  L-ascorbate" while §5.1 presents them.
- Fix (annotate the affected report sections) — see file-correction note below. Verified numerically:
  every `affinity_score` value equals `confidence × enrichment` to CSV rounding.

### C2. "L-ascorbate enriches 2–3× more" claim contradicts its own data  [FIXED-in-record]
- Table6 L/D ratios: Site 23 = 1.125×, Site 5 = 1.333×, Site 7 = 1.875× (max); REVERSED (L<D) at 11 sites:
  3,8,13,14,15,16,19,21,22,26,34.
- VALIDATION_REPORT.md:42's own example ("Site 5 — L-ascorbate 2.545 vs D-ascorbate 1.909") is a 1.33× ratio, not 2–3×.
- Contradicted by F2 (L/D overlap 0.931, p=0.601) and by the pharmacophore "D-ascorbate binds 3× weaker" (no ratio exceeds 1.88×).
- Also "enrichment strongest at Site 23" is false: L-ASC max enrichment in Table6 is 3.66 at **Site 2**; Site 23 = 1.35.

### C3. Falsification battery miscounted / PASS overstated  [FIXED-in-record]
- Battery is F1–F9 (9 tests). F2 (stereochemistry) FAIL, F9 (ensemble-docking discrimination) FAIL → 7/9, not "7-8 PASS".
- `test_results.json` (Sep 4) says 25/25 PASS; it pre-dates the falsification run and is a 0.6 s smoke suite.

### C4. Two authoritative docs disagree on the top site  [FIXED-in-record]
- `VALIDATION_REPORT.md:18` = Site 35 (weighted 0.683); `Table10_final_ranking_v2.csv` = Site 23 (0.664, HIGH_CONFIDENCE).
- Notable: `VALIDATION_REPORT.md:84` lists top-ranked Site 35 as **LOW_CONFIDENCE**.

### C5. "Held-out" split is contaminated  [FIXED-in-record]
- Split is by SITES, not compounds; held-out sites (1,4,6,11,12,15,24,28,29,36,39,41-44,46) appear in Table6
  and Table10. `sar_summary.txt` top-5 enrichment pairs include held-out sites 1 and 29.
- README's "held-out compounds test the model" claim is unsupported.

### C6. Docking contradiction at Site 23  [FIXED-in-record]
- `Table10` site 23: docking_n_models=0, docking_level=WEAK; `VALIDATION_REPORT_SAR.md:12-18` claims
  "docked 7 actives to Site 23, r=0.437, 7/7 bind". Also the r=0.437 sign is described inconsistently
  (positive vs "opposite direction").

### C7. Boltz-2 "affinity" SAR was run on an α10 MONOMER, not the α9α10 pentamer  [FIXED]
- `cofolding_study/generate_boltz2_yaml.py` used a SINGLE `RECEPTOR_SEQUENCE` = α10 (MGLRSH) as one chain.
  All 180 `boltz2_inputs/site*/` YAMLs were a single α10 chain + ligands (verified: model PDBs have only
  chain A, 3493 ATOM lines). So the "real Boltz-2 affinity" feed (r=0.20, mechanistic SAR) was computed on
  the wrong architecture.
- Only 83/180 affinity JSONs exist (partial completion).
- FIXED: generator now emits 2×α9 + 3×α10; 180 YAMLs regenerated (verified chain A–E pattern).
  NOTE: OLD `results/boltz2/` outputs are INVALID and must be re-run.

### C8. Boltz-2 "site" inputs were identical across sites  [FIXED]
- site21/site23/site5 YAMLs differed ONLY by directory/filename (same chain + same ligands).
  No per-site Boltz-2 affinity can be claimed from these runs. `boltz2_full_sar_results.csv` has no site column.
  FIXED in generator output (identical-input problem remains true for the regenerated files because the
  design has no per-site pocket constraint — that is a modeling decision, now documented).

### C9. F5 (receptor-state artifact) is only a PARTIAL pass  [NEW]
- Within F5, site 5 (p=0.040) and site 35 (p=0.012) DO show seed dependence while the verdict claims
  "no seed dependence detected". Sites 23/34/21/7 are seed-independent. Treat F5 as partial.

---

## MAJOR (provenance / reproducibility)

### M1. Site-ID collision / renumbering across files  [FIXED-in-record]
- `sar_summary.txt` ranks "site 1 composite=0.5066" (rank value; real site = 23), etc. `Table3`/`Table10`
  use real IDs. Cross-file joins unsafe without a mapping. `sar_top10_by_composite_score.csv` now carries
  an explicit `site_id` column (site 23) distinct from the `rank` column — the collision was rank-vs-ID,
  not a literal "site 1" header.

### M2. Record-count mismatch  [FIXED-in-record]
- Actual SAR records = 94 (22 L-ASC + 17 D-ASC + 29 O-ETHYL + 23 ACETATE + 3 RYANODINE), not 95
  (`VALIDATION_REPORT.md` lines 11, 20, 218, 251 all say 95; Table6 has 94 data rows).

### M3. CSV column misalignment  [FIXED]
- `sar_top10_by_composite_score.csv`: 11 header columns vs 12 data columns (extra unlabeled token,
  e.g. `L-ASC,L-ascorbate,active,1.353,...`). Column header corrected to include `ligand_activity_class`.

### M4. XAI computed on the wrong (overfit) model  [FIXED-in-record]
- `xai_results.json` permutation importances use base_auroc=0.8758 (train), papered model = CV AUROC 0.64.
  Every feature p-value non-significant (docking p=0.43, boltz2_probability p=0.45, boltz2_affinity p=0.28).
  No driver-feature claim is supported.

### M5. O-ETHYL ascorbate mislabeled `active`  [FIXED-in-record]
- Root cause: no stereochemistry flag → R/S penicillin/ascorbate enantiomer permutations counted as
  separate molecular species; O-ethyl ascorbate (a permutation artifact of L-ascorbate) was labeled "active".

### M6. "25/25 PASS" smoke suite not robust  [FIXED-in-record]
- 0.6 s total for all 25 phase-smoke tests; not a validation of the falsification results (F2/F9 FAIL).

---

## MINOR / consistency
- N1. Designed analogs include impossible chemistry: MOL_0007 (OCl), MOL_0008 (OBr), MOL_0006 (OF),
  MOL_0004 (OOH) in `Table11_prospective_compounds.csv`.
- N2. README MD row referenced nonexistent `data/allostery/md/` — fixed (now marked NOT AVAILABLE).
- N3. `REMAINING_JOBS.md` stale: "ach 0/2 predicted", "MD 2/4 done", "8-test battery" — all superseded
  (ach run dirs exist but 0 predictions (cache bug, now fixed); MD 0/4 done with 0 DCD; battery is 9 tests).
- N4. "8-test battery, 7/8 PASS" claim was in `REMAINING_JOBS.md` §31 Phase 27 (now corrected to
  F1–F9, F2/F9 FAIL, F5 partial). `BINDING_SITE_VALIDATION_REPORT.txt` has no such miscount; it is
  consistent with 9 tests.
- N5. `pharmacophore_features.csv` holds 5 features (not 51): "3/5 essential / 4/5" wording matches 5 features.

---

## STRUCK-OR-CORRECTED PREVIOUS FINDINGS (my own errors, fixed in this pass)
- ~~data/ symlinks broken / empty~~ → WRONG. `data/allostery -> /cluster/scratch/nbhatt04/allostery`
  exists since Sep 2 and resolves (468 AF3 CIFs, 90 Boltz-2 CIFs verified intact). README paths were valid.
- ~~slurm "Step 1: no state generation"~~ → WRONG. No such message exists anywhere. The actual cofolding
  failure was `FileExistsError: /cluster/home/nbhatt04/.boltz` (Boltz-2 cache.mkdir on a symlink path) in
  `step{1,2,3}_state_models.sh.` FIXED: CACHE now points to `/cluster/scratch/nbhatt04/.boltz` (real dir).
- ~~"α9/α10 sequences wrong" (implied corruption of runs)~~ → CORRECTED: the sequences in the ACTUAL run
  inputs were CORRECT (see VERIFIED-CORRECT). The alpha-9==alpha-10 bug existed ONLY in generator scripts.

---

## VERIFIED-CORRECT (do not flag as flaws)
- α9 sequence = MNWSH… (UniProt Q9UGM1, 479 aa); α10 = MGLRSH… (UniProt Q9GZZ6, 450 aa). Both confirmed
  in: all 346 `site_directed_cofolding/yaml_inputs/**/*.yaml` (A/B=α9, C/D/E=α10) and all 18 analyzed
  AF3 ternary model CIFs (entities 1-2 = α9 length 479, 3-5 = α10 length 450).
- `Table4_af3_results.csv` (7 runs, 175 rows) reads the correct-sequence AF3 outputs.
- F1, F3, F4, F8 pass; `phase09_specificity` 34/46 L>D, 38/46 L>acetate at data level.
- Boltz-2 site-directed cofolding YAML inputs (the 346) are correct-sequence heteropentamers.

## FIXES APPLIED 2026-09-10
1. `cofolding_study/generate_af3_ternary_inputs.py` — ALPHA9_SEQ/ALPHA10_SEQ corrected to verified
   Q9UGM1/Q9GZZ6 sequences; regenerated all 34 `af3_inputs/*_data.json` (verified heteromer).
2. `cofolding_study/generate_boltz2_yaml.py` — now builds 2×α9 + 3×α10 (A,B,C,D,E); regenerated all 180
   `boltz2_inputs/site*/` YAMLs (verified pattern MNWSH,MNWSH,MGLRSH,MGLRSH,MGLRSH). OLD results invalid;
   needs re-run (GPU).
3. `step{1,2,3}_state_models.sh` — CACHE=/cluster/scratch/nbhatt04/.boltz (fixes cache FileExistsError
   that produced 0 prediction outputs on the Sep-9 run).
4. README.md MD row corrected (no MD data exists).
5. Report/record corrections (95→94, 2–3× retraction, F5 partial, site-35 vs site-23, affinity
   proxy/circularity annotation, Boltz-2 monomer-α10 invalidation note appended to VALIDATION_REPORT_SAR.md,
   REMAINING_JOBS.md superseded notes: MD 0/4, ach re-run needed, AF3=468 CIFs, Boltz-2=90 CIFs,
   Phase-27 8-test→F1–F9). `BINDING_SITE_VALIDATION_REPORT.txt` verified consistent (no 8-test count there).
6. `sar_top10_by_composite_score.csv` realigned: 12-column header
   `rank,site_id,site_composite_score,ligand,ligand_name,ligand_activity_class,enrichment,fisher_p,
   mean_confidence,n_contacts,recurrent_contacts,site_evidence_level`; all 10 rows normalized to 12
   columns (2 short RYANODINE rows got a `ligand_name` token inserted). Backup: `/tmp/opencode/sar_top10.backup.csv`.

## STILL OPEN / REQUIRING DATA RE-RUN (not "fixed")
- Boltz-2 affinity re-run on the corrected α9α10 heteropentamer (180 YAMLs, ~4 h GPU) — old results invalid.
- Co-folding re-run (346 YAMLs via `sbatch master_submit.sh`) — blocked until user go-ahead (approx 43 h, 4 GPUs).
- Site-ID map across files; F2/F9 remain failing; affinity-matrix proxy still used in §5.1 (annotated).
- LBM repo push + flip-private (needs GitHub PAT).