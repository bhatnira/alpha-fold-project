# AGENTS.md — α9α10 nAChR PAM binding-site validation (lean_pipeline)

Read this note fully at the start of every session. Full detail: `FALLACIES_AUDIT.md` (same directory).

## Standing warnings (negative findings, treated as OPEN / NOT gains)
- "Boltz-2 affinity" in `VALIDATION_REPORT.md` §5.1 is a PROXY (`enrichment × confidence`); 67/94 rows have
  0 real Boltz-2 models; r²=0.87 is circular. Report also self-contradicts ("not yet computed"). NOT a gain.
- "L-ascorbate 2–3× enrichment" is false (max ratio 1.88× at site 7, 1.33× at site 5, 1.125× at site 23;
  reversed at 11 sites) and contradicts F2. L-ASC max enrichment is at Site 2 (3.66), not Site 23.
- Falsification battery is F1–F9 (9 tests): F2 and F9 FAIL (7/9). F5 is only PARTIAL (sites 5 p=0.040 and
  35 p=0.012 show seed dependence). `test_results.json` "25/25 PASS" is a 0.6 s smoke suite predating this.
- Top site contradiction: `VALIDATION_REPORT.md` says Site 35 (itself labelled LOW_CONFIDENCE);
  `Table10_final_ranking_v2.csv` says Site 23 (0.664 HIGH_CONFIDENCE).
- "Held-out" set = sites, not compounds; held-out sites leak into Table6/Table10 (incl. top-5 pairs sites 1, 29).
- Site 23: Table10 docking_n_models=0/WEAK vs `VALIDATION_REPORT_SAR.md` "docked 7 actives r=0.437, 7/7 bind".
- Boltz-2 "affinity" SAR was run on an α10 MONOMER (single chain A), not the α9α10 pentamer (C7), and the
  three "sites" had byte-identical inputs (C8). Old `cofolding_study/results/boltz2/` numbers are INVALID.
- XAI permutation importances were on an overfit model (train AUROC 0.876 vs CV 0.64); all feature p-values n.s.
- O-ETHYL ascorbate mislabelled `active` (stereochemistry permutation artifact).
- Site IDs renumbered inconsistently across files (rank vs real ID). Never join site_id across files
  without a verified mapping.
- SAR records = 94 (not 95). 5 pharmacophore features (not 51).

## Sequences (VERIFIED CORRECT, do not "fix")
- α9 = MNWSH… (UniProt Q9UGM1, 479 aa). α10 = MGLRSH… (UniProt Q9GZZ6, 450 aa).
- Both confirmed in the 346 `site_directed_cofolding/yaml_inputs/**/*.yaml` (A/B=α9, C/D/E=α10) and all 18
  analyzed AF3 ternary CIFs. `Table4_af3_results.csv` uses these correct outputs.

## Known code bugs — FIXED 2026-09-10 (verify before re-running, don't regress)
- `generate_af3_ternary_inputs.py`: ALPHA9_SEQ had been set to the α10 sequence. Now correct + regenerated
  34 `af3_inputs/*_data.json` (verified heteromer).
- `generate_boltz2_yaml.py`: had emitted a single α10 chain. Now builds 2×α9+3×α10; all 180
  `boltz2_inputs/site*/` regenerated. Boltz-2 affinity results MUST be re-run (old ones invalid).
- `site_directed_cofolding/slurm_scripts/step{1,2,3}_state_models.sh`: CACHE was `/cluster/home/nbhatt04/.boltz`
  (symlink) → caused `FileExistsError` and 0 predictions. Now `CACHE=/cluster/scratch/nbhatt04/.boltz`.
- README.md MD row now states MD data does not exist. `data/allostery` symlink is valid (do not "repair").

## Still open / blocked
- F2, F9 fail; affinity-matrix proxy still in §5.1 (annotated only); site-ID map missing; MD never ran (0 DCD).
- Boltz-2 affinity RE-RUN launched 2026-09-10 (job 3525810, 18-task array, corrected α9α10 heteropentamer
  YAMLs, warm cache, module boltz/2.2.1-gpu). Old (invalid, α10-monomer) outputs archived to
  `cofolding_study/results/boltz2_invalid_monomer_20260910/`; fresh `results/boltz2/` now being filled.
- Cofolding resubmission launched 2026-09-10 after CACHE fix: step1=3525813 (6 state/ach YAMLs) →
  step2=3525814 (160 binary+ternary, dep on step1) → step3=3525815 (180 full_panel, dep on step2).
- LBM repo: commits `196a29f`/`be14e54` need push + flip-to-private; user chose SKIP (no GitHub PAT).
- GPU wall-clock ceiling is 2 days; each cofolding step requests ≤23 h, affinity tasks 4 h.
- Do NOT mark items fixed until re-verified. Do NOT treat the flaws above as accomplishments/gains.
- This file + FALLACIES_AUDIT.md are standing records; keep them accurate when work changes anything.

## Validation module (22_binding_site_validation) — 2026-09-13
- Scripts 01-14 fixed + verified end-to-end (ran on live partial panel; all pass, incl. a tuple-key
  JSON bug fixed in 09). Interim classification: **B. Structurally supported (SiteAF3 pending)**.
- SiteAF3 (04_siteaf3_submit.sh) RESOLVED for DB/weights 2026-09-13: cluster **already hosts a full
  shared AF3 install** — DBs `/cluster/tufts/biocontainers/datasets/alphafold3/20241219/public_databases`
  (~630 GB), weights same dir `/models/af3.bin.zst`. NO download / NO TB quota needed. Only remaining
  SiteAF3 blocker = the AF3 python env itself (see `SITEAF3_AF3_SETUP.md` in this dir). Use env vars
  SITEAF3_MODEL_DIR / SITEAF3_DB_DIR / SITEAF3_AF3_ENV / SITEAF3_OUTPUT_DIR when submitting.
  Default output = `/cluster/scratch/nbhatt04/siteaf3_results`. Do NOT use the old defaults
  `/cluster/home/nbhatt04/models` or `/cluster/home/nbhatt04/public_databases`.
  Build AF3 env in `/cluster/scratch/nbhatt04/conda`, NOT home.
- Final full-panel re-run auto-fires on `afterok:3633143` via `sbatch --dependency=afterok:3633143
  submit_validation.sh` (job 3645926 submitted 2026-09-13). Post-array, full_panel = 180 conditions;
  report tier-1 = 5 residues (alpha9.120/176/224/82/175), discordant-inactive = 23 cpds/126 conditions.