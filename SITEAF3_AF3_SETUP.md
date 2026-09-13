# SiteAF3 / AlphaFold3 Setup Guide (Tufts PAX) — READ LATER

Created: 2026-09-13. Companion to `24_binding_site_validation` (04_siteaf3_submit.sh)
and `AGENTS.md`. Purpose: how to run SiteAF3 on the Tufts cluster WITHOUT downloading
the ~630 GB AF3 databases or the ~1 GB weights from Google.

---

## 1. KEY FACT: Tufts already hosts a shared, full AlphaFold3 install

The cluster maintains a cluster-wide AF3 database + weights set (world-readable).
Your earlier AF3 jobs (the 468 models in `allostery/af3_outputs`, jobs af3_ascorbic /
af3v2 / af3cofold / af3_a9a10, ~Sep 1 2026) almost certainly used exactly these paths.

| Component | Shared path |
|---|---|
| AF3 databases (~630 GB) | `/cluster/tufts/biocontainers/datasets/alphafold3/20241219/public_databases` |
| AF3 model weights (`af3.bin.zst`, ~1 GB) | `/cluster/tufts/biocontainers/datasets/alphafold3/20241219/models/` |
| AF3 container (source) | `module load alphafold3/3.0.3` → `run_alphafold.py` (singularity wrapper) |
| AF2 shared DB (NOT AF3, not needed here) | `/cluster/tufts/biocontainers/datasets/alphafold/db_20231031/` |

Verified contents of the shared DB dir (`public_databases/`):
`bfd-first_non_consensus_sequences.fasta`, `mgy_clusters_2022_05.fa`,
`uniref90_2022_05.fa`, `uniprot_all_2021_04.fa`, `mmcif_files/`,
`pdb_seqres_2022_09_28.fasta`, `nt_rna_2023_02_23_*`, `rfam_14_9_*`,
`rnacentral_active_seq_id_90_cov_80_linclust.fasta` — i.e. the full standard AF3 set.

Tufts demo script (their own, worked):
`/cluster/tufts/biocontainers/tests/alphafold3/alphafold3_demo.sh`

```bash
module load alphafold3/3.0.1
DB="/cluster/tufts/biocontainers/datasets/alphafold3/20241219/public_databases"
run_alphafold.py --output_dir=OUT --json_path=IN.json \
  --db_dir="$DB" \
  --model_dir="/cluster/tufts/biocontainers/datasets/alphafold3/20241219/models/" \
  --max_template_date="2021-09-30" \
  --run_data_pipeline="true" --run_inference="true"
```

## 2. What this means for disk quota

- NO download of the ~630 GB databases (already on cluster).
- NO Google DeepMind form for weights (already on cluster; still respect the
  AlphaFold3 Model Parameters Terms of Use — per-user terms apply anyway).
- The ONLY thing you still need is a small AF3 **python environment** for SiteAF3
  (~5–10 GB conda, fits nowhere near pre-2026 quotas — free space first).

Current quota situation: home 3.3 GB free (lean_pipeline = 23 GB, effectively nothing
else to remove), scratch 17.4 GB free after deleting the two pushed-to-GitHub repos
(HAIL-Poly 6.7 GB, Low-data ligand-based 9.0 GB). Build the AF3 env in `/cluster/scratch/nbhatt04/conda`
and write all SiteAF3 output to `/cluster/scratch/nbhatt04/siteaf3_results` — NOT home.
If the env + 80-config MSA panel exceed ~15 GB, ask TTS Research Storage for a scratch bump.

## 3. What is still missing for SiteAF3 (and only this)

Tufts installs AF3 as a SINGULARITY CONTAINER WRAPPER. There is no AF3 python
environment on the cluster that SiteAF3 can `import` from. SiteAF3 therefore needs:

1. AF3 open-source source code (Apache-2.0, no license gate):
   `git clone https://github.com/google-deepmind/alphafold3.git`
   then `git checkout 7a4a2f7` (the version SiteAF3 was built against), OR patch the
   current version (rename `.cached_ccd()` → `.Ccd()` in
   `./src/embeddings/embed_cond.py` and `./src/diffusion/run_cond_Diff.py`).
2. Replace `<AF3_repo>/src/alphafold3/model/model.py` with the SiteAF3-patched
   `site_directed_cofolding/SiteAF3_patched_AF3_model.py` (already present in
   `/cluster/home/nbhatt04/lean_pipeline/site_directed_cofolding/SiteAFsite_directed_cofolding/SiteAF3_patched_AF3_model.py`).
3. Build the conda env:
   ```bash
   conda activate <new_env>
   conda env update --file <AF3_repo>/environment.yml
   ```
4. Run SiteAF3 pointing at SHARED cluster paths — do NOT use home-dir defaults:

   ```bash
   export SITEAF3_AF3_ENV=/cluster/scratch/nbhatt04/conda/envs/siteaf3   # build in scratch, not home
   export SITEAF3_MODEL_DIR=/cluster/tufts/biocontainers/datasets/alphafold3/20241219/models
   export SITEAF3_DB_DIR=/cluster/tufts/biocontainers/datasets/alphafold3/20241219/public_databases
   # OPTIONAL: where outputs + MSA intermediates go (default already the scratch dir below)
   export SITEAF3_OUTPUT_DIR=/cluster/scratch/nbhatt04/siteaf3_results
   cd /cluster/home/nbhatt04/lean_pipeline/22_binding_site_validation
   bash 04_siteaf3_submit.sh            # or: sbatch submit_validation.sh (already sets these when env vars exist)
   ```

`04_siteaf3_submit.sh` overrides ({SITEAF3_MODEL_DIR:-/cluster/home/nbhatt04/models})
were the blocker; the shared paths above fix it. Defaults in
`phase2_siteaf3_validation/scripts/{setup,submit}_siteaf3.sh` still point at
`/cluster/home/nbhatt04/models` + `/cluster/home/nbhatt04/public_databases` — update
those too if reusing that older module.

## 4. Quick verification commands

```bash
ls /cluster/tufts/biocontainers/datasets/alphafold3/20241219/models/          # af3.bin.zst
ls /cluster/tufts/biocontainers/datasets/alphafold3/20241219/public_databases/ # the 9 DB files
module avail alphafold3                                                        # 3.0.0..3.0.3
cat /cluster/tufts/biocontainers/tests/alphafold3/alphafold3_demo.sh           # working demo
quota -s                                                                       # before building env
```

## 5. Research Technology contacts (reaching out later)

- Guides: https://rtguides.it.tufts.edu/bio/apps/alphafold.html (this is AF2 docs;
  AF3 demo lives at `/cluster/tufts/biocontainers/tests/alphafold3/`)
- OnDemand: https://ondemand.pax.tufts.edu → Bioinformatics Apps. NOTE: the OnDemand
  "AlphaFold" app is AlphaFold **2** with the shared `db_20231031` DB — NOT AF3.
- Support email: tts-research@tufts.edu