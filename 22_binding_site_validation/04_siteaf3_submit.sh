#!/bin/bash
# =============================================================================
# 4. SiteAF3 SLURM submission (site-conditioned candidate-site validation)
# =============================================================================
# Launches SiteAF3 for the candidate-site x stoichiometry x compound panel.
# Gated on completion of the sdcofold_step3 Boltz-2 array (--dependency).
#
# PREREQUISITES (must be verified before the first run):
#   * SiteAF3 python env (SiteAF3_env.yml) with AlphaFold3 3.0.3 source patched
#     per SiteAF3 README (AF3_code/model.py) -- currently NOT installed.
#   * AlphaFold3 model weights dir + MSA database dir.
# Tufts hosts shared AF3 DBs + weights (see SITEAF3_AF3_SETUP.md):
#   MODEL=/cluster/tufts/biocontainers/datasets/alphafold3/20241219/models
#   DB=/cluster/tufts/biocontainers/datasets/alphafold3/20241219/public_databases
# Override paths via env vars:
#   SITEAF3_AF3_ENV, SITEAF3_MODEL_DIR, SITEAF3_DB_DIR, SITEAF3_OUTPUT_DIR
# =============================================================================

#SBATCH --job-name=siteaf3_validation
#SBATCH --partition=gpu
#SBATCH --gres=gpu:a100:1
#SBATCH --constraint=a100-80G
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=23:00:00
#SBATCH --output=/cluster/home/nbhatt04/lean_pipeline/22_binding_site_validation/logs/siteaf3_%j_%a.out
#SBATCH --error=/cluster/home/nbhatt04/lean_pipeline/22_binding_site_validation/logs/siteaf3_%j_%a.err
#SBATCH --array=0-19%4

set -euo pipefail

PIPELINE=/cluster/home/nbhatt04/lean_pipeline
PKG=$PIPELINE/22_binding_site_validation
INPUT_DIR=$PKG/siteaf3_inputs
# Output defaults to /cluster/scratch: the 80-config panel + MSA intermediates
# do not fit the home quota. Override via SITEAF3_OUTPUT_DIR if needed.
RESULTS=${SITEAF3_OUTPUT_DIR:-/cluster/scratch/nbhatt04/siteaf3_results}
SITEAF3_REPO=$PIPELINE/site_directed_cofolding/SiteAF3
LOGS=$PKG/logs
mkdir -p "$RESULTS" "$LOGS"

# ---- solve python env for AlphaFold3-based SiteAF3 ----
AF3_ENV=${SITEAF3_AF3_ENV:-}
if [ -n "$AF3_ENV" ] && [ -x "$AF3_ENV/bin/python" ]; then
    PY="$AF3_ENV/bin/python"
elif PY=$(command -v python3); then
    PY="$PY"
else
    echo "ERROR: no AF3 python env found. Set SITEAF3_AF3_ENV." >&2
    exit 1
fi

MODEL_DIR=${SITEAF3_MODEL_DIR:-/cluster/home/nbhatt04/models}
DB_DIR=${SITEAF3_DB_DIR:-/cluster/home/nbhatt04/public_databases}
[ -d "$MODEL_DIR" ] || echo "WARN: model weights dir '$MODEL_DIR' missing"
[ -d "$DB_DIR" ] || echo "WARN: MSA DB dir '$DB_DIR' missing"

echo "Job $SLURM_JOB_ID array $SLURM_ARRAY_TASK_ID ($(hostname))"
nvidia-smi -L | head -2

ALL_JSONS=($(find "$INPUT_DIR" -name '*.json' -type f | sort))
TOTAL=${#ALL_JSONS[@]}
BATCH_SIZE=4
START_IDX=$((SLURM_ARRAY_TASK_ID * BATCH_SIZE))
END_IDX=$((START_IDX + BATCH_SIZE))
[ "$END_IDX" -gt "$TOTAL" ] && END_IDX=$TOTAL

for ((i=START_IDX; i<END_IDX; i++)); do
    json="${ALL_JSONS[$i]}"
    name=$(basename "$json" .json)
    outdir="$RESULTS/$name"
    if [ -f "$outdir"/seed_*/sample_*/model.cif ]; then
        echo "  SKIP $name (done)"
        continue
    fi
    mkdir -p "$outdir"
    echo "[$(date +%H:%M:%S)] ($((i+1))/$TOTAL) $name"
    "$PY" "$SITEAF3_REPO/run_SiteAF3.py" \
        --config_file "$json" \
        --output_dir "$outdir" \
        --model_weights_dir "$MODEL_DIR" \
        --db_dir "$DB_DIR" \
        --receptor_type protein \
        --ligand_type small_molecule \
        --use_pocket_masked_af3_msa_for_embedding \
        --use_pocket_diffusion_for_prediction \
        --verbose 2>&1 | tee "$outdir/log.txt" \
        || echo "  WARNING: $name failed"
done

echo "done $(date)"