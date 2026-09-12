#!/bin/bash
# =============================================================================
# SiteAF3 Submission Script
# =============================================================================
# Runs SiteAF3 for 5 sites × 8 compounds × 2 stoichiometries = 80 predictions
# Uses SLURM array with 4 GPUs
# =============================================================================

#SBATCH --job-name=siteaf3
#SBATCH --partition=gpu
#SBATCH --gres=gpu:a100:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=23:00:00
#SBATCH --output=/cluster/home/nbhatt04/lean_pipeline/phase2_siteaf3_validation/logs/siteaf3_%j_%a.out
#SBATCH --error=/cluster/home/nbhatt04/lean_pipeline/phase2_siteaf3_validation/logs/siteaf3_%j_%a.err
#SBATCH --array=0-19%4

set -e

module purge 2>/dev/null || true
module load alphafold3/3.0.3 2>/dev/null || true

PIPELINE=/cluster/home/nbhatt04/lean_pipeline
PHASE2=$PIPELINE/phase2_siteaf3_validation
INPUT_DIR=$PHASE2/siteaf3_inputs
RESULTS=$PHASE2/results/siteaf3
SITEAF3=$PHASE2/SiteAF3
MODEL_DIR=/cluster/home/nbhatt04/models
DB_DIR=/cluster/home/nbhatt04/public_databases
LOGS=$PHASE2/logs
mkdir -p $RESULTS $LOGS

echo "=========================================="
echo "SiteAF3 Predictions"
echo "Started: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Array ID: $SLURM_ARRAY_TASK_ID"
echo "Node: $(hostname)"
nvidia-smi -L | head -2
echo "=========================================="

# Collect all JSON configs
ALL_JSONS=($(find "$INPUT_DIR" -name "*.json" -type f | sort))
TOTAL=${#ALL_JSONS[@]}
echo "Total JSON configs: $TOTAL"

# Each array task processes a batch of 4
BATCH_SIZE=4
START_IDX=$((SLURM_ARRAY_TASK_ID * BATCH_SIZE))
END_IDX=$((START_IDX + BATCH_SIZE))
[ $END_IDX -gt $TOTAL ] && END_IDX=$TOTAL

echo "Processing indices $START_IDX to $((END_IDX - 1))"

for ((i=START_IDX; i<END_IDX; i++)); do
    json="${ALL_JSONS[$i]}"
    name=$(basename "$json" .json)
    OUTDIR="$RESULTS/$name"

    # Skip if done
    if [ -d "$OUTDIR" ] && [ -f "$OUTDIR/"*/confidences.json ]; then
        echo "  SKIP: $name (already done)"
        continue
    fi

    echo ""
    echo "[$(date +%H:%M:%S)] ($((i+1))/$TOTAL) Running: $name"
    mkdir -p "$OUTDIR"

    python3 "$SITEAF3/run_SiteAF3.py" \
        --config_file "$json" \
        --output_dir "$OUTDIR" \
        --model_weights_dir "$MODEL_DIR" \
        --db_dir "$DB_DIR" \
        --receptor_type protein \
        --ligand_type small_molecule \
        --hotspot_cutoff 8.0 \
        --pocket_cutoff 10.0 \
        --use_pocket_masked_af3_msa_for_embedding \
        --use_pocket_diffusion_for_prediction \
        --verbose \
        2>&1 | tee "$OUTDIR/log.txt" || echo "  WARNING: $name failed"

    echo "  Completed: $(date +%H:%M:%S)"
done

echo ""
echo "=========================================="
echo "SiteAF3 batch complete: $(date)"
echo "=========================================="
