#!/bin/bash
# =============================================================================
# Boltz-2 Affinity Predictions for All 30 Compounds × 3 Sites × 2 Complex Types
# =============================================================================
# Tests Hypothesis 1 & 2: Binding site sharing and SAR correlation
#
# Total jobs: 30 compounds × 3 sites × 2 complex types = 180 YAML files
# Strategy: Run in batches of 10 YAML files per SLURM job
# =============================================================================

#SBATCH --job-name=boltz2_sar
#SBATCH --partition=gpu
#SBATCH --gres=gpu:a100:1
#SBATCH --constraint=a100-80G
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=/cluster/home/nbhatt04/lean_pipeline/cofolding_study/slurm_scripts/boltz2_sar_%A_%a.out
#SBATCH --error=/cluster/home/nbhatt04/lean_pipeline/cofolding_study/slurm_scripts/boltz2_sar_%A_%a.err
#SBATCH --array=0-17%9

set -e

module purge 2>/dev/null || true
module load boltz/2.2.1-gpu
export SINGULARITY_BIND=/cluster/tufts,/cluster/scratch
export APPTAINER_BIND=/cluster/tufts,/cluster/scratch
BOLTZ="boltz"
BASE_DIR="/cluster/home/nbhatt04/lean_pipeline/cofolding_study"
INPUT_DIR="$BASE_DIR/boltz2_inputs"
RESULTS_DIR="$BASE_DIR/results/boltz2"
CACHE=/cluster/scratch/nbhatt04/.boltz
mkdir -p "$RESULTS_DIR"

echo "=========================================="
echo "Boltz-2 SAR Affinity Predictions"
echo "Started: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Array ID: $SLURM_ARRAY_TASK_ID"
echo "Node: $(hostname)"
echo "GPU: $CUDA_VISIBLE_DEVICES"
echo "=========================================="

# Collect all YAML files
ALL_YAMLS=($(find "$INPUT_DIR" -name "*.yaml" -type f | sort))
TOTAL=${#ALL_YAMLS[@]}

echo "Total YAML files: $TOTAL"

# Process batch for this array task
BATCH_SIZE=10
START_IDX=$((SLURM_ARRAY_TASK_ID * BATCH_SIZE))
END_IDX=$((START_IDX + BATCH_SIZE))
if [ $END_IDX -gt $TOTAL ]; then
    END_IDX=$TOTAL
fi

echo "Processing YAML files $START_IDX to $((END_IDX - 1))"

for ((i=START_IDX; i<END_IDX; i++)); do
    YAML="${ALL_YAMLS[$i]}"
    BASENAME=$(basename "$YAML" .yaml)
    SITE_DIR=$(basename $(dirname "$YAML"))
    OUTDIR="$RESULTS_DIR/$SITE_DIR/$BASENAME"
    
    echo ""
    echo "[$(date +%H:%M:%S)] Processing: $SITE_DIR/$BASENAME"
    
    if [ -d "$OUTDIR" ] && find "$OUTDIR" -name "affinity_*.json" 2>/dev/null | head -1 | grep -q .; then
        echo "  Already completed, skipping..."
        continue
    fi
    
    mkdir -p "$OUTDIR"
    
    # Run Boltz-2 predict with affinity head
    $BOLTZ predict "$YAML" \
        --out_dir "$OUTDIR" \
        --cache "$CACHE" \
        --model boltz2 \
        --accelerator gpu \
        --devices 1 \
        --diffusion_samples 1 \
        --output_format pdb \
        --override 2>&1 | tee "$OUTDIR/log.txt"
    
    echo "  Completed: $(date +%H:%M:%S)"
done

echo ""
echo "=========================================="
echo "Batch completed: $(date)"
echo "Results in: $RESULTS_DIR"
echo "=========================================="
