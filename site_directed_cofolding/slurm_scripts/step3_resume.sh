#!/bin/bash
# =============================================================================
# Step 3 resume: Full 30-Compound Panel (skips already-done outputs)
# 30 compounds x 3 sites x 2 stoichiometries = 180 YAML inputs
# Lower concurrency (%2) + higher memory (64G) to avoid node memory pressure
# =============================================================================

#SBATCH --job-name=sdcofold_step3
#SBATCH --partition=gpu
#SBATCH --gres=gpu:a100:1
#SBATCH --constraint=a100-80G
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=23:00:00
#SBATCH --output=/cluster/home/nbhatt04/lean_pipeline/site_directed_cofolding/logs/step3_%j_%a.out
#SBATCH --error=/cluster/home/nbhatt04/lean_pipeline/site_directed_cofolding/logs/step3_%j_%a.err
#SBATCH --array=0-35%2

set -e

module purge 2>/dev/null || true
module load boltz/2.2.1-gpu
export SINGULARITY_BIND=/cluster/tufts,/cluster/scratch
export APPTAINER_BIND=/cluster/tufts,/cluster/scratch

PIPELINE=/cluster/home/nbhatt04/lean_pipeline
BASE=$PIPELINE/site_directed_cofolding
INPUT_DIR=$BASE/yaml_inputs
RESULTS=$BASE/results
CACHE=/cluster/scratch/nbhatt04/.boltz

echo "=========================================="
echo "Step 3 RESUME: Full 30-Compound Panel"
echo "Started: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Array ID: $SLURM_ARRAY_TASK_ID"
echo "Node: $(hostname)"
nvidia-smi -L | head -2
echo "=========================================="

# Collect all full_panel YAMLs
ALL_YAMLS=($(find $INPUT_DIR/full_panel -name "*.yaml" -type f | sort))
TOTAL=${#ALL_YAMLS[@]}
echo "Total YAMLs: $TOTAL"

# Each array task processes a batch
BATCH_SIZE=5
START_IDX=$((SLURM_ARRAY_TASK_ID * BATCH_SIZE))
END_IDX=$((START_IDX + BATCH_SIZE))
[ $END_IDX -gt $TOTAL ] && END_IDX=$TOTAL

echo "Processing indices $START_IDX to $((END_IDX - 1))"

for ((i=START_IDX; i<END_IDX; i++)); do
    yaml="${ALL_YAMLS[$i]}"
    name=$(basename "$yaml" .yaml)
    OUTDIR="$RESULTS/full_panel/$name"

    # Skip if done
    if [ -d "$OUTDIR" ] && find "$OUTDIR" -name "confidence_*.json" 2>/dev/null | head -1 | grep -q .; then
        echo "  SKIP: $name (already done)"
        continue
    fi

    echo ""
    echo "[$(date +%H:%M:%S)] ($((i+1))/$TOTAL) Predicting: $name"
    mkdir -p "$OUTDIR"

    boltz predict "$yaml" \
        --out_dir "$OUTDIR" \
        --cache "$CACHE" \
        --model boltz2 \
        --accelerator gpu \
        --devices 1 \
        --recycling_steps 3 \
        --sampling_steps 200 \
        --diffusion_samples 3 \
        --output_format pdb \
        --num_workers 1 \
        --preprocessing-threads 1 \
        2>&1 | tee "$OUTDIR/log.txt" || echo "  WARNING: $name failed"

    echo "  Completed: $(date +%H:%M:%S)"
done

echo ""
echo "=========================================="
echo "Step 3 RESUME batch complete: $(date)"
echo "=========================================="