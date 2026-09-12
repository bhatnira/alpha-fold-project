#!/bin/bash
# =============================================================================
# Step 2: Site-Directed Binary + Ternary Cofolding
# =============================================================================
# 5 candidate sites x 8 compounds x 2 stoichiometries x 2 complex types
# = 160 YAML inputs
# Uses SLURM array with 4 concurrent GPU tasks
# =============================================================================

#SBATCH --job-name=sdcofold_step2
#SBATCH --partition=gpu
#SBATCH --gres=gpu:a100:1
#SBATCH --constraint=a100-80G
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=23:00:00
#SBATCH --begin=now+1minute
#SBATCH --output=/cluster/home/nbhatt04/lean_pipeline/site_directed_cofolding/logs/step2_%j_%a.out
#SBATCH --error=/cluster/home/nbhatt04/lean_pipeline/site_directed_cofolding/logs/step2_%j_%a.err
#SBATCH --array=0-39%4

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
echo "Step 2: Site-Directed Binary + Ternary"
echo "Started: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Array ID: $SLURM_ARRAY_TASK_ID"
echo "Node: $(hostname)"
nvidia-smi -L | head -2
echo "=========================================="

# Collect all binary + ternary YAMLs
ALL_YAMLS=($(find $INPUT_DIR/binary $INPUT_DIR/ternary -name "*.yaml" -type f | sort))
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
    
    # Determine output subdirectory
    if [[ "$yaml" == *"binary"* ]]; then
        out_subdir="binary"
    else
        out_subdir="ternary"
    fi
    
    OUTDIR="$RESULTS/$out_subdir/$name"
    
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
echo "Step 2 batch complete: $(date)"
echo "=========================================="
