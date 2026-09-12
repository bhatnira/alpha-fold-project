#!/bin/bash
# =============================================================================
# Step 1: State Models + ACh-Bound States
# =============================================================================
# Receptor-only predictions for open/desensitized states
# Receptor+ACh predictions for resting state
# 6 YAMLs total, ~2-3 GPU-hours
# =============================================================================

#SBATCH --job-name=sdcofold_step1
#SBATCH --partition=gpu
#SBATCH --gres=gpu:a100:1
#SBATCH --constraint=a100-80G
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --begin=now+1minute
#SBATCH --output=/cluster/home/nbhatt04/lean_pipeline/site_directed_cofolding/logs/step1_%j_%a.out
#SBATCH --error=/cluster/home/nbhatt04/lean_pipeline/site_directed_cofolding/logs/step1_%j_%a.err

set -e

module purge 2>/dev/null || true
module load boltz/2.2.1-gpu
export SINGULARITY_BIND=/cluster/tufts,/cluster/scratch
export APPTAINER_BIND=/cluster/tufts,/cluster/scratch

BOLTZ="boltz"
PIPELINE=/cluster/home/nbhatt04/lean_pipeline
BASE=$PIPELINE/site_directed_cofolding
INPUT_DIR=$BASE/yaml_inputs
RESULTS=$BASE/results
CACHE=/cluster/scratch/nbhatt04/.boltz
LOGDIR=$BASE/logs
mkdir -p $LOGDIR

echo "=========================================="
echo "Step 1: State Models + ACh States"
echo "Started: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
nvidia-smi -L | head -2
echo "=========================================="

# Collect all YAMLs for step 1
YAMLS=($(find $INPUT_DIR/state_models $INPUT_DIR/ach_states -name "*.yaml" -type f | sort))
TOTAL=${#YAMLS[@]}
echo "Total YAMLs: $TOTAL"

for yaml in "${YAMLS[@]}"; do
    name=$(basename "$yaml" .yaml)
    
    # Determine output subdirectory
    if [[ "$yaml" == *"state_models"* ]]; then
        out_subdir="state_models"
    else
        out_subdir="ach_states"
    fi
    
    OUTDIR="$RESULTS/$out_subdir/$name"
    
    # Skip if done
    if [ -d "$OUTDIR" ] && find "$OUTDIR" -name "confidence_*.json" 2>/dev/null | head -1 | grep -q .; then
        echo "  SKIP: $name (already done)"
        continue
    fi
    
    echo ""
    echo "[$(date +%H:%M:%S)] Predicting: $name"
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
echo "Step 1 complete: $(date)"
echo "=========================================="
