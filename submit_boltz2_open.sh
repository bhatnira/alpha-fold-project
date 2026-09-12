#!/bin/bash
# Run Boltz-2 predictions for open state templates
#SBATCH --job-name=boltz2_open
#SBATCH --partition=gpu
#SBATCH --gres=gpu:a100:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=boltz2_open_%j.out
#SBATCH --error=boltz2_open_%j.err

set -e

BOLTZ="/cluster/home/nbhatt04/.venvs/pymol/bin/boltz"
OUTDIR="/cluster/home/nbhatt04/lean_pipeline/02_receptor_ensemble/open_states/predictions"
mkdir -p "$OUTDIR"

echo "=========================================="
echo "Boltz-2 Open State Predictions"
echo "Started: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "GPU: $CUDA_VISIBLE_DEVICES"
echo "=========================================="

# Prediction 1: Open state (apo)
echo ""
echo "[$(date +%H:%M:%S)] Prediction 1: Open state (apo)"
$BOLTZ predict /cluster/home/nbhatt04/lean_pipeline/02_receptor_ensemble/open_states/boltz2_2to3_open.yaml \
    --out_dir "$OUTDIR/open_state" \
    --model boltz2 \
    --devices 1 \
    --diffusion_samples 1 \
    --output_format pdb \
    --override

# Prediction 2: Ternary complex (open state + PAM)
echo ""
echo "[$(date +%H:%M:%S)] Prediction 2: Ternary complex (open state + PAM)"
$BOLTZ predict /cluster/home/nbhatt04/lean_pipeline/02_receptor_ensemble/open_states/boltz2_2to3_ternary.yaml \
    --out_dir "$OUTDIR/ternary_complex" \
    --model boltz2 \
    --devices 1 \
    --diffusion_samples 1 \
    --output_format pdb \
    --override

echo ""
echo "[$(date +%H:%M:%S)] Predictions complete"
echo "Results in: $OUTDIR/"
echo ""
ls -la "$OUTDIR/"
echo ""
echo "Completed: $(date)"
