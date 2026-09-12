#!/bin/bash
#SBATCH --job-name=a9a10_full_pipeline
#SBATCH --partition=gpu
#SBATCH --gres=gpu:a100:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=/cluster/home/nbhatt04/lean_pipeline/logs/full_pipeline_%j.out
#SBATCH --error=/cluster/home/nbhatt04/lean_pipeline/logs/full_pipeline_%j.err
#SBATCH --array=1-6

set -e

PIPELINE_DIR="/cluster/home/nbhatt04/lean_pipeline"
PYTHON="/usr/bin/python3"
SCRATCH="/cluster/scratch/nbhatt04/allostery"
LOGDIR="$PIPELINE_DIR/logs"
mkdir -p "$LOGDIR" "$PIPELINE_DIR/figures"

echo "=========================================="
echo "Alpha9alpha10 Full Pipeline"
echo "Started: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Array ID: $SLURM_ARRAY_TASK_ID"
echo "Node: $(hostname)"
echo "GPU: $CUDA_VISIBLE_DEVICES"
echo "=========================================="

# Ensure packages are available
export PYTHONPATH="$HOME/.local/lib/python3.9/site-packages:$PYTHONPATH"

cd "$PIPELINE_DIR"

case $SLURM_ARRAY_TASK_ID in
    1)
        echo ""
        echo "[$(date +%H:%M:%S)] TASK 1: Full Docking (all 30 compounds)"
        $PYTHON scripts/full_docking.py
        ;;
    2)
        echo ""
        echo "[$(date +%H:%M:%S)] TASK 2: Interaction Fingerprints (all compounds)"
        $PYTHON scripts/full_ifp.py
        ;;
    3)
        echo ""
        echo "[$(date +%H:%M:%S)] TASK 3: SAR Validation (all 30 compounds)"
        $PYTHON scripts/full_sar_validation.py
        ;;
    4)
        echo ""
        echo "[$(date +%H:%M:%S)] TASK 4: Boltz-2 Analysis (existing outputs)"
        $PYTHON scripts/full_boltz2_analysis.py
        ;;
    5)
        echo ""
        echo "[$(date +%H:%M:%S)] TASK 5: Convergence & Evidence Matrix"
        $PYTHON scripts/full_convergence.py
        ;;
    6)
        echo ""
        echo "[$(date +%H:%M:%S)] TASK 6: Figure Generation"
        $PYTHON scripts/generate_figures.py
        ;;
esac

echo ""
echo "=========================================="
echo "Task Complete: $(date)"
echo "=========================================="
