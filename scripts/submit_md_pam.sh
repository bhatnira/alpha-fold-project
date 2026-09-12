#!/bin/bash
#SBATCH --job-name=md_pam
#SBATCH --partition=gpu
#SBATCH --gres=gpu:a100:4
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
#SBATCH --exclude=pax007
#SBATCH --output=/cluster/home/nbhatt04/lean_pipeline/logs/md_pam_%j.log
#SBATCH --error=/cluster/home/nbhatt04/lean_pipeline/logs/md_pam_%j.err

set -e

export OPENMM_PLUGIN_DIR=/cluster/home/nbhatt04/.conda/envs/md_cuda/lib/plugins
PYTHON=/cluster/home/nbhatt04/.conda/envs/md_cuda/bin/python

SIM=${SIM:-ascorbate}

echo "=========================================="
echo "PAM Ligand MD: $SIM"
echo "Started: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
nvidia-smi -L
echo "=========================================="

mkdir -p /cluster/home/nbhatt04/lean_pipeline/md_pam
cd /cluster/home/nbhatt04/lean_pipeline/md_pam

$PYTHON /cluster/home/nbhatt04/lean_pipeline/scripts/md_pam_runner.py --ligand $SIM --n_gpus 4

echo ""
echo "=========================================="
echo "Complete: $(date)"
echo "=========================================="
