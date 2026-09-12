#!/bin/bash
#SBATCH --job-name=a9a10_md_gpu
#SBATCH --partition=gpu
#SBATCH --gres=gpu:a100:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --exclude=pax007
#SBATCH --output=/cluster/home/nbhatt04/lean_pipeline/logs/md_gpu_%j.log
#SBATCH --error=/cluster/home/nbhatt04/lean_pipeline/logs/md_gpu_%j.err

set -e

export OPENMM_PLUGIN_DIR=/cluster/home/nbhatt04/.conda/envs/md_cuda/lib/plugins
PYTHON=/cluster/home/nbhatt04/.conda/envs/md_cuda/bin/python

echo "=========================================="
echo "OpenMM GPU MD - alpha9/alpha10 interface"
echo "Started: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
nvidia-smi -L | head -2
echo "=========================================="

mkdir -p /cluster/home/nbhatt04/lean_pipeline/md_gpu
cd /cluster/home/nbhatt04/lean_pipeline/md_gpu

$PYTHON /cluster/home/nbhatt04/lean_pipeline/scripts/md_gpu_openmm.py

echo ""
echo "=========================================="
echo "MD complete: $(date)"
echo "=========================================="
