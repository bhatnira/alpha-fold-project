#!/bin/bash
#SBATCH --job-name=a9a10_md100
#SBATCH --partition=gpu
#SBATCH --gres=gpu:a100:4
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=2-00:00:00
#SBATCH --exclude=pax007
#SBATCH --output=/cluster/home/nbhatt04/lean_pipeline/md_100ns/logs/md_%j.log
#SBATCH --error=/cluster/home/nbhatt04/lean_pipeline/md_100ns/logs/md_%j.err

# Usage: sbatch --export=SIM=apo_2to3 submit_md_100ns.sh
#        sbatch --export=SIM=ach_2to3 submit_md_100ns.sh
#        sbatch --export=SIM=apo_3to2 submit_md_100ns.sh

set -e

export OPENMM_PLUGIN_DIR=/cluster/home/nbhatt04/.conda/envs/md_cuda/lib/plugins
PYTHON=/cluster/home/nbhatt04/.conda/envs/md_cuda/bin/python
PIPELINE=/cluster/home/nbhatt04/lean_pipeline

SIM=${SIM:-apo_2to3}

echo "=========================================="
echo "100ns GPU MD Simulation"
echo "Simulation: $SIM"
echo "Started: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
nvidia-smi -L
echo "=========================================="

$PYTHON $PIPELINE/scripts/md_100ns_runner.py --sim $SIM

echo ""
echo "=========================================="
echo "Complete: $(date)"
echo "=========================================="
