#!/bin/bash
#SBATCH --job-name=boltz2_affinity
#SBATCH --partition=gpu
#SBATCH --gres=gpu:a100:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=boltz2_%j.out
#SBATCH --error=boltz2_%j.err

echo "Job started on $(hostname) at $(date)"
echo "SLURM job ID: $SLURM_JOB_ID"
echo "GPU: $CUDA_VISIBLE_DEVICES"

module purge 2>/dev/null
module load cuda 2>/dev/null

/cluster/home/nbhatt04/.venvs/pymol/bin/python \
    /cluster/home/nbhatt04/lean_pipeline/phase08_boltz_affinity/boltz2_affinity.py

echo "Job finished at $(date)"
