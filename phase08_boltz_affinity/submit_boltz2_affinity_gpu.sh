#!/bin/bash
#SBATCH --job-name=boltz2_affinity
#SBATCH --partition=gpu
#SBATCH --gres=gpu:a100:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=23:00:00
#SBATCH --exclude=pax007
#SBATCH --output=boltz2_affinity_%j.log

set -e

echo "=== Job started on $(hostname) at $(date) ==="
echo "SLURM job ID: $SLURM_JOB_ID"
echo "GPU devices: $CUDA_VISIBLE_DEVICES"

module purge 2>/dev/null || true
module load boltz/2.2.1-gpu

cd /cluster/home/nbhatt04/lean_pipeline/phase08_boltz_affinity

INPUT_DIR=/cluster/home/nbhatt04/lean_pipeline/phase08_boltz_affinity/yaml_inputs
OUT_DIR=/cluster/home/nbhatt04/lean_pipeline/phase08_boltz_affinity/boltz2_outputs
CACHE=/cluster/home/nbhatt04/.boltz

# Sanity checks on the compute node
echo "=== checks ==="
ls "$CACHE/boltz2_conf.ckpt" "$CACHE/boltz2_aff.ckpt" "$CACHE/mols/ALA.pkl"
nvidia-smi -L | head -2

echo ""
echo "=== Running Boltz-2 affinity prediction (12 inputs) ==="
boltz predict "$INPUT_DIR" \
    --out_dir "$OUT_DIR" \
    --cache "$CACHE" \
    --model boltz2 \
    --accelerator gpu \
    --devices 1 \
    --recycling_steps 3 \
    --sampling_steps 200 \
    --diffusion_samples 1 \
    --sampling_steps_affinity 200 \
    --diffusion_samples_affinity 5 \
    --output_format pdb \
    --num_workers 1 \
    --preprocessing-threads 1
RC=$?
echo "boltz predict exit code: $RC"

echo ""
echo "=== Job finished at $(date) ==="
exit $RC