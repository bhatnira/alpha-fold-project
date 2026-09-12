#!/bin/bash
#SBATCH --job-name=a9a10_gap2
#SBATCH --partition=gpu
#SBATCH --gres=gpu:a100:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=23:00:00
#SBATCH --exclude=pax007
#SBATCH --output=/cluster/home/nbhatt04/lean_pipeline/logs/gap2_%j.log

set -e

PIPELINE=/cluster/home/nbhatt04/lean_pipeline
BOLTZ_CACHE=/cluster/home/nbhatt04/.boltz

module purge 2>/dev/null || true
module load boltz/2.2.1-gpu

echo "=========================================="
echo "Gap-Fill Pipeline v2 (2-chain YAMLs)"
echo "Started: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
nvidia-smi -L | head -2
echo "=========================================="

# PHASE 6: ACh Occupancy Models (2-chain + ACh/PAM)
echo ""
echo "[$(date +%H:%M:%S)] PHASE 6: ACh Occupancy Modeling"
mkdir -p $PIPELINE/05_ach_occupancy/yaml_inputs_v2/predictions
for yaml_file in $PIPELINE/05_ach_occupancy/yaml_inputs_v2/*.yaml; do
    name=$(basename "$yaml_file" .yaml)
    echo "  Predicting: $name"
    boltz predict "$yaml_file" \
        --out_dir $PIPELINE/05_ach_occupancy/yaml_inputs_v2/predictions \
        --cache $BOLTZ_CACHE \
        --model boltz2 \
        --accelerator gpu \
        --devices 1 \
        --recycling_steps 3 \
        --sampling_steps 200 \
        --diffusion_samples 2 \
        --output_format pdb \
        --num_workers 1 \
        --preprocessing-threads 1 || echo "  WARNING: $name failed, continuing..."
done
echo "  Phase 6 complete: $(date)"

# PHASE 10: Site-Directed AF3 (2-chain + ligand)
echo ""
echo "[$(date +%H:%M:%S)] PHASE 10: Site-Directed AF3"
mkdir -p $PIPELINE/08_site_directed_af3/yaml_inputs_v2/predictions
for yaml_file in $PIPELINE/08_site_directed_af3/yaml_inputs_v2/*.yaml; do
    name=$(basename "$yaml_file" .yaml)
    echo "  Predicting: $name"
    boltz predict "$yaml_file" \
        --out_dir $PIPELINE/08_site_directed_af3/yaml_inputs_v2/predictions \
        --cache $BOLTZ_CACHE \
        --model boltz2 \
        --accelerator gpu \
        --devices 1 \
        --recycling_steps 3 \
        --sampling_steps 200 \
        --diffusion_samples 3 \
        --output_format pdb \
        --num_workers 1 \
        --preprocessing-threads 1 || echo "  WARNING: $name failed, continuing..."
done
echo "  Phase 10 complete: $(date)"

echo ""
echo "=========================================="
echo "All GPU jobs complete: $(date)"
echo "=========================================="
