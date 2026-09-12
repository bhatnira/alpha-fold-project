#!/bin/bash
#SBATCH --job-name=a9a10_gap_fill
#SBATCH --partition=gpu
#SBATCH --gres=gpu:a100:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=46:00:00
#SBATCH --exclude=pax007
#SBATCH --output=/cluster/home/nbhatt04/lean_pipeline/logs/gap_fill_%j.log

set -e

PIPELINE=/cluster/home/nbhatt04/lean_pipeline
BOLTZ_CACHE=/cluster/home/nbhatt04/.boltz

module purge 2>/dev/null || true
module load boltz/2.2.1-gpu

echo "=========================================="
echo "Gap-Fill Pipeline"
echo "Started: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
nvidia-smi -L | head -2
echo "=========================================="

# PHASE 3-4: Open and Desensitized States
echo ""
echo "[$(date +%H:%M:%S)] PHASE 3-4: State Model Generation"
mkdir -p $PIPELINE/02_receptor_ensemble/state_models/predictions
for yaml_file in $PIPELINE/02_receptor_ensemble/state_models/yaml_inputs/*.yaml; do
    name=$(basename "$yaml_file" .yaml)
    echo "  Predicting: $name"
    boltz predict "$yaml_file" \
        --out_dir $PIPELINE/02_receptor_ensemble/state_models/predictions \
        --cache $BOLTZ_CACHE \
        --model boltz2 \
        --accelerator gpu \
        --devices 1 \
        --recycling_steps 3 \
        --sampling_steps 200 \
        --diffusion_samples 2 \
        --output_format pdb \
        --num_workers 1 \
        --preprocessing-threads 1
done
echo "  Phase 3-4 complete: $(date)"

# PHASE 6: ACh Occupancy Models
echo ""
echo "[$(date +%H:%M:%S)] PHASE 6: ACh Occupancy Modeling"
mkdir -p $PIPELINE/05_ach_occupancy/yaml_inputs/predictions
for yaml_file in $PIPELINE/05_ach_occupancy/yaml_inputs/*.yaml; do
    name=$(basename "$yaml_file" .yaml)
    echo "  Predicting: $name"
    boltz predict "$yaml_file" \
        --out_dir $PIPELINE/05_ach_occupancy/yaml_inputs/predictions \
        --cache $BOLTZ_CACHE \
        --model boltz2 \
        --accelerator gpu \
        --devices 1 \
        --recycling_steps 3 \
        --sampling_steps 200 \
        --diffusion_samples 2 \
        --output_format pdb \
        --num_workers 1 \
        --preprocessing-threads 1
done
echo "  Phase 6 complete: $(date)"

# PHASE 10: Site-Directed AF3
echo ""
echo "[$(date +%H:%M:%S)] PHASE 10: Site-Directed AF3"
mkdir -p $PIPELINE/08_site_directed_af3/yaml_inputs/predictions
for yaml_file in $PIPELINE/08_site_directed_af3/yaml_inputs/sd_*.yaml; do
    name=$(basename "$yaml_file" .yaml)
    echo "  Predicting: $name"
    boltz predict "$yaml_file" \
        --out_dir $PIPELINE/08_site_directed_af3/yaml_inputs/predictions \
        --cache $BOLTZ_CACHE \
        --model boltz2 \
        --accelerator gpu \
        --devices 1 \
        --recycling_steps 3 \
        --sampling_steps 200 \
        --diffusion_samples 3 \
        --output_format pdb \
        --num_workers 1 \
        --preprocessing-threads 1
done
echo "  Phase 10 complete: $(date)"

echo ""
echo "=========================================="
echo "All GPU jobs complete: $(date)"
echo "=========================================="
