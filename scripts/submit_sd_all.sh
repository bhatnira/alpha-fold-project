#!/bin/bash
#SBATCH --job-name=boltz2_sd_all
#SBATCH --partition=gpu
#SBATCH --gres=gpu:a100:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=23:00:00
#SBATCH --exclude=pax007
#SBATCH --output=/cluster/home/nbhatt04/lean_pipeline/logs/boltz2_sd_all_%j.log
#SBATCH --error=/cluster/home/nbhatt04/lean_pipeline/logs/boltz2_sd_all_%j.err

set -e

module purge 2>/dev/null || true
module load boltz/2.2.1-gpu

PIPELINE=/cluster/home/nbhatt04/lean_pipeline
INPUT_DIR=$PIPELINE/08_site_directed_af3/yaml_inputs_v2
OUT_DIR=$PIPELINE/08_site_directed_af3/yaml_inputs_v2/predictions
CACHE=/cluster/home/nbhatt04/.boltz

echo "=========================================="
echo "Boltz-2 Site-Directed AF3 (All Sites)"
echo "Started: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
nvidia-smi -L | head -2
echo "=========================================="

cd $INPUT_DIR

for yaml_file in sd_site21_*.yaml sd_site25_*.yaml sd_site34_*.yaml; do
    name=$(basename "$yaml_file" .yaml)
    # Skip if already predicted
    if [ -d "$OUT_DIR/boltz_results_${name}" ]; then
        echo "  SKIP: $name (already done)"
        continue
    fi
    echo ""
    echo "=== Predicting: $name ==="
    boltz predict "$yaml_file" \
        --out_dir $OUT_DIR \
        --cache $CACHE \
        --model boltz2 \
        --accelerator gpu \
        --devices 1 \
        --recycling_steps 3 \
        --sampling_steps 200 \
        --diffusion_samples 3 \
        --output_format pdb \
        --num_workers 1 \
        --preprocessing-threads 1 || echo "  WARNING: $name failed"
done

echo ""
echo "=========================================="
echo "Site-directed AF3 complete: $(date)"
echo "=========================================="
