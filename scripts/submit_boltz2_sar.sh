#!/bin/bash
#SBATCH --job-name=boltz2_sar
#SBATCH --partition=gpu
#SBATCH --gres=gpu:a100:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=23:00:00
#SBATCH --exclude=pax007
#SBATCH --output=/cluster/home/nbhatt04/lean_pipeline/logs/boltz2_sar_%j.log
#SBATCH --error=/cluster/home/nbhatt04/lean_pipeline/logs/boltz2_sar_%j.err

set -e

module purge 2>/dev/null || true
module load boltz/2.2.1-gpu

PIPELINE=/cluster/home/nbhatt04/lean_pipeline
INPUT_DIR=$PIPELINE/phase08_boltz_affinity/yaml_inputs_sar
OUT_DIR=$PIPELINE/phase08_boltz_affinity/boltz2_sar_outputs
CACHE=/cluster/home/nbhatt04/.boltz

mkdir -p $OUT_DIR

echo "=========================================="
echo "Boltz-2 SAR Affinity Prediction"
echo "Started: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
nvidia-smi -L | head -2
echo "=========================================="

cd $INPUT_DIR

for yaml_file in sar_cpd*.yaml; do
    name=$(basename "$yaml_file" .yaml)
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
        --diffusion_samples 1 \
        --sampling_steps_affinity 200 \
        --diffusion_samples_affinity 5 \
        --output_format pdb \
        --num_workers 1 \
        --preprocessing-threads 1 || echo "  WARNING: $name failed"
done

echo ""
echo "=== Collecting results ==="
cd $PIPELINE
python3 scripts/collect_sar_results.py

echo ""
echo "=========================================="
echo "SAR Affinity complete: $(date)"
echo "=========================================="
