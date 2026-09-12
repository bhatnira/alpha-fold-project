#!/bin/bash
#SBATCH --job-name=a9a10_sdaf3
#SBATCH --partition=gpu
#SBATCH --gres=gpu:a100:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=23:00:00
#SBATCH --exclude=pax007
#SBATCH --output=/cluster/home/nbhatt04/lean_pipeline/logs/sitedirected_%j.log

set -e

module purge 2>/dev/null || true
module load boltz/2.2.1-gpu

cd /cluster/home/nbhatt04/lean_pipeline/08_site_directed_af3/yaml_inputs
mkdir -p predictions

echo "=== Site-directed AF3 started: $(date) ==="
echo "Node: $(hostname)"
nvidia-smi -L | head -2

for yaml_file in /cluster/home/nbhatt04/lean_pipeline/08_site_directed_af3/yaml_inputs/sd_*.yaml; do
    name=$(basename "$yaml_file" .yaml)
    echo ""
    echo "=== Predicting: $name ==="
    boltz predict "$yaml_file" \
        --out_dir predictions \
        --cache /cluster/home/nbhatt04/.boltz \
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

echo ""
echo "=== Site-directed AF3 complete: $(date) ==="
