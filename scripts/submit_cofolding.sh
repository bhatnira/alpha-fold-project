#!/bin/bash
#SBATCH --job-name=boltz2_cofold
#SBATCH --partition=gpu
#SBATCH --gres=gpu:a100:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=23:00:00
#SBATCH --exclude=pax007
#SBATCH --output=/cluster/home/nbhatt04/lean_pipeline/logs/boltz2_cofold_%j.log
#SBATCH --error=/cluster/home/nbhatt04/lean_pipeline/logs/boltz2_cofold_%j.err

set -e

module purge 2>/dev/null || true
module load boltz/2.2.1-gpu

PIPELINE=/cluster/home/nbhatt04/lean_pipeline
INPUT_BASE=$PIPELINE/cofolding_study/boltz2_inputs
OUTPUT_BASE=$PIPELINE/cofolding_study/results/boltz2
CACHE=/cluster/home/nbhatt04/.boltz

echo "=========================================="
echo "Boltz-2 Cofolding Study (All Sites)"
echo "Started: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
nvidia-smi -L | head -2
echo "=========================================="

for site_dir in $INPUT_BASE/site*; do
    site=$(basename $site_dir)
    echo ""
    echo "=== Site: $site ==="
    
    for yaml_file in $site_dir/*.yaml; do
        name=$(basename "$yaml_file" .yaml)
        
        # Skip if already done
        result_dir="$OUTPUT_BASE/$site/$name/boltz_results_$name"
        if [ -d "$result_dir" ] && [ -f "$result_dir/predictions/"*"/"*.pdb ] 2>/dev/null; then
            echo "  SKIP: $name (already done)"
            continue
        fi
        
        echo "  Predicting: $name"
        mkdir -p "$OUTPUT_BASE/$site/$name"
        boltz predict "$yaml_file" \
            --out_dir "$OUTPUT_BASE/$site/$name" \
            --cache $CACHE \
            --model boltz2 \
            --accelerator gpu \
            --devices 1 \
            --recycling_steps 3 \
            --sampling_steps 200 \
            --diffusion_samples 2 \
            --output_format pdb \
            --num_workers 1 \
            --preprocessing-threads 1 2>&1 | tee "$OUTPUT_BASE/$site/$name/log.txt" || echo "  WARNING: $name failed"
    done
done

echo ""
echo "=== Collecting PDB files ==="
find $OUTPUT_BASE -name "*.pdb" 2>/dev/null | wc -l
echo "PDB files generated"

echo ""
echo "=========================================="
echo "Cofolding complete: $(date)"
echo "=========================================="
