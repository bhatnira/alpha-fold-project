#!/bin/bash
# =============================================================================
# AF3 Cofolding: Ternary Complexes (Receptor + ACh + PAM)
# =============================================================================
# Tests whether cofolding changes in presence of ACh
# Tests binding site sharing across compound series
#
# Panel: 9 conditions × 2 stoichiometries × 5 seeds = 90 AF3 jobs
# =============================================================================

#SBATCH --job-name=af3_ternary
#SBATCH --partition=general
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=af3_ternary_%A_%a.out
#SBATCH --error=af3_ternary_%A_%a.err
#SBATCH --array=0-89%20

set -e

PYTHON="/cluster/home/nbhatt04/.venvs/pymol/bin/python"
BASE_DIR="/cluster/home/nbhatt04/lean_pipeline/cofolding_study"
INPUT_DIR="$BASE_DIR/af3_inputs"
RESULTS_DIR="$BASE_DIR/results/af3"
LOG_DIR="$BASE_DIR/results/af3_logs"
mkdir -p "$RESULTS_DIR" "$LOG_DIR"

echo "=========================================="
echo "AF3 Ternary Complex Cofolding"
echo "Started: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Array ID: $SLURM_ARRAY_TASK_ID"
echo "Node: $(hostname)"
echo "=========================================="

# Collect all input JSONs
ALL_INPUTS=($(find "$INPUT_DIR" -name "*_data.json" -type f | sort))
TOTAL=${#ALL_INPUTS[@]}

echo "Total AF3 input files: $TOTAL"

# Process single job for this array task
if [ $SLURM_ARRAY_TASK_ID -lt $TOTAL ]; then
    INPUT="${ALL_INPUTS[$SLURM_ARRAY_TASK_ID]}"
    BASENAME=$(basename "$INPUT" _data.json)
    OUTDIR="$RESULTS_DIR/$BASENAME"
    LOGFILE="$LOG_DIR/${BASENAME}.log"
    
    echo ""
    echo "[$(date +%H:%M:%S)] Processing: $BASENAME"
    echo "  Input: $INPUT"
    echo "  Output: $OUTDIR"
    
    mkdir -p "$OUTDIR"
    
    # Note: AF3 submission depends on available infrastructure
    # Options: ColabFold API, local AF3, or web interface
    # This script creates the submission command
    
    echo "  AF3 input ready at: $INPUT"
    echo "  To submit via ColabFold API, run:"
    echo "    $PYTHON -c \"from colabfold.batch import run_jobs; ...\""
    echo "  Or submit via web interface at: https://alphafoldserver.com"
    
    # Copy input to results for tracking
    cp "$INPUT" "$OUTDIR/"
    
    echo "  Input copied to results directory"
    echo "  Completed: $(date +%H:%M:%S)"
else
    echo "Array task $SLURM_ARRAY_TASK_ID exceeds total inputs ($TOTAL)"
fi

echo ""
echo "=========================================="
echo "Job completed: $(date)"
echo "=========================================="
