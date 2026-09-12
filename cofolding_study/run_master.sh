#!/bin/bash
# =============================================================================
# Master Submission Script: Cofolding Study
# =============================================================================
# Generates inputs, submits jobs, and tracks progress
# =============================================================================

set -e

PYTHON="/cluster/home/nbhatt04/.venvs/pymol/bin/python"
BASE_DIR="/cluster/home/nbhatt04/lean_pipeline/cofolding_study"
SCRIPTS_DIR="$BASE_DIR/slurm_scripts"

echo "=========================================="
echo "Cofolding Study Master Script"
echo "Started: $(date)"
echo "=========================================="

# Step 1: Generate AF3 inputs
echo ""
echo "[Step 1] Generating AF3 ternary complex inputs..."
$PYTHON "$BASE_DIR/generate_af3_ternary_inputs.py"
echo "  AF3 inputs generated in: $BASE_DIR/af3_inputs/"

# Step 2: Generate Boltz-2 YAML inputs
echo ""
echo "[Step 2] Generating Boltz-2 affinity YAML inputs..."
$PYTHON "$BASE_DIR/generate_boltz2_yaml.py"
echo "  Boltz-2 inputs generated in: $BASE_DIR/boltz2_inputs/"

# Step 3: Submit Boltz-2 jobs
echo ""
echo "[Step 3] Submitting Boltz-2 SAR affinity predictions..."
cd "$SCRIPTS_DIR"
BOLTZ_JOB_ID=$(sbatch --parsable run_boltz2_sar.sh)
echo "  Boltz-2 job submitted: $BOLTZ_JOB_ID"

# Step 4: Submit AF3 jobs (if infrastructure available)
echo ""
echo "[Step 4] AF3 ternary complex inputs ready"
echo "  AF3 inputs are in: $BASE_DIR/af3_inputs/"
echo "  Submit via ColabFold API or web interface at: https://alphafoldserver.com"
echo "  Or run: sbatch $SCRIPTS_DIR/run_af3_ternary.sh"

# Step 5: Create analysis tracking
echo ""
echo "[Step 5] Creating analysis tracking..."
$PYTHON -c "
import json
from pathlib import Path

tracking = {
    'status': 'submitted',
    'boltz2_job_id': '$BOLTZ_JOB_ID',
    'af3_status': 'inputs_generated',
    'next_steps': [
        'Monitor Boltz-2 jobs: squeue -u \$USER',
        'Submit AF3 jobs via ColabFold',
        'Run analysis when results available',
    ]
}

tracking_path = Path('$BASE_DIR/tracking.json')
with open(tracking_path, 'w') as f:
    json.dump(tracking, f, indent=2)
print(f'  Tracking: {tracking_path}')
"

echo ""
echo "=========================================="
echo "Setup Complete"
echo "=========================================="
echo ""
echo "Next steps:"
echo "  1. Monitor Boltz-2 jobs: squeue -u \$USER"
echo "  2. Submit AF3 jobs via ColabFold API"
echo "  3. Run analysis: python $BASE_DIR/analysis/interaction_fingerprint_analysis.py"
echo "  4. Run SAR correlation: python $BASE_DIR/analysis/sar_correlation_analysis.py"
echo ""
echo "Completed: $(date)"
