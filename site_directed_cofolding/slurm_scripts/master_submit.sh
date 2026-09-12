#!/bin/bash
# =============================================================================
# MASTER SUBMIT: Site-Directed Cofolding Study
# =============================================================================
# Orchestrates all 3 steps with dependency chaining
# Run this script to submit everything
# =============================================================================

set -e

SCRIPTS_DIR="/cluster/home/nbhatt04/lean_pipeline/site_directed_cofolding/slurm_scripts"
LOGS_DIR="/cluster/home/nbhatt04/lean_pipeline/site_directed_cofolding/logs"
mkdir -p "$LOGS_DIR"

echo "=========================================="
echo "Site-Directed Cofolding: Master Submit"
echo "Started: $(date)"
echo "=========================================="

# Step 1: State models + ACh states (6 YAMLs, ~2-3h)
echo ""
echo "[Step 1] Submitting state models + ACh states..."
STEP1_ID=$(sbatch --parsable "$SCRIPTS_DIR/step1_state_models.sh")
echo "  Job ID: $STEP1_ID"

# Step 2: Site-directed binary + ternary (160 YAMLs, ~20-24h)
# Depend on step 1 completing
echo ""
echo "[Step 2] Submitting site-directed binary + ternary..."
STEP2_ID=$(sbatch --parsable --depend=afterok:$STEP1_ID "$SCRIPTS_DIR/step2_site_directed.sh")
echo "  Job ID: $STEP2_ID"

# Step 3: Full 30-compound panel (180 YAMLs, ~20-24h)
# Depend on step 2 completing
echo ""
echo "[Step 3] Submitting full 30-compound panel..."
STEP3_ID=$(sbatch --parsable --depend=afterok:$STEP2_ID "$SCRIPTS_DIR/step3_full_panel.sh")
echo "  Job ID: $STEP3_ID"

echo ""
echo "=========================================="
echo "All jobs submitted"
echo "=========================================="
echo ""
echo "Job chain:"
echo "  Step 1 (state models):     $STEP1_ID"
echo "  Step 2 (site-directed):    $STEP2_ID  (depends on $STEP1_ID)"
echo "  Step 3 (full panel):       $STEP3_ID  (depends on $STEP2_ID)"
echo ""
echo "Monitor with:"
echo "  squeue -u \$USER"
echo "  squeue --jobs=$STEP1_ID,$STEP2_ID,$STEP3_ID"
echo ""
echo "Cancel with:"
echo "  scancel $STEP1_ID $STEP2_ID $STEP3_ID"
echo ""
echo "Completed: $(date)"
