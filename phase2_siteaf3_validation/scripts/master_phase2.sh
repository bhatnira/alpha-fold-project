#!/bin/bash
# =============================================================================
# PHASE 2: SiteAF3 Validation — Master Orchestration
# =============================================================================
# Runs AFTER current jobs complete (3332533, 3327732, 3327733)
# Submit this script when ready.
# =============================================================================

set -e

PIPELINE=/cluster/home/nbhatt04/lean_pipeline
PHASE2=$PIPELINE/phase2_siteaf3_validation
SCRIPTS=$PHASE2/scripts
LOGS=$PHASE2/logs
mkdir -p $LOGS

echo "=========================================="
echo "Phase 2: SiteAF3 Validation"
echo "Started: $(date)"
echo "=========================================="

# Step 1: Analyze current Boltz-2 cofolding results
echo ""
echo "[Step 1] Analyzing Boltz-2 cofolding results..."
python3 $SCRIPTS/analyze_cofolding_results.py
echo "  Analysis complete"

# Step 2: Generate SiteAF3 hotspot/pocket PDBs
echo ""
echo "[Step 2] Generating hotspot/pocket PDBs..."
python3 $SCRIPTS/generate_hotspot_pocket.py
echo "  Hotspot/pocket PDBs generated"

# Step 3: Generate SiteAF3 JSON configs
echo ""
echo "[Step 3] Generating SiteAF3 JSON configs..."
python3 $SCRIPTS/generate_siteaf3_configs.py
echo "  JSON configs generated"

# Step 4: Submit SiteAF3 runs (if SiteAF3 is installed)
if command -v run_SiteAF3.py &> /dev/null || [ -f "$PHASE2/siteaf3_setup_done.flag" ]; then
    echo ""
    echo "[Step 4] Submitting SiteAF3 runs..."
    sbatch $SCRIPTS/submit_siteaf3.sh
    echo "  SiteAF3 jobs submitted"
else
    echo ""
    echo "[Step 4] SiteAF3 not yet configured — skipping"
    echo "  Run: python3 $SCRIPTS/setup_siteaf3.py first"
fi

echo ""
echo "=========================================="
echo "Phase 2 setup complete: $(date)"
echo "=========================================="
