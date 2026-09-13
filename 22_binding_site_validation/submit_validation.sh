#!/bin/bash
# =============================================================================
# submit_validation.sh  -- submit the validation module gated on the
# sdcofold_step3 array (3633143).
#
#   STEP 1 (optional, runs after array): execute the CPU analysis pipeline.
#   STEP 2 (GPU): SiteAF3 array (80 configs, %4).
#
# Usage:
#   sbatch submit_validation.sh            # CPU pipeline + SiteAF3, dep 3633143
#   sbatch --dependency=afterok:3633143 submit_validation.sh
#
# NOTE: SiteAF3 requires an AlphaFold3 python env + model weights + MSA DBs.
# Set SITEAF3_AF3_ENV / SITEAF3_MODEL_DIR / SITEAF3_DB_DIR if defaults are wrong.
# Output goes to /cluster/scratch by default (SITEAF3_OUTPUT_DIR override).
# Shared Tufts AF3 data: see SITEAF3_AF3_SETUP.md.
# =============================================================================

#SBATCH --job-name=bs_validation
#SBATCH --output=/cluster/home/nbhatt04/lean_pipeline/22_binding_site_validation/logs/validation_%j.out
#SBATCH --error=/cluster/home/nbhatt04/lean_pipeline/22_binding_site_validation/logs/validation_%j.err
#SBATCH --time=02:00:00

set -euo pipefail
PIPELINE=/cluster/home/nbhatt04/lean_pipeline
PKG=$PIPELINE/22_binding_site_validation
mkdir -p "$PKG/logs"

echo "[submit_validation] depending on job 3633143 (sdcofold_step3)"
echo "[submit_validation] $(date) start"

bash "$PKG/run_all.sh"

echo "[submit_validation] CPU analysis done $(date)"
echo "[submit_validation] launching SiteAF3 (SBFLAKE) ..."

if [ -n "${SITEAF3_AF3_ENV:-}" ] || [ -n "${SITEAF3_MODEL_DIR:-}" ]; then
    S4=$(sbatch --parsable "$PKG/04_siteaf3_submit.sh")
    echo "[submit_validation] submitted SiteAF3 array j/$S4"
else
    echo "[submit_validation] WARN: SiteAF3 AF3 env / weights not configured;"
    echo "                     run 04_siteaf3_submit.sh manually after setup."
fi

echo "[submit_validation] done $(date)"