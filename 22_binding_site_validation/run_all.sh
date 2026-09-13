#!/bin/bash
# =============================================================================
# run_all.sh  -- run the full validation module (scripts 01-14)
# =============================================================================
# Run AFTER the sdcofold_step3 Boltz-2 array (3633143) completes.
# SiteAF3 predictions are gated separately (submit_validation.sh) and their
# outputs fold into scripts 06/10/13/14 on re-run.
# =============================================================================

set -euo pipefail

PY=/cluster/scratch/nbhatt04/conda/envs/polymer_ai/bin/python
PKG=/cluster/home/nbhatt04/lean_pipeline/22_binding_site_validation
cd "$PKG"

echo "[run_all] $(date) start"
"$PY" 01_existing_evidence_audit.py
"$PY" 02_af3_boltz2_ensemble_analysis.py
"$PY" 03_siteaf3_prepare.py
"$PY" 05_boltz2_consensus.py
"$PY" 06_three_method_convergence.py
"$PY" 07_medchem_qc.py
"$PY" 08_pocket_plasticity.py
"$PY" 09_sar_activity_cliff_validation.py
"$PY" 10_residue_consensus.py
"$PY" 11_mutagenesis_plan.py
"$PY" 12_md_decision.py
"$PY" 13_evidence_matrix.py
"$PY" 14_generate_report.py
echo "[run_all] $(date) done -> outputs/14_validation_report.md"