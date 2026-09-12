#!/bin/bash
# Master pipeline runner for alpha9alpha10 nAChR allosteric-site discovery
# Runs all 22 phases as background SLURM jobs
# Uses A100 GPUs for Boltz-2 predictions

#SBATCH --job-name=lean_pipeline
#SBATCH --partition=gpu
#SBATCH --gres=gpu:a100:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=lean_pipeline_%j.out
#SBATCH --error=lean_pipeline_%j.err

set -e

PIPELINE_DIR="/cluster/home/nbhatt04/lean_pipeline"
PYTHON="/cluster/home/nbhatt04/.venvs/pymol/bin/python"

echo "=========================================="
echo "Alpha9alpha10 Pipeline Runner"
echo "Started: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "GPU: $CUDA_VISIBLE_DEVICES"
echo "=========================================="

cd "$PIPELINE_DIR"

# Phase I - SAR Dataset
echo ""
echo "[$(date +%H:%M:%S)] Phase I: SAR Dataset"
$PYTHON phase01_sar/build_sar.py

# Phase II - Receptor Ensemble
echo ""
echo "[$(date +%H:%M:%S)] Phase II: Receptor Ensemble"
$PYTHON phase02_receptor/inventory_receptors.py

# Phase III - Pocket Discovery
echo ""
echo "[$(date +%H:%M:%S)] Phase III: Pocket Discovery"
if [ -f "phase03_pockets/run_fpocket.sh" ]; then
    bash phase03_pockets/run_fpocket.sh
fi

# Phase IV - Convergence
echo ""
echo "[$(date +%H:%M:%S)] Phase IV: AF3 + Boltz-2 Convergence"
$PYTHON phase04_convergence/map_convergence.py

# Phase V - Docking
echo ""
echo "[$(date +%H:%M:%S)] Phase V: Ensemble Docking"
$PYTHON phase05_docking/run_docking.py

# Phase VI - Fingerprints
echo ""
echo "[$(date +%H:%M:%S)] Phase VI: Interaction Fingerprints"
$PYTHON phase06_fingerprints/build_ifp.py

# Phase VII - SAR Model
echo ""
echo "[$(date +%H:%M:%S)] Phase VII: SAR Model"
$PYTHON phase07_sar_model/derive_sar_rules.py

# Phase VIII - Boltz-2 Affinity (GPU)
echo ""
echo "[$(date +%H:%M:%S)] Phase VIII: Boltz-2 Affinity (GPU)"
$PYTHON phase08_boltz_affinity/boltz2_affinity.py

# Phase IX - Specificity
echo ""
echo "[$(date +%H:%M:%S)] Phase IX: Specificity"
$PYTHON phase09_specificity/run_specificity.py

# Phase X - Electrostatics
echo ""
echo "[$(date +%H:%M:%S)] Phase X: Electrostatics"
$PYTHON phase10_electrostatics/compute_electrostatics.py

# Phase XI - Robustness
echo ""
echo "[$(date +%H:%M:%S)] Phase XI: Robustness"
$PYTHON phase11_robustness/compute_robustness.py

# Phase XII - Subtype
echo ""
echo "[$(date +%H:%M:%S)] Phase XII: Subtype Comparison"
$PYTHON phase12_subtype/compare_subtypes.py

# Phase XIII - Ryanodine
echo ""
echo "[$(date +%H:%M:%S)] Phase XIII: Ryanodine"
$PYTHON phase13_ryanodine/analyze_ryanodine.py

# Phase XIV - Evidence
echo ""
echo "[$(date +%H:%M:%S)] Phase XIV: Evidence Matrix"
$PYTHON phase14_evidence/build_evidence_matrix.py

# Phase XV - Pharmacophore
echo ""
echo "[$(date +%H:%M:%S)] Phase XV: Pharmacophore"
$PYTHON phase15_pharmacophore/build_pharmacophore.py

# Phase XVI - Molecule Generation
echo ""
echo "[$(date +%H:%M:%S)] Phase XVI: Molecule Generation"
$PYTHON phase16_molgen/generate_molecules.py

# Phase XVII - Cascade Filter
echo ""
echo "[$(date +%H:%M:%S)] Phase XVII: Cascade Filter"
$PYTHON phase17_cascade/filter_library.py

# Phase XVIII - Boltz-2 Design (GPU)
echo ""
echo "[$(date +%H:%M:%S)] Phase XVIII: Boltz-2 Design (GPU)"
$PYTHON phase18_boltz_design/run_boltz2.py

# Phase XIX - Selectivity
echo ""
echo "[$(date +%H:%M:%S)] Phase XIX: Selectivity"
$PYTHON phase19_selectivity/compute_selectivity.py

# Phase XX - Ranking
echo ""
echo "[$(date +%H:%M:%S)] Phase XX: Ranking"
$PYTHON phase20_ranking/compute_ranking.py

# Phase XXI - Validation
echo ""
echo "[$(date +%H:%M:%S)] Phase XXI: Validation"
$PYTHON phase21_validation/design_validation.py

# Phase XXII - Closed Loop
echo ""
echo "[$(date +%H:%M:%S)] Phase XXII: Closed Loop"
$PYTHON phase22_closed_loop/create_loop_template.py

echo ""
echo "=========================================="
echo "Pipeline Complete: $(date)"
echo "=========================================="
