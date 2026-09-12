#!/bin/bash
# =============================================================================
# SiteAF3 Setup Script
# =============================================================================
# Installs SiteAF3 dependencies and verifies AF3 model weights
# Run once before submitting SiteAF3 jobs
# =============================================================================

set -e

PIPELINE=/cluster/home/nbhatt04/lean_pipeline
PHASE2=$PIPELINE/phase2_siteaf3_validation
SITEAF3=$PHASE2/SiteAF3

echo "=========================================="
echo "SiteAF3 Setup"
echo "Started: $(date)"
echo "=========================================="

# Check AF3 model weights
echo ""
echo "[1] Checking AF3 model weights..."
MODEL_DIR="/cluster/home/nbhatt04/models"
if [ -d "$MODEL_DIR" ]; then
    echo "  Model directory: $MODEL_DIR"
    ls "$MODEL_DIR"/*.pt 2>/dev/null | head -3
    echo "  OK"
else
    echo "  WARNING: Model directory not found at $MODEL_DIR"
    echo "  You may need to download AF3 weights"
fi

# Check AF3 databases
echo ""
echo "[2] Checking AF3 databases..."
DB_DIR="/cluster/home/nbhatt04/public_databases"
if [ -d "$DB_DIR" ]; then
    echo "  Database directory: $DB_DIR"
    ls "$DB_DIR"/*.fasta 2>/dev/null | head -3
    echo "  OK"
else
    echo "  WARNING: Database directory not found at $DB_DIR"
fi

# Check SiteAF3 installation
echo ""
echo "[3] Checking SiteAF3..."
cd "$SITEAF3"
python3 -c "import src.embeddings.embed_cond; print('  embeddings OK')" 2>/dev/null || echo "  embeddings: NOT INSTALLED"
python3 -c "import src.diffusion.run_cond_Diff; print('  diffusion OK')" 2>/dev/null || echo "  diffusion: NOT INSTALLED"
python3 -c "from alphafold3.model import model; print('  alphafold3 OK')" 2>/dev/null || echo "  alphafold3: NOT INSTALLED"
python3 -c "from Bio import PDB; print('  biopython OK')" 2>/dev/null || echo "  biopython: NOT INSTALLED"

# Check generated inputs
echo ""
echo "[4] Checking generated inputs..."
n_hotspot=$(find "$PHASE2/hotspot_pocket_pdbs" -name "*.pdb" | wc -l)
n_json=$(find "$PHASE2/siteaf3_inputs" -name "*.json" | wc -l)
echo "  Hotspot/pocket PDBs: $n_hotspot"
echo "  SiteAF3 JSON configs: $n_json"

echo ""
echo "=========================================="
echo "Setup check complete: $(date)"
echo "=========================================="
