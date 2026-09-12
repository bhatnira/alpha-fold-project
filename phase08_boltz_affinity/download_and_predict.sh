#!/bin/bash
#SBATCH --job-name=boltz2_download_predict
#SBATCH --partition=gpu
#SBATCH --gres=gpu:a100:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=boltz2_%j.out
#SBATCH --error=boltz2_%j.err

set -e

echo "=== Job started on $(hostname) at $(date) ==="
echo "SLURM job ID: $SLURM_JOB_ID"
echo "GPU: $CUDA_VISIBLE_DEVICES"

module purge 2>/dev/null || true
module load cuda 2>/dev/null || true

PYTHON=/cluster/home/nbhatt04/.venvs/pymol/bin/python

echo ""
echo "=== Step 1: Verify CCD data and model weights ==="
MOL_DIR=~/.boltz/mols
MOL_TAR=~/.boltz/mols.tar
AFF_CKPT=~/.boltz/boltz2_aff.ckpt
CONF_CKPT=~/.boltz/boltz2_conf.ckpt

# Check if mols dir has ALA.pkl (canonical amino acid)
if [ -f "$MOL_DIR/ALA.pkl" ]; then
    echo "CCD data OK (ALA.pkl found)"
else
    echo "CCD data missing or incomplete, re-downloading..."
    rm -rf "$MOL_DIR" "$MOL_TAR"
    
    echo "Downloading mols.tar (~2.3GB)..."
    wget -q --show-progress -O "$MOL_TAR" \
        "https://huggingface.co/boltz-community/boltz-2/resolve/main/mols.tar" || \
    curl -L -o "$MOL_TAR" \
        "https://huggingface.co/boltz-community/boltz-2/resolve/main/mols.tar"
    
    echo "Extracting mols.tar..."
    tar xf "$MOL_TAR" -C ~/.boltz/
    
    echo "Verifying extraction..."
    ls "$MOL_DIR/ALA.pkl" && echo "CCD extraction OK"
fi

# Check affinity model
FILE_SIZE=$(stat -c%s "$AFF_CKPT" 2>/dev/null || echo "0")
if [ "$FILE_SIZE" -gt 1000000 ]; then
    echo "Affinity model OK ($FILE_SIZE bytes)"
else
    echo "Affinity model missing or corrupt ($FILE_SIZE bytes), re-downloading..."
    rm -f "$AFF_CKPT"
    
    echo "Downloading boltz2_aff.ckpt..."
    wget -q --show-progress -O "$AFF_CKPT" \
        "https://model-gateway.boltz.bio/boltz2_aff.ckpt" || \
    wget -q --show-progress -O "$AFF_CKPT" \
        "https://huggingface.co/boltz-community/boltz-2/resolve/main/boltz2_aff.ckpt"
    
    FILE_SIZE=$(stat -c%s "$AFF_CKPT" 2>/dev/null || echo "0")
    echo "Affinity model downloaded: $FILE_SIZE bytes"
fi

echo ""
echo "=== Step 2: Quick validation test ==="
$PYTHON -c "
from boltz.data.mol import load_canonicals
from pathlib import Path
mol_dir = Path.home() / '.boltz' / 'mols'
ccd = load_canonicals(str(mol_dir))
print(f'Loaded {len(ccd)} canonical molecules')
print('ALA loaded:', 'ALA' in ccd)
"

echo ""
echo "=== Step 3: Run Boltz-2 affinity predictions ==="
cd /cluster/home/nbhatt04/lean_pipeline/phase08_boltz_affinity
$PYTHON boltz2_affinity.py

echo ""
echo "=== Job finished at $(date) ==="
