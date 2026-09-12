#!/bin/bash
#SBATCH --job-name=openmm_test
#SBATCH --partition=gpu
#SBATCH --gres=gpu:a100:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:05:00
#SBATCH --output=/cluster/home/nbhatt04/lean_pipeline/logs/openmm_test_%j.log
#SBATCH --error=/cluster/home/nbhatt04/lean_pipeline/logs/openmm_test_%j.err

set -e

PYTHON=/cluster/home/nbhatt04/.conda/envs/md_cuda/bin/python

echo "=== OpenMM CUDA Test ==="
echo "Node: $(hostname)"
nvidia-smi -L

$PYTHON -c "
import openmm
import os

os.environ['OPENMM_PLUGIN_DIR'] = '/cluster/home/nbhatt04/.conda/envs/md_cuda/lib/python3.11/site-packages/openmm/plugins'

print('OpenMM version:', openmm.__version__)
print('Platforms:', [openmm.Platform.getPlatform(i).getName() for i in range(openmm.Platform.getNumPlatforms())])

try:
    platform = openmm.Platform.getPlatformByName('CUDA')
    print('CUDA platform: OK')
    print('CUDA device:', platform.getPropertyValue('DeviceIndex'))
except Exception as e:
    print('CUDA error:', e)
"

echo "=== Test complete ==="
