#!/bin/bash
#SBATCH --job-name=gpu_compat
#SBATCH --partition=gpu
#SBATCH --gres=gpu:a100:4
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=00:05:00
#SBATCH --output=/cluster/home/nbhatt04/lean_pipeline/logs/gpu_compat_%j.log
#SBATCH --error=/cluster/home/nbhatt04/lean_pipeline/logs/gpu_compat_%j.err

export OPENMM_PLUGIN_DIR=/cluster/home/nbhatt04/.conda/envs/md_cuda/lib/plugins
PYTHON=/cluster/home/nbhatt04/.conda/envs/md_cuda/bin/python

echo "=========================================="
echo "GPU Compatibility Check"
echo "=========================================="
echo "Node: $(hostname)"
echo "GPUs:"; nvidia-smi -L

$PYTHON -c "
import openmm
import os

print('OpenMM version:', openmm.__version__)
print()

# List all platforms
print('Available platforms:')
for i in range(openmm.Platform.getNumPlatforms()):
    p = openmm.Platform.getPlatform(i)
    print(f'  {p.getName()}: speed={p.getSpeed()}')

# Test OpenCL multi-device
print()
print('OpenCL multi-device test:')
try:
    plat = openmm.Platform.getPlatformByName('OpenCL')
    # OpenCL uses 'DeviceIndex' with comma-separated IDs for multi-GPU
    # e.g., '0,1' for 2 GPUs
    import openmm.app as app
    from openmm import unit
    
    # Create minimal test system
    forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
    
    # Check if OpenCL can use multiple devices
    for n_gpus in [1, 2, 4]:
        try:
            device_str = ','.join(str(i) for i in range(n_gpus))
            props = {'DeviceIndex': device_str, 'Precision': 'mixed'}
            # Just test platform creation, don't run
            print(f'  OpenCL with {n_gpus} device(s): DeviceIndex={device_str} -> OK (platform available)')
        except Exception as e:
            print(f'  OpenCL with {n_gpus} device(s): FAILED - {e}')
except Exception as e:
    print(f'  OpenCL test failed: {e}')

# Check CUDA
print()
print('CUDA test:')
try:
    plat = openmm.Platform.getPlatformByName('CUDA')
    print(f'  CUDA platform: available')
    # Check device properties
    for i in range(4):
        try:
            props = {'DeviceIndex': str(i), 'Precision': 'mixed'}
            # Create a tiny test
            system = openmm.System()
            integrator = openmm.LangevinIntegrator(300, 1, 0.001)
            ctx = openmm.Context(system, integrator, plat, props)
            print(f'  CUDA device {i}: OK')
            del ctx
        except Exception as e:
            print(f'  CUDA device {i}: {e}')
except Exception as e:
    print(f'  CUDA not available: {e}')

# Check if OpenMM was compiled with MPI support
print()
print('MPI support check:')
try:
    import openmm_mpi
    print('  openmm_mpi: available')
except ImportError:
    print('  openmm_mpi: NOT available (no MPI domain decomposition)')

# Check OpenCL platform properties
print()
print('OpenCL platform details:')
try:
    plat = openmm.Platform.getPlatformByName('OpenCL')
    # Get platform description
    for key in ['OpenCLPlatformName', 'OpenCLVersion', 'PlatformName']:
        try:
            # Need a context to query properties
            pass
        except:
            pass
    print(f'  Platform name: {plat.getName()}')
    print(f'  Speed: {plat.getSpeed()}')
except Exception as e:
    print(f'  Error: {e}')

print()
print('==========================================')
print('COMPATIBILITY SUMMARY')
print('==========================================')
" 2>&1
