#!/usr/bin/env python3
"""
GPU-accelerated MD using OpenMM with CUDA.
Runs 5ns MD on alpha9/alpha10 ECD interface for cryptic pocket analysis.
"""
import os, sys, time, re
import openmm
from openmm import app, unit
from openmm.app import PDBFile, Simulation, DCDReporter, StateDataReporter
from pdbfixer import PDBFixer

WORKDIR = "/cluster/home/nbhatt04/lean_pipeline/md_gpu"
os.makedirs(WORKDIR, exist_ok=True)

SRC_PDB = "/cluster/home/nbhatt04/lean_pipeline/05_ach_occupancy/yaml_inputs_v2/predictions/boltz_results_ach_only_2to3/predictions/ach_only_2to3/ach_only_2to3_model_0.pdb"
CLEAN_PDB = os.path.join(WORKDIR, "clean_protein.pdb")
FIXED_PDB = os.path.join(WORKDIR, "fixed_protein.pdb")

# Strip HETATM and CONECT lines
print("Cleaning PDB (keeping ATOM only)...")
with open(SRC_PDB) as fin, open(CLEAN_PDB, "w") as fout:
    for line in fin:
        if line.startswith("ATOM"):
            fout.write(line)

# Fix structure with PDBFixer (add terminal caps, missing atoms, hydrogens)
print("Fixing structure with PDBFixer...")
fixer = PDBFixer(filename=CLEAN_PDB)
fixer.findMissingResidues()
fixer.findMissingAtoms()
fixer.addMissingAtoms()
fixer.addMissingHydrogens(7.4)
with open(FIXED_PDB, "w") as f:
    PDBFile.writeFile(fixer.topology, fixer.positions, f)
print(f"  Fixed structure written to: {FIXED_PDB}")

FORCE_FIELD = "amber14-all.xml"
WATER_MODEL = "amber14/tip3pfb.xml"
TEMPERATURE = 310 * unit.kelvin
FRICTION = 1.0 / unit.picoseconds
TIMESTEP = 2 * unit.femtoseconds

print("=" * 60)
print("OpenMM GPU MD - alpha9/alpha10 interface")
print("=" * 60)

# Step 1: Load fixed structure
print("\n[1/8] Loading fixed structure...")
pdb = PDBFile(FIXED_PDB)
print(f"  Atoms: {pdb.topology.getNumAtoms()}")
print(f"  Chains: {pdb.topology.getNumChains()}")

# Step 2: Create force field
print("\n[2/8] Creating force field...")
forcefield = app.ForceField(FORCE_FIELD, WATER_MODEL)

# Solvate directly (hydrogens already added by PDBFixer)
print("  Solvating...")
modeller = app.Modeller(pdb.topology, pdb.positions)
modeller.addSolvent(forcefield, model='tip3p', padding=1.0*unit.nanometers, ionicStrength=0.15*unit.molar)
print(f"  Atoms after solvation: {modeller.topology.getNumAtoms()}")

# Step 3: Create system
print("\n[3/8] Creating system (PME, 1nm cutoff)...")
system = forcefield.createSystem(
    modeller.topology,
    nonbondedMethod=app.PME,
    nonbondedCutoff=1.0 * unit.nanometers,
    constraints=app.HBonds,
    ewaldErrorTolerance=0.0005,
)

# Add barostat for NPT
system.addForce(openmm.MonteCarloBarostat(1.0 * unit.atmospheres, TEMPERATURE, 25))

# Step 4: Create simulation
print("\n[4/8] Setting up simulation...")
integrator = openmm.LangevinMiddleIntegrator(TEMPERATURE, FRICTION, TIMESTEP)

# Try CUDA first, fall back to OpenCL
for pname in ["CUDA", "OpenCL", "CPU"]:
    try:
        platform = openmm.Platform.getPlatformByName(pname)
        properties = {}
        if pname == "CUDA":
            properties = {"DeviceIndex": "0", "Precision": "mixed"}
        elif pname == "OpenCL":
            properties = {"DeviceIndex": "0", "Precision": "mixed"}
        simulation = Simulation(modeller.topology, system, integrator, platform, properties)
        print(f"  Platform: {pname}")
        break
    except Exception as e:
        print(f"  {pname} failed: {e}")
        continue
simulation.context.setPositions(modeller.positions)
simulation.context.setVelocitiesToTemperature(TEMPERATURE)

print(f"  Platform: {simulation.context.getPlatform().getName()}")
print(f"  Device: {properties.get('DeviceIndex', 'N/A')}")

# Step 5: Energy Minimization
print("\n[5/8] Energy minimization...")
t0 = time.time()
simulation.minimizeEnergy(maxIterations=5000, tolerance=10 * unit.kilojoule_per_mole / unit.nanometer)
state = simulation.context.getState(getEnergy=True)
pe = state.getPotentialEnergy()
print(f"  PE after minimization: {pe}")
print(f"  Time: {time.time() - t0:.1f}s")

# Step 6: NVT Equilibration (100ps)
print("\n[6/8] NVT equilibration (100ps)...")
simulation.context.setVelocitiesToTemperature(TEMPERATURE)
simulation.step(50000)  # 50000 * 2fs = 100ps
print("  NVT complete")

# Step 7: NPT Equilibration (100ps)
print("\n[7/8] NPT equilibration (100ps)...")
simulation.step(50000)  # 100ps
print("  NPT complete")

# Step 8: Production MD (5ns)
print("\n[8/8] Production MD (5ns = 2.5M steps)...")
simulation.reporters.append(DCDReporter(os.path.join(WORKDIR, "production.dcd"), 2500))
simulation.reporters.append(StateDataReporter(
    os.path.join(WORKDIR, "production.log"), 2500,
    step=True, potentialEnergy=True, temperature=True, density=True, speed=True
))

t0 = time.time()
TOTAL_STEPS = 2500000  # 5ns
REPORT_EVERY = 25000   # log every 50ps

for i in range(0, TOTAL_STEPS, REPORT_EVERY):
    simulation.step(REPORT_EVERY)
    elapsed = time.time() - t0
    ns_done = (i + REPORT_EVERY) * 2e-6  # 2fs timestep
    rate = ns_done / (elapsed / 3600) if elapsed > 0 else 0
    print(f"  {ns_done:.1f}ns / 5.0ns ({100*ns_done/5:.0f}%) - {rate:.1f} ns/hr - {elapsed:.0f}s elapsed")

total_time = time.time() - t0
print(f"\nProduction complete: {total_time:.0f}s ({total_time/3600:.2f}h)")
print(f"Performance: {5/(total_time/3600):.1f} ns/hr")

# Save final structure
state = simulation.context.getState(getPositions=True)
PDBFile.writeFile(modeller.topology, state.getPositions(), open(os.path.join(WORKDIR, "final.pdb"), "w"))

print(f"\nResults saved to: {WORKDIR}")
print(f"  production.dcd  - trajectory")
print(f"  production.log  - energy/temperature")
print(f"  final.pdb       - final structure")
