#!/usr/bin/env python3
"""
100ns GPU MD runner for alpha9alpha10 nAChR.
Supports: apo_2to3, ach_2to3, apo_3to2
"""
import os, sys, time, argparse
import openmm
from openmm import app, unit
from openmm.app import PDBFile, Simulation, DCDReporter, StateDataReporter
from pdbfixer import PDBFixer

WORKDIR = "/cluster/home/nbhatt04/lean_pipeline/md_100ns"
os.makedirs(WORKDIR, exist_ok=True)
os.makedirs(os.path.join(WORKDIR, "logs"), exist_ok=True)

PDB_SOURCES = {
    "apo_2to3": "/cluster/home/nbhatt04/lean_pipeline/05_ach_occupancy/yaml_inputs_v2/predictions/boltz_results_ach_only_2to3/predictions/ach_only_2to3/ach_only_2to3_model_0.pdb",
    "ach_2to3": "/cluster/home/nbhatt04/lean_pipeline/05_ach_occupancy/yaml_inputs_v2/predictions/boltz_results_ach_only_2to3/predictions/ach_only_2to3/ach_only_2to3_model_0.pdb",
    "apo_3to2": "/cluster/home/nbhatt04/lean_pipeline/05_ach_occupancy/yaml_inputs_v2/predictions/boltz_results_ach_only_3to2/predictions/ach_only_3to2/ach_only_3to2_model_0.pdb",
}

FORCE_FIELD = "amber14-all.xml"
WATER_MODEL = "amber14/tip3pfb.xml"
TEMPERATURE = 310 * unit.kelvin
FRICTION = 1.0 / unit.picoseconds
TIMESTEP = 2 * unit.femtoseconds
TOTAL_NS = 100
TOTAL_STEPS = TOTAL_NS * 500000  # 100ns at 2fs
REPORT_EVERY = 50000  # 100ps

parser = argparse.ArgumentParser()
parser.add_argument("--sim", type=str, default="apo_2to3", choices=list(PDB_SOURCES.keys()))
args = parser.parse_args()

sim_name = args.sim
sim_dir = os.path.join(WORKDIR, sim_name)
os.makedirs(sim_dir, exist_ok=True)

print("=" * 60)
print(f"100ns GPU MD: {sim_name}")
print("=" * 60)

# Step 1: Load and fix structure
print("\n[1/6] Loading and fixing structure...")
src_pdb = PDB_SOURCES[sim_name]
clean_pdb = os.path.join(sim_dir, "clean.pdb")

# Strip HETATM
with open(src_pdb) as fin, open(clean_pdb, "w") as fout:
    for line in fin:
        if line.startswith("ATOM"):
            fout.write(line)

fixer = PDBFixer(filename=clean_pdb)
fixer.findMissingResidues()
fixer.findMissingAtoms()
fixer.addMissingAtoms()
fixer.addMissingHydrogens(7.4)

print(f"  Atoms: {fixer.topology.getNumAtoms()}")

# Step 2: Force field and solvation
print("\n[2/6] Creating system...")
forcefield = app.ForceField(FORCE_FIELD, WATER_MODEL)
modeller = app.Modeller(fixer.topology, fixer.positions)
modeller.addSolvent(forcefield, model='tip3p', padding=1.0*unit.nanometers, ionicStrength=0.15*unit.molar)
print(f"  Atoms after solvation: {modeller.topology.getNumAtoms()}")

# Step 3: Create system
system = forcefield.createSystem(
    modeller.topology,
    nonbondedMethod=app.PME,
    nonbondedCutoff=1.0 * unit.nanometers,
    constraints=app.HBonds,
    ewaldErrorTolerance=0.0005,
)
system.addForce(openmm.MonteCarloBarostat(1.0 * unit.atmospheres, TEMPERATURE, 25))

# Step 4: Create simulation
print("\n[3/6] Setting up simulation...")
integrator = openmm.LangevinMiddleIntegrator(TEMPERATURE, FRICTION, TIMESTEP)
platform = openmm.Platform.getPlatformByName("OpenCL")
properties = {"DeviceIndex": "0,1,2,3", "Precision": "mixed"}
simulation = Simulation(modeller.topology, system, integrator, platform, properties)
simulation.context.setPositions(modeller.positions)
simulation.context.setVelocitiesToTemperature(TEMPERATURE)
print(f"  Platform: {simulation.context.getPlatform().getName()}")
print(f"  Devices: {properties['DeviceIndex']}")

# Step 5: Energy minimization
print("\n[4/6] Energy minimization...")
t0 = time.time()
simulation.minimizeEnergy(maxIterations=5000, tolerance=10*unit.kilojoule_per_mole/unit.nanometer)
print(f"  Minimized in {time.time()-t0:.1f}s")

# Step 6: Equilibration (NVT 100ps + NPT 100ps)
print("\n[5/6] Equilibration (200ps)...")
simulation.context.setVelocitiesToTemperature(TEMPERATURE)
simulation.step(50000)  # NVT 100ps
simulation.step(50000)  # NPT 100ps
print("  Equilibration complete")

# Step 7: Production (100ns)
print(f"\n[6/6] Production MD ({TOTAL_NS}ns)...")
simulation.reporters.append(DCDReporter(os.path.join(sim_dir, "production.dcd"), REPORT_EVERY))
simulation.reporters.append(StateDataReporter(
    os.path.join(sim_dir, "production.log"), REPORT_EVERY,
    step=True, potentialEnergy=True, temperature=True, density=True, speed=True
))

t0 = time.time()
for i in range(0, TOTAL_STEPS, REPORT_EVERY):
    simulation.step(REPORT_EVERY)
    elapsed = time.time() - t0
    ns_done = (i + REPORT_EVERY) * 2e-6
    rate = ns_done / (elapsed / 86400) if elapsed > 0 else 0
    eta_h = (TOTAL_NS - ns_done) / (rate / 24) if rate > 0 else 0
    print(f"  {ns_done:.1f}/{TOTAL_NS}ns ({100*ns_done/TOTAL_NS:.0f}%) - {rate:.1f} ns/day - ETA {eta_h:.1f}h")

total_time = time.time() - t0
print(f"\nComplete: {total_time/3600:.2f}h ({TOTAL_NS/(total_time/86400):.1f} ns/day)")

# Save final structure
state = simulation.context.getState(getPositions=True)
PDBFile.writeFile(modeller.topology, state.getPositions(), open(os.path.join(sim_dir, "final.pdb"), "w"))

# Summary
print(f"\nResults: {sim_dir}")
print(f"  production.dcd  - trajectory")
print(f"  production.log  - energy/temperature")
print(f"  final.pdb       - final structure")
