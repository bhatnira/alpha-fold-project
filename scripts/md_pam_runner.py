#!/usr/bin/env python3
"""
Prepare and run 100ns MD for PAM ligands at Site 23.
Ligands: ascorbic acid, O-ethyl ascorbate, ryanodine
"""
import os, sys, time, argparse
import openmm
from openmm import app, unit
from openmm.app import PDBFile, Simulation, DCDReporter, StateDataReporter
from pdbfixer import PDBFixer
import numpy as np

WORKDIR = "/cluster/home/nbhatt04/lean_pipeline/md_pam"
os.makedirs(WORKDIR, exist_ok=True)

LIGANDS = {
    "ascorbate": "OC[C@H](O)[C@H]1OC(=O)C(O)=C1O",
    "oethyl_ascorbate": "CCOC(=O)C(=O)OC[C@@H](O)[C@@H]1OC(=O)C(O)=C1O",
    "ryanodine": "CC1(C)CC(C)(C)C(=O)OC2CC(C)CCC3C(C)CCC4(O)C5CC6OC(O)C(C)C6OC5C43C21",
}

# Site 23 center (from existing Boltz-2 predictions, chain B alpha10 interface)
SITE23_CENTER = np.array([24.564, 17.239, -19.641])  # B:145 ASP CA

FORCE_FIELD = "amber14-all.xml"
WATER_MODEL = "amber14/tip3pfb.xml"
TEMPERATURE = 310 * unit.kelvin
FRICTION = 1.0 / unit.picoseconds
TIMESTEP = 2 * unit.femtoseconds
TOTAL_NS = 100
TOTAL_STEPS = TOTAL_NS * 500000
REPORT_EVERY = 50000  # 100ps

parser = argparse.ArgumentParser()
parser.add_argument("--ligand", type=str, default="ascorbate", choices=list(LIGANDS.keys()))
parser.add_argument("--n_gpus", type=int, default=4)
args = parser.parse_args()

lig_name = args.ligand
lig_smiles = LIGANDS[lig_name]
sim_dir = os.path.join(WORKDIR, lig_name)
os.makedirs(sim_dir, exist_ok=True)

print("=" * 60)
print(f"100ns GPU MD: {lig_name} at Site 23")
print("=" * 60)

# Load receptor structure (from Boltz-2 prediction with ligand)
src_pdb = "/cluster/home/nbhatt04/lean_pipeline/08_site_directed_af3/yaml_inputs_v2/predictions/boltz_results_sd_site23_ascorbate/predictions/sd_site23_ascorbate/sd_site23_ascorbate_model_0.pdb"

print("\n[1/6] Loading receptor structure...")
clean_pdb = os.path.join(sim_dir, "clean.pdb")
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

# Force field and solvation
print("\n[2/6] Creating system...")
forcefield = app.ForceField(FORCE_FIELD, WATER_MODEL)
modeller = app.Modeller(fixer.topology, fixer.positions)
modeller.addSolvent(forcefield, model='tip3p', padding=1.0*unit.nanometers, ionicStrength=0.15*unit.molar)
print(f"  Atoms after solvation: {modeller.topology.getNumAtoms()}")

# Create system
system = forcefield.createSystem(
    modeller.topology,
    nonbondedMethod=app.PME,
    nonbondedCutoff=1.0 * unit.nanometers,
    constraints=app.HBonds,
    ewaldErrorTolerance=0.0005,
)
system.addForce(openmm.MonteCarloBarostat(1.0 * unit.atmospheres, TEMPERATURE, 25))

# Create simulation
print("\n[3/6] Setting up simulation...")
integrator = openmm.LangevinMiddleIntegrator(TEMPERATURE, FRICTION, TIMESTEP)
platform = openmm.Platform.getPlatformByName("OpenCL")
gpu_ids = ",".join(str(i) for i in range(args.n_gpus))
properties = {"DeviceIndex": gpu_ids, "Precision": "mixed"}
simulation = Simulation(modeller.topology, system, integrator, platform, properties)
simulation.context.setPositions(modeller.positions)
simulation.context.setVelocitiesToTemperature(TEMPERATURE)
print(f"  Platform: OpenCL, GPUs: {gpu_ids}")

# Energy minimization
print("\n[4/6] Energy minimization...")
t0 = time.time()
simulation.minimizeEnergy(maxIterations=5000, tolerance=10*unit.kilojoule_per_mole/unit.nanometer)
print(f"  Minimized in {time.time()-t0:.1f}s")

# Equilibration
print("\n[5/6] Equilibration (200ps)...")
simulation.context.setVelocitiesToTemperature(TEMPERATURE)
simulation.step(50000)  # NVT 100ps
simulation.step(50000)  # NPT 100ps
print("  Equilibration complete")

# Production
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

print(f"\nResults: {sim_dir}")
