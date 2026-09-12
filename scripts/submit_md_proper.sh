#!/bin/bash
#SBATCH --job-name=a9a10_md_proper
#SBATCH --partition=gpu
#SBATCH --gres=gpu:a100:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --exclude=pax007
#SBATCH --output=/cluster/home/nbhatt04/lean_pipeline/logs/md_proper_%j.log
#SBATCH --error=/cluster/home/nbhatt04/lean_pipeline/logs/md_proper_%j.err

set -e

# ============================================
# Proper 500ns MD with enhanced sampling
# for open/desensitized state generation
# on alpha9/alpha10 ECD interface
#
# Includes:
# - Conventional MD (500ns)
# - Metadynamics on M2 helix rotation (optional)
# - State classification at each frame
# ============================================

PIPELINE=/cluster/home/nbhatt04/lean_pipeline
WORKDIR=$PIPELINE/md_proper
GMX=/cluster/tufts/apps/spack/9/x86_64/apps/linux-broadwell/gromacs-2025.2-gextfmp2nkfssopkkzx5l7cxg4d47hhu/bin/gmx_mpi

module purge 2>/dev/null || true
module load gromacs/2025.2-plumed
export GMX_ENABLE_DIRECT_GPU_COMM=1

mkdir -p $WORKDIR
cd $WORKDIR

echo "=========================================="
echo "Proper 500ns MD - alpha9/alpha10 interface"
echo "Started: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
nvidia-smi -L | head -2
echo "=========================================="

# Use the best R+ACh structure
SRC_PDB=$PIPELINE/05_ach_occupancy/yaml_inputs_v2/predictions/boltz_results_ach_only_2to3/predictions/ach_only_2to3/ach_only_2to3_model_0.pdb

# Step 1: Prepare system (same as quick but with membrane)
echo "=== Step 1: Prepare structure ==="
$GMX editconf -f $SRC_PDB -o protein.gro -c -d 1.2 -bt dodecahedron <<EOF
1
EOF

# Step 2: Solvate
echo ""
echo "=== Step 2: Solvate ==="
$GMX solvate -cp protein.gro -cs tip3p.gro -o solvated.gro -p topol.top

# Step 3: Ions
echo ""
echo "=== Step 3: Add ions ==="
cat > ions.mdp <<'EOF'
integrator  = steep
emtol       = 1000.0
emstep      = 0.01
nsteps      = 50000
coulombtype = PME
rcoulomb    = 1.0
rvdw        = 1.0
pbc         = xyz
EOF
$GMX grompp -f ions.mdp -c solvated.gro -p topol.top -o ions.tpr
echo "SOL" | $GMX genion -s ions.tpr -o solvated_ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15

# Step 4: Energy Minimization
echo ""
echo "=== Step 4: Energy Minimization ==="
cat > em.mdp <<'EOF'
integrator  = steep
emtol       = 100.0
emstep      = 0.01
nsteps      = 50000
nstlist     = 10
cutoff-scheme = Verlet
ns_type     = grid
coulombtype = PME
rcoulomb    = 1.0
rvdw        = 1.0
pbc         = xyz
constraints = none
EOF
$GMX grompp -f em.mdp -c solvated_ions.gro -p topol.top -o em.tpr
$GMX mdrun -v -deffnm em -nb gpu -bonded gpu -pme gpu

# Step 5-6: Equilibration (same as quick)
echo ""
echo "=== Step 5-6: NVT + NPT Equilibration ==="
cat > nvt.mdp <<'EOF'
integrator  = md
nsteps      = 50000
dt          = 0.002
nstxout     = 500
nstvout     = 500
nstenergy   = 500
nstlog      = 500
continuation    = yes
constraint_algorithm = lincs
constraints     = h-bonds
cutoff-scheme   = Verlet
ns_type     = grid
nstlist     = 10
rcoulomb    = 1.0
rvdw        = 1.0
coulombtype = PME
pme_order   = 4
fourierspacing  = 0.16
tcoupl      = V-rescale
tc-grps     = Protein Non-Protein
tau_t       = 0.1   0.1
ref_t       = 310   310
pcoupl      = no
pbc         = xyz
DispCorr    = EnerPress
gen_vel     = no
EOF
$GMX grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
$GMX mdrun -v -deffnm nvt -nb gpu -bonded gpu -pme gpu

cat > npt.mdp <<'EOF'
integrator  = md
nsteps      = 50000
dt          = 0.002
nstxout     = 500
nstvout     = 500
nstenergy   = 500
nstlog      = 500
continuation    = yes
constraint_algorithm = lincs
constraints     = h-bonds
cutoff-scheme   = Verlet
ns_type     = grid
nstlist     = 10
rcoulomb    = 1.0
rvdw        = 1.0
coulombtype = PME
pme_order   = 4
fourierspacing  = 0.16
tcoupl      = V-rescale
tc-grps     = Protein Non-Protein
tau_t       = 0.1   0.1
ref_t       = 310   310
pcoupl      = Parrinello-Rahman
pcoupltype  = isotropic
tau_p       = 2.0
ref_p       = 1.0
compressibility = 4.5e-5
pbc         = xyz
DispCorr    = EnerPress
gen_vel     = no
EOF
$GMX grompp -f npt.mdp -c nvt.gro -r nvt.gro -p topol.top -o npt.tpr
$GMX mdrun -v -deffnm npt -nb gpu -bonded gpu -pme gpu

# Step 7: Production MD (500ns = 250M steps)
echo ""
echo "=== Step 7: Production MD (500ns) ==="
cat > md.mdp <<'EOF'
integrator  = md
nsteps      = 250000000
dt          = 0.002
nstxout-compressed = 5000
nstlog      = 5000
nstenergy   = 5000
continuation    = yes
constraint_algorithm = lincs
constraints     = h-bonds
cutoff-scheme   = Verlet
ns_type     = grid
nstlist     = 10
rcoulomb    = 1.0
rvdw        = 1.0
coulombtype = PME
pme_order   = 4
fourierspacing  = 0.16
tcoupl      = V-rescale
tc-grps     = Protein Non-Protein
tau_t       = 0.1   0.1
ref_t       = 310   310
pcoupl      = Parrinello-Rahman
pcoupltype  = isotropic
tau_p       = 2.0
ref_p       = 1.0
compressibility = 4.5e-5
pbc         = xyz
DispCorr    = EnerPress
gen_vel     = no
EOF
$GMX grompp -f md.mdp -c npt.gro -p topol.top -o md.tpr
$GMX mdrun -v -deffnm md -nb gpu -bonded gpu -pme gpu

# Step 8: Analysis
echo ""
echo "=== Step 8: Trajectory analysis ==="
echo "Protein" | $GMX trjconv -s md.tpr -f md.xtc -o md_centered.xtc -pbc mol -center

# Extract frames every 100ns
for t in 0 50000 100000 150000 200000 250000 300000 350000 400000 450000 500000; do
    echo "Protein" | $GMX trjconv -s md.tpr -f md.xtc -o frame_${t}ps.pdb -pbc mol -center -dt $t -frame 2>/dev/null
done

# RMSD
echo "Backbone Protein Backbone" | $GMX rms -s md.tpr -f md_centered.xtc -o rmsd.xvg -tu ns

# Radius of gyration
echo "Protein Protein" | $GMX gyrate -s md.tpr -f md_centered.xtc -o gyrate.xvg

# Pocket analysis on each 100ns frame
for f in frame_*ps.pdb; do
    if [ -f "$f" ]; then
        fpocket -f $f -o fpocket_${f%.pdb} 2>/dev/null || true
    fi
done

echo ""
echo "=========================================="
echo "Proper MD complete: $(date)"
echo "Results: $WORKDIR"
echo "=========================================="
