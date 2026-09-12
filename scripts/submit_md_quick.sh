#!/bin/bash
#SBATCH --job-name=a9a10_md_quick
#SBATCH --partition=gpu
#SBATCH --gres=gpu:a100:1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=10:00:00
#SBATCH --exclude=pax007
#SBATCH --ntasks=1
#SBATCH --output=/cluster/home/nbhatt04/lean_pipeline/logs/md_quick_%j.log
#SBATCH --error=/cluster/home/nbhatt04/lean_pipeline/logs/md_quick_%j.err

set -e

PIPELINE=/cluster/home/nbhatt04/lean_pipeline
WORKDIR=$PIPELINE/md_quick
GMX=/cluster/tufts/apps/spack/9/x86_64/apps/linux-broadwell/gromacs-2025.2-gextfmp2nkfssopkkzx5l7cxg4d47hhu/bin/gmx_mpi

module purge 2>/dev/null || true
module load gromacs/2025.2

mkdir -p $WORKDIR
cd $WORKDIR

echo "=========================================="
echo "Quick 5ns MD - alpha9/alpha10 interface (CPU)"
echo "Started: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "=========================================="

SRC_PDB=$PIPELINE/05_ach_occupancy/yaml_inputs_v2/predictions/boltz_results_ach_only_2to3/predictions/ach_only_2to3/ach_only_2to3_model_0.pdb

echo ""
echo "=== Step 1: Clean PDB ==="
grep "^ATOM" $SRC_PDB > clean_protein.pdb

echo ""
echo "=== Step 2: Generate topology ==="
echo "1" | $GMX pdb2gmx -f clean_protein.pdb -o protein.gro -p topol.top -water tip3p -ff amber99sb-ildn

echo ""
echo "=== Step 3: Create box ==="
$GMX editconf -f protein.gro -o boxed.gro -c -d 1.2 -bt dodecahedron

echo ""
echo "=== Step 4: Solvate ==="
$GMX solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top

echo ""
echo "=== Step 5: Add ions ==="
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
$GMX grompp -f ions.mdp -c solvated.gro -p topol.top -o ions.tpr -maxwarn 2
echo "SOL" | $GMX genion -s ions.tpr -o solvated_ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15

echo ""
echo "=== Step 6: Energy Minimization ==="
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
$GMX grompp -f em.mdp -c solvated_ions.gro -p topol.top -o em.tpr -maxwarn 2
$GMX mdrun -v -deffnm em -nb cpu -pme cpu

echo ""
echo "=== Step 7: NVT Equilibration (100ps) ==="
cat > nvt.mdp <<'EOF'
integrator  = md
nsteps      = 50000
dt          = 0.002
nstxout-compressed = 1000
nstlog      = 1000
nstenergy   = 1000
continuation    = no
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
DispCorr    = EnerPres
gen_vel     = yes
gen_seed    = -1
EOF
$GMX grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn 2
$GMX mdrun -v -deffnm nvt -nb cpu -pme cpu

echo ""
echo "=== Step 8: NPT Equilibration (100ps) ==="
cat > npt.mdp <<'EOF'
integrator  = md
nsteps      = 50000
dt          = 0.002
nstxout-compressed = 1000
nstlog      = 1000
nstenergy   = 1000
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
DispCorr    = EnerPres
gen_vel     = no
EOF
$GMX grompp -f npt.mdp -c nvt.gro -r nvt.gro -p topol.top -o npt.tpr -maxwarn 2
$GMX mdrun -v -deffnm npt -nb cpu -pme cpu

echo ""
echo "=== Step 9: Production MD (5ns = 2.5M steps) ==="
cat > md.mdp <<'EOF'
integrator  = md
nsteps      = 2500000
dt          = 0.002
nstxout-compressed = 2500
nstlog      = 2500
nstenergy   = 2500
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
DispCorr    = EnerPres
gen_vel     = no
EOF
$GMX grompp -f md.mdp -c npt.gro -p topol.top -o md.tpr -maxwarn 2
$GMX mdrun -v -deffnm md -nb cpu -pme cpu

echo ""
echo "=== Step 10: Analysis ==="
echo "Protein" | $GMX trjconv -s md.tpr -f md.xtc -o md_centered.xtc -pbc mol -center

for t in 0 500000 1000000 1500000 2000000 2500000; do
    ns=$((t / 1000000))
    echo "Protein" | $GMX trjconv -s md.tpr -f md.xtc -o frame_${ns}ns.pdb -pbc mol -center -dt $t -frame 2>/dev/null || true
done

echo "Backbone Protein Backbone" | $GMX rms -s md.tpr -f md_centered.xtc -o rmsd.xvg -tu ns
echo "Backbone Backbone" | $GMX rmsf -s md.tpr -f md_centered.xtc -o rmsf.xvg -res
echo "Protein Protein" | $GMX gyrate -s md.tpr -f md_centered.xtc -o gyrate.xvg

echo ""
echo "=== Step 11: Pocket analysis ==="
which fpocket >/dev/null 2>&1 && for f in frame_*ns.pdb; do
    [ -f "$f" ] && fpocket -f $f -o fpocket_${f%.pdb} 2>/dev/null || true
done || echo "fpocket not in PATH"

echo ""
echo "=========================================="
echo "Quick MD complete: $(date)"
echo "Results: $WORKDIR"
echo "=========================================="
