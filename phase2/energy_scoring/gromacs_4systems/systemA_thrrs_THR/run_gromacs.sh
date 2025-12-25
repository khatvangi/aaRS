#!/bin/bash
# GROMACS simulation script for systemA_thrrs_THR
# Run from: /storage/kiran-stuff/aaRS/phase2/energy_scoring/gromacs_4systems/systemA_thrrs_THR

set -e

GMX="gmx"  # or gmx_mpi

echo "=== Setting up systemA_thrrs_THR ==="

# 1. Generate protein topology
$GMX pdb2gmx -f protein.pdb -o protein_processed.gro -water tip3p -ff amber99sb-ildn -ignh

# 2. Create box
$GMX editconf -f protein_processed.gro -o boxed.gro -c -d 1.2 -bt cubic

# 3. Solvate
$GMX solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top

# 4. Add ions (neutralize)
$GMX grompp -f em.mdp -c solvated.gro -p topol.top -o ions.tpr -maxwarn 5
echo "SOL" | $GMX genion -s ions.tpr -o solvated_ions.gro -p topol.top -pname NA -nname CL -neutral

# 5. Energy minimization
$GMX grompp -f em.mdp -c solvated_ions.gro -p topol.top -o em.tpr -maxwarn 5
$GMX mdrun -v -deffnm em

# 6. NVT equilibration
$GMX grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn 5
$GMX mdrun -v -deffnm nvt

# 7. NPT equilibration
$GMX grompp -f npt.mdp -c nvt.gro -r nvt.gro -p topol.top -o npt.tpr -maxwarn 5
$GMX mdrun -v -deffnm npt

# 8. Production MD
$GMX grompp -f md.mdp -c npt.gro -p topol.top -o md.tpr -maxwarn 5
$GMX mdrun -v -deffnm md

echo "=== systemA_thrrs_THR complete ==="
