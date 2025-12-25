#!/usr/bin/env python3
"""
Setup all 8 GROMACS MD simulations for AF3/Vina pose validation.

4 paired comparisons:
1. Modern ThrRS +Zn: THR vs SER
2. Modern ThrRS -Zn: THR vs SER
3. Ancestral ProRS Edit: PRO vs ALA
4. Ancestral ThrRS Cat (no Zn): THR vs ILE

Uses CHARMM36m force field with TIP3P water.
Amino acid ligands parameterized by same force field (no GAFF).
"""

import subprocess
import shutil
from pathlib import Path

# Paths
VINA_DIR = Path("/storage/kiran-stuff/aaRS/phase2/energy_scoring/vina_182jobs")
WORK_DIR = Path("/storage/kiran-stuff/aaRS/phase2/energy_scoring/gromacs_md_validation")
MDP_DIR = WORK_DIR / "mdp"
GMX = "/usr/local/gromacs/bin/gmx_mpi"

# 8 simulations: (system_name, vina_job_prefix, ligand, has_zn)
SYSTEMS = [
    # Pair 1: Modern ThrRS +Zn
    ("modern_thrrs_zn_THR", "modern_thrrs_ecoli_zn", "THR", True),
    ("modern_thrrs_zn_SER", "modern_thrrs_ecoli_zn", "SER", True),
    # Pair 2: Modern ThrRS -Zn
    ("modern_thrrs_noZn_THR", "modern_thrrs", "THR", False),
    ("modern_thrrs_noZn_SER", "modern_thrrs", "SER", False),
    # Pair 3: Ancestral ProRS Edit
    ("anc_prors_edit_PRO", "anc_prors_edit", "PRO", False),
    ("anc_prors_edit_ALA", "anc_prors_edit", "ALA", False),
    # Pair 4: Ancestral ThrRS Cat (no Zn)
    ("anc_thrrs_cat_THR", "anc_thrrs_cat", "THR", False),
    ("anc_thrrs_cat_ILE", "anc_thrrs_cat", "ILE", False),
]


def run_cmd(cmd, cwd=None, check=True):
    """Run shell command"""
    print(f"  $ {' '.join(cmd)}")
    result = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    if check and result.returncode != 0:
        print(f"STDERR: {result.stderr}")
        raise RuntimeError(f"Command failed: {' '.join(cmd)}")
    return result


def prepare_input_pdb(vina_job_prefix, ligand, sys_dir, has_zn):
    """
    Prepare input PDB with protein as chain A and ligand as chain L.
    For Zn systems, keep Zn as chain Z.
    """
    vina_job = f"{vina_job_prefix}_{ligand}_{ligand}"
    complex_pdb = VINA_DIR / vina_job / "complex.pdb"

    if not complex_pdb.exists():
        raise FileNotFoundError(f"Missing: {complex_pdb}")

    protein_lines = []
    ligand_lines = []
    zn_lines = []

    with open(complex_pdb) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                chain = line[21]
                res = line[17:20].strip()

                if res == "ZN":
                    if has_zn:
                        # Keep Zn, assign to chain Z
                        new_line = line[:21] + "Z" + line[22:]
                        zn_lines.append(new_line)
                elif chain == "A":
                    protein_lines.append(line)
                elif chain == "B":
                    # Ligand - assign to chain L
                    new_line = line[:21] + "L" + line[22:]
                    ligand_lines.append(new_line)

    # Write combined PDB
    output_pdb = sys_dir / "input.pdb"
    with open(output_pdb, "w") as f:
        f.writelines(protein_lines)
        f.write("TER\n")
        f.writelines(ligand_lines)
        f.write("TER\n")
        if zn_lines:
            f.writelines(zn_lines)
            f.write("TER\n")
        f.write("END\n")

    # Also write separate files for pdb2gmx
    with open(sys_dir / "protein.pdb", "w") as f:
        f.writelines(protein_lines)
        f.write("END\n")

    with open(sys_dir / "ligand.pdb", "w") as f:
        f.writelines(ligand_lines)
        f.write("END\n")

    if zn_lines:
        with open(sys_dir / "zn.pdb", "w") as f:
            f.writelines(zn_lines)
            f.write("END\n")

    return len(protein_lines), len(ligand_lines), len(zn_lines)


def setup_system(name, vina_prefix, ligand, has_zn):
    """Setup a single system"""
    print(f"\n{'='*60}")
    print(f"Setting up: {name}")
    print(f"{'='*60}")

    sys_dir = WORK_DIR / "systems" / name
    sys_dir.mkdir(parents=True, exist_ok=True)

    # 1. Prepare input PDB
    n_prot, n_lig, n_zn = prepare_input_pdb(vina_prefix, ligand, sys_dir, has_zn)
    print(f"  Protein atoms: {n_prot}")
    print(f"  Ligand atoms: {n_lig}")
    print(f"  Zn atoms: {n_zn}")

    # 2. Copy MDP files
    for mdp in ["minim.mdp", "nvt.mdp", "npt.mdp", "md.mdp"]:
        shutil.copy(MDP_DIR / mdp, sys_dir / mdp)

    # 3. Create run script
    create_run_script(sys_dir, name, has_zn)

    print(f"  Created: {sys_dir}")
    return sys_dir


def create_run_script(sys_dir, name, has_zn):
    """Create the GROMACS run script for this system"""

    zn_section = ""
    if has_zn:
        zn_section = """
# Handle Zn - add to topology manually and create position restraints
echo "Adding Zn ion to topology..."
# Zn will be added as ion during genion step or manually
ZN_POSRES="-DPOSRES_ZN"
"""
    else:
        zn_section = 'ZN_POSRES=""'

    script = f'''#!/bin/bash
#
# GROMACS MD run script for {name}
# Uses CHARMM36m + TIP3P
#

set -e
cd "$(dirname "$0")"

GMX="/usr/local/gromacs/bin/gmx_mpi"
NCORES=32

echo "============================================================"
echo "System: {name}"
echo "============================================================"

# 1. Process protein with pdb2gmx (CHARMM36m)
echo "Step 1: Processing protein..."
$GMX pdb2gmx -f protein.pdb -o protein_processed.gro -water tip3p -ff charmm36m -ignh <<EOF
1
EOF

# 2. Process ligand (amino acid - same force field)
echo "Step 2: Processing ligand..."
# For amino acid ligands, we can use pdb2gmx with termini
$GMX pdb2gmx -f ligand.pdb -o ligand_processed.gro -water tip3p -ff charmm36m -ignh -ter <<EOF
1
1
1
EOF

# Rename ligand topology
mv topol.top ligand.top
mv posre.itp posre_ligand.itp

# 3. Combine protein and ligand
echo "Step 3: Combining structures..."
# Use editconf to combine (manual merge needed)
cat protein_processed.gro > complex.gro
# Append ligand coordinates (skip header/footer)
tail -n +3 ligand_processed.gro | head -n -1 >> complex_temp.gro

# Create combined topology
cat > topol.top << 'TOPEOF'
; Combined topology for {name}
#include "charmm36m.ff/forcefield.itp"
#include "charmm36m.ff/tip3p.itp"

; Include protein topology
#include "topol_Protein_chain_A.itp"

; Position restraints for protein
#ifdef POSRES
#include "posre_Protein_chain_A.itp"
#endif

; Include ligand topology
#include "topol_Protein_chain_L.itp"

; Position restraints for ligand
#ifdef POSRES_LIG
#include "posre_ligand.itp"
#endif

#include "charmm36m.ff/ions.itp"

[ system ]
{name}

[ molecules ]
Protein_chain_A 1
Protein_chain_L 1
TOPEOF

{zn_section}

# 4. Create simulation box
echo "Step 4: Creating box..."
$GMX editconf -f complex.gro -o boxed.gro -c -d 1.2 -bt dodecahedron

# 5. Solvate
echo "Step 5: Solvating..."
$GMX solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top

# 6. Add ions
echo "Step 6: Adding ions..."
$GMX grompp -f minim.mdp -c solvated.gro -p topol.top -o ions.tpr -maxwarn 5
echo "SOL" | $GMX genion -s ions.tpr -o ionized.gro -p topol.top -pname NA -nname CL -neutral

# 7. Energy minimization
echo "Step 7: Energy minimization..."
$GMX grompp -f minim.mdp -c ionized.gro -p topol.top -o em.tpr -maxwarn 5
$GMX mdrun -v -deffnm em -ntmpi 1 -ntomp $NCORES

# 8. NVT equilibration
echo "Step 8: NVT equilibration (250 ps)..."
$GMX grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn 5
$GMX mdrun -v -deffnm nvt -ntmpi 1 -ntomp $NCORES

# 9. NPT equilibration
echo "Step 9: NPT equilibration (500 ps)..."
$GMX grompp -f npt.mdp -c nvt.gro -r nvt.gro -p topol.top -o npt.tpr -maxwarn 5
$GMX mdrun -v -deffnm npt -ntmpi 1 -ntomp $NCORES

# 10. Production MD
echo "Step 10: Production MD (5 ns)..."
$GMX grompp -f md.mdp -c npt.gro -p topol.top -o md.tpr -maxwarn 5
$GMX mdrun -v -deffnm md -ntmpi 1 -ntomp $NCORES -nb gpu -pme gpu -bonded gpu

echo "============================================================"
echo "DONE: {name}"
echo "============================================================"
'''

    script_path = sys_dir / "run_md.sh"
    with open(script_path, "w") as f:
        f.write(script)

    script_path.chmod(0o755)


def main():
    print("=" * 60)
    print("GROMACS MD Validation Setup")
    print("8 simulations (4 paired comparisons)")
    print("=" * 60)

    # Create directories
    (WORK_DIR / "systems").mkdir(exist_ok=True)

    # Setup each system
    for name, vina_prefix, ligand, has_zn in SYSTEMS:
        try:
            setup_system(name, vina_prefix, ligand, has_zn)
        except Exception as e:
            print(f"  ERROR: {e}")

    print("\n" + "=" * 60)
    print("SETUP COMPLETE")
    print("=" * 60)
    print(f"\nSystems created in: {WORK_DIR / 'systems'}")
    print("\nTo run a simulation:")
    print("  cd systems/<name>")
    print("  ./run_md.sh")


if __name__ == "__main__":
    main()
