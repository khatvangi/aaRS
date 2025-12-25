#!/usr/bin/env python3
"""
setup_gromacs.py

Set up 4 GROMACS systems for MM-PBSA/GBSA energy calculations:
  A: modern ThrRS+Zn + THR (cognate)
  B: modern ThrRS+Zn + SER (near-cognate)
  C: modern ThrRS+Zn + ILE (non-cognate)
  D: ProRS editing domain + PRO (cognate for editing)

Steps:
1. Convert CIF to PDB
2. Separate protein, ligand, Zn
3. Generate GROMACS topology with pdb2gmx
4. Add ligand parameters
5. Create simulation mdp files
"""

import os
import shutil
from pathlib import Path

BASE_DIR = Path("/storage/kiran-stuff/aaRS/phase2/energy_scoring/gromacs_4systems")

# Source CIF files
SYSTEMS = {
    "systemA_thrrs_THR": {
        "cif": "/storage/kiran-stuff/aaRS/phase2/af3_modern_matrix/thrrs_ecoli_zn_jobs/af3_output/modern_thrrs_ecoli_zn_THR/modern_thrrs_ecoli_zn_THR_model.cif",
        "ligand": "THR",
        "has_zn": True,
        "description": "Modern ThrRS + Zn + THR (cognate)"
    },
    "systemB_thrrs_SER": {
        "cif": "/storage/kiran-stuff/aaRS/phase2/af3_modern_matrix/thrrs_ecoli_zn_jobs/af3_output/modern_thrrs_ecoli_zn_SER/modern_thrrs_ecoli_zn_SER_model.cif",
        "ligand": "SER",
        "has_zn": True,
        "description": "Modern ThrRS + Zn + SER (near-cognate)"
    },
    "systemC_thrrs_ILE": {
        "cif": "/storage/kiran-stuff/aaRS/phase2/af3_modern_matrix/thrrs_ecoli_zn_jobs/af3_output/modern_thrrs_ecoli_zn_ILE/modern_thrrs_ecoli_zn_ILE_model.cif",
        "ligand": "ILE",
        "has_zn": True,
        "description": "Modern ThrRS + Zn + ILE (non-cognate)"
    },
    "systemD_prors_edit": {
        "cif": "/storage/kiran-stuff/aaRS/phase2/batch_local_domains/af3_output/anc_prors_edit_PRO/anc_prors_edit_PRO_model.cif",
        "ligand": "PRO",
        "has_zn": False,
        "description": "Ancestral ProRS editing domain + PRO"
    }
}


def cif_to_pdb(cif_path, pdb_path):
    """Convert CIF to PDB using gemmi or biopython"""
    try:
        import gemmi
        structure = gemmi.read_structure(str(cif_path))
        structure.write_pdb(str(pdb_path))
        return True
    except ImportError:
        # Fallback to biopython
        from Bio.PDB import MMCIFParser, PDBIO
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure("struct", str(cif_path))
        io = PDBIO()
        io.set_structure(structure)
        io.save(str(pdb_path))
        return True


def separate_components(pdb_path, output_dir, ligand_code, has_zn):
    """Separate protein, ligand, and optionally Zn from PDB"""

    protein_lines = []
    ligand_lines = []
    zn_lines = []

    with open(pdb_path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                resname = line[17:20].strip()

                if resname == ligand_code:
                    ligand_lines.append(line)
                elif resname == "ZN":
                    zn_lines.append(line)
                elif resname in ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU",
                                "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE",
                                "PRO", "SER", "THR", "TRP", "TYR", "VAL"]:
                    # Standard amino acid - part of protein
                    protein_lines.append(line)
                else:
                    # Other residues go to protein
                    protein_lines.append(line)

    # Write protein
    with open(output_dir / "protein.pdb", "w") as f:
        f.writelines(protein_lines)
        f.write("END\n")

    # Write ligand
    with open(output_dir / "ligand.pdb", "w") as f:
        f.writelines(ligand_lines)
        f.write("END\n")

    # Write Zn if present
    if has_zn and zn_lines:
        with open(output_dir / "zn.pdb", "w") as f:
            f.writelines(zn_lines)
            f.write("END\n")

    return len(protein_lines), len(ligand_lines), len(zn_lines)


def write_mdp_files(output_dir):
    """Write GROMACS mdp parameter files"""

    # Energy minimization
    em_mdp = """
; Energy minimization parameters
integrator  = steep
emtol       = 1000.0
emstep      = 0.01
nsteps      = 50000

; Neighbor searching
nstlist     = 10
cutoff-scheme = Verlet
ns_type     = grid
rlist       = 1.2

; Electrostatics
coulombtype = PME
rcoulomb    = 1.2

; VdW
vdwtype     = Cut-off
rvdw        = 1.2

; Constraints
constraints = none
"""

    # NVT equilibration
    nvt_mdp = """
; NVT equilibration
integrator  = md
dt          = 0.002
nsteps      = 50000  ; 100 ps

; Output
nstxout     = 5000
nstvout     = 5000
nstenergy   = 500
nstlog      = 500

; Neighbor searching
nstlist     = 10
cutoff-scheme = Verlet
ns_type     = grid
rlist       = 1.2

; Electrostatics
coulombtype = PME
rcoulomb    = 1.2

; VdW
vdwtype     = Cut-off
rvdw        = 1.2

; Temperature coupling
tcoupl      = V-rescale
tc-grps     = System
tau_t       = 0.1
ref_t       = 300

; Pressure coupling (off for NVT)
pcoupl      = no

; Constraints
constraints = h-bonds
constraint_algorithm = lincs

; Velocity generation
gen_vel     = yes
gen_temp    = 300
gen_seed    = -1
"""

    # NPT equilibration
    npt_mdp = """
; NPT equilibration
integrator  = md
dt          = 0.002
nsteps      = 50000  ; 100 ps

; Output
nstxout     = 5000
nstvout     = 5000
nstenergy   = 500
nstlog      = 500

; Neighbor searching
nstlist     = 10
cutoff-scheme = Verlet
ns_type     = grid
rlist       = 1.2

; Electrostatics
coulombtype = PME
rcoulomb    = 1.2

; VdW
vdwtype     = Cut-off
rvdw        = 1.2

; Temperature coupling
tcoupl      = V-rescale
tc-grps     = System
tau_t       = 0.1
ref_t       = 300

; Pressure coupling
pcoupl      = Parrinello-Rahman
pcoupltype  = isotropic
tau_p       = 2.0
ref_p       = 1.0
compressibility = 4.5e-5

; Constraints
constraints = h-bonds
constraint_algorithm = lincs

; Velocity generation
gen_vel     = no
"""

    # Production MD
    md_mdp = """
; Production MD
integrator  = md
dt          = 0.002
nsteps      = 500000  ; 1 ns

; Output
nstxout     = 5000
nstvout     = 5000
nstxtcout   = 1000   ; Write xtc every 2 ps
nstenergy   = 1000
nstlog      = 1000

; Neighbor searching
nstlist     = 10
cutoff-scheme = Verlet
ns_type     = grid
rlist       = 1.2

; Electrostatics
coulombtype = PME
rcoulomb    = 1.2

; VdW
vdwtype     = Cut-off
rvdw        = 1.2

; Temperature coupling
tcoupl      = V-rescale
tc-grps     = System
tau_t       = 0.1
ref_t       = 300

; Pressure coupling
pcoupl      = Parrinello-Rahman
pcoupltype  = isotropic
tau_p       = 2.0
ref_p       = 1.0
compressibility = 4.5e-5

; Constraints
constraints = h-bonds
constraint_algorithm = lincs

; Velocity generation
gen_vel     = no
"""

    with open(output_dir / "em.mdp", "w") as f:
        f.write(em_mdp.strip())

    with open(output_dir / "nvt.mdp", "w") as f:
        f.write(nvt_mdp.strip())

    with open(output_dir / "npt.mdp", "w") as f:
        f.write(npt_mdp.strip())

    with open(output_dir / "md.mdp", "w") as f:
        f.write(md_mdp.strip())


def write_run_script(output_dir, system_name, has_zn):
    """Write bash script to run GROMACS simulation"""

    script = f"""#!/bin/bash
# GROMACS simulation script for {system_name}
# Run from: {output_dir}

set -e

GMX="gmx"  # or gmx_mpi

echo "=== Setting up {system_name} ==="

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

echo "=== {system_name} complete ==="
"""

    with open(output_dir / "run_gromacs.sh", "w") as f:
        f.write(script)

    os.chmod(output_dir / "run_gromacs.sh", 0o755)


def main():
    print("="*60)
    print("Setting up 4 GROMACS systems")
    print("="*60)

    for system_name, config in SYSTEMS.items():
        print(f"\n--- {system_name}: {config['description']} ---")

        output_dir = BASE_DIR / system_name
        output_dir.mkdir(exist_ok=True)

        # Check source CIF exists
        cif_path = Path(config["cif"])
        if not cif_path.exists():
            print(f"  ERROR: CIF not found: {cif_path}")
            continue

        # Copy source CIF
        shutil.copy(cif_path, output_dir / "source.cif")
        print(f"  Copied: {cif_path.name}")

        # Convert CIF to PDB
        pdb_path = output_dir / "complex.pdb"
        try:
            cif_to_pdb(cif_path, pdb_path)
            print(f"  Converted to PDB")
        except Exception as e:
            print(f"  ERROR converting CIF: {e}")
            continue

        # Separate components
        n_prot, n_lig, n_zn = separate_components(
            pdb_path, output_dir,
            config["ligand"], config["has_zn"]
        )
        print(f"  Separated: protein={n_prot} atoms, ligand={n_lig} atoms, Zn={n_zn}")

        # Write mdp files
        write_mdp_files(output_dir)
        print(f"  Created: em.mdp, nvt.mdp, npt.mdp, md.mdp")

        # Write run script
        write_run_script(output_dir, system_name, config["has_zn"])
        print(f"  Created: run_gromacs.sh")

        # Write info file
        with open(output_dir / "INFO.txt", "w") as f:
            f.write(f"System: {system_name}\n")
            f.write(f"Description: {config['description']}\n")
            f.write(f"Source CIF: {config['cif']}\n")
            f.write(f"Ligand: {config['ligand']}\n")
            f.write(f"Has Zn: {config['has_zn']}\n")

    print("\n" + "="*60)
    print("SETUP COMPLETE")
    print("="*60)
    print(f"\nFiles created in: {BASE_DIR}")
    print("\nTo run a system:")
    print("  cd systemA_thrrs_THR")
    print("  bash run_gromacs.sh")
    print("\nNote: You may need to manually add ligand parameters")
    print("      using ACPYPE or similar for non-standard residues.")


if __name__ == "__main__":
    main()
