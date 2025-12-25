#!/usr/bin/env python3
"""
pymol_ligand_overlay.py

PyMOL script to overlay all 20 amino acid ligands on modern_prors.
Aligns proteins, shows ligands superimposed to visualize binding pocket.

Run with: pymol -cq pymol_ligand_overlay.py
Or interactively: run pymol_ligand_overlay.py
"""

from pymol import cmd, util
import os

# Base path
BASE = "/storage/kiran-stuff/aaRS/phase2/af3_modern_matrix"

# All 20 amino acids
AMINO_ACIDS = [
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"
]

# Color scheme - cognate (PRO) in green, others by property
COLORS = {
    # Cognate - bright green
    "PRO": "green",
    # Hydrophobic - blues
    "ALA": "lightblue", "VAL": "skyblue", "ILE": "slate", "LEU": "deepblue",
    "MET": "marine", "PHE": "density", "TRP": "purpleblue", "TYR": "violet",
    # Polar - yellows/oranges
    "SER": "yellow", "THR": "orange", "ASN": "lightorange", "GLN": "brightorange",
    "CYS": "sulfur",
    # Charged - reds/magentas
    "ASP": "red", "GLU": "firebrick",
    "LYS": "magenta", "ARG": "hotpink", "HIS": "pink",
    # Small
    "GLY": "white",
}

def load_structures():
    """Load all 20 structures"""
    for aa in AMINO_ACIDS:
        cif_path = f"{BASE}/modern_prors_{aa}/af3_output/modern_prors_{aa}/modern_prors_{aa}_model.cif"
        if os.path.exists(cif_path):
            cmd.load(cif_path, f"prors_{aa}")
            print(f"Loaded: {aa}")
        else:
            print(f"NOT FOUND: {aa}")

def align_structures():
    """Align all structures to PRO (cognate) reference"""
    ref = "prors_PRO"

    for aa in AMINO_ACIDS:
        if aa == "PRO":
            continue
        obj = f"prors_{aa}"
        if obj in cmd.get_names():
            # Align by protein chain A
            cmd.align(f"{obj} and chain A", f"{ref} and chain A")
            print(f"Aligned: {aa}")

def setup_visualization():
    """Set up the visualization"""

    # Hide everything first
    cmd.hide("everything")

    # Show protein from reference (PRO) as cartoon
    cmd.show("cartoon", "prors_PRO and chain A")
    cmd.color("gray80", "prors_PRO and chain A")

    # Show all ligands as sticks
    for aa in AMINO_ACIDS:
        obj = f"prors_{aa}"
        if obj in cmd.get_names():
            # Select ligand (chain B typically)
            cmd.show("sticks", f"{obj} and chain B")
            cmd.color(COLORS.get(aa, "white"), f"{obj} and chain B")

            # Make cognate PRO thicker and prominent
            if aa == "PRO":
                cmd.set("stick_radius", 0.4, f"{obj} and chain B")
            else:
                cmd.set("stick_radius", 0.2, f"{obj} and chain B")

    # Create selection for all ligands
    cmd.select("all_ligands", " or ".join([f"(prors_{aa} and chain B)" for aa in AMINO_ACIDS]))

    # Zoom to ligand binding site
    cmd.zoom("all_ligands", buffer=10)

    # Nice rendering settings
    cmd.set("ray_shadows", 0)
    cmd.set("ambient", 0.4)
    cmd.set("specular", 0.2)
    cmd.bg_color("white")

    # Label the cognate
    cmd.label("prors_PRO and chain B and name CA", '"PRO"')

    print("\nVisualization ready!")
    print("- Protein: gray cartoon (from PRO structure)")
    print("- Ligands: colored sticks (PRO = green, thick)")
    print("- Zoom: ligand binding site")

def save_images():
    """Save PNG images"""
    outdir = "/storage/kiran-stuff/aaRS/phase3_balanced"

    # Front view
    cmd.ray(2400, 2400)
    cmd.png(f"{outdir}/FIG_ligand_overlay_front.png", dpi=300)
    print(f"Saved: {outdir}/FIG_ligand_overlay_front.png")

    # Rotate 90 degrees
    cmd.turn("y", 90)
    cmd.ray(2400, 2400)
    cmd.png(f"{outdir}/FIG_ligand_overlay_side.png", dpi=300)
    print(f"Saved: {outdir}/FIG_ligand_overlay_side.png")

    # Top view
    cmd.turn("y", -90)
    cmd.turn("x", 90)
    cmd.ray(2400, 2400)
    cmd.png(f"{outdir}/FIG_ligand_overlay_top.png", dpi=300)
    print(f"Saved: {outdir}/FIG_ligand_overlay_top.png")

def main():
    """Main function"""
    print("="*50)
    print("Loading modern_prors structures (20 ligands)")
    print("="*50)

    cmd.reinitialize()

    load_structures()
    align_structures()
    setup_visualization()
    save_images()

    # Save session
    session_path = "/storage/kiran-stuff/aaRS/phase3_balanced/ligand_overlay.pse"
    cmd.save(session_path)
    print(f"\nSession saved: {session_path}")
    print("\nDone! Open ligand_overlay.pse in PyMOL to view interactively.")

# Run if called directly
if __name__ == "__main__" or __name__ == "pymol":
    main()
