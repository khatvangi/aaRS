# Generate pocket close-ups for AF3 complexes
# Usage:
# pymol -cq pymol_scripts/af3_pocket_views.pml -- \
#   phase2/af3_output_full/fulllength_deep_pro/fulllength_deep_pro_model.cif LUCA-Pro \
#   figures/af3_LUCA_pro_pocket.png

import sys
from pymol import cmd

pdb_path = sys.argv[1]
label = sys.argv[2]
out_png = sys.argv[3]

ligand_sel = "(resn PRI+PYR) and hetatm"

cmd.reinitialize()
cmd.load(pdb_path, "prot")

cmd.select("lig", ligand_sel)
cmd.select("prot_chains", "all and not lig and not resn HOH+WAT")

cmd.hide("everything")
cmd.show("cartoon", "prot_chains")
cmd.color("gray80", "prot_chains")

cmd.show("sticks", "lig")
cmd.color("marine", "lig")

cmd.select("pocket", "byres (prot_chains within 5 of lig)")
cmd.show("sticks", "pocket")
cmd.color("tv_orange", "pocket")

cmd.dist("hb_lig", "(pocket and (elem N+O))", "(lig and (elem N+O))", 3.5, 0)
cmd.color("yellow", "hb_lig")

cmd.bg_color("white")
cmd.set("cartoon_sampling", 10)
cmd.set("ray_opaque_background", 0)
cmd.zoom("lig", 10)

cmd.ray(1600, 1200)
cmd.png(out_png, dpi=300)
cmd.quit()
