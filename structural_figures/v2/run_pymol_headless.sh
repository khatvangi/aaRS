#!/bin/bash
# Run PyMOL in proper headless mode

export PYMOL_PATH=/home/kiran/miniforge3/bin/pymol

# Run PyMOL with -c (command line), -q (quiet), -k (keep going)
$PYMOL_PATH -cq /storage/kiran-stuff/aaRS/structural_figures/v2/pymol_analysis.pml > /storage/kiran-stuff/aaRS/structural_figures/v2/pymol_output.log 2>&1

echo "PyMOL analysis complete. Check pymol_output.log for details."
