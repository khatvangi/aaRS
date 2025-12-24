#!/bin/bash
# Script to regenerate all ProRS manuscript figures
# Usage: bash regenerate_all_figures.sh

set -e  # Exit on error

echo "========================================="
echo "Regenerating ProRS Manuscript Figures"
echo "========================================="
echo ""

PYMOL="/home/kiran/miniforge3/bin/pymol"
PYTHON="python3"

# Check if PyMOL is available
if ! command -v $PYMOL &> /dev/null; then
    echo "Error: PyMOL not found at $PYMOL"
    echo "Please check PyMOL installation"
    exit 1
fi

echo "Using PyMOL: $PYMOL"
echo "Using Python: $PYTHON"
echo ""

# Figure 1: LUCA Promiscuity
echo "[1/4] Generating Figure 1: LUCA ProRS Promiscuity Overlay..."
$PYMOL -c generate_figure1_luca_promiscuity.pml
if [ $? -eq 0 ]; then
    echo "      ✓ Figure 1 complete"
else
    echo "      ✗ Figure 1 failed"
    exit 1
fi
echo ""

# Figure 2: Editing Domain
echo "[2/4] Generating Figure 2: Editing Domain Inverted Specificity..."
$PYMOL -c generate_figure2_editing_inverted.pml
if [ $? -eq 0 ]; then
    echo "      ✓ Figure 2 complete"
else
    echo "      ✗ Figure 2 failed"
    exit 1
fi
echo ""

# Figure 3: LUCA vs Modern
echo "[3/4] Generating Figure 3: LUCA vs Modern Comparison..."
$PYMOL -c generate_figure3_luca_vs_modern.pml
if [ $? -eq 0 ]; then
    echo "      ✓ Figure 3 complete"
else
    echo "      ✗ Figure 3 failed"
    exit 1
fi
echo ""

# Figure 5: ipTM Bar Charts
echo "[4/4] Generating Figure 5: ipTM Bar Charts..."
$PYTHON generate_figure5_iptm_bars.py
if [ $? -eq 0 ]; then
    echo "      ✓ Figure 5 complete"
else
    echo "      ✗ Figure 5 failed"
    exit 1
fi
echo ""

echo "========================================="
echo "All figures generated successfully!"
echo "========================================="
echo ""
echo "Generated files:"
ls -lh figure*.png figure*.pdf 2>/dev/null | awk '{print "  ", $9, "("$5")"}'
echo ""
echo "View summary: cat FIGURE_GENERATION_SUMMARY.txt"
echo "View documentation: less ADDITIONAL_FIGURES_README.md"
echo ""
