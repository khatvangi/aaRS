#!/bin/bash
# Script to view all generated figures
# Usage: bash view_all_figures.sh

echo "=== ProRS Manuscript Figures ==="
echo ""
echo "Location: /storage/kiran-stuff/aaRS/"
echo ""

# Check if figures exist
figures=(
    "figure1_luca_promiscuity_clean.png"
    "figure2_editing_inverted_clean.png"
    "figure3_luca_vs_modern_clean.png"
    "figure5_iptm_bars.png"
    "figure5_iptm_comprehensive.png"
)

echo "Available Figures:"
echo "-------------------"
for fig in "${figures[@]}"; do
    if [ -f "$fig" ]; then
        size=$(du -h "$fig" | cut -f1)
        dims=$(identify -format "%wx%h" "$fig" 2>/dev/null || echo "N/A")
        echo "✓ $fig ($size, $dims)"
    else
        echo "✗ $fig (NOT FOUND)"
    fi
done

echo ""
echo "PyMOL Session Files (for interactive viewing/editing):"
echo "-------------------------------------------------------"
sessions=(
    "figure1_session.pse"
    "figure2_session.pse"
    "figure3_session.pse"
)

for sess in "${sessions[@]}"; do
    if [ -f "$sess" ]; then
        size=$(du -h "$sess" | cut -f1)
        echo "✓ $sess ($size)"
    else
        echo "✗ $sess (NOT FOUND)"
    fi
done

echo ""
echo "To view a figure:"
echo "  - PNG: eog figure1_luca_promiscuity_clean.png"
echo "  - PDF: evince figure5_iptm_bars.pdf"
echo "  - PyMOL session: /home/kiran/miniforge3/bin/pymol figure1_session.pse"
echo ""

echo "To regenerate all figures:"
echo "  bash regenerate_all_figures.sh"
echo ""

# Optional: Open all PNG files in image viewer
read -p "Open all PNG figures in image viewer? (y/n): " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    echo "Opening figures..."
    for fig in "${figures[@]}"; do
        if [ -f "$fig" ]; then
            eog "$fig" &
        fi
    done
fi
