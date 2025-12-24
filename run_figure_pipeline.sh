#!/bin/bash
# Master script to generate all manuscript figures
# Run from: /storage/kiran-stuff/aaRS

set -e  # Exit on error

echo "=========================================="
echo "MANUSCRIPT FIGURE GENERATION PIPELINE"
echo "=========================================="

# Check we're in the right directory
if [ ! -d "phase1" ] || [ ! -d "phase2" ]; then
    echo "ERROR: Please run this from /storage/kiran-stuff/aaRS"
    exit 1
fi

# Create output directory
mkdir -p manuscript_figures/structures

echo ""
echo "Step 1: Checking Python dependencies..."
python3 -c "import matplotlib, seaborn, pandas, numpy, Bio" 2>/dev/null
if [ $? -ne 0 ]; then
    echo "Installing required packages..."
    pip install --user biopython matplotlib seaborn pandas numpy ete3
fi
echo "✓ Dependencies OK"

echo ""
echo "Step 2: Extracting data from project..."
python3 extract_manuscript_data.py
if [ $? -ne 0 ]; then
    echo "ERROR: Data extraction failed"
    exit 1
fi
echo "✓ Data extraction complete"

echo ""
echo "Step 3: Generating figures..."
python3 generate_all_figures.py
if [ $? -ne 0 ]; then
    echo "ERROR: Figure generation failed"
    exit 1
fi
echo "✓ Figures generated"

echo ""
echo "Step 4: Generating PyMOL scripts..."
python3 generate_pymol_script.py
echo "✓ PyMOL scripts ready"

echo ""
echo "=========================================="
echo "✓ PIPELINE COMPLETE"
echo "=========================================="
echo ""
echo "Generated files in manuscript_figures/:"
ls -lh manuscript_figures/*.pdf manuscript_figures/*.png 2>/dev/null | awk '{print "  "$9" ("$5")"}'

echo ""
echo "=========================================="
echo "NEXT STEPS"
echo "=========================================="
echo ""
echo "1. Review PNG previews:"
echo "   - Figure1_phylogeny_domains.png"
echo "   - Figure2_af3_results.png"
echo "   - Figure3_domain_evolution.png"
echo ""
echo "2. Render structures in PyMOL:"
echo "   pymol -c manuscript_figures/render_structures.pml"
echo "   OR follow manual guide in manuscript_figures/PyMOL_manual_guide.md"
echo ""
echo "3. For publication (600 DPI):"
echo "   Edit generate_all_figures.py and change:"
echo "   plt.savefig(..., dpi=600)"
echo ""
echo "4. Upload figures to manuscript!"
echo ""
echo "=========================================="

# Summary statistics
echo ""
echo "Quick Stats:"
echo "  - PDF figures: $(ls manuscript_figures/*.pdf 2>/dev/null | wc -l)"
echo "  - PNG previews: $(ls manuscript_figures/*.png 2>/dev/null | wc -l)"
if [ -f "manuscript_figures/extracted_data.json" ]; then
    echo "  - Data file: manuscript_figures/extracted_data.json ($(wc -c < manuscript_figures/extracted_data.json) bytes)"
fi

echo ""
echo "All done! 🎉"
