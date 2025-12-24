#!/bin/bash
# Run fpocket for accurate binding pocket volume calculations
# This will provide publication-quality volume measurements

echo "=========================================="
echo "Running fpocket for pocket volume analysis"
echo "=========================================="

# Check if fpocket is installed
if ! command -v fpocket &> /dev/null; then
    echo "fpocket not found. Installing..."
    conda install -c conda-forge fpocket -y
fi

# Create output directory
mkdir -p /storage/kiran-stuff/aaRS/structural_figures/v2/fpocket_results
cd /storage/kiran-stuff/aaRS/structural_figures/v2/fpocket_results

echo ""
echo "Running fpocket on LUCA ProRS pocket..."
fpocket -f ../luca_pocket.pdb

echo ""
echo "Running fpocket on Modern ProRS pocket..."
fpocket -f ../modern_pocket.pdb

echo ""
echo "=========================================="
echo "fpocket analysis complete!"
echo "=========================================="
echo ""
echo "Results saved to:"
echo "  - luca_pocket_out/"
echo "  - modern_pocket_out/"
echo ""
echo "Key files to check:"
echo "  - *_info.txt: Summary with pocket volumes"
echo "  - *_pockets.pqr: Pocket coordinates and volumes"
echo ""
echo "Extract volumes from *_info.txt files:"
grep -A 20 "Pocket 1" luca_pocket_out/*_info.txt 2>/dev/null || echo "LUCA results pending"
grep -A 20 "Pocket 1" modern_pocket_out/*_info.txt 2>/dev/null || echo "Modern results pending"
echo ""
echo "Use these volumes to update Figure 5 and manuscript text!"
