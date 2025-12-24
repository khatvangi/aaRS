#!/bin/bash
AF3_CMD="/storage/kiran-stuff/alphafold3/af3_wrapper_fixperms.sh"

echo "================================================================================"
echo "VALIDATION SUITE (No sudo - fixed permissions)"
echo "================================================================================"

cd inputs/af3_jsons_validation

for json in *.json; do
    name=$(basename "$json" .json)
    
    # Skip already completed
    if [ -d "../../outputs/$name" ]; then
        echo "⏭️  Skipping $name (already done)"
        continue
    fi
    
    echo ""
    echo "Processing: $name ($(date))"
    
    $AF3_CMD --json_path="$json" 2>&1 | tee ../../logs/${name}.log
    
    # Move output (no sudo)
    if [ -d "af3_output" ]; then
        mkdir -p ../../outputs/$name
        cp -r af3_output/* ../../outputs/$name/
        rm -rf af3_output
        echo "✓ $name complete"
    fi
done

cd ../..
echo "================================================================================"
echo "ALL VALIDATION TESTS COMPLETE"
echo "================================================================================"
