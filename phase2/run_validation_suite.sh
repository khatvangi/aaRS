#!/bin/bash
set -e

AF3_CMD="/storage/kiran-stuff/alphafold3/af3_wrapper.sh"

echo "================================================================================"
echo "VALIDATION SUITE: Modern Controls + Editing + Reverse + Negative"
echo "================================================================================"

START_TIME=$(date +%s)

cd inputs/af3_jsons_validation

for json in *.json; do
    name=$(basename "$json" .json)
    echo ""
    echo "--------------------------------------------------------------------------------"
    echo "Processing: $name"
    echo "Started: $(date)"
    echo "--------------------------------------------------------------------------------"
    
    if $AF3_CMD --json_path="$json" 2>&1 | tee ../../logs/${name}.log; then
        echo "✓ AF3 completed for $name"
    else
        echo "✗ AF3 failed for $name"
        continue
    fi
    
    # Move output
    if [ -d "af3_output" ]; then
        mkdir -p ../../outputs/$name
        sudo cp -r af3_output/* ../../outputs/$name/ 2>/dev/null || true
        sudo chown -R kiran:kiran ../../outputs/$name 2>/dev/null || true
        sudo rm -rf af3_output
        echo "✓ Output saved to: outputs/$name"
    fi
    
    echo "Completed: $(date)"
done

cd ../..

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
HOURS=$((ELAPSED / 3600))
MINUTES=$(((ELAPSED % 3600) / 60))

echo ""
echo "================================================================================"
echo "VALIDATION SUITE COMPLETE"
echo "================================================================================"
echo "Total runtime: ${HOURS}h ${MINUTES}m"
echo "Outputs saved in: outputs/"
