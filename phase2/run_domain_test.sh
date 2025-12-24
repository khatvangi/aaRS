#!/bin/bash
set -e

AF3_CMD="/storage/kiran-stuff/alphafold3/af3_wrapper.sh"

echo "================================================================================"
echo "PHASE 2: CATALYTIC DOMAIN TESTING (aa 200-700)"
echo "Testing ancestral promiscuity with memory-efficient domain models"
echo "================================================================================"

START_TIME=$(date +%s)

for json_file in inputs/af3_jsons_domain/*.json; do
    name=$(basename "$json_file" .json)
    echo ""
    echo "--------------------------------------------------------------------------------"
    echo "Processing: $name"
    echo "Started: $(date)"
    echo "--------------------------------------------------------------------------------"
    
    cd inputs/af3_jsons_domain
    
    if $AF3_CMD --json_path=$(basename "$json_file") 2>&1 | tee ../../logs/${name}.log; then
        echo "✓ AF3 completed for $name"
    else
        echo "✗ AF3 failed for $name"
        cd ../..
        continue
    fi
    
    cd ../..
    
    if [ -d "inputs/af3_jsons_domain/af3_output" ]; then
        sudo mkdir -p outputs/$name && sudo chown -R kiran:kiran outputs/
        sudo mv inputs/af3_jsons_domain/af3_output/* outputs/$name/
        echo "✓ Output saved to: outputs/$name"
    fi
    
    echo "Completed: $(date)"
done

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
HOURS=$((ELAPSED / 3600))
MINUTES=$(((ELAPSED % 3600) / 60))

echo ""
echo "================================================================================"
echo "PHASE 2 DOMAIN TESTING COMPLETE"
echo "================================================================================"
echo "Total runtime: ${HOURS}h ${MINUTES}m"
echo "Outputs saved in: outputs/"
