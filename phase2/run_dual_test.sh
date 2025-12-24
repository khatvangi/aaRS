#!/bin/bash
set -e

# Use full path to af3
AF3_CMD="/storage/kiran-stuff/alphafold3/af3_wrapper.sh"

echo "================================================================================"
echo "PHASE 2: DUAL ANCESTOR AF3 TESTING"
echo "Testing: Shallow (eukaryotic) vs Deep (LUCA) ancestral aaRS-tRNA complexes"
echo "================================================================================"

START_TIME=$(date +%s)

for json_file in inputs/af3_jsons/*.json; do
    name=$(basename "$json_file" .json)
    echo ""
    echo "--------------------------------------------------------------------------------"
    echo "Processing: $name"
    echo "Started: $(date)"
    echo "--------------------------------------------------------------------------------"
    
    # Run AF3 from the json directory
    cd inputs/af3_jsons
    
    # Run with full path
    if $AF3_CMD --json_path=$(basename "$json_file") 2>&1; then
        echo "✓ AF3 completed for $name"
    else
        echo "✗ AF3 failed for $name"
        cd ../..
        continue
    fi
    
    cd ../..
    
    # Move output
    if [ -d "inputs/af3_jsons/af3_output" ]; then
        mkdir -p outputs/$name
        mv inputs/af3_jsons/af3_output/* outputs/$name/ 2>/dev/null || true
        echo "✓ Output saved to: outputs/$name"
        
        # Quick quality check
        if [ -f "outputs/$name/*/summary_confidences_0.json" ]; then
            echo "  Confidence scores:"
            python3 << PYEOF
import json
import glob
files = glob.glob("outputs/$name/*/summary_confidences_0.json")
if files:
    with open(files[0]) as f:
        data = json.load(f)
        print(f"    ipTM:  {data.get('iptm', 'N/A')}")
        print(f"    pLDDT: {data.get('ptm', 'N/A')}")
PYEOF
        fi
    fi
    
    echo "Completed: $(date)"
done

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
HOURS=$((ELAPSED / 3600))
MINUTES=$(((ELAPSED % 3600) / 60))

echo ""
echo "================================================================================"
echo "PHASE 2 COMPLETE"
echo "================================================================================"
echo "Total runtime: ${HOURS}h ${MINUTES}m"
echo "Outputs saved in: outputs/"

