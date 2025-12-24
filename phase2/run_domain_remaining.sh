#!/bin/bash
set -e

AF3_CMD="/storage/kiran-stuff/alphafold3/af3_wrapper.sh"

echo "================================================================================"
echo "RUNNING REMAINING 3 DOMAIN TESTS"
echo "================================================================================"

# Jobs to run (skip deep_domain_pro which is done)
JOBS=(
    "deep_domain_thr"
    "shallow_domain_pro"
    "shallow_domain_thr"
)

for name in "${JOBS[@]}"; do
    echo ""
    echo "--------------------------------------------------------------------------------"
    echo "Processing: $name"
    echo "Started: $(date)"
    echo "--------------------------------------------------------------------------------"
    
    cd inputs/af3_jsons_domain
    
    if $AF3_CMD --json_path="${name}.json" 2>&1 | tee ../../logs/${name}.log; then
        echo "✓ AF3 completed for $name"
    else
        echo "✗ AF3 failed for $name"
        cd ../..
        continue
    fi
    
    cd ../..
    
    # Move output with proper permissions
    if [ -d "inputs/af3_jsons_domain/af3_output" ]; then
        mkdir -p outputs/$name
        sudo cp -r inputs/af3_jsons_domain/af3_output/* outputs/$name/
        sudo chown -R kiran:kiran outputs/$name
        sudo rm -rf inputs/af3_jsons_domain/af3_output
        echo "✓ Output saved to: outputs/$name"
    fi
    
    echo "Completed: $(date)"
done

echo ""
echo "================================================================================"
echo "ALL REMAINING JOBS COMPLETE"
echo "================================================================================"
