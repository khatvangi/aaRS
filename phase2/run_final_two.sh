#!/bin/bash
# Run the final 2 jobs without set -e

AF3_CMD="/storage/kiran-stuff/alphafold3/af3_wrapper.sh"

echo "================================================================================"
echo "RUNNING FINAL 2 JOBS"
echo "================================================================================"

for name in "shallow_domain_pro" "shallow_domain_thr"; do
    echo ""
    echo "--------------------------------------------------------------------------------"
    echo "Processing: $name ($(date))"
    echo "--------------------------------------------------------------------------------"
    
    cd inputs/af3_jsons_domain
    $AF3_CMD --json_path="${name}.json" 2>&1 | tee ../../logs/${name}.log
    cd ../..
    
    # Move output
    if [ -d "inputs/af3_jsons_domain/af3_output" ]; then
        mkdir -p outputs/$name
        sudo cp -r inputs/af3_jsons_domain/af3_output/* outputs/$name/ || true
        sudo chown -R kiran:kiran outputs/$name || true
        sudo rm -rf inputs/af3_jsons_domain/af3_output || true
        echo "✓ Output saved"
    fi
    
    echo "Completed: $(date)"
done

echo "ALL JOBS COMPLETE!"
