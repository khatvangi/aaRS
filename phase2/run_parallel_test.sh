#!/bin/bash
set -e

echo "================================================================================"
echo "PHASE 2: PARALLEL DUAL-GPU AF3 TESTING"
echo "Running 2 jobs simultaneously (GPU 0 + GPU 1)"
echo "================================================================================"

START_TIME=$(date +%s)

# Define the 4 test cases
TESTS=(
    "shallow_ancestral_pro"
    "shallow_ancestral_thr"
    "deep_ancestral_pro"
    "deep_ancestral_thr"
)

# Function to run a single AF3 job
run_job() {
    local name=$1
    local gpu=$2
    local wrapper=$3
    
    echo ""
    echo "--------------------------------------------------------------------------------"
    echo "Processing: $name (GPU $gpu)"
    echo "Started: $(date)"
    echo "--------------------------------------------------------------------------------"
    
    cd inputs/af3_jsons
    
    if $wrapper --json_path="${name}.json" 2>&1 | tee ../../logs/${name}_gpu${gpu}.log; then
        echo "✓ AF3 completed for $name"
    else
        echo "✗ AF3 failed for $name"
        cd ../..
        return 1
    fi
    
    cd ../..
    
    # Move output
    if [ -d "inputs/af3_jsons/af3_output" ]; then
        mkdir -p outputs/$name
        mv inputs/af3_jsons/af3_output/* outputs/$name/ 2>/dev/null || true
        echo "✓ Output saved to: outputs/$name"
    fi
    
    echo "Completed: $(date)"
}

# Run jobs in pairs (GPU 0 + GPU 1)
for i in 0 2; do
    if [ $i -lt ${#TESTS[@]} ]; then
        # Launch first job on GPU 0
        run_job "${TESTS[$i]}" 0 "/storage/kiran-stuff/alphafold3/af3_wrapper_gpu0.sh" &
        PID0=$!
        
        # Launch second job on GPU 1 (if exists)
        if [ $((i+1)) -lt ${#TESTS[@]} ]; then
            run_job "${TESTS[$((i+1))]}" 1 "/storage/kiran-stuff/alphafold3/af3_wrapper_gpu1.sh" &
            PID1=$!
            
            # Wait for both to complete
            wait $PID0
            wait $PID1
        else
            wait $PID0
        fi
    fi
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
