#!/bin/bash
set -e

AF3_CMD="/storage/kiran-stuff/alphafold3/af3_wrapper.sh"

echo "================================================================================"
echo "TESTING: Single job (shallow_ancestral_pro) to verify no OOM"
echo "================================================================================"

json_file="inputs/af3_jsons/shallow_ancestral_pro.json"
name="shallow_ancestral_pro"

echo "Processing: $name"
echo "Started: $(date)"

cd inputs/af3_jsons

if $AF3_CMD --json_path=$(basename "$json_file") 2>&1 | tee ../../logs/test_single.log; then
    echo "✓ AF3 completed successfully!"
else
    echo "✗ AF3 failed"
    exit 1
fi

cd ../..

if [ -d "inputs/af3_jsons/af3_output" ]; then
    mkdir -p outputs/$name
    mv inputs/af3_jsons/af3_output/* outputs/$name/
    echo "✓ Output saved"
fi

echo "Test complete: $(date)"
