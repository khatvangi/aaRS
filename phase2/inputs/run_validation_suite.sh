#!/bin/bash
set -e

AF3_CMD="/storage/kiran-stuff/alphafold3/af3_wrapper.sh"

echo "================================================================================"
echo "VALIDATION SUITE: Modern Controls + Editing Domain + Reverse Tests"
echo "================================================================================"

cd inputs/af3_jsons_validation

for json in *.json; do
    name=$(basename "$json" .json)
    echo ""
    echo "Processing: $name ($(date))"
    
    $AF3_CMD --json_path="$json" 2>&1 | tee ../../logs/${name}.log
    
    if [ -d "af3_output" ]; then
        mkdir -p ../../outputs/$name
        sudo cp -r af3_output/* ../../outputs/$name/
        sudo chown -R kiran:kiran ../../outputs/$name
        sudo rm -rf af3_output
        echo "✓ $name complete"
    fi
done

cd ../..
echo "================================================================================"
echo "VALIDATION SUITE COMPLETE"
echo "================================================================================"
