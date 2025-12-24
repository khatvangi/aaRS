#!/bin/bash
set -e  # Exit immediately if any command fails

# Configuration
WRAPPER="/storage/kiran-stuff/alphafold3/af3_wrapper.sh"

# Define the specific list of 6 validation jobs
JOBS=(
    "modern_ecoli_cat_trp.json"
    "modern_ecoli_cat_phe.json"
    "modern_ecoli_full_pro.json"
    "modern_ecoli_full_thr.json"
    "modern_human_full_thr.json"
    "modern_human_full_pro.json"
)

echo "Starting batch processing of ${#JOBS[@]} validation jobs..."
echo "Wrapper: $WRAPPER"
echo "------------------------------------------------"

for json_file in "${JOBS[@]}"; do
    if [ -f "$json_file" ]; then
        echo ">>> [$(date '+%H:%M:%S')] Starting AF3 run for: $json_file"
        
        # Run the wrapper
        $WRAPPER --json_path="$json_file"
        
        echo ">>> [$(date '+%H:%M:%S')] Completed: $json_file"
        echo "------------------------------------------------"
    else
        echo "ERROR: Input file '$json_file' not found in current directory!"
        echo "Please generate it first using the python script."
        exit 1
    fi
done

echo "All validation jobs completed successfully."
