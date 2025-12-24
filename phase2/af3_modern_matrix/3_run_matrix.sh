#!/bin/bash

# Path to your existing wrapper
WRAPPER="/storage/kiran-stuff/alphafold3/af3_wrapper.sh"

# Save the current location
BASE_DIR=$(pwd)

echo "Starting Batch AF3 Run (Full MSA Mode)..."

# Loop through all folders starting with "modern_"
for RUN_DIR in modern_*; do
    if [ -d "$RUN_DIR" ]; then
        echo "------------------------------------------------"
        echo "Processing: $RUN_DIR"
        
        # 1. Enter the run directory
        cd "$RUN_DIR" || exit

        # 2. Identify the JSON file
        JSON_FILE="${RUN_DIR}.json"

        # 3. Run the Wrapper
        # Since we couldn't copy the MSAs, this will run Jackhmmer every time.
        # This is robust and guarantees correct results.
        bash "$WRAPPER" --json_path="$JSON_FILE"

        # 4. Return to base directory
        cd "$BASE_DIR" || exit
    fi
done

echo "------------------------------------------------"
echo "All jobs complete."
