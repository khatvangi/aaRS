#!/bin/bash

# =================================================================
# CONFIGURATION
# =================================================================
# Path to your FIXED GPU 1 wrapper
WRAPPER="/storage/kiran-stuff/alphafold3/af3_gpu1.sh"

# Directory where outputs will be saved (relative to current folder)
OUTPUT_DIR="af3_output"

echo "================================================================="
echo "Starting Local Domain Matrix Run (60 jobs)"
echo "Target GPU: 1 (via af3_gpu1.sh)"
echo "Input Directory: $(pwd)"
echo "================================================================="

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Loop through all JSON files in the current directory
for JSON_FILE in *.json; do
    # Get the run name without extension (e.g., anc_prors_cat_ALA)
    RUN_NAME="${JSON_FILE%.json}"
    
    # ---------------------------------------------------------
    # 1. Check if already complete
    # ---------------------------------------------------------
    # We look for the summary_confidences.json which indicates a successful finish
    SUMMARY_FILE="$OUTPUT_DIR/$RUN_NAME/${RUN_NAME}_summary_confidences.json"
    
    if [ -f "$SUMMARY_FILE" ]; then
        echo "✅ Skipping $RUN_NAME (Already complete)"
        continue
    fi

    # ---------------------------------------------------------
    # 2. Run AlphaFold 3
    # ---------------------------------------------------------
    echo "------------------------------------------------"
    echo "🚀 Processing: $RUN_NAME"
    echo "   Time: $(date)"
    
    # Execute the wrapper
    # We pass the JSON path. The wrapper handles mounting and output location.
    bash "$WRAPPER" --json_path="$JSON_FILE"
    
    # Check exit status
    if [ $? -eq 0 ]; then
        echo "✅ Finished: $RUN_NAME"
    else
        echo "❌ FAILED: $RUN_NAME"
    fi
    
    # Optional: Brief pause to let GPU memory clear
    sleep 2
done

echo "================================================================="
echo "All jobs in this batch have been processed."
echo "================================================================="
