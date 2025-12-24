#!/bin/bash

# =================================================================
# CONFIGURATION
# =================================================================
# Path to your standard wrapper (uses GPU 0 by default, or adjust if needed)
# If you want to use GPU 1, change this to: /storage/kiran-stuff/alphafold3/af3_gpu1.sh
WRAPPER="/storage/kiran-stuff/alphafold3/af3_wrapper.sh"

# Output directory naming
# AF3 wrapper usually creates an 'af3_output' folder in the current dir
OUTPUT_BASE="af3_output"

echo "================================================================="
echo "Starting Modern ThrRS + Zinc Matrix (20 jobs)"
echo "Wrapper: $WRAPPER"
echo "================================================================="

# Loop through all the Zinc-containing JSONs
for JSON_FILE in modern_thrrs_ecoli_zn_*.json; do
    # Strip extension to get run name
    RUN_NAME="${JSON_FILE%.json}"
    
    # ---------------------------------------------------------
    # 1. Check if already complete
    # ---------------------------------------------------------
    # Look for the summary file to avoid re-running completed jobs
    SUMMARY_FILE="$OUTPUT_BASE/$RUN_NAME/${RUN_NAME}_summary_confidences.json"
    
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
    
    # Run the full pipeline (MSA + Model)
    # We do NOT use --run_data_pipeline=false because we need fresh MSAs 
    # for each ligand-specific data.json
    bash "$WRAPPER" --json_path="$JSON_FILE"
    
    # Check exit status
    if [ $? -eq 0 ]; then
        echo "✅ Finished: $RUN_NAME"
    else
        echo "❌ FAILED: $RUN_NAME"
    fi
    
done

echo "================================================================="
echo "Batch Run Complete."
