#!/bin/bash

# Path to your matrix folder
BASE_DIR="af3_modern_matrix"

# Function to copy MSAs
copy_msas() {
    local ENZYME_PREFIX=$1  # e.g., "modern_prors" or "modern_thrrs"
    local TEMPLATE_NAME="${ENZYME_PREFIX}_ALA"
    
    # 1. Locate the Source MSA Folder
    # It is usually: [RunFolder]/af3_output/[JobName]/msas
    SOURCE_PATH="$BASE_DIR/$TEMPLATE_NAME/af3_output/$TEMPLATE_NAME/msas"

    echo "------------------------------------------------"
    echo "Processing $ENZYME_PREFIX..."
    echo "Looking for source MSAs at: $SOURCE_PATH"

    if [ ! -d "$SOURCE_PATH" ]; then
        echo "ERROR: Source MSA folder not found!"
        echo "Please verify that the '$TEMPLATE_NAME' run has finished successfully."
        return
    fi

    # 2. Distribute to all other variants
    for DIR in "$BASE_DIR/${ENZYME_PREFIX}_"*; do
        RUN_NAME=$(basename "$DIR")

        # Skip the template itself
        if [ "$RUN_NAME" == "$TEMPLATE_NAME" ]; then
            continue
        fi

        # 3. Create the Destination Structure
        # We must mimic AF3's expected structure:
        # [RunFolder]/af3_output/[RunName]/msas
        DEST_DIR="$DIR/af3_output/$RUN_NAME/msas"
        
        mkdir -p "$DEST_DIR"

        # 4. Copy the MSA files
        cp "$SOURCE_PATH/"* "$DEST_DIR/"
        
        echo "  -> Pre-seeded MSAs for: $RUN_NAME"
    done
}

# Execute for both enzymes
copy_msas "modern_prors"
copy_msas "modern_thrrs"

echo "------------------------------------------------"
echo "Done. All runs are primed with MSAs."
