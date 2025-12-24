#!/bin/bash
set -e
echo "=== Aligning tRNAs with MAFFT (simple mode) ==="

align_one_trna() {
    local TRNA=$1
    local INPUT="data/raw/tRNA_${TRNA}_all.fasta"
    local OUTPUT_FASTA="data/interim/tRNA_${TRNA}_aligned.fasta"
    
    echo "  Aligning tRNA-$TRNA..."
    
    # Use simpler MAFFT mode (less memory)
    mafft --auto \
          --thread 4 \
          --quiet \
          $INPUT > $OUTPUT_FASTA 2>/dev/null || \
    mafft --retree 2 \
          --thread 4 \
          --quiet \
          $INPUT > $OUTPUT_FASTA
    
    # Stats
    COUNT=$(grep -c "^>" $OUTPUT_FASTA || echo 0)
    if [ $COUNT -gt 0 ]; then
        echo "✓ tRNA-${TRNA}: $COUNT sequences aligned"
    else
        echo "✗ tRNA-${TRNA}: alignment failed"
    fi
}

export -f align_one_trna

# Sequential to avoid OOM
align_one_trna Pro
align_one_trna Thr
align_one_trna Ser
align_one_trna Val

echo "✓ tRNA alignment complete"
