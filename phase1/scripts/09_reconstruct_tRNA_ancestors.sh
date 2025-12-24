#!/bin/bash
set -e
echo "=== tRNA ancestral reconstruction ==="

reconstruct_one_trna() {
    local TRNA=$1
    local INPUT="data/interim/tRNA_${TRNA}_aligned.fasta"
    local PREFIX="results/tRNA_${TRNA}"
    
    echo "  Reconstructing tRNA-$TRNA..."
    
    if [ ! -f "$INPUT" ] || [ ! -s "$INPUT" ]; then
        echo "  ERROR: Alignment file missing or empty: $INPUT"
        return 1
    fi
    
    # Combined tree + ASR
    iqtree3 -s $INPUT \
            -st DNA \
            -m MFP \
            -bb 1000 \
            -asr \
            -nt 16 \
            --prefix ${PREFIX} \
            --quiet || echo "  WARNING: IQ-TREE failed for $TRNA"
    
    echo "✓ tRNA-${TRNA} complete"
}

export -f reconstruct_one_trna

reconstruct_one_trna Pro &
reconstruct_one_trna Thr &
reconstruct_one_trna Ser &
reconstruct_one_trna Val &

wait

echo "✓ tRNA ancestral reconstruction complete"
