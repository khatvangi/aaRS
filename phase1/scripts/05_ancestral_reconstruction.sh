#!/bin/bash
set -e
echo "=== Ancestral reconstruction with IQ-TREE 3 ==="

reconstruct_one_aars() {
    local AARS=$1
    local INPUT="data/interim/${AARS}_aligned.fasta"
    local PREFIX="results/${AARS}"
    
    echo "  Reconstructing $AARS..."
    
    # Combined: tree + ASR in one run
    iqtree3 -s $INPUT \
            -st AA \
            -m MFP \
            -bb 1000 \
            -alrt 1000 \
            -asr \
            -nt 16 \
            --prefix ${PREFIX} \
            --quiet
    
    echo "✓ ${AARS} reconstruction complete"
}

export -f reconstruct_one_aars

reconstruct_one_aars ProRS &
reconstruct_one_aars ThrRS &
reconstruct_one_aars SerRS &
reconstruct_one_aars ValRS &

wait

echo "✓ All aaRS reconstructions complete"
