#!/bin/bash
# 04_structural_alignment.sh
# High-quality alignment with MAFFT L-INS-i

set -e

echo "=== Building structural alignments with MAFFT ==="

# Function to align one aaRS
align_one_aars() {
    local AARS=$1
    local INPUT="data/interim/${AARS}_catalytic.fasta"
    local OUTPUT="data/interim/${AARS}_aligned.fasta"
    
    echo "  Aligning $AARS with L-INS-i (16 threads)..."
    
    mafft --localpair \
          --maxiterate 1000 \
          --thread 16 \
          --quiet \
          $INPUT > $OUTPUT
    
    # Calculate alignment statistics
    python3 - <<EOF
from Bio import AlignIO
aln = AlignIO.read("$OUTPUT", "fasta")
print(f"✓ ${AARS}: {len(aln)} sequences, {aln.get_alignment_length()} columns")
EOF
}

export -f align_one_aars

# Parallel alignment (4 aaRS × 16 threads each = 64 threads total)
parallel -j4 align_one_aars ::: ProRS ThrRS SerRS ValRS

echo "✓ Alignment complete"
