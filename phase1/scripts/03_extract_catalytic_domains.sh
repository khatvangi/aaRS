#!/bin/bash
set -e
echo "=== Preparing aaRS sequences (using full-length) ==="

mkdir -p data/interim

for AARS in ProRS ThrRS SerRS ValRS; do
    INPUT="data/raw/${AARS}_filtered.fasta"
    OUTPUT="data/interim/${AARS}_catalytic.fasta"
    
    echo "  Processing $AARS..."
    
    # Just copy full-length sequences - skip domain extraction
    cp "$INPUT" "$OUTPUT"
    
    COUNT=$(grep -c "^>" "$OUTPUT")
    echo "  ✓ $AARS: $COUNT sequences"
done

echo "✓ Sequence preparation complete"
