#!/bin/bash
set -e

echo "=== Setting up Phase 2 directory ==="

mkdir -p phase2/{inputs,outputs,scripts,logs,results}
cd phase2

# Copy ancestral sequences from Phase 1
cp ../phase1/results/Anc-ProThrRS.fasta inputs/
cp ../phase1/results/Anc-tRNA-ProThr.fasta inputs/

echo "✓ Phase 2 structure created"
tree -L 2
