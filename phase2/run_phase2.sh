#!/bin/bash
set -e

echo "=== Phase 2: AlphaFold3 Modeling Pipeline ==="

# Step 1: Prepare inputs
echo -e "\n[Step 1/3] Preparing AF3 input JSONs..."
python3 scripts/01_prepare_af3_inputs.py || { echo "✗ Step 1 failed"; exit 1; }

# Step 2: Run AF3 (requires manual setup)
echo -e "\n[Step 2/3] Running AlphaFold3..."
echo "⚠ NOTE: You must configure scripts/02_run_af3.sh with your AF3 installation"
echo "⚠ Press Enter to continue (will run placeholder), or Ctrl+C to exit and configure AF3"
read
bash scripts/02_run_af3.sh || { echo "✗ Step 2 failed"; exit 1; }

# Step 3: Parse results
echo -e "\n[Step 3/3] Parsing AF3 outputs..."
python3 scripts/03_parse_af3_results.py || { echo "✗ Step 3 failed"; exit 1; }

echo -e "\n✓✓✓ Phase 2 pipeline complete ✓✓✓"
