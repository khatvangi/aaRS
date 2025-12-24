#!/usr/bin/env bash
set -euo pipefail

ROOT="/storage/kiran-stuff/aaRS/phase2"
SCRIPTS="$ROOT/scripts/alpharank_like"
RESULTS="$ROOT/results"

mkdir -p "$RESULTS"

echo "[1/4] Building manifest..."
python "$SCRIPTS/build_manifest.py"

echo "[2/4] Extracting AF3 metrics..."
python "$SCRIPTS/extract_af3_metrics.py"

echo "[3/4] Computing pocket IoU..."
python "$SCRIPTS/pocket_iou.py"

echo "[4/4] Summarizing competitions..."
python "$SCRIPTS/summarize_competitions.py"

echo "Done. Results in $RESULTS"
