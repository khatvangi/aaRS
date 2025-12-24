#!/bin/bash

# AlphaFold3 Docker Wrapper for compute capability 7.x GPUs
IMAGE="alphafold3-local"
if ! docker image inspect $IMAGE >/dev/null 2>&1; then
    IMAGE="cford38/alphafold3"
fi

echo "Using Docker image: $IMAGE"
echo "Configured for compute capability 7.x GPUs (using XLA attention)"

# Set XLA flags for compute capability 7.x GPUs
export XLA_FLAGS="--xla_disable_hlo_passes=custom-kernel-fusion-rewriter"

docker run --rm \
  --gpus all \
  -e XLA_FLAGS="$XLA_FLAGS" \
  -v /storage/kiran-stuff/alphafold3_models:/models \
  -v /storage/kiran-stuff/alphafold_databases:/databases \
  -v "$(pwd)":/workspace \
  -w /workspace \
  $IMAGE \
  python /app/alphafold/run_alphafold.py \
  --model_dir=/models \
  --db_dir=/databases \
  --output_dir=/workspace/af3_output \
  --flash_attention_implementation=xla \
  "$@"
