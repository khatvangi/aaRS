#!/bin/bash
# AlphaFold3 wrapper for boron, prefer GPU 1

IMAGE="alphafold3-local"
JSON_PATH="$1"

if [ -z "$JSON_PATH" ]; then
  echo "Usage: $0 /storage/.../job.json"
  exit 1
fi

echo "Running AF3 on boron GPU 1 for $(basename "$JSON_PATH")"
echo "Using Docker image: $IMAGE"

docker run --rm \
  --gpus all \
  -e XLA_FLAGS="--xla_disable_hlo_passes=custom-kernel-fusion-rewriter" \
  -v /storage:/storage \
  -v /storage/kiran-stuff/alphafold3_models:/models \
  -v /storage/kiran-stuff/alphafold_databases:/databases \
  -v "$(pwd)":/workspace \
  -w /workspace \
  "$IMAGE" \
  python /app/alphafold/run_alphafold.py \
    --gpu_device=1 \
    --model_dir=/models \
    --db_dir=/databases \
    --output_dir=/workspace/af3_output \
    --flash_attention_implementation=xla \
    --num_recycles=1 \
    --num_diffusion_samples=1 \
    --buckets=256,512,768,1024,1280,1536,2048,2560 \
    --json_path="$JSON_PATH"

