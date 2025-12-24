#!/bin/bash
# Run MMseqs2-based AlphaFold 3 for modern_thrrs_ALA

source ~/miniforge3/etc/profile.d/conda.sh
conda activate af3_mmseqs2

export XLA_FLAGS="--xla_gpu_enable_triton_gemm=false"
export XLA_PYTHON_CLIENT_PREALLOCATE=true
export XLA_CLIENT_MEM_FRACTION=0.95
export MMSEQS_TMP_DIR=/storage/kiran-stuff/mmseqs2_tmp
export PATH=$PATH:/storage/kiran-stuff/alphafold/colabfold_batch/colabfold-conda/bin

echo "========================================="
echo "MMseqs2-based AlphaFold 3"
echo "========================================="
echo "Protein: modern_thrrs_ALA (ThrRS)"
echo "Sequence length: 642 aa"
echo "Ligand: ALA"
echo "Database searches: Sequential (4 databases)"
echo ""
echo "Comparison:"
echo "  Jackhmmer MSA: 20 minutes (1200 seconds)"
echo "  MMseqs2 MSA: ? (measuring now...)"
echo ""

START_TIME=$(date +%s)

~/miniforge3/envs/af3_mmseqs2/bin/python \
  /storage/kiran-stuff/AF3_mmseqs2/run_alphafold.py \
  --json_path=modern_thrrs_ALA.json \
  --model_dir=/storage/kiran-stuff/alphafold3_models \
  --db_dir=/storage/kiran-stuff/alphafold_databases \
  --output_dir=mmseqs2_output \
  --buckets='256,512,768,1024' \
  --flash_attention_implementation=xla \
  --msa_tool=mmseqs2 \
  --run_data_pipeline=true \
  --run_inference=true 2>&1 | tee mmseqs2_run.log

EXIT_CODE=${PIPESTATUS[0]}

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
MINUTES=$((ELAPSED / 60))
SECONDS=$((ELAPSED % 60))

echo ""
echo "========================================="
echo "Run completed!"
echo "Exit code: $EXIT_CODE"
echo "Total time: ${MINUTES}m ${SECONDS}s (${ELAPSED} seconds)"
echo "========================================="

# Extract MSA generation time from log
MSA_TIME=$(grep "Getting protein MSAs took" mmseqs2_run.log | awk '{print $6}')
if [ ! -z "$MSA_TIME" ]; then
    MSA_MINUTES=$(echo "scale=1; $MSA_TIME / 60" | bc)
    echo ""
    echo "Performance Comparison:"
    echo "  Jackhmmer MSA: 1200 seconds (20.0 minutes)"
    echo "  MMseqs2 MSA:   ${MSA_TIME} seconds (${MSA_MINUTES} minutes)"
    
    SPEEDUP=$(echo "scale=1; 1200 / $MSA_TIME" | bc)
    echo "  Speedup: ${SPEEDUP}x faster"
fi

exit $EXIT_CODE
