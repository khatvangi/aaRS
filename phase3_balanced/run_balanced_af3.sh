#!/bin/bash
#
# run_balanced_af3.sh
#
# Runs AF3 jobs from manifests/jobs.csv
# - Max 2 concurrent jobs (one per GPU)
# - Alternates GPU by seed%2
# - Logs to outputs/af3/logs/
# - Skips jobs with existing output dirs
#
# PLACEHOLDER: Replace AF3_COMMAND with your actual AF3 runner command

set -e

BASE_DIR="$(cd "$(dirname "$0")" && pwd)"
MANIFEST="$BASE_DIR/manifests/jobs.csv"
OUTPUT_DIR="$BASE_DIR/outputs/af3"
LOG_DIR="$OUTPUT_DIR/logs"

# PLACEHOLDER - Replace with your actual AF3 command
# Example: AF3_COMMAND="python /path/to/run_alphafold.py"
AF3_COMMAND="PLACEHOLDER_AF3_COMMAND"

if [[ "$AF3_COMMAND" == "PLACEHOLDER"* ]]; then
    echo "ERROR: AF3_COMMAND is still a placeholder"
    echo "Please edit this script and set the correct AF3 command"
    exit 1
fi

mkdir -p "$LOG_DIR"

# Track running jobs (PIDs)
declare -A GPU_PIDS
GPU_PIDS[0]=""
GPU_PIDS[1]=""

wait_for_gpu() {
    local gpu=$1
    if [[ -n "${GPU_PIDS[$gpu]}" ]]; then
        wait "${GPU_PIDS[$gpu]}" 2>/dev/null || true
        GPU_PIDS[$gpu]=""
    fi
}

run_job() {
    local job_name=$1
    local input_json=$2
    local seed=$3

    local gpu=$((seed % 2))
    local job_output="$OUTPUT_DIR/$job_name"
    local log_file="$LOG_DIR/${job_name}.log"

    # Skip if output exists
    if [[ -d "$job_output" ]]; then
        echo "SKIP: $job_name (output exists)"
        return
    fi

    # Wait for GPU to be free
    wait_for_gpu $gpu

    echo "RUN: $job_name on GPU $gpu"

    (
        export CUDA_VISIBLE_DEVICES=$gpu
        $AF3_COMMAND "$input_json" "$job_output" > "$log_file" 2>&1
    ) &

    GPU_PIDS[$gpu]=$!
}

# Read manifest and run jobs
tail -n +2 "$MANIFEST" | while IFS=, read -r condition ligand seed job_name input_json include_zn; do
    run_job "$job_name" "$input_json" "$seed"
done

# Wait for remaining jobs
wait_for_gpu 0
wait_for_gpu 1

echo "All jobs completed"
