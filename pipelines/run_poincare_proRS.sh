#!/usr/bin/env bash

# Pipeline: ProRS MSA -> PoincaréMSA 2D embedding
# - Uses existing PoincareMSA clone at /storage/kiran-stuff/PoincareMSA
# - Input MSA: aaRS/phase1/data/interim/ProRS_aligned.fasta (already aligned)
# - Outputs: aaRS/poincare_runs/proRS/poincare_knn${KNN}_gamma${GAMMA}/PM*.csv + PNG

set -euo pipefail

# Paths (override via env vars if needed)
POINCARE_REPO="${POINCARE_REPO:-/storage/kiran-stuff/PoincareMSA}"
MSA_SRC="${MSA_SRC:-/storage/kiran-stuff/aaRS/phase1/data/interim/ProRS_aligned.fasta}"
RUN_ROOT="${RUN_ROOT:-/storage/kiran-stuff/aaRS/poincare_runs/proRS}"

# Hyperparameters
GAP_TH="${GAP_TH:-0.9}"
KNN="${KNN:-5}"
GAMMA="${GAMMA:-2}"
BATCHSIZE="${BATCHSIZE:-4}"
EPOCHS="${EPOCHS:-1000}"
SEED="${SEED:-42}"

CONDA_ENV="${CONDA_ENV:-env_poincare}"
# Reduce OpenMP / MKL noise and avoid SHM issues in restricted envs
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
export MKL_NUM_THREADS=${MKL_NUM_THREADS:-1}
export OMP_THREAD_LIMIT=${OMP_THREAD_LIMIT:-1}
export OMP_WAIT_POLICY=${OMP_WAIT_POLICY:-PASSIVE}
export JOBLIB_TEMP_FOLDER=${JOBLIB_TEMP_FOLDER:-/tmp}
export MKL_THREADING_LAYER=${MKL_THREADING_LAYER:-GNU}
export KMP_AFFINITY=${KMP_AFFINITY:-disabled}
export KMP_INIT_AT_FORK=${KMP_INIT_AT_FORK:-FALSE}

if [[ ! -d "${POINCARE_REPO}" ]]; then
  echo "PoincareMSA repo not found at ${POINCARE_REPO}" >&2
  exit 1
fi

if [[ ! -f "${MSA_SRC}" ]]; then
  echo "Input MSA not found at ${MSA_SRC}" >&2
  exit 1
fi

mkdir -p "${RUN_ROOT}"
MSA_MFASTA="${RUN_ROOT}/proRS_full.mfasta"

# Ensure MSA has .mfasta extension expected by the prep script
cp "${MSA_SRC}" "${MSA_MFASTA}"

# Activate conda env (set MKL vars to avoid nounset issues in activate scripts)
export MKL_INTERFACE_LAYER=${MKL_INTERFACE_LAYER:-LP64}
if command -v conda >/dev/null 2>&1; then
  eval "$(conda shell.bash hook)"
  conda activate "${CONDA_ENV}"
else
  echo "conda not found; please load it before running this script." >&2
  exit 1
fi

echo "=== Running PoincaréMSA for ProRS ==="
echo "Repo       : ${POINCARE_REPO}"
echo "MSA        : ${MSA_MFASTA}"
echo "Run root   : ${RUN_ROOT}"
echo "Gap thr    : ${GAP_TH}"
echo "knn/gamma  : ${KNN}/${GAMMA}"
echo "epochs/bsz : ${EPOCHS}/${BATCHSIZE}"

PREP_OUT="${RUN_ROOT}/prep_gap${GAP_TH}"
mkdir -p "${PREP_OUT}"

# Step 1: MSA -> per-sequence PSSMs
(
  cd "${POINCARE_REPO}"
  bash scripts/prepare_data/create_projection.sh \
    scripts/prepare_data \
    "${MSA_MFASTA}" \
    "${PREP_OUT}" \
    proRS \
    "${GAP_TH}"
)

# Step 2: Build Poincaré embedding (2D)
OUTPUT_DIR="${RUN_ROOT}/poincare_knn${KNN}_gamma${GAMMA}/"
mkdir -p "${OUTPUT_DIR}"
(
  cd "${POINCARE_REPO}/scripts/build_poincare_map"
  python main.py \
    --input_path "${PREP_OUT}/fasta${GAP_TH}/" \
    --output_path "${OUTPUT_DIR}" \
    --knn "${KNN}" \
    --gamma "${GAMMA}" \
    --batchsize "${BATCHSIZE}" \
    --epochs "${EPOCHS}" \
    --plot False \
    --debugplot 0 \
    --seed "${SEED}"
)

echo "=== Done. Outputs under ${OUTPUT_DIR} (CSV + PNG). ==="
