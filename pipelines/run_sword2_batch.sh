#!/usr/bin/env bash

# Batch runner: convert AF3 mmCIFs -> PDB and run SWORD2 domain/PU decomposition.
# - Assumes SWORD2 is cloned at /storage/kiran-stuff/SWORD2
# - Default inputs: three ProRS AF3 models (edit list below as needed)
# - Outputs land in aaRS/sword2_results/runs/

set -euo pipefail

SWORD_DIR="${SWORD_DIR:-/storage/kiran-stuff/SWORD2}"
OUT_ROOT="${OUT_ROOT:-/storage/kiran-stuff/aaRS/sword2_results}"
CONDA_ENV="${CONDA_ENV:-sword2}"

# List of structures to process ("path:CHAIN"). Edit to suit.
STRUCTURES=(
  "/storage/kiran-stuff/aaRS/phase2/outputs/deep_editing_pro/deep_editing_pro/deep_editing_pro_model.cif:A"
  "/storage/kiran-stuff/aaRS/phase2/outputs/deep_domain_pro/deep_domain_pro_model.cif:A"
  "/storage/kiran-stuff/aaRS/phase2/outputs/modern_prours_pro/modern_prours_pro/modern_prours_pro_model.cif:A"
)

if [[ ! -d "${SWORD_DIR}" ]]; then
  echo "SWORD2 repo not found at ${SWORD_DIR}" >&2
  exit 1
fi

mkdir -p "${OUT_ROOT}"
PDB_CACHE="${OUT_ROOT}/pdb_inputs"
RUN_DIR="${OUT_ROOT}/runs"
mkdir -p "${PDB_CACHE}" "${RUN_DIR}"

# Activate conda env
if command -v conda >/dev/null 2>&1; then
  eval "$(conda shell.bash hook)"
  conda activate "${CONDA_ENV}"
else
  echo "conda not found; please load it before running this script." >&2
  exit 1
fi

convert_cif_to_pdb() {
  local cif_path="$1"
  local pdb_path="$2"
  local chain_id="$3"
  python - "$cif_path" "$pdb_path" "$chain_id" <<'PY'
import sys
from pathlib import Path
from Bio.PDB import MMCIFParser, PDBIO, Select

cif_path = Path(sys.argv[1])
pdb_path = Path(sys.argv[2])
chain_id = sys.argv[3]

if not cif_path.exists():
    sys.exit(f"Missing input file: {cif_path}")

parser = MMCIFParser(QUIET=True)
structure = parser.get_structure("model", cif_path)

class ChainSelect(Select):
    def accept_chain(self, chain):
        return (chain_id == "") or (chain.id == chain_id)

io = PDBIO()
io.set_structure(structure)
io.save(pdb_path.as_posix(), select=ChainSelect())
PY
}

echo "=== Running SWORD2 batch ==="
echo "Repo   : ${SWORD_DIR}"
echo "Outputs: ${RUN_DIR}"

cd "${SWORD_DIR}"

for entry in "${STRUCTURES[@]}"; do
  IFS=":" read -r cif_path chain_id <<<"${entry}"
  base="$(basename "${cif_path}")"
  stub="${base%.*}"
  pdb_path="${PDB_CACHE}/${stub}_chain${chain_id}.pdb"

  echo "-> Converting ${cif_path} (chain ${chain_id})"
  convert_cif_to_pdb "${cif_path}" "${pdb_path}" "${chain_id}"

  echo "-> SWORD2 on ${pdb_path}"
  ./SWORD2.py \
    -i "${pdb_path}" \
    -c "${chain_id}" \
    -o "${RUN_DIR}" \
    --disable-energies \
    --disable-plots
done

echo "=== Done. Summaries in ${RUN_DIR}/*/SWORD2_summary.json ==="
