#!/bin/bash
#
# GROMACS system preparation script
# Usage: ./prepare_system.sh <system_name> <vina_job_dir>
#
# Example: ./prepare_system.sh modern_thrrs_zn_THR modern_thrrs_ecoli_zn_THR_THR
#

set -e

NAME=$1
VINA_JOB=$2

if [ -z "$NAME" ] || [ -z "$VINA_JOB" ]; then
    echo "Usage: $0 <system_name> <vina_job_dir>"
    exit 1
fi

GMX="/usr/local/gromacs/bin/gmx_mpi"
VINA_DIR="/storage/kiran-stuff/aaRS/phase2/energy_scoring/vina_182jobs"
WORK_DIR="/storage/kiran-stuff/aaRS/phase2/energy_scoring/gromacs_md_validation"
MDP_DIR="$WORK_DIR/mdp"

# Create system directory
SYS_DIR="$WORK_DIR/systems/$NAME"
mkdir -p "$SYS_DIR"
cd "$SYS_DIR"

echo "============================================================"
echo "Preparing: $NAME"
echo "From: $VINA_JOB"
echo "============================================================"

# Copy input structure
cp "$VINA_DIR/$VINA_JOB/complex.pdb" ./complex_input.pdb

# Check for Zn
HAS_ZN=$(grep -c "^HETATM.*ZN" complex_input.pdb || true)
echo "Zn atoms found: $HAS_ZN"

# Extract protein (chain A) - clean for pdb2gmx
echo "Extracting protein (chain A)..."
grep "^ATOM" complex_input.pdb | grep " A " > protein_raw.pdb || true
# Also get HETATM that are part of protein (like modified residues)
grep "^HETATM" complex_input.pdb | grep " A " >> protein_raw.pdb 2>/dev/null || true
echo "END" >> protein_raw.pdb

# Extract ligand (chain B)
echo "Extracting ligand (chain B)..."
grep "^HETATM" complex_input.pdb | grep " B " > ligand_raw.pdb
echo "END" >> ligand_raw.pdb

# Extract Zn if present (chain C typically)
if [ "$HAS_ZN" -gt 0 ]; then
    echo "Extracting Zn..."
    grep "ZN" complex_input.pdb | head -1 > zn.pdb
    echo "END" >> zn.pdb
fi

# Count atoms
echo "Protein atoms: $(grep -c "^ATOM\|^HETATM" protein_raw.pdb)"
echo "Ligand atoms: $(grep -c "^HETATM" ligand_raw.pdb)"

# Copy MDP files
cp "$MDP_DIR"/*.mdp ./

echo ""
echo "Files created in: $SYS_DIR"
echo "  - protein_raw.pdb"
echo "  - ligand_raw.pdb"
if [ "$HAS_ZN" -gt 0 ]; then
    echo "  - zn.pdb"
fi
echo ""
echo "Next: Run pdb2gmx manually to generate topologies"
