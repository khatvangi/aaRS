#!/bin/bash
# Test AM1-BCC on THR after pre-adding hydrogens with obabel

cd ligands/THR

echo "=== Step 1: Add hydrogens with obabel ==="
obabel ligand.pdb -O ligand_H.pdb -h -p 7.0

echo ""
echo "=== Step 2: Check hydrogen count ==="
grep "^HETATM\|^ATOM" ligand_H.pdb | wc -l

echo ""
echo "=== Step 3: Run antechamber with AM1-BCC on hydrogenated structure ==="
/home/kiran/miniforge3/envs/amber_mmgbsa/bin/antechamber \
  -i ligand_H.pdb \
  -fi pdb \
  -o lig_am1bcc.mol2 \
  -fo mol2 \
  -at gaff2 \
  -c bcc \
  -nc 0 \
  -s 2 \
  -rn LIG

echo ""
if [ -f "lig_am1bcc.mol2" ]; then
    echo "=== SUCCESS! ==="
    echo "Generated lig_am1bcc.mol2"
    head -20 lig_am1bcc.mol2
else
    echo "=== FAILED ==="
    echo "sqm.out:"
    tail -20 sqm.out
fi
