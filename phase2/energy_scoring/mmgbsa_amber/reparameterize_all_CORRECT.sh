#!/bin/bash
# CORRECT ligand parameterization: obabel (add H) → antechamber (GAFF2 + Gasteiger)
# CRITICAL: Must add H BEFORE antechamber, not with -j flag

set -e

mkdir -p ligands_fixed

AA_LIST="ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL"

echo "========================================================================"
echo "CORRECT Ligand Parameterization"
echo "Method: obabel (add H at pH 7) → antechamber (GAFF2 + Gasteiger)"
echo "========================================================================"
echo ""

for aa in $AA_LIST; do
    echo "=== $aa ==="

    mkdir -p ligands_fixed/$aa

    # Determine net charge at pH ~7
    nc=0
    if [ "$aa" = "ARG" ] || [ "$aa" = "LYS" ]; then
        nc=1
        echo "  Net charge: +1 (basic)"
    elif [ "$aa" = "ASP" ] || [ "$aa" = "GLU" ]; then
        nc=-1
        echo "  Net charge: -1 (acidic)"
    else
        nc=0
        echo "  Net charge: 0 (neutral/zwitterion)"
    fi

    # Step 1: Add hydrogens with obabel
    echo "  1. Adding hydrogens with obabel..."
    obabel ligands/$aa/ligand.pdb -O ligands_fixed/$aa/lig_H.pdb -h -p 7.0 2>&1 | grep -v "Warning" || true

    if [ ! -f "ligands_fixed/$aa/lig_H.pdb" ]; then
        echo "  ✗ FAILED: obabel did not create lig_H.pdb"
        continue
    fi

    n_atoms_H=$(grep -c "^HETATM\|^ATOM" ligands_fixed/$aa/lig_H.pdb || echo "0")
    echo "  2. Hydrogenated: $n_atoms_H atoms"

    # Step 2: Parameterize hydrogenated structure with GAFF2 + Gasteiger
    echo "  3. Running antechamber (GAFF2 + Gasteiger)..."
    cd ligands_fixed/$aa

    /home/kiran/miniforge3/envs/amber_mmgbsa/bin/antechamber \
        -i lig_H.pdb \
        -fi pdb \
        -o lig.mol2 \
        -fo mol2 \
        -at gaff2 \
        -c gas \
        -nc $nc \
        -rn LIG \
        -s 2 \
        > antechamber.log 2>&1

    if [ $? -ne 0 ]; then
        echo "  ✗ FAILED: antechamber failed"
        tail -10 antechamber.log
        cd ../..
        continue
    fi

    # Step 3: parmchk2
    echo "  4. Running parmchk2..."
    /home/kiran/miniforge3/envs/amber_mmgbsa/bin/parmchk2 \
        -i lig.mol2 \
        -f mol2 \
        -o lig.frcmod \
        > parmchk2.log 2>&1

    if [ $? -ne 0 ]; then
        echo "  ✗ FAILED: parmchk2 failed"
        cd ../..
        continue
    fi

    cd ../..

    # Verify outputs and count H
    if [ ! -f "ligands_fixed/$aa/lig.mol2" ] || [ ! -f "ligands_fixed/$aa/lig.frcmod" ]; then
        echo "  ✗ FAILED: Output files missing"
        continue
    fi

    # Count atoms and H
    n_atoms=$(grep -A 1000 "@<TRIPOS>ATOM" ligands_fixed/$aa/lig.mol2 | grep -B 1000 "@<TRIPOS>BOND" | grep -v "@<TRIPOS>" | wc -l)
    n_H=$(grep -A 1000 "@<TRIPOS>ATOM" ligands_fixed/$aa/lig.mol2 | grep -B 1000 "@<TRIPOS>BOND" | grep -E " h[0-9]* | H[0-9]* " | wc -l)

    echo "  ✓ SUCCESS: $n_atoms total atoms, $n_H hydrogens"
    echo ""
done

echo "========================================================================"
echo "Summary:"
echo "========================================================================"
echo ""

for aa in $AA_LIST; do
    if [ -f "ligands_fixed/$aa/lig.mol2" ]; then
        n_atoms=$(grep -A 1000 "@<TRIPOS>ATOM" ligands_fixed/$aa/lig.mol2 | grep -B 1000 "@<TRIPOS>BOND" | grep -v "@<TRIPOS>" | wc -l)
        n_H=$(grep -A 1000 "@<TRIPOS>ATOM" ligands_fixed/$aa/lig.mol2 | grep -B 1000 "@<TRIPOS>BOND" | grep -E " h[0-9] | H[0-9] " | wc -l)
        printf "  ✓ %-3s: %2d atoms (%2d H)\n" "$aa" "$n_atoms" "$n_H"
    else
        echo "  ✗ $aa: FAILED"
    fi
done

echo ""
echo "========================================================================"
echo "Next: Update run_mmpbsa_batch.py to use ligands_fixed/"
echo "========================================================================"
