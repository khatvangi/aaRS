#!/usr/bin/env python3
# 12_quality_control.py
"""Validate Phase 1 outputs"""

import sys
from pathlib import Path
from Bio import SeqIO, AlignIO
import json

def check_file_exists(filepath, description):
    """Check if file exists and is non-empty"""
    path = Path(filepath)
    if not path.exists():
        print(f"✗ FAIL: {description} - file not found: {filepath}")
        return False
    
    if path.stat().st_size == 0:
        print(f"✗ FAIL: {description} - file is empty: {filepath}")
        return False
    
    print(f"✓ PASS: {description}")
    return True

def validate_fasta_sequences(filepath, min_length=100):
    """Validate FASTA sequences"""
    try:
        records = list(SeqIO.parse(filepath, "fasta"))
        
        if len(records) == 0:
            print(f"  ✗ No sequences in {filepath}")
            return False
        
        for rec in records:
            if len(rec.seq) < min_length:
                print(f"  ✗ Sequence {rec.id} too short: {len(rec.seq)} < {min_length}")
                return False
        
        print(f"  ✓ {len(records)} sequences, lengths: {min([len(r.seq) for r in records])}-{max([len(r.seq) for r in records])}")
        return True
    
    except Exception as e:
        print(f"  ✗ Error parsing FASTA: {e}")
        return False

def validate_alignment(filepath):
    """Validate MSA"""
    try:
        aln = AlignIO.read(filepath, "fasta")
        n_seqs = len(aln)
        aln_len = aln.get_alignment_length()
        
        # Check for all-gap columns
        gap_cols = 0
        for i in range(aln_len):
            col = aln[:, i]
            if col.count("-") == n_seqs:
                gap_cols += 1
        
        gap_pct = 100 * gap_cols / aln_len
        
        if gap_pct > 50:
            print(f"  ✗ Too many gap columns: {gap_pct:.1f}%")
            return False
        
        print(f"  ✓ {n_seqs} sequences × {aln_len} columns, {gap_pct:.1f}% gap columns")
        return True
    
    except Exception as e:
        print(f"  ✗ Error parsing alignment: {e}")
        return False

def validate_tree(filepath):
    """Validate phylogenetic tree"""
    try:
        from ete3 import Tree
        tree = Tree(filepath, format=1)
        
        n_leaves = len(tree.get_leaves())
        n_internal = len([n for n in tree.traverse() if not n.is_leaf()])
        
        print(f"  ✓ {n_leaves} leaves, {n_internal} internal nodes")
        return True
    
    except Exception as e:
        print(f"  ✗ Error parsing tree: {e}")
        return False

def main():
    print("=== Phase 1 Quality Control ===\n")
    
    all_passed = True
    
    # 1. Check aaRS outputs
    print("1. aaRS Sequences and Alignments")
    for aars in ["ProRS", "ThrRS", "SerRS", "ValRS"]:
        # Raw sequences
        if check_file_exists(f"data/raw/{aars}_filtered.fasta", f"{aars} raw sequences"):
            validate_fasta_sequences(f"data/raw/{aars}_filtered.fasta", min_length=200)
        else:
            all_passed = False
        
        # Catalytic domains
        if check_file_exists(f"data/interim/{aars}_catalytic.fasta", f"{aars} catalytic domains"):
            validate_fasta_sequences(f"data/interim/{aars}_catalytic.fasta", min_length=100)
        else:
            all_passed = False
        
        # Alignments
        if check_file_exists(f"data/interim/{aars}_aligned.fasta", f"{aars} alignment"):
            validate_alignment(f"data/interim/{aars}_aligned.fasta")
        else:
            all_passed = False
        
        # Trees
        if check_file_exists(f"results/{aars}_tree.treefile", f"{aars} phylogeny"):
            validate_tree(f"results/{aars}_tree.treefile")
        else:
            all_passed = False
    
    print("\n2. tRNA Sequences and Alignments")
    for trna in ["Pro", "Thr", "Ser", "Val"]:
        # Raw sequences
        if check_file_exists(f"data/raw/tRNA_{trna}_all.fasta", f"tRNA-{trna} raw sequences"):
            validate_fasta_sequences(f"data/raw/tRNA_{trna}_all.fasta", min_length=60)
        else:
            all_passed = False
        
        # Alignments
        if check_file_exists(f"data/interim/tRNA_{trna}_aligned.fasta", f"tRNA-{trna} alignment"):
            validate_alignment(f"data/interim/tRNA_{trna}_aligned.fasta")
        else:
            all_passed = False
        
        # Trees
        if check_file_exists(f"results/tRNA_{trna}_tree.treefile", f"tRNA-{trna} phylogeny"):
            validate_tree(f"results/tRNA_{trna}_tree.treefile")
        else:
            all_passed = False
    
    print("\n3. Final Ancestral Sequences")
    
    # Ancestral aaRS
    if check_file_exists("results/Anc-ProThrRS.fasta", "Ancestral ProThrRS"):
        if validate_fasta_sequences("results/Anc-ProThrRS.fasta", min_length=200):
            rec = SeqIO.read("results/Anc-ProThrRS.fasta", "fasta")
            print(f"  ✓ Anc-ProThrRS: {len(rec.seq)} residues")
        else:
            all_passed = False
    else:
        all_passed = False
    
    # Ancestral tRNA
    if check_file_exists("results/Anc-tRNA-ProThr.fasta", "Ancestral tRNA-ProThr"):
        if validate_fasta_sequences("results/Anc-tRNA-ProThr.fasta", min_length=60):
            rec = SeqIO.read("results/Anc-tRNA-ProThr.fasta", "fasta")
            print(f"  ✓ Anc-tRNA-ProThr: {len(rec.seq)} nucleotides")
            
            # Check tRNA structure
            if 70 <= len(rec.seq) <= 90:
                print("  ✓ Length consistent with tRNA (70-90 nt)")
            else:
                print(f"  ⚠ WARNING: Unusual tRNA length: {len(rec.seq)} nt")
        else:
            all_passed = False
    else:
        all_passed = False
    
    # Final verdict
    print("\n" + "="*50)
    if all_passed:
        print("✓ ALL CHECKS PASSED - Phase 1 complete and valid")
        sys.exit(0)
    else:
        print("✗ SOME CHECKS FAILED - Review errors above")
        sys.exit(1)

if __name__ == "__main__":
    main()
