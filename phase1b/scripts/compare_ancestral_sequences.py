#!/usr/bin/env python3
"""Compare shallow (eukaryotic) vs deep (LUCA) ancestral sequences"""

from Bio import SeqIO, pairwise2
from Bio.Align import substitution_matrices
import numpy as np

print("="*70)
print("COMPARING ANCESTRAL RECONSTRUCTIONS")
print("="*70)

# Load both ancestors
shallow = SeqIO.read("../phase1/results/Anc-ProThrRS.fasta", "fasta")
deep = SeqIO.read("results/Anc-ProThrRS-LUCA.fasta", "fasta")

print(f"\nShallow (Eukaryotic-biased): {len(shallow.seq)} aa")
print(f"Deep (Cross-domain):         {len(deep.seq)} aa")

# Global alignment
alignments = pairwise2.align.globalxx(str(shallow.seq), str(deep.seq))
best_aln = alignments[0]

identity = sum(1 for a, b in zip(best_aln.seqA, best_aln.seqB) if a == b)
total = len(best_aln.seqA)
pct_id = 100 * identity / total

print(f"\nPairwise identity: {pct_id:.1f}%")

if pct_id > 70:
    print("  ✓ Very similar - shallow ancestor is reliable")
elif pct_id > 50:
    print("  ⚠ Moderately similar - some differences")
else:
    print("  ✗ Very different - shallow ancestor was biased")

# Amino acid composition comparison
def get_composition(seq):
    total = len(seq)
    return {
        'hydrophobic': 100 * sum(seq.count(aa) for aa in 'AVILMFYW') / total,
        'charged': 100 * sum(seq.count(aa) for aa in 'DEKR') / total,
        'polar': 100 * sum(seq.count(aa) for aa in 'STNQ') / total,
        'special': 100 * sum(seq.count(aa) for aa in 'CGP') / total,
    }

shallow_comp = get_composition(str(shallow.seq))
deep_comp = get_composition(str(deep.seq))

print("\nAmino acid composition:")
print(f"  Property        Shallow   Deep    Diff")
print(f"  {'='*45}")
for prop in ['hydrophobic', 'charged', 'polar', 'special']:
    diff = deep_comp[prop] - shallow_comp[prop]
    sign = '+' if diff > 0 else ''
    print(f"  {prop:12s}    {shallow_comp[prop]:5.1f}%   {deep_comp[prop]:5.1f}%   {sign}{diff:4.1f}%")

print("\n" + "="*70)
print("RECOMMENDATION FOR PHASE 2")
print("="*70)

if pct_id < 60:
    print("✓ USE BOTH ancestors in Phase 2:")
    print("  1. Shallow (eukaryotic) - tests recent evolution")
    print("  2. Deep (LUCA) - tests deep ancestral promiscuity")
    print("\nThis gives TWO independent tests of your hypothesis!")
else:
    print("✓ Ancestors are similar - use Deep (LUCA) version")
    print("  This has better phylogenetic support")

