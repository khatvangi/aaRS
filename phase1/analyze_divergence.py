#!/usr/bin/env python3
"""Deep dive: Is 7-17% identity actually too divergent?"""

from Bio import SeqIO, AlignIO
import numpy as np
from collections import Counter

print("="*70)
print("DIVERGENCE ANALYSIS: Is This Actually A Problem?")
print("="*70)

# Load ProRS alignment
aln = AlignIO.read("data/interim/ProRS_aligned.fasta", "fasta")

# Extract species/taxonomy from sequence IDs
print("\n[1] Taxonomic Composition")
print("-"*70)

species_distribution = {"Bacteria": 0, "Archaea": 0, "Eukaryota": 0, "Unknown": 0}

for rec in aln:
    desc = rec.description.lower()
    # Simple heuristic classification
    if any(x in desc for x in ['ecoli', 'coli', 'bacillus', 'subtilis', 'mycobacterium']):
        species_distribution["Bacteria"] += 1
    elif any(x in desc for x in ['methanocaldococcus', 'pyrococcus', 'sulfolobus', 'archaeo']):
        species_distribution["Archaea"] += 1
    elif any(x in desc for x in ['human', 'mouse', 'homo', 'mus', 'drosophila', 'cerevisiae', 'sapiens']):
        species_distribution["Eukaryota"] += 1
    else:
        species_distribution["Unknown"] += 1

print("Taxonomic breakdown:")
for domain, count in species_distribution.items():
    pct = 100 * count / len(aln)
    print(f"  {domain}: {count} ({pct:.0f}%)")

# Check if bacteria-dominated
if species_distribution["Bacteria"] > 0.7 * len(aln):
    print("\n⚠ BACTERIA-DOMINATED SAMPLE")
    print("  This should have HIGHER identity (~40-60%)")
    print("  Low identity suggests alignment issues or very divergent paralogs")
elif sum([species_distribution["Bacteria"], species_distribution["Archaea"]]) < 5:
    print("\n⚠ PROKARYOTE-POOR SAMPLE")
    print("  Missing deep evolutionary diversity")
else:
    print("\n✓ BALANCED CROSS-DOMAIN SAMPLE")
    print("  Low identity (7-17%) is EXPECTED for deep phylogeny")

# [2] Check alignment quality
print("\n[2] Alignment Quality Assessment")
print("-"*70)

aln_array = np.array([list(str(rec.seq)) for rec in aln])
n_seqs, aln_len = aln_array.shape

# Calculate conservation per column
conservation_scores = []
for col_idx in range(aln_len):
    col = aln_array[:, col_idx]
    # Remove gaps
    col_no_gaps = col[col != '-']
    if len(col_no_gaps) > 0:
        # Most common residue frequency
        most_common_count = Counter(col_no_gaps).most_common(1)[0][1]
        conservation = most_common_count / len(col_no_gaps)
        conservation_scores.append(conservation)

mean_conservation = np.mean(conservation_scores)
print(f"Mean column conservation: {mean_conservation:.3f}")

# Classify columns
highly_conserved = sum(1 for x in conservation_scores if x > 0.8)
moderately_conserved = sum(1 for x in conservation_scores if 0.5 < x <= 0.8)
variable = sum(1 for x in conservation_scores if x <= 0.5)

print(f"Column types:")
print(f"  Highly conserved (>80%): {highly_conserved} ({100*highly_conserved/len(conservation_scores):.1f}%)")
print(f"  Moderately conserved (50-80%): {moderately_conserved} ({100*moderately_conserved/len(conservation_scores):.1f}%)")
print(f"  Variable (<50%): {variable} ({100*variable/len(conservation_scores):.1f}%)")

if highly_conserved > 0.2 * len(conservation_scores):
    print("\n✓ GOOD: >20% highly conserved columns")
    print("  Sufficient signal for ancestral reconstruction")
elif highly_conserved > 0.1 * len(conservation_scores):
    print("\n⚠ MARGINAL: 10-20% highly conserved columns")
    print("  Ancestral reconstruction possible but uncertain")
else:
    print("\n✗ POOR: <10% highly conserved columns")
    print("  Alignment may be unreliable")

# [3] Check for alignment artifacts
print("\n[3] Alignment Artifact Detection")
print("-"*70)

gap_content = []
for rec in aln:
    gaps = str(rec.seq).count('-')
    gap_pct = 100 * gaps / len(rec.seq)
    gap_content.append(gap_pct)

mean_gap = np.mean(gap_content)
max_gap = np.max(gap_content)

print(f"Gap content:")
print(f"  Mean: {mean_gap:.1f}%")
print(f"  Max: {max_gap:.1f}%")

if mean_gap < 30:
    print("  ✓ Reasonable gap content")
elif mean_gap < 50:
    print("  ⚠ High gap content")
else:
    print("  ✗ Very high gap content (possible alignment failure)")

# [4] Literature comparison
print("\n[4] Literature Comparison")
print("-"*70)
print("What does the literature say about aaRS divergence?")
print()
print("From Douglas et al. (2023) 'aaRS phylogeny':")
print("  - aaRS are among the MOST CONSERVED proteins")
print("  - Cross-domain alignments typically show 20-40% identity")
print("  - Below 20% indicates either:")
print("    a) Very deep evolutionary time (pre-LUCA)")
print("    b) Alignment artifacts")
print("    c) Inclusion of paralogs")
print()

if mean_conservation > 0.3:
    print("Our data: CONSISTENT with deep phylogeny")
    print("  Column conservation (~30%) suggests real homology")
elif mean_conservation > 0.2:
    print("Our data: MARGINAL but likely valid")
    print("  Some signal present but noisy")
else:
    print("Our data: QUESTIONABLE")
    print("  May need to check for paralogs or alignment errors")

# [5] Final verdict
print("\n" + "="*70)
print("VERDICT: Is Low Pairwise Identity A Problem?")
print("="*70)

if highly_conserved > 0.15 * len(conservation_scores) and mean_conservation > 0.25:
    print("✓ NO - Low *pairwise* identity is expected for cross-domain phylogeny")
    print("  What matters is *column conservation*, which is adequate")
    print("  Ancestral reconstruction is valid, though uncertain")
elif species_distribution["Bacteria"] > 0.8 * len(aln):
    print("⚠ MAYBE - Sample is bacteria-dominated, identity should be higher")
    print("  Consider adding more diverse sequences")
else:
    print("✗ YES - Both pairwise identity AND column conservation are low")
    print("  Alignment quality is suspect")
    print("  Recommend re-aligning or filtering sequences")

