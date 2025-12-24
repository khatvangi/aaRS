#!/usr/bin/env python3
"""Combine prokaryotic + eukaryotic sequences and align"""

from Bio import SeqIO
from pathlib import Path
import subprocess

def combine_sequences(aars):
    """Combine prokaryotic and eukaryotic sequences"""
    
    # Load prokaryotic
    prok_file = f"data/raw/{aars}_prokaryotes.fasta"
    prok_seqs = list(SeqIO.parse(prok_file, "fasta"))
    
    # Load eukaryotic from Phase 1
    euk_file = f"../phase1/data/raw/{aars}_filtered.fasta"
    euk_seqs = []
    if Path(euk_file).exists():
        euk_seqs = list(SeqIO.parse(euk_file, "fasta"))
    
    print(f"\n{aars}:")
    print(f"  Prokaryotes: {len(prok_seqs)}")
    print(f"  Eukaryotes: {len(euk_seqs)}")
    
    # Combine
    all_seqs = prok_seqs + euk_seqs
    
    # Remove duplicates by ID
    seen_ids = set()
    unique_seqs = []
    for seq in all_seqs:
        seq_id = seq.id.split('|')[1] if '|' in seq.id else seq.id
        if seq_id not in seen_ids:
            seen_ids.add(seq_id)
            unique_seqs.append(seq)
    
    print(f"  Combined: {len(unique_seqs)} unique sequences")
    
    # Check domain distribution
    bacteria = sum(1 for s in unique_seqs if any(x in s.description.lower() for x in ['thermus', 'escherichia', 'bacillus', 'thermotoga', 'aquifex']))
    archaea = sum(1 for s in unique_seqs if any(x in s.description.lower() for x in ['methanocaldococcus', 'pyrococcus', 'sulfolobus', 'archaeo']))
    eukarya = len(unique_seqs) - bacteria - archaea
    
    print(f"  Distribution: Bacteria={bacteria}, Archaea={archaea}, Eukarya={eukarya}")
    
    if bacteria < 5 or archaea < 3:
        print(f"  ⚠ Warning: Insufficient prokaryotic diversity")
    else:
        print(f"  ✓ Good cross-domain sampling")
    
    # Save combined
    output_file = f"data/raw/{aars}_combined.fasta"
    SeqIO.write(unique_seqs, output_file, "fasta")
    
    return len(unique_seqs)

def align_sequences(aars):
    """Align with MAFFT"""
    
    input_file = f"data/raw/{aars}_combined.fasta"
    output_file = f"data/interim/{aars}_aligned.fasta"
    
    print(f"\n  Aligning {aars} with MAFFT...")
    
    # Use L-INS-i (accurate but slower)
    cmd = [
        "mafft",
        "--localpair",
        "--maxiterate", "1000",
        "--thread", "16",
        "--quiet",
        input_file
    ]
    
    with open(output_file, 'w') as f:
        subprocess.run(cmd, stdout=f, check=True)
    
    print(f"  ✓ Alignment saved to {output_file}")

def main():
    print("="*70)
    print("COMBINING PROKARYOTIC + EUKARYOTIC SEQUENCES")
    print("="*70)
    
    Path("data/interim").mkdir(parents=True, exist_ok=True)
    
    for aars in ["ProRS", "ThrRS"]:
        count = combine_sequences(aars)
        if count > 0:
            align_sequences(aars)
    
    print("\n" + "="*70)
    print("✓ COMBINATION AND ALIGNMENT COMPLETE")
    print("="*70)

if __name__ == "__main__":
    main()

