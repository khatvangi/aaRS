#!/usr/bin/env python3
"""Collect ProRS/ThrRS from bacteria + archaea for deep phylogeny"""

import requests
from Bio import SeqIO, Entrez
import time
from pathlib import Path

Entrez.email = "your.email@example.com"

# Curated list: thermophiles + model organisms (high-quality genomes)
TARGET_ORGANISMS = {
    "bacteria": [
        # Thermophiles (closer to ancestral state)
        "Thermus thermophilus",
        "Thermotoga maritima",
        "Aquifex aeolicus",
        "Deinococcus radiodurans",
        
        # Deeply branching
        "Chloroflexus aurantiacus",
        "Thermomicrobium roseum",
        
        # Model organisms (well-annotated)
        "Escherichia coli",
        "Bacillus subtilis",
        "Pseudomonas aeruginosa",
        "Mycobacterium tuberculosis",
        
        # Diverse phyla
        "Synechocystis sp. PCC 6803",  # Cyanobacteria
        "Geobacter sulfurreducens",    # Deltaproteobacteria
        "Streptomyces coelicolor",     # Actinobacteria
        "Chlorobium tepidum",          # Green sulfur bacteria
        "Bacteroides fragilis",        # Bacteroidetes
    ],
    
    "archaea": [
        # Euryarchaeota
        "Methanocaldococcus jannaschii",
        "Methanothermobacter thermautotrophicus",
        "Pyrococcus furiosus",
        "Thermococcus kodakarensis",
        "Archaeoglobus fulgidus",
        
        # Crenarchaeota
        "Sulfolobus acidocaldarius",
        "Pyrobaculum aerophilum",
        "Thermoproteus neutrophilus",
        
        # Other archaea
        "Halobacterium salinarum",
        "Thermoplasma acidophilum",
    ]
}

def search_and_fetch_aars(organism, aars_name):
    """Search NCBI for specific aaRS from organism"""
    
    # Search query
    if aars_name == "ProRS":
        query = f'"{organism}"[Organism] AND (prolyl-tRNA synthetase[Protein Name] OR ProS[Gene Name])'
    elif aars_name == "ThrRS":
        query = f'"{organism}"[Organism] AND (threonyl-tRNA synthetase[Protein Name] OR ThrS[Gene Name])'
    else:
        return None
    
    try:
        # Search
        handle = Entrez.esearch(db="protein", term=query, retmax=5)
        record = Entrez.read(handle)
        handle.close()
        
        if not record['IdList']:
            return None
        
        # Fetch first result (most relevant)
        protein_id = record['IdList'][0]
        handle = Entrez.efetch(db="protein", id=protein_id, rettype="fasta", retmode="text")
        seq_record = SeqIO.read(handle, "fasta")
        handle.close()
        
        # Add organism to description
        seq_record.description = f"{seq_record.description} [{organism}]"
        
        return seq_record
        
    except Exception as e:
        print(f"  Error fetching {organism} {aars_name}: {e}")
        return None

def collect_all_sequences():
    """Collect ProRS and ThrRS from all target organisms"""
    
    Path("data/raw").mkdir(parents=True, exist_ok=True)
    
    for aars in ["ProRS", "ThrRS"]:
        print(f"\n{'='*70}")
        print(f"Collecting {aars}")
        print(f"{'='*70}\n")
        
        sequences = []
        
        for domain, organisms in TARGET_ORGANISMS.items():
            print(f"\n{domain.upper()}:")
            for org in organisms:
                print(f"  {org}...", end=" ", flush=True)
                
                seq = search_and_fetch_aars(org, aars)
                if seq:
                    sequences.append(seq)
                    print(f"✓ ({len(seq.seq)} aa)")
                else:
                    print("✗ not found")
                
                time.sleep(0.5)  # Rate limiting
        
        # Save
        output_file = f"data/raw/{aars}_prokaryotes.fasta"
        SeqIO.write(sequences, output_file, "fasta")
        
        print(f"\n✓ Saved {len(sequences)} {aars} sequences to {output_file}")
        
        # Count by domain
        bacteria_count = sum(1 for s in sequences if any(b in s.description for b in TARGET_ORGANISMS["bacteria"]))
        archaea_count = sum(1 for s in sequences if any(a in s.description for a in TARGET_ORGANISMS["archaea"]))
        
        print(f"  Bacteria: {bacteria_count}")
        print(f"  Archaea: {archaea_count}")

if __name__ == "__main__":
    print("="*70)
    print("PROKARYOTIC aaRS COLLECTION")
    print("Target: Deep phylogeny (bacteria + archaea)")
    print("="*70)
    
    collect_all_sequences()
    
    print("\n" + "="*70)
    print("COLLECTION COMPLETE")
    print("="*70)
    print("\nNext steps:")
    print("  1. Add eukaryotic sequences from Phase 1")
    print("  2. Align and reconstruct deep ancestral node")

