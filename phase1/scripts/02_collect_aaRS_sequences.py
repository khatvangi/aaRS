#!/usr/bin/env python3
# 02_collect_aaRS_sequences.py
"""Parallel collection of aaRS sequences from UniProt"""

import requests
import time
import json
from pathlib import Path
from multiprocessing import Pool
from Bio import SeqIO
from tqdm import tqdm

def query_uniprot(aars_type, keywords, max_results=1000):
    """Query UniProt REST API for aaRS sequences"""
    # Build query string
    keyword_query = " OR ".join([f'protein_name:"{kw}"' for kw in keywords])
    query = f"({keyword_query}) AND reviewed:true"
    
    url = "https://rest.uniprot.org/uniprotkb/stream"
    params = {
        "query": query,
        "format": "fasta",
        "size": max_results,
    }
    
    print(f"  Querying UniProt for {aars_type}...")
    try:
        response = requests.get(url, params=params, timeout=120)
        response.raise_for_status()
        return response.text
    except Exception as e:
        print(f"  ERROR: {e}")
        return None

def filter_by_species(fasta_text, target_species_file):
    """Filter sequences to target species only"""
    target_species = set()
    with open(target_species_file) as f:
        target_species = {line.strip() for line in f}
    
    from io import StringIO
    filtered = []
    
    for record in SeqIO.parse(StringIO(fasta_text), "fasta"):
        # Extract organism from header (OS=...)
        if "OS=" in record.description:
            org_start = record.description.index("OS=") + 3
            org_end = record.description.find(" OX=", org_start)
            if org_end == -1:
                org_end = record.description.find(" GN=", org_start)
            if org_end == -1:
                org_end = len(record.description)
            
            organism = record.description[org_start:org_end].strip()
            
            if organism in target_species:
                filtered.append(record)
    
    return filtered

def collect_one_aars(args):
    """Worker function for parallel collection"""
    aars_type, info = args
    
    output_raw = Path(f"data/raw/{aars_type}_raw.fasta")
    output_filtered = Path(f"data/raw/{aars_type}_filtered.fasta")
    
    # Skip if already done
    if output_filtered.exists():
        n_seqs = len(list(SeqIO.parse(output_filtered, "fasta")))
        print(f"✓ {aars_type}: {n_seqs} sequences (cached)")
        return aars_type, n_seqs
    
    # Query UniProt
    fasta_text = query_uniprot(aars_type, info["uniprot_keywords"])
    if not fasta_text:
        return aars_type, 0
    
    # Save raw
    output_raw.parent.mkdir(parents=True, exist_ok=True)
    with open(output_raw, "w") as f:
        f.write(fasta_text)
    
    # Filter to target species
    filtered = filter_by_species(fasta_text, "data/species_list.txt")
    SeqIO.write(filtered, output_filtered, "fasta")
    
    print(f"✓ {aars_type}: {len(filtered)} sequences")
    return aars_type, len(filtered)

def main():
    # Load config
    with open("data/config.json") as f:
        config = json.load(f)
    
    # Parallel collection (4 aaRS types)
    print("=== Collecting aaRS sequences from UniProt ===")
    
    with Pool(processes=4) as pool:
        results = pool.map(collect_one_aars, config["aars"].items())
    
    # Summary
    print("\n=== Collection Summary ===")
    for aars, count in results:
        print(f"  {aars}: {count} sequences")
    
    # Save checkpoint
    with open("checkpoints/02_collection.json", "w") as f:
        json.dump(dict(results), f)

if __name__ == "__main__":
    main()
