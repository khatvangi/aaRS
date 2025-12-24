#!/usr/bin/env python3
"""Extract ancestral node sequences from IQ-TREE output"""

import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pathlib import Path

def parse_iqtree_asr(state_file):
    """Parse IQ-TREE ancestral state file (.state or .iqtree)"""
    ancestors = {}
    current_node = None
    current_seq = []
    
    if not Path(state_file).exists():
        print(f"WARNING: {state_file} not found, trying alternative...")
        # Try without suffix
        base = str(state_file).replace('.state', '')
        for ext in ['.state', '.iqtree', '.ancestral']:
            alt = f"{base}{ext}"
            if Path(alt).exists():
                state_file = alt
                print(f"Found: {alt}")
                break
        else:
            raise FileNotFoundError(f"Could not find ASR file for {state_file}")
    
    with open(state_file) as f:
        for line in f:
            line = line.strip()
            
            # Node header
            if line.startswith("Node"):
                if current_node and current_seq:
                    ancestors[current_node] = "".join(current_seq)
                
                current_node = line.split()[0]
                current_seq = []
            
            # Sequence line (position AA prob)
            elif line and current_node:
                parts = line.split()
                if len(parts) >= 2:
                    pos = parts[0]
                    aa = parts[1]
                    if aa.isalpha() and len(aa) == 1:
                        current_seq.append(aa)
        
        # Last node
        if current_node and current_seq:
            ancestors[current_node] = "".join(current_seq)
    
    return ancestors

def extract_anc_prothrs():
    """Extract ancestral ProRS/ThrRS sequence"""
    print("=== Extracting Anc-ProThrRS ===")
    
    # Find the actual ASR file
    pro_file = None
    for candidate in ['results/ProRS.state', 'results/ProRS_asr.state', 'results/ProRS.iqtree']:
        if Path(candidate).exists():
            pro_file = candidate
            break
    
    if not pro_file:
        print("ERROR: No ProRS ASR file found")
        print("Available files:")
        for f in Path('results').glob('ProRS*'):
            print(f"  {f}")
        return
    
    pro_ancestors = parse_iqtree_asr(pro_file)
    
    if not pro_ancestors:
        print(f"ERROR: No ancestors found in {pro_file}")
        return
    
    # Take root node (usually Node1 or similar)
    root_nodes = [n for n in pro_ancestors.keys() if 'Node' in n]
    if not root_nodes:
        print(f"ERROR: No Node entries found. Available keys: {list(pro_ancestors.keys())[:5]}")
        return
    
    # Sort to get earliest node
    root_nodes.sort()
    mrca_node = root_nodes[0]
    
    anc_seq = pro_ancestors[mrca_node]
    
    record = SeqRecord(
        Seq(anc_seq),
        id="Anc-ProThrRS",
        description=f"Ancestral ProRS/ThrRS from {mrca_node}",
    )
    
    SeqIO.write([record], "results/Anc-ProThrRS.fasta", "fasta")
    print(f"✓ Anc-ProThrRS: {len(anc_seq)} residues (from {mrca_node})")

def main():
    extract_anc_prothrs()

if __name__ == "__main__":
    main()
