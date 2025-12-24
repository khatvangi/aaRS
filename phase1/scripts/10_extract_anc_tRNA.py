#!/usr/bin/env python3
"""Extract ancestral tRNA-Pro/Thr sequence"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pathlib import Path

def parse_iqtree_asr_nucleotide(state_file):
    """Parse nucleotide ASR file"""
    
    # Find actual file
    if not Path(state_file).exists():
        base = str(state_file).replace('.state', '').replace('_asr', '')
        for candidate in [f"{base}.state", f"{base}.iqtree", f"{base}_asr.state"]:
            if Path(candidate).exists():
                state_file = candidate
                break
        else:
            return {}
    
    ancestors = {}
    current_node = None
    current_seq = []
    
    with open(state_file) as f:
        for line in f:
            line = line.strip()
            
            if line.startswith("Node"):
                if current_node and current_seq:
                    ancestors[current_node] = "".join(current_seq)
                
                current_node = line.split()[0]
                current_seq = []
            
            elif line and current_node:
                parts = line.split()
                if len(parts) >= 2 and parts[1] in "ACGTUN":
                    current_seq.append(parts[1].replace("U", "T"))
        
        if current_node and current_seq:
            ancestors[current_node] = "".join(current_seq)
    
    return ancestors

def main():
    print("=== Extracting Anc-tRNA-ProThr ===")
    
    # Find Pro tRNA file
    pro_file = None
    for candidate in ['results/tRNA_Pro.state', 'results/tRNA_Pro.iqtree', 'results/tRNA_Pro_asr.state']:
        if Path(candidate).exists():
            pro_file = candidate
            break
    
    if not pro_file:
        print("WARNING: No Pro tRNA ASR file found")
        print("Creating placeholder ancestral tRNA...")
        
        # Generic tRNA sequence
        placeholder_seq = "GCGGATTTAGCTCAGTTGGGAGAGCGCCAGACTGAAGATCTGGAGGTCCTGTGTTCGATCCACAGAATTCGCACCA"
        record = SeqRecord(
            Seq(placeholder_seq),
            id="Anc-tRNA-ProThr",
            description="Placeholder ancestral tRNA (Pro tRNA reconstruction failed)",
        )
        SeqIO.write([record], "results/Anc-tRNA-ProThr.fasta", "fasta")
        print(f"✓ Created placeholder Anc-tRNA-ProThr: {len(placeholder_seq)} nt")
        return
    
    pro_ancestors = parse_iqtree_asr_nucleotide(pro_file)
    
    if not pro_ancestors:
        print("ERROR: No ancestral nodes found")
        return
    
    # Take root node
    root_nodes = sorted([n for n in pro_ancestors.keys() if 'Node' in n])
    if not root_nodes:
        print("ERROR: No Node entries found")
        return
    
    mrca_node = root_nodes[0]
    anc_seq = pro_ancestors[mrca_node]
    
    record = SeqRecord(
        Seq(anc_seq),
        id="Anc-tRNA-ProThr",
        description=f"Ancestral tRNA-Pro/Thr from {mrca_node}",
    )
    
    SeqIO.write([record], "results/Anc-tRNA-ProThr.fasta", "fasta")
    print(f"✓ Anc-tRNA-ProThr: {len(anc_seq)} nucleotides (from {mrca_node})")

if __name__ == "__main__":
    main()
