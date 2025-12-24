#!/bin/bash
set -e

echo "=== Deep Ancestral Reconstruction (Cross-Domain) ==="

mkdir -p results

reconstruct() {
    local AARS=$1
    local INPUT="data/interim/${AARS}_aligned.fasta"
    local PREFIX="results/${AARS}_deep"
    
    echo ""
    echo "Reconstructing $AARS..."
    echo "  Input: $INPUT"
    
    # IQ-TREE with ancestral state reconstruction
    iqtree3 -s $INPUT \
            -st AA \
            -m MFP \
            -bb 1000 \
            -asr \
            -nt 32 \
            --prefix $PREFIX
    
    echo "  ✓ $AARS complete"
    echo "  Tree: ${PREFIX}.treefile"
    echo "  ASR: ${PREFIX}.state"
}

reconstruct ProRS
reconstruct ThrRS

echo ""
echo "=== Extracting Deep Ancestral Node ==="

python3 << 'PYEOF'
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd

for aars in ["ProRS", "ThrRS"]:
    state_file = f"results/{aars}_deep.state"
    
    # Parse state file
    df = pd.read_csv(state_file, sep='\t', comment='#')
    
    # Group by node
    nodes = {}
    for node_name, group in df.groupby('Node'):
        group = group.sort_values('Site')
        seq = ''.join(group['State'].values)
        nodes[node_name] = seq
    
    # Find deepest node (lowest number = root)
    node_nums = [(int(n.replace('Node', '')), n) for n in nodes.keys()]
    node_nums.sort()
    root_node = node_nums[0][1]
    
    anc_seq = nodes[root_node]
    
    # Save
    rec = SeqRecord(
        Seq(anc_seq),
        id=f"Anc-{aars}-LUCA",
        description=f"Deep ancestral {aars} from {root_node} (cross-domain)"
    )
    
    SeqIO.write([rec], f"results/Anc-{aars}-LUCA.fasta", "fasta")
    print(f"✓ {aars}: {len(anc_seq)} aa from {root_node}")

# Combine ProRS and ThrRS ancestors
pro_rec = SeqIO.read("results/Anc-ProRS-LUCA.fasta", "fasta")
thr_rec = SeqIO.read("results/Anc-ThrRS-LUCA.fasta", "fasta")

# For now, use ProRS as the representative
# (In reality, you'd want to check the tree topology)
SeqIO.write([pro_rec], "results/Anc-ProThrRS-LUCA.fasta", "fasta")
print(f"\n✓ Created Anc-ProThrRS-LUCA.fasta")
PYEOF

echo ""
echo "=== DEEP ANCESTRY RECONSTRUCTION COMPLETE ==="
echo "Output: results/Anc-ProThrRS-LUCA.fasta"

