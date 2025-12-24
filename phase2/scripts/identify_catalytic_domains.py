#!/usr/bin/env python3
"""
Identify catalytic domain boundaries for ProRS/ThrRS.
Based on Pfam/InterPro annotations for Class II aaRS.
"""

def read_fasta_clean(fasta_path):
    with open(fasta_path) as f:
        lines = f.readlines()
    seq = ''.join(line.strip() for line in lines if not line.startswith('>'))
    # Remove gaps
    import re
    return re.sub(r'[-.\*X]', '', seq)

# Load sequences
shallow = read_fasta_clean('../phase1/results/Anc-ProThrRS.fasta')
deep = read_fasta_clean('../phase1b/results/Anc-ProThrRS-LUCA.fasta')

print("="*60)
print("ANCESTRAL PROTEIN SEQUENCES")
print("="*60)
print(f"Shallow ancestor: {len(shallow)} aa")
print(f"Deep ancestor:    {len(deep)} aa")

print("\n" + "="*60)
print("CLASS II aaRS DOMAIN ARCHITECTURE (TYPICAL)")
print("="*60)
print("ProRS/ThrRS are Class IIa aminoacyl-tRNA synthetases")
print()
print("Typical domain organization:")
print("  1. N-terminal domain      (~150-250 aa) - tRNA anticodon binding")
print("  2. Catalytic domain       (~400-500 aa) - AMINO ACID BINDING + ATP")
print("  3. C-terminal domain      (~100-200 aa) - tRNA acceptor stem")
print("  4. Editing domain (some)  (~100-150 aa) - Proofreading")
print()
print("For promiscuity testing, we need: CATALYTIC DOMAIN (domain 2)")
print("  - Contains motifs 1, 2, 3 (signature sequences)")
print("  - Binds amino acid + ATP")
print("  - This is where Pro/Thr discrimination occurs")

print("\n" + "="*60)
print("RECOMMENDED DOMAIN EXTRACTION")
print("="*60)
print("Strategy: Extract residues ~200-700 (centered on catalytic core)")
print("  - Avoids N-terminal anticodon domain")
print("  - Includes full catalytic domain + editing domain")
print("  - Size: ~500 aa (well within GPU memory)")
print()
print("Alternative (conservative): Residues 150-800 (~650 aa)")
print()

# Show first 100 aa to look for obvious domain boundaries
print("First 100 aa of shallow ancestor:")
print(shallow[:100])
print()
print("Region 200-300 (likely catalytic core start):")
print(shallow[200:300])

print("\n" + "="*60)
print("DECISION REQUIRED")
print("="*60)
print("Option A (Conservative): Extract aa 150-800  (~650 aa)")
print("Option B (Targeted):     Extract aa 200-700  (~500 aa)")
print("Option C (Minimal):      Extract aa 250-650  (~400 aa)")
print()
print("Recommendation: Start with Option B (200-700)")
print("  - Includes full catalytic + editing domains")
print("  - Safe margin for domain boundaries")
print("  - Size: ~500 aa = NO OOM issues")
