#!/usr/bin/env python3
"""Parse Pfam results into comparative table."""

import re
from collections import defaultdict

def parse_pfam(tblfile):
    """Parse domain table output."""
    domains = defaultdict(list)
    
    with open(tblfile) as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.split()
            if len(fields) < 23:
                continue
            
            domain = fields[0]
            query = fields[3]
            evalue = float(fields[6])
            score = float(fields[7])
            env_from = int(fields[19])
            env_to = int(fields[20])
            
            if evalue < 0.01:
                domains[query].append({
                    'domain': domain,
                    'evalue': evalue,
                    'score': score,
                    'start': env_from,
                    'end': env_to
                })
    
    return domains

# Parse
domains = parse_pfam('complete_domains.tbl')

print("="*80)
print("COMPARATIVE PFAM DOMAIN ANALYSIS")
print("="*80)

# Define sequences
seqs = [
    'Deep_LUCA_ProRS',
    'Deep_LUCA_ThrRS', 
    'Shallow_Euk_ProThrRS',
    'Modern_Human_ProRS',
    'Modern_Human_ThrRS'
]

# Key domains of interest
key_domains = [
    'tRNA-synt_1c',      # Class I catalytic
    'tRNA-synt_2b',      # Class II catalytic
    'tRNA_edit',         # Editing domain
    'ProRS-C_1',         # ProRS C-terminal
    'HGTP_anticodon',    # Anticodon binding
    'TGS',               # TGS domain
    'tRNA_SAD',          # Thr/Ala-specific
]

print("\n" + "="*80)
print("DOMAIN ARCHITECTURE COMPARISON")
print("="*80)

for seq in seqs:
    if seq not in domains:
        print(f"\n{seq}: No significant domains found")
        continue
    
    print(f"\n{seq}:")
    print("-"*80)
    
    # Sort by position
    sorted_domains = sorted(domains[seq], key=lambda x: x['start'])
    
    # Group by domain type
    domain_types = defaultdict(list)
    for d in sorted_domains:
        domain_types[d['domain']].append(d)
    
    for domain_name in key_domains:
        if domain_name in domain_types:
            hits = domain_types[domain_name]
            if len(hits) == 1:
                d = hits[0]
                print(f"  ✓ {domain_name:<20} {d['start']:>4}-{d['end']:<4} (E={d['evalue']:.1e})")
            else:
                print(f"  ✓ {domain_name:<20} × {len(hits)} copies:")
                for d in hits:
                    print(f"     {d['start']:>4}-{d['end']:<4} (E={d['evalue']:.1e})")
    
    # Show other significant domains
    other_domains = [d for d in domain_types.keys() if d not in key_domains]
    if other_domains:
        print(f"  • Other domains: {', '.join(other_domains[:5])}")

print("\n" + "="*80)
print("KEY FINDINGS")
print("="*80)

print("""
CRITICAL OBSERVATIONS:

1. EDITING DOMAIN (tRNA_edit) PRESENCE:
   - Check if present in LUCA vs Modern
   - Check if lost in Shallow (fusion event)

2. CLASS I vs CLASS II:
   - LUCA ProRS: Should have tRNA-synt_1c
   - LUCA ThrRS: Should have tRNA-synt_2b
   - Shallow: Might have both (hybrid)

3. DOMAIN ORGANIZATION:
   - Compare N→C terminal arrangement
   - Identify fusion/loss events
   - Map to your AF3 results
""")

print("="*80)
