#!/usr/bin/env python3
import json
import glob

jobs = ["deep_domain_pro", "deep_domain_thr"]

print("="*70)
print("COMPARISON: PRO vs THR BINDING (Deep/LUCA Ancestor)")
print("="*70)

for job in jobs:
    summary = glob.glob(f"outputs/{job}/*_summary_confidences.json")[0]
    with open(summary) as f:
        data = json.load(f)
    
    overall_iptm = data['iptm']
    protein_ligand = data['chain_pair_iptm'][0][2]  # A-C interaction
    
    ligand = "PRO (cognate)" if "pro" in job else "THR (test)"
    
    print(f"\n{job}:")
    print(f"  Ligand: {ligand}")
    print(f"  Overall ipTM: {overall_iptm:.3f}")
    print(f"  Protein-Ligand ipTM: {protein_ligand:.3f}")
    if protein_ligand > 0.6:
        print(f"  ✅ BINDING DETECTED")
    else:
        print(f"  ⚠️  WEAK/NO BINDING")

print("\n" + "="*70)
print("HYPOTHESIS TEST")
print("="*70)
delta = abs(data['chain_pair_iptm'][0][2] - 0.75)  # Compare to PRO
print(f"ΔipTM (THR vs PRO): Will calculate after THR analysis")
print("\nIf THR ipTM ≈ PRO ipTM → PROMISCUITY DETECTED ✅")
print("If THR ipTM << PRO ipTM → SPECIFIC BINDING (no promiscuity)")
print("="*70)
