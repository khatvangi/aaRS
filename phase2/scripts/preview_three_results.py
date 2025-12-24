#!/usr/bin/env python3
import json, glob

jobs = {
    'deep_pro': 'deep_domain_pro',
    'deep_thr': 'deep_domain_thr', 
    'shallow_pro': 'shallow_domain_pro'
}

print("="*70)
print("PREVIEW: 3/4 JOBS COMPLETE")
print("="*70)

for key, job in jobs.items():
    try:
        summary = glob.glob(f"outputs/{job}/*_summary_confidences.json")[0]
        with open(summary) as f:
            data = json.load(f)
        iptm = data['chain_pair_iptm'][0][2]
        ancestor = "Deep" if "deep" in key else "Shallow"
        ligand = "PRO" if "pro" in key else "THR"
        print(f"{ancestor:8} + {ligand}: {iptm:.3f}")
    except:
        print(f"{key}: Not found")

print("\nWaiting for: shallow_domain_thr...")
print("="*70)
