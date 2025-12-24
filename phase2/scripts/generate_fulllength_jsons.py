#!/usr/bin/env python3
import json, re
from pathlib import Path

def clean_fasta(path, is_rna=False):
    with open(path) as f:
        seq = ''.join(l.strip() for l in f if not l.startswith('>'))
    seq = re.sub(r'[-.\*X]', '', seq)
    if is_rna:
        seq = seq.replace('T', 'U')
    return seq

shallow = clean_fasta('../phase1/results/Anc-ProThrRS.fasta')
deep = clean_fasta('../phase1b/results/Anc-ProThrRS-LUCA.fasta')
rna = clean_fasta('../phase1/results/Anc-tRNA-ProThr.fasta', True)

print(f"Shallow: {len(shallow)} aa")
print(f"Deep: {len(deep)} aa")

ligands = {'PRO': 'PRO', 'THR': 'THR'}
tests = [
    {'name': 'fulllength_shallow_pro', 'protein': shallow, 'ligand': 'PRO'},
    {'name': 'fulllength_shallow_thr', 'protein': shallow, 'ligand': 'THR'},
    {'name': 'fulllength_deep_pro', 'protein': deep, 'ligand': 'PRO'},
    {'name': 'fulllength_deep_thr', 'protein': deep, 'ligand': 'THR'}
]

Path('inputs/af3_jsons_fulllength').mkdir(exist_ok=True)

for t in tests:
    data = {
        "name": t['name'],
        "dialect": "alphafold3",
        "version": 1,
        "sequences": [
            {"protein": {"id": ["A"], "sequence": t['protein']}},
            {"rna": {"id": ["B"], "sequence": rna}},
            {"ligand": {"id": ["C"], "ccdCodes": [ligands[t['ligand']]]}}
        ],
        "modelSeeds": [1]
    }
    
    path = Path(f"inputs/af3_jsons_fulllength/{t['name']}.json")
    with open(path, 'w') as f:
        json.dump(data, f, indent=2)
    print(f"✓ {t['name']}: {len(t['protein'])} aa")

print("\n⚠️  WARNING: These are 1908-2037 aa proteins")
print("Expected GPU memory: 20-25 GB (may OOM)")
print("Estimated time: 60-90 min per job")
