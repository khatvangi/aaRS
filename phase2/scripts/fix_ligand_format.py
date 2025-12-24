#!/usr/bin/env python3
"""Fix ligand format: Use CCD codes only (no SMILES)."""

import json
from pathlib import Path

# Read existing JSONs and fix ligand format
json_dir = Path('inputs/af3_jsons')

for json_file in json_dir.glob('*.json'):
    print(f"Fixing: {json_file.name}")
    
    with open(json_file) as f:
        data = json.load(f)
    
    # Find and fix the ligand entry
    for seq in data['sequences']:
        if 'ligand' in seq:
            ligand = seq['ligand']
            
            # Keep only CCD code, remove SMILES
            if 'smiles' in ligand:
                del ligand['smiles']
                print(f"  Removed SMILES, kept CCD: {ligand['ccdCodes']}")
    
    # Write back
    with open(json_file, 'w') as f:
        json.dump(data, f, indent=2)
    
    print(f"  ✓ Fixed")

print("\n✓ All ligand formats corrected!")
print("Launch with: nohup bash run_dual_test.sh > logs/phase2_fixed.log 2>&1 &")
