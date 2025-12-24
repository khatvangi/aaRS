#!/usr/bin/env python3
"""
Quick inspection of AF3 output JSON files
"""

import json
from pathlib import Path
import sys

def inspect_json(json_path):
    """Print structure of JSON file"""
    print(f"\n{'='*80}")
    print(f"File: {json_path.name}")
    print('='*80)
    
    with open(json_path) as f:
        data = json.load(f)
    
    def print_dict(d, indent=0):
        """Recursively print dict structure"""
        for key, value in d.items():
            if isinstance(value, dict):
                print("  " * indent + f"{key}: {{dict with {len(value)} keys}}")
                if indent < 2:  # Only go 2 levels deep
                    print_dict(value, indent + 1)
            elif isinstance(value, list):
                print("  " * indent + f"{key}: [list with {len(value)} items]")
                if len(value) > 0 and indent < 1:
                    print("  " * (indent+1) + f"First item type: {type(value[0])}")
                    if isinstance(value[0], (int, float, str)):
                        print("  " * (indent+1) + f"Sample values: {value[:3]}")
            else:
                print("  " * indent + f"{key}: {value}")
    
    print_dict(data)

def main(base_dir):
    """Inspect all JSON files in a model directory"""
    base_dir = Path(base_dir)
    
    # Look at the first model
    first_model = list(base_dir.glob('fulllength_*'))[0]
    print(f"\nInspecting model: {first_model.name}")
    
    # Check sample-0
    sample0_dir = first_model / 'seed-1_sample-0'
    summary_json = list(sample0_dir.glob('*_summary_confidences.json'))[0]
    conf_json = list(sample0_dir.glob('*_confidences.json'))[0]
    
    inspect_json(summary_json)
    inspect_json(conf_json)

if __name__ == '__main__':
    if len(sys.argv) > 1:
        base_dir = sys.argv[1]
    else:
        base_dir = '/storage/kiran-stuff/aaRS/phase2/af3_output_full'
    
    main(base_dir)
