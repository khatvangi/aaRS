#!/usr/bin/env python3
"""Parse AlphaFold3 outputs and extract binding metrics"""

import json
from pathlib import Path
import pandas as pd

def parse_af3_output(output_dir):
    """Extract pDockQ, ipTM, pLDDT from AF3 output"""
    
    # AF3 output structure (adjust based on actual format):
    # - ranking_debug.json contains confidence scores
    # - result_model_*.pdb contains structures
    
    metrics = {
        "pDockQ": None,
        "ipTM": None,
        "pLDDT": None,
        "status": "not_run"
    }
    
    ranking_file = output_dir / "ranking_debug.json"
    
    if ranking_file.exists():
        try:
            with open(ranking_file) as f:
                data = json.load(f)
            
            # Extract metrics (adjust keys based on actual AF3 output)
            if "iptm" in data:
                metrics["ipTM"] = data["iptm"]
            if "plddt" in data:
                metrics["pLDDT"] = data["plddt"]
            # pDockQ might need to be calculated from structure
            
            metrics["status"] = "success"
        except Exception as e:
            metrics["status"] = f"error: {e}"
    else:
        metrics["status"] = "output_missing"
    
    return metrics

# Parse all results
results = []

for case_dir in Path("outputs").iterdir():
    if case_dir.is_dir():
        name = case_dir.name
        metrics = parse_af3_output(case_dir)
        
        results.append({
            "case": name,
            **metrics
        })

# Save results table
df = pd.DataFrame(results)
df.to_csv("results/af3_metrics.csv", index=False)

print("=== AlphaFold3 Results Summary ===\n")
print(df.to_string(index=False))

# Check hypothesis
print("\n=== Hypothesis Test ===")

ancestral_cases = df[df['case'].str.contains('ancestral')]
modern_cases = df[df['case'].str.contains('modern')]

if len(ancestral_cases) > 0 and all(ancestral_cases['status'] == 'success'):
    print("Ancestral promiscuity test:")
    for _, row in ancestral_cases.iterrows():
        binding = "BINDS" if row.get('ipTM', 0) > 0.5 else "NO BINDING"
        print(f"  {row['case']}: ipTM={row.get('ipTM', 'N/A'):.3f} → {binding}")
else:
    print("⚠ Ancestral cases not complete")

if len(modern_cases) > 0 and all(modern_cases['status'] == 'success'):
    print("\nModern specificity test:")
    for _, row in modern_cases.iterrows():
        binding = "BINDS" if row.get('ipTM', 0) > 0.5 else "NO BINDING"
        expected = "EXPECTED" if "cognate" in row['case'] else "UNEXPECTED"
        print(f"  {row['case']}: ipTM={row.get('ipTM', 'N/A'):.3f} → {binding} ({expected})")
else:
    print("⚠ Modern cases not complete")

df.to_csv("results/af3_metrics.csv", index=False)
print(f"\n✓ Results saved to results/af3_metrics.csv")
