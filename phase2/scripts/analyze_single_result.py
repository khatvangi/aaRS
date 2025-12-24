#!/usr/bin/env python3
"""Quick analysis of AF3 output to extract confidence metrics."""

import json
import glob
from pathlib import Path

def analyze_output(output_dir):
    """Analyze a single AF3 output directory."""
    
    print("="*70)
    print(f"ANALYZING: {output_dir}")
    print("="*70)
    
    # Find the summary confidence file
    summary_files = glob.glob(f"{output_dir}/*_summary_confidences.json")
    
    if not summary_files:
        print("❌ No summary confidence file found")
        return
    
    summary_file = summary_files[0]
    print(f"\n📄 Reading: {Path(summary_file).name}")
    
    with open(summary_file) as f:
        data = json.load(f)
    
    # Extract key metrics
    print("\n" + "="*70)
    print("KEY CONFIDENCE METRICS")
    print("="*70)
    
    # Overall metrics
    if 'iptm' in data:
        iptm = data['iptm']
        print(f"\n🎯 ipTM (Interface pTM):  {iptm:.3f}")
        if iptm > 0.5:
            print("   ✅ GOOD - Complex likely forms (ipTM > 0.5)")
        else:
            print("   ⚠️  LOW - Weak or no binding predicted")
    
    if 'ptm' in data:
        ptm = data['ptm']
        print(f"\n📊 pTM (predicted TM-score): {ptm:.3f}")
        if ptm > 0.5:
            print("   ✅ GOOD - High confidence prediction")
        else:
            print("   ⚠️  MEDIUM - Moderate confidence")
    
    if 'ranking_score' in data:
        ranking = data['ranking_score']
        print(f"\n⭐ Ranking Score: {ranking:.3f}")
    
    # Chain-level metrics (handle both dict and list formats)
    if 'chain_iptm' in data:
        print("\n" + "-"*70)
        print("CHAIN-LEVEL ipTM (Interaction Confidence)")
        print("-"*70)
        chain_iptm = data['chain_iptm']
        if isinstance(chain_iptm, dict):
            for chain_pair, score in chain_iptm.items():
                print(f"  {chain_pair}: {score:.3f}")
        elif isinstance(chain_iptm, list):
            print(f"  Chain ipTM matrix: {chain_iptm}")
    
    if 'chain_ptm' in data:
        print("\n" + "-"*70)
        print("CHAIN-LEVEL pTM (Structure Confidence)")
        print("-"*70)
        chain_ptm = data['chain_ptm']
        if isinstance(chain_ptm, dict):
            for chain, score in chain_ptm.items():
                print(f"  Chain {chain}: {score:.3f}")
        elif isinstance(chain_ptm, list):
            print(f"  Chain pTM: {chain_ptm}")
    
    # Check for multiple samples
    sample_dirs = glob.glob(f"{output_dir}/seed-*_sample-*")
    if sample_dirs:
        print("\n" + "="*70)
        print(f"MULTIPLE SAMPLES: {len(sample_dirs)} predictions generated")
        print("="*70)
        
        sample_scores = []
        for sample_dir in sorted(sample_dirs):
            sample_conf = glob.glob(f"{sample_dir}/*_summary_confidences.json")
            if sample_conf:
                with open(sample_conf[0]) as f:
                    sample_data = json.load(f)
                    sample_iptm = sample_data.get('iptm', 0)
                    sample_ptm = sample_data.get('ptm', 0)
                    sample_scores.append((Path(sample_dir).name, sample_iptm, sample_ptm))
        
        if sample_scores:
            print("\nConfidence scores across samples:")
            for sample_name, iptm, ptm in sample_scores:
                status = "✅" if iptm > 0.5 else "⚠️"
                print(f"  {status} {sample_name}: ipTM = {iptm:.3f}, pTM = {ptm:.3f}")
            
            best_sample = max(sample_scores, key=lambda x: x[1])
            print(f"\n🏆 Best ipTM: {best_sample[0]} (ipTM = {best_sample[1]:.3f})")
    
    # Show full JSON for debugging
    print("\n" + "="*70)
    print("FULL CONFIDENCE DATA (for debugging)")
    print("="*70)
    print(json.dumps(data, indent=2))
    
    print("\n" + "="*70)
    print("⚠️  ISSUE DETECTED")
    print("="*70)
    print("\nControl test (PRO cognate) shows LOW ipTM (0.310)")
    print("This suggests one of:")
    print("  1. Catalytic domain alone insufficient (needs other domains)")
    print("  2. tRNA not properly positioned near binding site")
    print("  3. Ligand (PRO) too small for stable complex prediction")
    print("  4. Domain boundaries need adjustment")
    print("\nRecommendation: Wait for all 4 jobs to compare patterns")
    print("="*70)

# Analyze the completed output
analyze_output("outputs/deep_domain_pro")
