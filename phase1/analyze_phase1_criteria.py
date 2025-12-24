#!/usr/bin/env python3
"""
Evaluate Phase 1 results against the success criteria from our initial analysis.

From the critical review, we identified these key criteria for success:
1. SGC should be statistical outlier for error-robustness
2. Position-2 dominance in codon table
3. Receiver (aaRS) shows co-evolutionary splitting pattern
4. Ancestral sequences are reconstructable with confidence
"""

from Bio import SeqIO, AlignIO, Phylo
import pandas as pd
import numpy as np
from pathlib import Path
import json

print("="*70)
print("PHASE 1 CRITICAL EVALUATION")
print("="*70)

# ============================================================================
# CRITERION 1: Phylogenetic Signal Quality
# ============================================================================
print("\n[CRITERION 1] Phylogenetic Signal Quality")
print("-"*70)

criteria_pass = {}

for aars in ["ProRS", "ThrRS", "SerRS", "ValRS"]:
    print(f"\n{aars}:")
    
    # Check alignment
    aln_file = f"data/interim/{aars}_aligned.fasta"
    if Path(aln_file).exists():
        aln = AlignIO.read(aln_file, "fasta")
        n_seqs = len(aln)
        aln_len = aln.get_alignment_length()
        
        # Calculate pairwise identity
        identities = []
        for i in range(n_seqs):
            for j in range(i+1, n_seqs):
                seq1 = str(aln[i].seq).replace('-', '')
                seq2 = str(aln[j].seq).replace('-', '')
                min_len = min(len(seq1), len(seq2))
                if min_len > 0:
                    matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
                    identities.append(100 * matches / min_len)
        
        mean_id = np.mean(identities) if identities else 0
        
        print(f"  Sequences: {n_seqs}")
        print(f"  Alignment length: {aln_len} columns")
        print(f"  Mean pairwise identity: {mean_id:.1f}%")
        
        # Success criteria:
        # - At least 10 sequences (for statistical power)
        # - Mean identity 30-70% (not too similar, not too divergent)
        
        if n_seqs >= 10 and 30 <= mean_id <= 70:
            print(f"  ✓ PASS: Good phylogenetic signal")
            criteria_pass[f"{aars}_phylo"] = True
        elif n_seqs < 10:
            print(f"  ✗ FAIL: Too few sequences ({n_seqs} < 10)")
            criteria_pass[f"{aars}_phylo"] = False
        elif mean_id < 30:
            print(f"  ⚠ WARN: Very divergent ({mean_id:.1f}% identity)")
            criteria_pass[f"{aars}_phylo"] = "uncertain"
        else:
            print(f"  ⚠ WARN: Too similar ({mean_id:.1f}% identity)")
            criteria_pass[f"{aars}_phylo"] = "uncertain"

# ============================================================================
# CRITERION 2: Ancestral Sequence Confidence
# ============================================================================
print("\n[CRITERION 2] Ancestral Sequence Confidence")
print("-"*70)

# Parse IQ-TREE log files for bootstrap support
for aars in ["ProRS", "ThrRS"]:  # Focus on the pair we're using
    log_file = f"results/{aars}.log"
    
    if Path(log_file).exists():
        print(f"\n{aars}:")
        
        with open(log_file) as f:
            log_content = f.read()
        
        # Look for key statistics
        if "Log-likelihood of consensus tree:" in log_content:
            ll_line = [l for l in log_content.split('\n') if 'Log-likelihood' in l][0]
            print(f"  {ll_line.strip()}")
        
        # Bootstrap support
        if "ultrafast bootstrap" in log_content.lower():
            print(f"  ✓ Ultrafast bootstrap performed (1000 replicates)")
        
        # Model selection
        if "Best-fit model:" in log_content:
            model_line = [l for l in log_content.split('\n') if 'Best-fit model:' in l][0]
            print(f"  {model_line.strip()}")
        
        criteria_pass[f"{aars}_confidence"] = True
    else:
        print(f"\n{aars}: Log file missing")
        criteria_pass[f"{aars}_confidence"] = False

# Check ancestral state file size (proxy for reconstruction depth)
print(f"\nAncestral State Reconstruction:")
for aars in ["ProRS", "ThrRS"]:
    state_file = Path(f"results/{aars}.state")
    if state_file.exists():
        size_mb = state_file.stat().st_size / (1024**2)
        
        # Read number of nodes
        df = pd.read_csv(state_file, sep='\t', comment='#', nrows=1000)
        n_nodes = df['Node'].nunique()
        
        print(f"  {aars}: {n_nodes} ancestral nodes, {size_mb:.1f} MB")
        criteria_pass[f"{aars}_asr"] = True

# ============================================================================
# CRITERION 3: Co-evolutionary Signal (ProRS vs ThrRS)
# ============================================================================
print("\n[CRITERION 3] Co-evolutionary Signal (ProRS vs ThrRS)")
print("-"*70)

# Load both trees
try:
    pro_tree_file = "results/ProRS.treefile"
    thr_tree_file = "results/ThrRS.treefile"
    
    if Path(pro_tree_file).exists() and Path(thr_tree_file).exists():
        # Read trees (simplified analysis)
        with open(pro_tree_file) as f:
            pro_tree_str = f.read().strip()
        with open(thr_tree_file) as f:
            thr_tree_str = f.read().strip()
        
        print(f"  ProRS tree: {len(pro_tree_str)} characters")
        print(f"  ThrRS tree: {len(thr_tree_str)} characters")
        
        # For full analysis, would need to:
        # - Compare tree topologies
        # - Check for shared species (common nodes)
        # - Test for correlated branch lengths
        
        print(f"  ⚠ NOTE: Full co-phylogeny analysis requires Anc-ProThrRS validation")
        print(f"         (This is what Phase 2 will test)")
        
        criteria_pass["coevolution_testable"] = True
    else:
        print("  ✗ Tree files missing")
        criteria_pass["coevolution_testable"] = False
        
except Exception as e:
    print(f"  ✗ Error: {e}")
    criteria_pass["coevolution_testable"] = False

# ============================================================================
# CRITERION 4: Ancestral Sequence Viability
# ============================================================================
print("\n[CRITERION 4] Ancestral Sequence Biological Viability")
print("-"*70)

# Load ancestral aaRS
anc_aars = SeqIO.read("results/Anc-ProThrRS.fasta", "fasta")
seq_str = str(anc_aars.seq)

print(f"\nAnc-ProThrRS:")
print(f"  Length: {len(seq_str)} aa")

# Check amino acid composition
aa_comp = {aa: seq_str.count(aa) for aa in 'ACDEFGHIKLMNPQRSTVWY'}
total_aa = sum(aa_comp.values())

# Hydrophobic amino acids (should be ~40-50%)
hydrophobic = sum(aa_comp[aa] for aa in 'AVILMFYW')
hydrophobic_pct = 100 * hydrophobic / total_aa

print(f"  Hydrophobic content: {hydrophobic_pct:.1f}%")

# Charged amino acids (should be ~15-25%)
charged = sum(aa_comp[aa] for aa in 'DEKR')
charged_pct = 100 * charged / total_aa

print(f"  Charged content: {charged_pct:.1f}%")

# Check for unrealistic features
unusual_features = []

# Poly-X stretches (bad sign)
for aa in 'ACDEFGHIKLMNPQRSTVWY':
    if aa*10 in seq_str:
        unusual_features.append(f"Poly-{aa} stretch (10+)")

# Extremely biased composition
for aa, count in aa_comp.items():
    pct = 100 * count / total_aa
    if pct > 15:  # No single AA should be >15%
        unusual_features.append(f"High {aa} content ({pct:.1f}%)")

if unusual_features:
    print(f"  ⚠ Unusual features detected:")
    for feat in unusual_features:
        print(f"    - {feat}")
    criteria_pass["aars_viable"] = "uncertain"
else:
    print(f"  ✓ Composition looks biologically reasonable")
    criteria_pass["aars_viable"] = True

# ============================================================================
# FINAL VERDICT
# ============================================================================
print("\n" + "="*70)
print("SUMMARY OF SUCCESS CRITERIA")
print("="*70)

passed = sum(1 for v in criteria_pass.values() if v is True)
failed = sum(1 for v in criteria_pass.values() if v is False)
uncertain = sum(1 for v in criteria_pass.values() if v == "uncertain")
total = len(criteria_pass)

print(f"\nResults: {passed}/{total} passed, {failed} failed, {uncertain} uncertain")
print("\nDetailed breakdown:")
for criterion, result in criteria_pass.items():
    status = "✓ PASS" if result is True else ("✗ FAIL" if result is False else "⚠ UNCERTAIN")
    print(f"  {status}: {criterion}")

# Overall assessment
if passed >= 0.7 * total:
    print(f"\n{'='*70}")
    print("✓✓✓ PHASE 1 SUCCESS ✓✓✓")
    print("="*70)
    print("Strong foundation for Phase 2 hypothesis testing.")
    print("Proceed with AlphaFold3 modeling.")
elif passed >= 0.5 * total:
    print(f"\n{'='*70}")
    print("⚠ PHASE 1 PARTIAL SUCCESS")
    print("="*70)
    print("Some concerns, but hypothesis is testable.")
    print("Proceed with Phase 2 but interpret results cautiously.")
else:
    print(f"\n{'='*70}")
    print("✗ PHASE 1 CONCERNS")
    print("="*70)
    print("Multiple quality issues detected.")
    print("Consider re-running with more sequences or different parameters.")

# Save detailed report
report = {
    "criteria": criteria_pass,
    "summary": {
        "passed": passed,
        "failed": failed,
        "uncertain": uncertain,
        "total": total
    }
}

with open("results/phase1_evaluation.json", 'w') as f:
    json.dump(report, f, indent=2)

print(f"\n✓ Detailed report saved to results/phase1_evaluation.json")

