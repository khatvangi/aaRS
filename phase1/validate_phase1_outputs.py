#!/usr/bin/env python3
"""Comprehensive validation of Phase 1 outputs"""

from Bio import SeqIO, AlignIO
from pathlib import Path
import sys

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'

def test_pass(msg):
    print(f"{bcolors.OKGREEN}✓ PASS{bcolors.ENDC}: {msg}")
    return True

def test_fail(msg):
    print(f"{bcolors.FAIL}✗ FAIL{bcolors.ENDC}: {msg}")
    return False

def test_warn(msg):
    print(f"{bcolors.WARNING}⚠ WARN{bcolors.ENDC}: {msg}")
    return True

def check_protein_sequence(seq, name):
    """Validate protein sequence composition"""
    valid_aa = set('ACDEFGHIKLMNPQRSTVWY')
    ambiguous = set('XBZ')
    gap = set('-.')
    
    seq_set = set(str(seq).upper())
    
    # Check for invalid characters
    invalid = seq_set - valid_aa - ambiguous - gap
    if invalid:
        return test_fail(f"{name}: Invalid amino acids found: {invalid}")
    
    # Check composition
    aa_count = sum(1 for x in str(seq) if x in valid_aa)
    ambig_count = sum(1 for x in str(seq) if x in ambiguous)
    gap_count = sum(1 for x in str(seq) if x in gap)
    
    if aa_count == 0:
        return test_fail(f"{name}: No valid amino acids found")
    
    # Calculate percentages
    total = len(seq)
    aa_pct = 100 * aa_count / total
    ambig_pct = 100 * ambig_count / total
    gap_pct = 100 * gap_count / total
    
    print(f"  {name}: {aa_count}/{total} valid aa ({aa_pct:.1f}%)")
    
    if ambig_pct > 10:
        test_warn(f"{name}: {ambig_pct:.1f}% ambiguous residues")
    
    if gap_pct > 50:
        return test_fail(f"{name}: {gap_pct:.1f}% gaps (too high)")
    
    return test_pass(f"{name}: Valid protein sequence")

def check_nucleotide_sequence(seq, name):
    """Validate nucleotide sequence composition"""
    valid_nt = set('ACGTU')
    ambiguous = set('NRYWSMKBDHV')
    gap = set('-.')
    
    seq_set = set(str(seq).upper())
    
    # Check for invalid characters
    invalid = seq_set - valid_nt - ambiguous - gap
    if invalid:
        return test_fail(f"{name}: Invalid nucleotides found: {invalid}")
    
    # Check composition
    nt_count = sum(1 for x in str(seq).upper() if x in valid_nt)
    ambig_count = sum(1 for x in str(seq).upper() if x in ambiguous)
    gap_count = sum(1 for x in str(seq) if x in gap)
    
    if nt_count == 0:
        return test_fail(f"{name}: No valid nucleotides found")
    
    total = len(seq)
    nt_pct = 100 * nt_count / total
    
    print(f"  {name}: {nt_count}/{total} valid nt ({nt_pct:.1f}%)")
    
    if ambig_count > 0:
        test_warn(f"{name}: {ambig_count} ambiguous nucleotides")
    
    return test_pass(f"{name}: Valid nucleotide sequence")

def validate_ancestral_sequences():
    """Validate the main Phase 1 outputs"""
    print(f"\n{bcolors.HEADER}{'='*70}{bcolors.ENDC}")
    print(f"{bcolors.HEADER}PHASE 1 DATA INTEGRITY CHECK{bcolors.ENDC}")
    print(f"{bcolors.HEADER}{'='*70}{bcolors.ENDC}\n")
    
    all_pass = True
    
    # Test 1: Ancestral aaRS
    print(f"{bcolors.BOLD}Test 1: Ancestral aaRS (Anc-ProThrRS){bcolors.ENDC}")
    aars_file = Path("results/Anc-ProThrRS.fasta")
    
    if not aars_file.exists():
        all_pass &= test_fail("File not found")
    else:
        try:
            rec = SeqIO.read(aars_file, "fasta")
            print(f"  ID: {rec.id}")
            print(f"  Length: {len(rec.seq)} aa")
            
            # Length check
            if len(rec.seq) < 200:
                all_pass &= test_fail(f"Sequence too short ({len(rec.seq)} aa)")
            elif len(rec.seq) > 3000:
                all_pass &= test_fail(f"Sequence too long ({len(rec.seq)} aa)")
            else:
                test_pass(f"Length reasonable ({len(rec.seq)} aa)")
            
            # Composition check
            all_pass &= check_protein_sequence(rec.seq, "Anc-ProThrRS")
            
            # Check for poly-X regions (signs of errors)
            seq_str = str(rec.seq)
            for aa in 'ACDEFGHIKLMNPQRSTVWY':
                if aa * 10 in seq_str:
                    all_pass &= test_warn(f"Poly-{aa} region detected (10+ consecutive)")
            
        except Exception as e:
            all_pass &= test_fail(f"Error reading file: {e}")
    
    # Test 2: Ancestral tRNA (gapped alignment)
    print(f"\n{bcolors.BOLD}Test 2: Ancestral tRNA (Anc-tRNA-ProThr){bcolors.ENDC}")
    trna_file = Path("results/Anc-tRNA-ProThr.fasta")
    
    if not trna_file.exists():
        all_pass &= test_fail("File not found")
    else:
        try:
            rec = SeqIO.read(trna_file, "fasta")
            print(f"  ID: {rec.id}")
            print(f"  Alignment length: {len(rec.seq)} positions")
            
            # Remove gaps
            ungapped = str(rec.seq).replace('-', '').replace('.', '')
            print(f"  Ungapped length: {len(ungapped)} nt")
            
            # tRNA length check
            if len(ungapped) < 60:
                all_pass &= test_fail(f"tRNA too short ({len(ungapped)} nt, expected ~76)")
            elif len(ungapped) > 100:
                all_pass &= test_fail(f"tRNA too long ({len(ungapped)} nt, expected ~76)")
            elif 70 <= len(ungapped) <= 90:
                test_pass(f"tRNA length correct ({len(ungapped)} nt)")
            else:
                test_warn(f"tRNA length unusual ({len(ungapped)} nt, typical is 76)")
            
            # Composition check
            all_pass &= check_nucleotide_sequence(ungapped, "Anc-tRNA-ProThr")
            
        except Exception as e:
            all_pass &= test_fail(f"Error reading file: {e}")
    
    # Test 3: Raw sequence collections
    print(f"\n{bcolors.BOLD}Test 3: Raw Sequence Collections{bcolors.ENDC}")
    
    for aars in ["ProRS", "ThrRS", "SerRS", "ValRS"]:
        raw_file = Path(f"data/raw/{aars}_filtered.fasta")
        if raw_file.exists():
            count = len(list(SeqIO.parse(raw_file, "fasta")))
            if count >= 10:
                test_pass(f"{aars}: {count} sequences")
            else:
                test_warn(f"{aars}: Only {count} sequences (low diversity)")
        else:
            all_pass &= test_fail(f"{aars}: Raw sequences missing")
    
    # Test 4: Alignments
    print(f"\n{bcolors.BOLD}Test 4: Multiple Sequence Alignments{bcolors.ENDC}")
    
    for aars in ["ProRS", "ThrRS", "SerRS", "ValRS"]:
        aln_file = Path(f"data/interim/{aars}_aligned.fasta")
        if aln_file.exists():
            try:
                aln = AlignIO.read(aln_file, "fasta")
                n_seqs = len(aln)
                aln_len = aln.get_alignment_length()
                
                # Check gap content
                gap_cols = 0
                for i in range(aln_len):
                    col = aln[:, i]
                    if col.count('-') == n_seqs:
                        gap_cols += 1
                
                gap_pct = 100 * gap_cols / aln_len
                
                if gap_pct > 80:
                    all_pass &= test_fail(f"{aars}: {gap_pct:.1f}% all-gap columns")
                else:
                    test_pass(f"{aars}: {n_seqs} seqs × {aln_len} cols ({gap_pct:.1f}% gaps)")
            except:
                all_pass &= test_fail(f"{aars}: Could not read alignment")
        else:
            all_pass &= test_fail(f"{aars}: Alignment missing")
    
    # Test 5: Phylogenetic trees
    print(f"\n{bcolors.BOLD}Test 5: Phylogenetic Trees{bcolors.ENDC}")
    
    for aars in ["ProRS", "ThrRS", "SerRS", "ValRS"]:
        tree_file = Path(f"results/{aars}.treefile")
        if tree_file.exists():
            with open(tree_file) as f:
                tree_str = f.read().strip()
            if tree_str.count('(') >= 5 and tree_str.count(';') == 1:
                test_pass(f"{aars}: Valid Newick tree")
            else:
                all_pass &= test_warn(f"{aars}: Tree format unusual")
        else:
            all_pass &= test_fail(f"{aars}: Tree file missing")
    
    # Test 6: ASR state files
    print(f"\n{bcolors.BOLD}Test 6: Ancestral State Reconstruction Files{bcolors.ENDC}")
    
    for aars in ["ProRS", "ThrRS", "SerRS", "ValRS"]:
        state_file = Path(f"results/{aars}.state")
        if state_file.exists():
            size_mb = state_file.stat().st_size / (1024 * 1024)
            if size_mb > 0.1:
                test_pass(f"{aars}: ASR file present ({size_mb:.1f} MB)")
            else:
                all_pass &= test_warn(f"{aars}: ASR file very small ({size_mb:.2f} MB)")
        else:
            all_pass &= test_fail(f"{aars}: ASR file missing")
    
    # Final summary
    print(f"\n{bcolors.HEADER}{'='*70}{bcolors.ENDC}")
    if all_pass:
        print(f"{bcolors.OKGREEN}{bcolors.BOLD}✓✓✓ ALL TESTS PASSED ✓✓✓{bcolors.ENDC}")
        print(f"{bcolors.OKGREEN}Phase 1 data integrity confirmed - ready for Phase 2{bcolors.ENDC}")
        return 0
    else:
        print(f"{bcolors.FAIL}{bcolors.BOLD}✗✗✗ SOME TESTS FAILED ✗✗✗{bcolors.ENDC}")
        print(f"{bcolors.FAIL}Review errors above before proceeding to Phase 2{bcolors.ENDC}")
        return 1

if __name__ == "__main__":
    sys.exit(validate_ancestral_sequences())
