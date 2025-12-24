#!/bin/bash
# RUN_PHASE1.sh
# Complete Phase 1 execution

set -e

echo "=== Phase 1: Ancestral Reconstruction ==="
echo "Estimated time: 2-3 days"
echo "Hardware: 64 cores, ~100GB RAM required"
echo ""

# 1. Setup environment (run once)
if [ ! -d "venv" ]; then
    echo "Setting up environment (first time only)..."
    bash scripts/00_setup_environment.sh
fi

# Activate environment
source venv/bin/activate
export PATH=$PWD/bin:$PATH

# 2. Run master pipeline
bash scripts/Phase1_master_pipeline.sh

# 3. Quality control
echo ""
echo "=== Running Quality Control ==="
python3 scripts/12_quality_control.py

if [ $? -eq 0 ]; then
    echo ""
    echo "✓✓✓ Phase 1 COMPLETE ✓✓✓"
    echo ""
    echo "Outputs ready for Phase 2 (AlphaFold3):"
    echo "  - results/Anc-ProThrRS.fasta"
    echo "  - results/Anc-tRNA-ProThr.fasta"
    echo ""
    echo "Next step: Run Phase 2 (AlphaFold3 modeling)"
else
    echo ""
    echo "✗✗✗ Phase 1 FAILED - Check logs in logs/ ✗✗✗"
    exit 1
fi
```

---

## **14. Directory Structure After Completion**
```
phase1/
├── data/
│   ├── config.json
│   ├── species_list.txt
│   ├── raw/
│   │   ├── ProRS_filtered.fasta (40 species)
│   │   ├── ThrRS_filtered.fasta
│   │   ├── SerRS_filtered.fasta
│   │   ├── ValRS_filtered.fasta
│   │   ├── tRNA_Pro_all.fasta (150+ sequences)
│   │   ├── tRNA_Thr_all.fasta
│   │   ├── tRNA_Ser_all.fasta
│   │   └── tRNA_Val_all.fasta
│   └── interim/
│       ├── ProRS_catalytic.fasta (domains only)
│       ├── ProRS_aligned.fasta (MSA)
│       ├── ... (similar for all aaRS/tRNA)
├── results/
│   ├── Anc-ProThrRS.fasta ⭐ (MAIN OUTPUT)
│   ├── Anc-tRNA-ProThr.fasta ⭐ (MAIN OUTPUT)
│   ├── ProRS_tree.treefile
│   ├── ProRS_asr.state (all ancestral nodes)
│   ├── ... (similar for all)
├── logs/
│   ├── master.log
│   ├── step_01.log
│   └── ... (all step logs)
├── checkpoints/
│   ├── phase1_progress.txt
│   └── phase1_summary.json
├── scripts/
│   └── (all 13 scripts above)
├── bin/ (installed tools)
└── venv/ (Python environment)
