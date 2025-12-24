#!/bin/bash
# Phase1_master_pipeline.sh
# Coordinates all Phase 1 steps with checkpointing

set -e

LOGDIR="logs"
CHECKPOINT="checkpoints/phase1_progress.txt"

mkdir -p $LOGDIR checkpoints

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a $LOGDIR/master.log
}

run_step() {
    local step_num=$1
    local step_name=$2
    local script=$3
    
    # Check if already completed
    if grep -q "^${step_num}$" $CHECKPOINT 2>/dev/null; then
        log "✓ Step $step_num ($step_name) already completed (skipping)"
        return 0
    fi
    
    log "=== Step $step_num: $step_name ==="
    
    # Run script with logging
    if [[ $script == *.sh ]]; then
        bash $script 2>&1 | tee $LOGDIR/step_${step_num}.log
    elif [[ $script == *.py ]]; then
        python3 $script 2>&1 | tee $LOGDIR/step_${step_num}.log
    fi
    
    # Check exit status
    if [ ${PIPESTATUS[0]} -eq 0 ]; then
        echo "$step_num" >> $CHECKPOINT
        log "✓ Step $step_num completed successfully"
    else
        log "✗ Step $step_num FAILED"
        exit 1
    fi
}

# Main pipeline
log "=== Phase 1: Ancestral Reconstruction Pipeline ==="
log "Hardware: 64 CPU cores, 2 RTX GPUs"
log "Start time: $(date)"

run_step "01" "Define targets" "scripts/01_define_targets.py"
run_step "02" "Collect aaRS sequences" "scripts/02_collect_aaRS_sequences.py"
run_step "03" "Extract catalytic domains" "scripts/03_extract_catalytic_domains.sh"
run_step "04" "Structural alignment" "scripts/04_structural_alignment.sh"
run_step "05" "Ancestral reconstruction (aaRS)" "scripts/05_ancestral_reconstruction.sh"
run_step "06" "Extract ancestral aaRS" "scripts/06_extract_ancestors.py"
run_step "07" "Collect tRNA sequences" "scripts/07_collect_tRNAs.py"
run_step "08" "Align tRNAs" "scripts/08_align_tRNAs.sh"
run_step "09" "Reconstruct tRNA ancestors" "scripts/09_reconstruct_tRNA_ancestors.sh"
run_step "10" "Extract ancestral tRNA" "scripts/10_extract_anc_tRNA.py"

log "=== Phase 1 Complete ==="
log "End time: $(date)"
log "Outputs:"
log "  - results/Anc-ProThrRS.fasta"
log "  - results/Anc-tRNA-ProThr.fasta"
log "  - results/*_tree.treefile (phylogenies)"
log "  - results/*_asr.state (all ancestral states)"

# Generate summary report
python3 - <<'EOF'
from Bio import SeqIO
import json

summary = {
    "phase": 1,
    "completed": True,
    "outputs": {}
}

# Check outputs
for fname in ["Anc-ProThrRS.fasta", "Anc-tRNA-ProThr.fasta"]:
    path = f"results/{fname}"
    try:
        rec = SeqIO.read(path, "fasta")
        summary["outputs"][fname] = {
            "length": len(rec.seq),
            "id": rec.id
        }
    except:
        summary["outputs"][fname] = "MISSING"

with open("checkpoints/phase1_summary.json", "w") as f:
    json.dump(summary, f, indent=2)

print("\n=== Summary ===")
print(json.dumps(summary, indent=2))
EOF
