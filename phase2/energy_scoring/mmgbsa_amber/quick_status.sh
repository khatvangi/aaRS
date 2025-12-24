#!/bin/bash
# Quick status check

echo "=========================================="
echo "Quick Status Check - $(date '+%H:%M:%S')"
echo "=========================================="
echo ""

# Main batch progress
completed=$(find jobs -name "FINAL_RESULTS_MMPBSA.dat" 2>/dev/null | wc -l)
minimized=$(find jobs -name "min.rst" 2>/dev/null | wc -l)
running=$(ps aux | grep -c "[s]ander.*min.in")

echo "Main Batch:"
echo "  MM/GBSA complete:     $completed/182 ($(( completed * 100 / 182 ))%)"
echo "  Minimizations done:   $minimized/182"
echo "  Sander processes:     $running"
echo ""

# Validation jobs
echo "Validation Jobs (ThrRS+Zn):"
val_ready=0
for job in anc_thrrs_cat_zn_THR anc_thrrs_cat_zn_SER modern_thrrs_ecoli_zn_THR modern_thrrs_ecoli_zn_SER modern_thrrs_ecoli_THR_zinc; do
    if [ -f "jobs/$job/min.rst" ]; then
        echo "  ✓ $job"
        ((val_ready++))
    else
        echo "  ⏳ $job"
    fi
done
echo ""
echo "  Ready: $val_ready/5"
echo ""

if [ $val_ready -eq 5 ]; then
    echo "✓ All validation jobs ready - check validation_monitor.log"
else
    echo "⏳ Waiting for $(( 5 - val_ready )) more validation jobs..."
fi

echo ""
echo "=========================================="
