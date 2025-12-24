#!/bin/bash
# Auto-run validation when ThrRS+Zn jobs complete

VALIDATION_JOBS=(
    "anc_thrrs_cat_zn_THR"
    "anc_thrrs_cat_zn_SER"
    "modern_thrrs_ecoli_zn_THR"
    "modern_thrrs_ecoli_zn_SER"
    "modern_thrrs_ecoli_THR_zinc"
)

echo "=========================================="
echo "Monitoring validation jobs..."
echo "=========================================="
echo ""

while true; do
    # Check if all validation jobs are complete
    ready=0
    total=${#VALIDATION_JOBS[@]}

    for job in "${VALIDATION_JOBS[@]}"; do
        if [ -f "jobs/$job/min.rst" ]; then
            ((ready++))
        fi
    done

    timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    echo "[$timestamp] Validation jobs ready: $ready/$total"

    # Check main batch progress
    completed=$(find jobs -name "FINAL_RESULTS_MMPBSA.dat" 2>/dev/null | wc -l)
    echo "[$timestamp] Main batch: $completed/182 complete"

    if [ $ready -eq $total ]; then
        echo ""
        echo "=========================================="
        echo "✓ All validation jobs ready!"
        echo "=========================================="
        echo ""
        echo "Running validation: Gasteiger vs AM1-BCC"
        echo ""

        /storage/kiran-stuff/blast_env/bin/python validate_gasteiger_vs_am1bcc.py

        exit_code=$?
        echo ""
        echo "=========================================="
        if [ $exit_code -eq 0 ]; then
            echo "✓ Validation completed successfully"
        else
            echo "✗ Validation failed (exit code: $exit_code)"
        fi
        echo "=========================================="

        break
    fi

    # Wait 2 minutes before checking again
    sleep 120
done
