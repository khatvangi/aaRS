#!/bin/bash
# Monitor MM/GBSA batch progress

echo "=========================================="
echo "MM/GBSA Batch Monitor"
echo "=========================================="
echo ""

# Count sander processes
NSANDER=$(ps aux | grep -c "[s]ander.*min.in")
echo "Running minimizations: $NSANDER"

# Count completed jobs
NCOMPLETED=$(find jobs -name "FINAL_RESULTS_MMPBSA.dat" 2>/dev/null | wc -l)
NTOTAL=$(wc -l < ../mmgbsa_real/manifest_af3only.csv)
NTOTAL=$((NTOTAL - 1))  # Subtract header
echo "Completed MM/GBSA: $NCOMPLETED/$NTOTAL"

# Count jobs with minimization done
NMIN_DONE=$(find jobs -name "min.rst" 2>/dev/null | wc -l)
echo "Minimizations complete: $NMIN_DONE/$NTOTAL"

# Check for any error patterns
if [ -f "mmgbsa_results.csv" ]; then
    echo ""
    echo "Status summary:"
    tail -n +2 mmgbsa_results.csv | awk -F',' '{print $3}' | sort | uniq -c | sort -rn
fi

# Show latest batch log
echo ""
echo "Latest progress:"
tail -5 batch.log

echo ""
echo "=========================================="
