# AF3Score AlphaRank-like results

- Manifest: `/storage/kiran-stuff/aaRS/phase2/results/manifest_runs.csv`
- Metrics: `/storage/kiran-stuff/aaRS/phase2/results/metrics_long.csv`
- Competition summary: `/storage/kiran-stuff/aaRS/phase2/results/competition_summary.csv`
- Table: `/storage/kiran-stuff/aaRS/phase2/results/TABLE1_AUTO.md`
- Plot: `/storage/kiran-stuff/aaRS/phase2/results/FIG_competition_barplot.png`

Regenerate:
```bash
python scripts/alpharank_like/build_manifest.py
python scripts/alpharank_like/extract_af3_metrics.py
python scripts/alpharank_like/pocket_iou.py
python scripts/alpharank_like/summarize_competitions.py
```