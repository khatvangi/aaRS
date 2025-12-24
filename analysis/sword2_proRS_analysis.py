"""
Parse SWORD2 summaries for ProRS models and emit a compact metrics table.

CLI:
python -m analysis.sword2_proRS_analysis \
    --runs-dir aaRS/sword2_results/runs \
    --regions-json aaRS/metadata/proRS_regions.json
"""
from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd

from .utils import count_overlaps, ensure_dir, read_regions_config

logging.basicConfig(level=logging.INFO, format="%(levelname)s:%(message)s")
logger = logging.getLogger(__name__)

STRUCTURE_PREFIXES = [
    "deep_editing_pro_model_chainA_A",
    "deep_domain_pro_model_chainA_A",
    "modern_prours_pro_model_chainA_A",
]


def newest_run_dir(runs_dir: Path, prefix: str) -> Optional[Path]:
    matches = sorted([p for p in runs_dir.glob(f"{prefix}*") if p.is_dir()])
    if not matches:
        return None
    return matches[-1]


def parse_sword2_summary(path: Path) -> Dict:
    with path.open() as fh:
        return json.load(fh)


def extract_metrics(summary: Dict, regions: Optional[Dict[str, List[tuple]]] = None) -> Dict:
    out: Dict[str, Optional[object]] = {}
    domains = summary.get("domains") or summary.get("Domains") or summary.get("Optimal partition", {}).get("Domains")
    pus = summary.get("PUs") or summary.get("Protein Units") or []
    solutions = summary.get("solutions") or summary.get("Solutions") or []

    # n_domains: try solutions first, then domains length
    if solutions:
        first_sol = solutions[0]
        out["n_domains"] = first_sol.get("n_domains") or first_sol.get("Nb. domains")
        out["n_solutions"] = len(solutions)
        out["best_solution_score"] = first_sol.get("score") or first_sol.get("Quality")
    else:
        out["n_domains"] = len(domains) if domains else None
        out["n_solutions"] = 0
        out["best_solution_score"] = None

    out["n_PUs"] = len(pus) if pus else None

    # Overlap counts
    if regions and pus:
        catalytic_regs = regions.get("catalytic_site", [])
        editing_regs = regions.get("editing_site", [])
        if catalytic_regs:
            out["n_PUs_overlapping_catalytic"] = count_overlaps(pus, catalytic_regs).n_overlapping
        if editing_regs:
            out["n_PUs_overlapping_editing"] = count_overlaps(pus, editing_regs).n_overlapping
    return out


def main():
    parser = argparse.ArgumentParser(description="SWORD2 ProRS analysis")
    parser.add_argument("--runs-dir", type=Path, required=True, help="Directory containing SWORD2 run subfolders")
    parser.add_argument("--regions-json", type=Path, default=None, help="Optional JSON defining residue regions")
    parser.add_argument("--output-csv", type=Path, default=Path("aaRS/sword2_results/proRS_sword2_summary.csv"))
    args = parser.parse_args()

    regions = read_regions_config(args.regions_json) if args.regions_json else None

    records: List[Dict] = []
    for prefix in STRUCTURE_PREFIXES:
        latest = newest_run_dir(args.runs_dir, prefix)
        if latest is None:
            logger.warning("No runs found for prefix %s", prefix)
            continue
        summary_path = latest / "SWORD2_summary.json"
        if not summary_path.exists():
            logger.warning("Missing SWORD2_summary.json in %s", latest)
            continue
        summary = parse_sword2_summary(summary_path)
        metrics = extract_metrics(summary, regions)
        metrics["structure_id"] = prefix
        metrics["run_dir"] = str(latest)
        records.append(metrics)
        logger.info("Parsed %s -> %s", prefix, summary_path)

    if not records:
        raise SystemExit("No SWORD2 summaries parsed.")

    df = pd.DataFrame.from_records(records)
    ensure_dir(args.output_csv.parent)
    df.to_csv(args.output_csv, index=False)
    logger.info("Wrote SWORD2 summary CSV: %s", args.output_csv)


if __name__ == "__main__":
    main()
