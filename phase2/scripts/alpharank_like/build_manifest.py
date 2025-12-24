#!/usr/bin/env python3
"""Build a manifest mapping AF3 input JSONs to output directories.

Outputs: results/manifest_runs.csv
Columns include input_json, input_base, output_dir, match_confidence, status,
and presence flags for key output files.
"""
import argparse
import csv
import json
import logging
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

PROJECT_ROOT = Path("/storage/kiran-stuff/aaRS/phase2")
INPUT_DIRS = [
    PROJECT_ROOT / "af3_inputs_modern_allAAs",
    PROJECT_ROOT / "af3_inputs_ancient",
    PROJECT_ROOT / "af3_modern_matrix",
    PROJECT_ROOT / "batch_cloud_fulllength",
]
OUTPUT_ROOTS = [
    PROJECT_ROOT / "af3_output_full",
    PROJECT_ROOT / "af3_modern_matrix",
    PROJECT_ROOT / "batch_cloud_fulllength",
    PROJECT_ROOT / "batch_local_domains",
    PROJECT_ROOT / "af3_gaps",
]
RESULTS_DIR = PROJECT_ROOT / "results"


def setup_logger(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format="%(levelname)s: %(message)s")


def find_inputs() -> List[Path]:
    inputs: List[Path] = []
    for d in INPUT_DIRS:
        if not d.exists():
            continue
        for path in d.rglob("*.json"):
            name = path.name
            parts = path.as_posix().split("/")
            if any(tok.endswith("_Results") or tok.startswith("seed-") for tok in parts):
                continue
            if name.endswith(("_data.json", "_confidences.json", "summary_confidences.json", "confidences.json")):
                continue
            inputs.append(path)
    return sorted(inputs)


@dataclass
class OutputRun:
    path: Path
    job_name: str
    files: List[str] = field(default_factory=list)
    flags: Dict[str, bool] = field(default_factory=dict)
    match_confidence: str = "exact"
    status: str = "ok"
    n_cif_replicates: int = 0


def choose_timestamp_dir(base_dir: Path) -> Tuple[Path, str]:
    subdirs = [p for p in base_dir.iterdir() if p.is_dir()]
    ts_dirs = [p for p in subdirs if p.name.isdigit() or p.name.replace("-", "").isdigit()]
    if ts_dirs:
        ts = sorted(ts_dirs, key=lambda p: p.name)[-1]
        return ts, "timestamp_exact"
    return base_dir, "exact"


def scan_outputs() -> List[OutputRun]:
    runs: List[OutputRun] = []
    for root in OUTPUT_ROOTS:
        if not root.exists():
            continue
        for dirpath, _, filenames in os.walk(root):
            if filenames:
                runs.append(OutputRun(Path(dirpath), Path(dirpath).name, filenames, {}))
    return runs


def find_output_for_input(input_base: str) -> Tuple[Optional[OutputRun], str]:
    for root in OUTPUT_ROOTS:
        base_dir = root / input_base
        if base_dir.exists() and base_dir.is_dir():
            chosen, conf = choose_timestamp_dir(base_dir)
            return OutputRun(path=chosen, job_name=input_base, match_confidence=conf), conf
    # fallback: scan for any dir named input_base under roots
    for root in OUTPUT_ROOTS:
        matches = list(root.rglob(input_base))
        if matches:
            return OutputRun(path=matches[-1], job_name=input_base, match_confidence="fallback"), "fallback"
    return None, "none"


def write_manifest(rows: List[Dict[str, str]], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0].keys()) if rows else [
        "input_json",
        "input_base",
        "output_dir",
        "match_confidence",
        "status",
        "job_name_out",
        "n_cif_replicates",
        "has_ranking_scores",
        "has_summary_confidences",
        "model_files",
    ]
    with out_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def build_manifest(verbose: bool = False) -> None:
    setup_logger(verbose)
    inputs = find_inputs()
    logging.info("Found %d input JSONs", len(inputs))

    rows: List[Dict[str, str]] = []
    for inp in inputs:
        base = inp.stem
        match, conf = find_output_for_input(base)
        if match is None or not match.path.exists():
            rows.append({
                "input_json": inp.as_posix(),
                "input_base": base,
                "output_dir": "",
                "match_confidence": conf,
                "status": "missing_output",
                "job_name_out": "",
                "n_cif_replicates": 0,
                "has_ranking_scores": 0,
                "has_summary_confidences": 0,
                "model_files": "",
            })
            continue
        out_dir = match.path
        rep_files = list(out_dir.rglob("seed-*_sample-*.cif"))
        has_rank = any(out_dir.rglob("ranking_scores.csv"))
        has_summary = any(out_dir.rglob("*summary_confidences.json"))
        model_files = rep_files
        status = "ok"
        match_conf = conf
        if not rep_files:
            model_files = list(out_dir.rglob("*_model.cif"))
            if model_files:
                status = "collapsed_output"
                rep_files = model_files
                match_conf = match_conf or "fallback"
        rows.append({
                "input_json": inp.as_posix(),
                "input_base": base,
                "output_dir": out_dir.as_posix(),
                "match_confidence": match_conf,
                "status": status,
                "job_name_out": base,
                "n_cif_replicates": len(rep_files),
                "has_ranking_scores": 1 if has_rank else 0,
                "has_summary_confidences": 1 if has_summary else 0,
                "model_files": ";".join(sorted([p.name for p in model_files])),
        })

    manifest_path = RESULTS_DIR / "manifest_runs.csv"
    write_manifest(rows, manifest_path)
    logging.info("Wrote manifest to %s", manifest_path)


def main() -> None:
    ap = argparse.ArgumentParser(description="Build AF3Score run manifest.")
    ap.add_argument("-v", "--verbose", action="store_true", help="Verbose logging")
    args = ap.parse_args()
    build_manifest(verbose=args.verbose)


if __name__ == "__main__":
    main()
