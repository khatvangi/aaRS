#!/usr/bin/env python3
"""Extract AF3 metrics from outputs using the manifest."""
import argparse
import logging
import os
import re
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd

PROJECT_ROOT = Path("/storage/kiran-stuff/aaRS/phase2")
RESULTS_DIR = PROJECT_ROOT / "results"
MANIFEST = RESULTS_DIR / "manifest_runs.csv"


def setup_logger(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format="%(levelname)s: %(message)s")


def parse_input_base(base: str) -> Dict[str, str]:
    base = base or ""
    info = {
        "run_type": "standard",
        "era": "unknown",
        "enzyme": "",
        "ligand": "",
        "ligand_1": "",
        "ligand_2": "",
    }
    if base.startswith("COMPETITION"):
        info["run_type"] = "competition"
        rest = base.replace("COMPETITION_", "")
        if rest.startswith("FULL_"):
            rest = rest.replace("FULL_", "", 1)
        if "_vs_" in rest:
            left, right = rest.split("_vs_", 1)
        else:
            left, right = rest, ""
        left_parts = left.split("_")
        if left_parts:
            info["era"] = left_parts[0].lower() if left_parts[0] in ("anc", "modern") else "unknown"
        if len(left_parts) > 1:
            info["enzyme"] = left_parts[1]
        if len(left_parts) > 2:
            info["ligand_1"] = left_parts[2]
        if right:
            info["ligand_2"] = right.split("_")[0]
    else:
        toks = base.split("_")
        if toks:
            info["era"] = toks[0].lower() if toks[0] in ("anc", "modern") else "unknown"
        if len(toks) > 1:
            info["enzyme"] = toks[1]
        if toks:
            info["ligand"] = toks[-1]
    return info


def extract_metrics(manifest_path: Path, out_path: Path, verbose: bool = False) -> None:
    setup_logger(verbose)
    man = pd.read_csv(manifest_path)
    rows: List[Dict] = []
    for _, row in man.iterrows():
        input_json = row["input_json"]
        job = row["input_base"]
        output_dir = row["output_dir"]
        status = row["status"]
        match_conf = row["match_confidence"]
        if status not in ("ok", "collapsed_output") or not isinstance(output_dir, str) or not output_dir:
            continue
        info = parse_input_base(os.path.basename(input_json).replace(".json", ""))
        out_dir = Path(output_dir)
        struct_files = list(out_dir.rglob("seed-*_sample-*.cif"))
        if not struct_files:
            struct_files = list(out_dir.rglob("*_model.cif"))
        rank_df = None
        rank_path = next(iter(out_dir.rglob("ranking_scores.csv")), None)
        if rank_path:
            try:
                rank_df = pd.read_csv(rank_path)
            except Exception as e:
                logging.warning("Failed to read ranking_scores %s: %s", rank_path, e)
        def score_for(seed, sample):
            if rank_df is None or not {"seed","sample","ranking_score"}.issubset(rank_df.columns):
                return None
            hit = rank_df[(rank_df["seed"]==seed) & (rank_df["sample"]==sample)]
            if not hit.empty:
                return float(hit["ranking_score"].iloc[0])
            return None
        for struct in struct_files:
            seed = -1
            sample = -1
            m = re.search(r"seed-(\d+)_sample-(\d+)", struct.name)
            if m:
                seed = int(m.group(1))
                sample = int(m.group(2))
            replicate_id = str(sample if sample >=0 else struct.name)
            row_dict = {
                "input_json": input_json,
                "job_name": job,
                "output_dir": output_dir,
                "structure_file": struct.as_posix(),
                "era": info.get("era"),
                "enzyme": info.get("enzyme"),
                "run_type": info.get("run_type"),
                "ligand_code": info.get("ligand"),
                "ligand_1": info.get("ligand_1"),
                "ligand_2": info.get("ligand_2"),
                "seed": seed,
                "sample": sample,
                "replicate_id": replicate_id,
                "ranking_score": score_for(seed, sample),
                "match_confidence": match_conf,
                "status": status,
            }
            rows.append(row_dict)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out_path, index=False)
    logging.info("Wrote metrics to %s", out_path)


def main() -> None:
    ap = argparse.ArgumentParser(description="Extract AF3 metrics from outputs.")
    ap.add_argument("-m", "--manifest", default=str(MANIFEST), help="Path to manifest_runs.csv")
    ap.add_argument("-o", "--out", default=str(RESULTS_DIR / "metrics_long.csv"), help="Output CSV path")
    ap.add_argument("-v", "--verbose", action="store_true", help="Verbose logging")
    args = ap.parse_args()
    extract_metrics(Path(args.manifest), Path(args.out), verbose=args.verbose)


if __name__ == "__main__":
    main()
