"""
Join Poincaré and SWORD2 outputs into a single snapshot CSV.
Paths are relative to project root; adjust via CLI if needed.
"""
from __future__ import annotations

import argparse
import logging
from pathlib import Path

import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(levelname)s:%(message)s")
logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(description="aaRS geometry snapshot")
    parser.add_argument(
        "--poincare-merged",
        type=Path,
        default=Path("aaRS/poincare_runs/proRS/merged_proRS_poincare_metadata.csv"),
    )
    parser.add_argument(
        "--sword2-summary",
        type=Path,
        default=Path("aaRS/sword2_results/proRS_sword2_summary.csv"),
    )
    parser.add_argument(
        "--poincare-stats",
        type=Path,
        default=Path("aaRS/summary/proRS_poincare_stats.json"),
        help="Optional stats JSON from final Poincaré analysis",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("aaRS/summary/proRS_geometry_metrics.csv"),
    )
    parser.add_argument(
        "--snapshot",
        type=Path,
        default=Path("aaRS/summary/proRS_geometry_snapshot.json"),
        help="Optional distilled JSON snapshot",
    )
    args = parser.parse_args()

    df_poincare = pd.read_csv(args.poincare_merged)
    df_sword2 = pd.read_csv(args.sword2_summary)

    # simple outer join: keep as many fields as possible
    # since identifiers differ (sequence vs structure), just concatenate with keys
    df_poincare["source"] = "poincare"
    df_sword2["source"] = "sword2"
    combined = pd.concat([df_poincare, df_sword2], ignore_index=True, sort=False)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(args.output, index=False)
    logger.info("Wrote combined snapshot CSV: %s", args.output)

    # Optional distilled snapshot
    snapshot = {}
    if args.poincare_stats.exists():
        try:
            import json

            stats_json = json.loads(args.poincare_stats.read_text())
            snapshot["LUCA_radius"] = stats_json.get("LUCA", {}).get("radius_mean")
            by_k = stats_json.get("by_kingdom", {})
            for k in ["Bacteria", "Archaea", "Eukarya"]:
                if k in by_k:
                    snapshot[f"median_{k.lower()}_radius"] = by_k[k].get("radius_median")
            snapshot["LUCA_vs_non_LUCA_p"] = stats_json.get("luca_vs_non_luca", {}).get("p_value")
            edit = stats_json.get("editing", {})
            if isinstance(edit, dict):
                snapshot["editing_vs_no_editing_p"] = edit.get("p_value")
            corr = stats_json.get("radius_promiscuity_correlation", {})
                # avoid KeyError on missing
            snapshot["radius_promiscuity_spearman_r"] = corr.get("spearman_rho")
            snapshot["radius_promiscuity_spearman_p"] = corr.get("spearman_p")
        except Exception as exc:
            logger.warning("Could not read poincare stats JSON: %s", exc)

    if args.snapshot:
        # add sword2 per-structure basics
        try:
            import json

            s2 = df_sword2.to_dict(orient="records")
            sword2_map = {}
            def int_or_none(val):
                try:
                    if pd.isna(val):
                        return None
                    return int(val)
                except Exception:
                    return None
            for rec in s2:
                name = rec.get("structure_id") or rec.get("structure")
                if not name:
                    continue
                sword2_map[name] = {
                    "n_domains": int_or_none(rec.get("n_domains")),
                    "n_PUs": int_or_none(rec.get("n_PUs")),
                    "n_solutions": int_or_none(rec.get("n_solutions")),
                }
            snapshot["sword2"] = sword2_map
            args.snapshot.parent.mkdir(parents=True, exist_ok=True)
            args.snapshot.write_text(json.dumps(snapshot, indent=2))
            logger.info("Wrote distilled snapshot JSON: %s", args.snapshot)
        except Exception as exc:
            logger.warning("Could not write snapshot JSON: %s", exc)


if __name__ == "__main__":
    main()
