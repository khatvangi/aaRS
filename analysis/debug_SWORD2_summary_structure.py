"""
Inspect SWORD2_summary.json structures for ProRS runs.
Finds latest run per known prefix and logs top-level keys and list lengths.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import List, Optional

PREFIXES = [
    "deep_editing_pro_model_chainA_A",
    "deep_domain_pro_model_chainA_A",
    "modern_prours_pro_model_chainA_A",
]


def newest_run_dir(runs_dir: Path, prefix: str) -> Optional[Path]:
    matches = sorted([p for p in runs_dir.glob(f"{prefix}*") if p.is_dir()])
    if not matches:
        return None
    return matches[-1]


def inspect_summary(path: Path) -> str:
    data = json.loads(path.read_text())
    lines: List[str] = []
    lines.append(f"File: {path}")
    lines.append("Top-level keys: " + ", ".join(data.keys()))
    for key in ["domains", "PUs", "solutions"]:
        if key in data and isinstance(data[key], list):
            lines.append(f"{key}: len={len(data[key])}")
            if data[key]:
                lines.append(f"  example: {data[key][0]}")
    if "Optimal partition" in data:
        opt = data["Optimal partition"]
        lines.append("Optimal partition keys: " + ", ".join(opt.keys()))
    return "\n".join(lines)


def main():
    runs_dir = Path("sword2_results/runs")
    out_path = Path("summary/proRS_SWORD2_debug.txt")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    logs: List[str] = []
    for prefix in PREFIXES:
        latest = newest_run_dir(runs_dir, prefix)
        if latest is None:
            logs.append(f"{prefix}: no run dir found")
            continue
        summary_path = latest / "SWORD2_summary.json"
        if not summary_path.exists():
            logs.append(f"{prefix}: missing SWORD2_summary.json in {latest}")
            continue
        logs.append(inspect_summary(summary_path))
        logs.append("")
    out_path.write_text("\n".join(logs))
    print("\n".join(logs))


if __name__ == "__main__":
    main()
