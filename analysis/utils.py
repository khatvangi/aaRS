"""
Shared helpers for locating files and reading simple configs.
Only uses stdlib, pandas, numpy.
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def find_first_csv(directory: Path, pattern: str = "*.csv") -> Path:
    """Return the first CSV that matches the pattern, raises FileNotFoundError if none."""
    for path in sorted(directory.glob(pattern)):
        return path
    raise FileNotFoundError(f"No CSV matching {pattern} found under {directory}")


def ensure_dir(path: Path) -> None:
    """Create a directory if it does not exist."""
    path.mkdir(parents=True, exist_ok=True)


def read_regions_config(path: Path) -> Optional[Dict[str, List[Tuple[str, int, int]]]]:
    """
    Read a JSON defining residue regions.
    Expected format: {\"catalytic_site\": [[\"A\", 50, 120], ...], \"editing_site\": [...], ...}
    Returns None if the file is missing.
    """
    if not path:
        return None
    if not path.exists():
        logger.warning("Regions JSON not found at %s; skipping overlaps.", path)
        return None
    with path.open() as fh:
        data = json.load(fh)
    parsed: Dict[str, List[Tuple[str, int, int]]] = {}
    for key, entries in data.items():
        parsed[key] = []
        for entry in entries:
            if len(entry) != 3:
                logger.warning("Skipping malformed region entry for %s: %s", key, entry)
                continue
            chain, start, end = entry
            parsed[key].append((str(chain), int(start), int(end)))
    return parsed


@dataclass
class OverlapResult:
    n_overlapping: int
    indices: List[int]


def count_overlaps(
    pus: Sequence[Dict],
    regions: Sequence[Tuple[str, int, int]],
    pu_chain_key: str = "chain",
    pu_start_key: str = "start",
    pu_end_key: str = "end",
) -> OverlapResult:
    """
    Count how many Protein Units overlap any provided region.
    Each PU entry is expected to have chain/start/end keys (robust to missing keys).
    """
    overlapping_indices: List[int] = []
    for idx, pu in enumerate(pus):
        try:
            chain = str(pu.get(pu_chain_key, "A"))
            start = int(pu.get(pu_start_key))
            end = int(pu.get(pu_end_key))
        except Exception:
            continue
        for reg_chain, reg_start, reg_end in regions:
            if chain != reg_chain:
                continue
            # overlap if intervals intersect
            if max(start, reg_start) <= min(end, reg_end):
                overlapping_indices.append(idx)
                break
    return OverlapResult(n_overlapping=len(overlapping_indices), indices=overlapping_indices)
