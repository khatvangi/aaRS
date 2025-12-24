#!/usr/bin/env python3
import json
from pathlib import Path

BASE = Path("phase2/af3_ligand_inputs")

# We only touch these 6 core files, not the fixed_inputs/ or old clean ones
TARGET_FILES = [
    "proRS_LUCA_Pro.json",
    "proRS_LUCA_Thr.json",
    "proRS_modern_Pro.json",
    "proRS_modern_Thr.json",
    "proRS_editing_Pro.json",
    "proRS_editing_Thr.json",
]


def infer_ligand_code(path: Path) -> str:
    name = path.name
    if "_Pro" in name:
        return "PRO"
    if "_Thr" in name:
        return "THR"
    raise ValueError(f"Cannot infer ligand code (PRO/THR) from filename: {name}")


def clean_and_fix_json(path: Path, out_path: Path):
    print(f"[INFO] Processing {path}")
    with path.open() as f:
        data = json.load(f)

    # We assume the top-level is a single AF3 job dict:
    # { name, dialect, version, sequences: [...], ... }
    if not isinstance(data, dict) or "sequences" not in data:
        raise ValueError(f"Unexpected JSON structure in {path}; expected dict with 'sequences'.")

    # Remove any root-level junk keys that AF3 doesn't like
    for bad_key in ("small_molecules", "bondedAtomPairs", "userCCD"):
        if bad_key in data:
            print(f"  - Removing root key '{bad_key}'")
            data.pop(bad_key, None)

    # Ensure dialect/version set correctly
    data.setdefault("dialect", "alphafold3")
    data.setdefault("version", 1)

    sequences = data["sequences"]
    if not isinstance(sequences, list):
        raise ValueError(f"'sequences' is not a list in {path}")

    # Drop any existing ligand entries
    new_sequences = []
    for entry in sequences:
        if not isinstance(entry, dict):
            raise ValueError(f"Non-dict entry in 'sequences' for {path}")
        if "ligand" in entry:
            print("  - Dropping existing ligand entry")
            continue
        new_sequences.append(entry)

    ligand_code = infer_ligand_code(path)

    # Append correct ligand block
    ligand_entry = {
        "ligand": {
            "id": "L",          # chain ID, must be uppercase letter
            "ccdCodes": [ligand_code],
        }
    }
    print(f"  - Adding ligand block with ccdCodes=['{ligand_code}']")
    new_sequences.append(ligand_entry)

    data["sequences"] = new_sequences

    # Write to *_clean.json
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as f:
        json.dump(data, f, indent=2)

    print(f"[OK] Wrote cleaned JSON to {out_path}")


def main():
    for fname in TARGET_FILES:
        in_path = BASE / fname
        if not in_path.exists():
            print(f"[WARN] Missing {in_path}, skipping")
            continue
        out_path = BASE / in_path.with_suffix("").name  # strip .json
        out_path = out_path.with_name(out_path.name + "_clean.json")
        out_path = out_path.with_suffix(".json")
        clean_and_fix_json(in_path, out_path)


if __name__ == "__main__":
    main()

