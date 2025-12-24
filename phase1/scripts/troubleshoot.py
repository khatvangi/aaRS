#!/usr/bin/env python3
# troubleshoot.py
"""Diagnose common Phase 1 issues"""

import subprocess
import sys
from pathlib import Path

def check_tool(cmd, version_flag="--version"):
    """Check if a tool is installed"""
    try:
        result = subprocess.run([cmd, version_flag], capture_output=True, text=True, timeout=5)
        print(f"✓ {cmd}: installed")
        return True
    except:
        print(f"✗ {cmd}: NOT FOUND")
        return False

def check_database(path):
    """Check if database exists (try multiple locations)"""
    search_paths = [
        Path(path),
        Path("data") / path,
        Path("data") / Path(path).name
    ]
    
    for p in search_paths:
        if p.exists():
            print(f"✓ Database: {p}")
            return str(p)
    
    print(f"✗ Database missing: {path}")
    return None

def main():
    print("=== Phase 1 Troubleshooting ===\n")
    
    print("1. Required Tools")
    tools_ok = all([
        check_tool("python3"),
        check_tool("iqtree3", "-h"),
        check_tool("mafft", "--version"),
        check_tool("hmmscan", "-h"),
        check_tool("cmalign", "-h"),
        check_tool("parallel", "--version"),
    ])
    
    print("\n2. Databases")
    db_paths = {}
    db_paths['pfam'] = check_database("Pfam-A.hmm")
    db_paths['pfam_h3m'] = check_database("Pfam-A.hmm.h3m")
    db_paths['rfam'] = check_database("Rfam.cm")
    db_paths['rfam_i1m'] = check_database("Rfam.cm.i1m")
    
    db_ok = all(db_paths.values())
    
    # Write database config
    if db_ok:
        import json
        config = {
            "pfam_hmm": str(db_paths['pfam']),
            "rfam_cm": str(db_paths['rfam'])
        }
        with open("checkpoints/database_paths.json", "w") as f:
            json.dump(config, f, indent=2)
        print("\n✓ Database paths saved to checkpoints/database_paths.json")
    
    print("\n3. Directory Structure")
    dirs_ok = all([
        Path(d).is_dir() for d in ["data", "results", "logs", "scripts", "checkpoints"]
    ])
    
    if dirs_ok:
        print("✓ All directories present")
    else:
        print("✗ Missing directories - run 00_setup_environment.sh")
    
    print("\n4. Python Packages")
    try:
        import Bio
        import requests
        import ete3
        print("✓ All Python packages installed")
        py_ok = True
    except ImportError as e:
        print(f"✗ Missing Python package: {e}")
        py_ok = False
    
    print("\n" + "="*50)
    if tools_ok and db_ok and dirs_ok and py_ok:
        print("✓ System ready for Phase 1")
        sys.exit(0)
    else:
        print("✗ System NOT ready - fix errors above")
        sys.exit(1)

if __name__ == "__main__":
    main()
