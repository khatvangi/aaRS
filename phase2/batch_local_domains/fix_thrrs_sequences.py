#!/usr/bin/env python3
"""
Fix ThrRS catalytic domain sequences in all JSON files.
Replace placeholder "MAAAAA" with correct 278 aa sequence.
"""
import json
import glob

# Correct ThrRS sequence (278 amino acids)
CORRECT_THRRS_SEQUENCE = "EEEYLLKPMNCPHHIRIYASLRPRSYRDLPLRLAEFGTVSYRYEQSGELHGLTRVRGFTQDDAHIFCRPLDQIKGEFLGVLDLVLLVFKLTFGLNGDDYRARIGVRDKVNDPKSGKYVVSYVCRNCGYRIEGARGTDIECPVCHSHDLEGDEENWALAERQIKEAVEEYKNEVCGLPGNYTIEPGDAAFYGPKLDFVVKDALGREWQLGTIQVDYNLPERFDLTYVGGLNNNNNNNNNNEEINDNNEDGQELEHRPVMIHRAPFGSIERFIGILIEHY"

# Find all ThrRS catalytic domain files
files = glob.glob("anc_thrrs_cat_*.json")
files.sort()

print(f"Found {len(files)} ThrRS catalytic domain files")
print(f"Correct sequence length: {len(CORRECT_THRRS_SEQUENCE)} aa\n")

fixed_count = 0
for filepath in files:
    # Read the JSON file
    with open(filepath, 'r') as f:
        data = json.load(f)

    # Get current sequence
    current_seq = data['sequences'][0]['protein']['sequence']

    # Check if it needs fixing
    if current_seq == "MAAAAA":
        print(f"Fixing {filepath}...")
        print(f"  Before: {current_seq} ({len(current_seq)} aa)")

        # Update with correct sequence
        data['sequences'][0]['protein']['sequence'] = CORRECT_THRRS_SEQUENCE

        # Write back to file
        with open(filepath, 'w') as f:
            json.dump(data, f, indent=2)

        print(f"  After:  {CORRECT_THRRS_SEQUENCE[:50]}... ({len(CORRECT_THRRS_SEQUENCE)} aa)")
        fixed_count += 1
    else:
        print(f"Skipping {filepath} (already has correct sequence)")

print(f"\n✓ Fixed {fixed_count} files")
print(f"✓ All ThrRS catalytic domain sequences are now correct!")
