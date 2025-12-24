#!/usr/bin/env python3
# 01_define_targets.py
"""Define taxonomic sampling and molecular targets"""

import json
from pathlib import Path

# Taxonomically diverse species for deep phylogenetic sampling
SPECIES_LIST = {
    # Bacteria (thermophiles prioritized for ancestral inference)
    "bacteria": [
        "Escherichia coli",
        "Bacillus subtilis",
        "Thermus thermophilus",
        "Thermotoga maritima",
        "Deinococcus radiodurans",
        "Aquifex aeolicus",
        "Synechocystis sp. PCC 6803",
        "Mycobacterium tuberculosis",
        "Pseudomonas aeruginosa",
        "Helicobacter pylori",
        "Campylobacter jejuni",
        "Chlamydia trachomatis",
        "Rickettsia prowazekii",
        "Borrelia burgdorferi",
        "Treponema pallidum",
    ],
    # Archaea (essential for rooting)
    "archaea": [
        "Methanocaldococcus jannaschii",
        "Pyrococcus furiosus",
        "Thermococcus kodakarensis",
        "Sulfolobus acidocaldarius",
        "Archaeoglobus fulgidus",
        "Halobacterium salinarum",
        "Methanothermobacter thermautotrophicus",
        "Pyrobaculum aerophilum",
        "Nanoarchaeum equitans",
        "Thermoplasma acidophilum",
    ],
    # Eukaryotes (crown group diversity)
    "eukaryota": [
        "Homo sapiens",
        "Mus musculus",
        "Danio rerio",
        "Drosophila melanogaster",
        "Caenorhabditis elegans",
        "Saccharomyces cerevisiae",
        "Schizosaccharomyces pombe",
        "Candida albicans",
        "Plasmodium falciparum",
        "Trypanosoma brucei",
        "Leishmania major",
        "Dictyostelium discoideum",
        "Arabidopsis thaliana",
        "Oryza sativa",
        "Chlamydomonas reinhardtii",
    ],
}

# aaRS targets with Pfam domains
AARS_TARGETS = {
    "ProRS": {
        "name": "Prolyl-tRNA synthetase",
        "class": "II",
        "pfam": "PF03129",  # HGTP_anticodon (Class II, Pro-specific)
        "uniprot_keywords": ["Prolyl-tRNA synthetase", "Proline--tRNA ligase"],
    },
    "ThrRS": {
        "name": "Threonyl-tRNA synthetase",
        "class": "II",
        "pfam": "PF03129",
        "uniprot_keywords": ["Threonyl-tRNA synthetase", "Threonine--tRNA ligase"],
    },
    "SerRS": {
        "name": "Seryl-tRNA synthetase",
        "class": "II",
        "pfam": "PF00587",  # tRNA-synt_2b
        "uniprot_keywords": ["Seryl-tRNA synthetase", "Serine--tRNA ligase"],
    },
    "ValRS": {
        "name": "Valyl-tRNA synthetase",
        "class": "I",
        "pfam": "PF00133",  # tRNA-synt_1
        "uniprot_keywords": ["Valyl-tRNA synthetase", "Valine--tRNA ligase"],
    },
}

# tRNA targets (anticodon families)
TRNA_TARGETS = {
    "Pro": {"anticodons": ["AGG", "CGG", "TGG"], "rfam": "RF00005"},
    "Thr": {"anticodons": ["AGT", "CGT", "TGT", "GGT"], "rfam": "RF00005"},
    "Ser": {"anticodons": ["AGA", "CGA", "TGA", "GCT", "ACT"], "rfam": "RF00005"},
    "Val": {"anticodons": ["AAC", "CAC", "TAC", "GAC"], "rfam": "RF00005"},
}

def save_config():
    """Save configuration to JSON"""
    config = {
        "species": SPECIES_LIST,
        "aars": AARS_TARGETS,
        "trna": TRNA_TARGETS,
        "n_species": sum(len(v) for v in SPECIES_LIST.values()),
        "n_cores": 64,
        "n_gpus": 2,
    }
    
    with open("data/config.json", "w") as f:
        json.dump(config, f, indent=2)
    
    # Create species flat list
    all_species = []
    for domain in SPECIES_LIST.values():
        all_species.extend(domain)
    
    with open("data/species_list.txt", "w") as f:
        f.write("\n".join(all_species))
    
    print(f"✓ Configuration saved: {len(all_species)} species, 4 aaRS, 4 tRNAs")

if __name__ == "__main__":
    Path("data").mkdir(exist_ok=True)
    save_config()
