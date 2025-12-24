# AF3 Results Analysis - Corrected Extraction
## chain_pair_iptm[0][1] = Protein-AA Binding Score

### Key Findings

## 1. MODERN E. coli ThrRS (with Zn cofactor, 401 aa)

| Amino Acid | AA_iptm | Interpretation |
|------------|---------|----------------|
| **THR (cognate)** | **0.97** | Excellent specificity |
| SER (near-cognate) | 0.95 | High - explains editing requirement |
| VAL | 0.90 | Moderate |
| CYS, ASN | 0.89 | Moderate |
| ALA | 0.88 | Moderate |
| Other AAs | 0.79-0.87 | Lower |

## 2. MODERN E. coli ProRS (572 aa, no Zn)

| Amino Acid | AA_iptm | Interpretation |
|------------|---------|----------------|
| **PRO (cognate)** | **0.95** | Excellent |
| ALA | 0.93 | HIGH - promiscuous! |
| VAL | 0.92 | HIGH - promiscuous! |
| LEU | 0.90 | HIGH - promiscuous! |
| ILE, CYS | 0.89 | HIGH - promiscuous! |
| SER, THR | 0.87-0.88 | High |
| Other AAs | 0.72-0.82 | Moderate-High |

**ProRS shows HIGH binding to many non-cognate AAs** - maintains ancestral promiscuity!

## 3. ANCIENT ThrRS Catalytic Domain (278 aa, with/without Zn)

| Amino Acid | AA_iptm (no Zn) | AA_iptm (with Zn) |
|------------|-----------------|-------------------|
| **THR** | 0.85 | 0.89 |
| ARG | 0.87 | 0.88 |
| ILE | 0.87 | 0.86 |
| ALA | 0.86 | 0.84 |
| Other AAs | 0.72-0.86 | 0.80-0.88 |

**Ancient ThrRS was MUCH MORE PROMISCUOUS** - THR not significantly higher than other AAs!

## 4. ANCIENT ProRS Catalytic Domain (500 aa)

| Amino Acid | AA_iptm | Interpretation |
|------------|---------|----------------|
| GLU | 0.89 | Highest! |
| ILE | 0.86 | High |
| PRO (cognate) | 0.85 | Not highest |
| LYS | 0.85 | Equal to PRO |
| Other AAs | 0.78-0.84 | Similar range |

**Ancient ProRS showed NO preferential binding to PRO**

---

## COMPETITION EXPERIMENTS

### Modern ThrRS: THR vs ILE (with Zn)
- THR: **0.96**
- ILE: 0.79
- **DISCRIMINATES** - Evolved specificity!

### Ancient ThrRS: THR vs ILE (with Zn)
- THR: **0.89**
- ILE: 0.88
- **NO DISCRIMINATION** - Promiscuous ancestor!

---

## KEY EVOLUTIONARY INSIGHT

| Metric | Ancient | Modern | Evolution |
|--------|---------|--------|-----------|
| ThrRS-THR | 0.85-0.89 | **0.97** | +0.08-0.12 INCREASED |
| ThrRS-non-cognate | 0.72-0.88 | 0.79-0.95 | Similar range |
| ProRS-PRO | 0.85 | **0.95** | +0.10 INCREASED |
| ProRS-non-cognate | 0.78-0.89 | 0.72-0.93 | STILL HIGH |

**ThrRS evolved ~10% higher specificity for cognate THR**
**ProRS evolved ~10% higher binding for PRO but maintained promiscuity for other AAs**

---

## Data Notes

- 187 unique predictions extracted
- Zinc runs: [0][1] = AA binding, [0][2] = Zn binding
- tRNA runs (has_rna=True): Very low AA_iptm (0.17-0.38) - different binding geometry
- Some 663 aa modern ThrRS runs show low scores (0.36-0.66) - may be different construct
