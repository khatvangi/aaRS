# CRITICAL FINDINGS: Physical Mechanisms Behind ipTM Scores

## Analysis Overview
Comprehensive ligand interaction analysis of **1,177 AF3 structures** across all enzyme conditions.

**Key metrics extracted:**
- Zn²⁺ coordination distances
- Protein-ligand contact counts (4Å cutoff)
- H-bond networks
- Active site positioning

---

## Finding 1: Zn Coordination is NOT Discriminatory

### Initial Hypothesis (WRONG)
- Zn²⁺ discriminates between amino acids via coordination chemistry
- THR coordinates → accepted (high ipTM)
- ILE cannot coordinate → rejected (low ipTM)

### Actual Discovery
**ALL 20 amino acids coordinate Zn in modern ThrRS!**

| Ligand | ipTM | Zn Distance | Coordinating? |
|--------|------|-------------|---------------|
| THR    | 0.97 | 2.14 Å     | ✓ YES         |
| SER    | 0.95 | 2.14 Å     | ✓ YES         |
| VAL    | 0.90 | 2.18 Å     | ✓ YES         |
| ILE    | 0.83 | 2.19 Å     | ✓ YES         |
| PRO    | 0.82 | 2.20 Å     | ✓ YES         |
| ARG    | 0.79 | 2.88 Å     | ✓ YES         |

**19/20 ligands < 3.0 Å** (only TRP at 3.18 Å is borderline)

**Conclusion:** Zn does NOT discriminate - it ORGANIZES the active site!

---

## Finding 2: Zn Coordination EVOLVED

### Ancestral ThrRS + Zn
- **0/20 ligands coordinate Zn**
- All ligands 17-32 Å from Zn ion
- Zn present but NOT integrated into binding pocket
- Ligands float in pocket with minimal organization

| Metric | Ancestral + Zn | Modern + Zn |
|--------|----------------|-------------|
| Zn coordination (< 3Å) | 0% (0/20) | 95% (19/20) |
| Mean Zn distance | 29.8 Å | 2.3 Å |
| THR Zn distance | 31.4 Å | 2.1 Å |

### Evolutionary Timeline

```
Ancestral no Zn:
  THR ipTM: 0.85
  Zn: absent

         ↓ Add Zn (improves structure remotely)

Ancestral + Zn:
  THR ipTM: 0.89 (↑ 5%)
  Zn: 31 Å away (not coordinating)
  Mechanism: Electrostatic effects? Pocket reshaping?

         ↓ Evolve active site geometry

Modern + Zn:
  THR ipTM: 0.97 (↑ 9%)
  Zn: 2.1 Å (COORDINATING)
  Mechanism: Zn integrated into pocket → geometry optimization
```

**Key Insight:** Adding Zn to ancestral enzyme helps BEFORE it's coordinating! Evolution then positions Zn into the active site for maximum effect.

---

## Finding 3: ipTM Correlates with Physical Contacts, Not Coordination

### Modern ThrRS Correlations with ipTM

| Parameter | Correlation (R) | Mechanism |
|-----------|----------------|-----------|
| **Zn distance** | -0.564 | Closer = better (but all coordinate) |
| **Contact atoms** | -0.366 | Fewer contacts = better fit? |
| **H-bonds** | +0.093 | Weak correlation |

### Top vs Bottom Binders

**High ipTM (THR, SER):**
- Zn: ~2.14 Å (coordinating)
- Contacts: 15-17 atoms
- H-bonds: 6
- Tight geometric fit

**Low ipTM (ARG, TRP):**
- Zn: 2.88-3.18 Å (still coordinating!)
- Contacts: 17-30 atoms (more, not less!)
- H-bonds: 3-5
- Bulkier, poorer geometric fit

**Conclusion:** Discrimination comes from **POCKET GEOMETRY**, not metal coordination chemistry.

---

## Finding 4: Editing Domain Uses Contact Networks for Selectivity

### ProRS Editing Domain (No Zn)

**Best structures per ligand:**

| Ligand | ipTM | Contacts (4Å) | H-bonds | Role |
|--------|------|---------------|---------|------|
| **THR** | **0.87** | **25** | **7** | **Error (binds best)** |
| **PRO** | **0.82** | **18** | **7** | **Cognate (binds worse!)** |
| ALA    | 0.80 | 5  | 0 | Error |
| VAL    | 0.78 | 23 | 6 | Error |
| LEU    | 0.77 | 25 | 6 | Error |

### Inverted Selectivity Mechanism

**THR binds 6% better than PRO** in editing domain

**Molecular basis:**
1. **THR makes 39% more contacts** (25 vs 18 atoms within 4Å)
2. Same H-bond count (7), but different geometry
3. THR -OH group creates additional interaction surface

**This is the double-sieve mechanism:**
- Catalytic site: PRO ranks #1 (cognate preferred)
- Editing site: THR ranks #1 (error preferred) ← **INVERTED!**

Errors that escape catalytic site are SELECTIVELY bound and hydrolyzed by editing domain.

---

## Finding 5: The Real Zinc Mechanism

### What Zn Actually Does

**NOT:** Discriminate via coordination chemistry
**NOT:** Create a filter that rejects hydrophobics

**YES:** Organize active site geometry
**YES:** Position ligands for optimal interactions
**YES:** Enable pocket evolution toward specificity

### Evidence

1. **All ligands coordinate** in modern enzyme → no inherent discrimination
2. **Ancestral enzyme:** Zn helps before coordination evolves → structural role
3. **Modern enzyme:** Tight correlation between Zn distance and ipTM (-0.564) → geometric optimization

### Revised Model

```
ANCESTRAL ThrRS:
┌─────────────────┐
│  Active Site    │
│                 │
│   [ligands]     │
│    floating     │
│                 │
│                 │
│  Zn²⁺ (distant) │ ← electrostatic effects
└─────────────────┘

MODERN ThrRS:
┌─────────────────┐
│  Active Site    │
│      Zn²⁺       │ ← central organizer
│    ╱  │  ╲      │
│  [ligand]       │ ← coordinated, positioned
│   optimized     │
│   geometry      │
└─────────────────┘
```

The Zn ion became a **STRUCTURAL SCAFFOLD** that enables precise pocket geometry, not a chemical discriminator.

---

## Finding 6: SER is Still Trapped (Despite Coordination)

### The "Zinc Trap" Still Exists

| Ligand | ipTM | Zn Dist | Contacts | H-bonds | % of THR |
|--------|------|---------|----------|---------|----------|
| **THR** | 0.97 | 2.14 Å | 17 | 6 | 100% |
| **SER** | 0.95 | 2.14 Å | 15 | 6 | 97.9% |

**SER coordinates Zn identically to THR** (both 2.14 Å)

**Discrimination is minimal:** Only 2.1% difference in ipTM

**This validates the need for EDITING DOMAIN!**

Modern ThrRS cannot distinguish SER from THR reliably, even with optimized Zn coordination. The editing domain provides additional quality control.

---

## Summary: The Complete Molecular Picture

### ProRS Strategy (Kinetic Solution)
1. **Catalytic domain:** Promiscuous (PRO #1, but ALA=98%, VAL=97%, THR=92%)
2. **Editing domain:** Selective for errors (THR #1 > PRO #2)
3. **Mechanism:** Contact networks and H-bonds, no metal
4. **Result:** Double-sieve kinetic proofreading

### ThrRS Strategy (Structural Solution)
1. **Evolution:** Integrate Zn into active site (distant → coordinating)
2. **Mechanism:** Zn organizes pocket geometry for tight fits
3. **Discrimination:** Pocket shape, not coordination chemistry
4. **Limitation:** SER still 98% of THR (zinc trap)
5. **Result:** Structural filter + minimal editing

### Key Insight
**Neither enzyme uses simple chemical discrimination!**

Both rely on **GEOMETRIC COMPLEMENTARITY:**
- ProRS editing: Shape complementarity to errors
- ThrRS pocket: Shape complementarity to cognate (with Zn scaffold)

The remarkable specificity we observe comes from **EVOLVED PROTEIN ARCHITECTURE**, not from inherent chemical properties of the substrates.

---

## Data Files Generated

1. **`figures/data/comprehensive_ligand_analysis.csv`**
   Main dataset: 1,187 structures × 20+ features

2. **`figures/data/comprehensive_ligand_analysis.json`**
   Detailed data including H-bond atom lists

3. **`figures/mechanisms/physical_mechanisms.png`**
   5-panel figure showing all mechanisms

4. **`figures/schematics/fig2a_double_sieve_mechanism.png`**
   Double-sieve schematic

5. **`figures/schematics/fig4a_zinc_filter_mechanism.png`**
   Zinc filter schematic (updated with new understanding)

---

## Implications for Manuscript

### Revise These Statements
1. ~~"Zn discriminates via coordination chemistry"~~
   → "Zn organizes active site geometry for selective binding"

2. ~~"Hydrophobic residues cannot coordinate Zn"~~
   → "All residues coordinate Zn; discrimination comes from pocket fit"

3. ~~"Zinc filter rejects ILE/VAL"~~
   → "Zn-organized pocket geometry disfavors bulky hydrophobics"

### Emphasize These Findings
1. **Zn coordination evolved** (0% ancestral → 95% modern)
2. **Zn helps before coordinating** (ancestral +Zn shows improvement)
3. **Editing domain has inverted selectivity** (THR > PRO)
4. **Geometric complementarity** is the universal mechanism

### Figure Updates Needed
1. Update zinc filter mechanism (panel 4a) to show organization vs discrimination
2. Add panel showing ancestral vs modern Zn positioning
3. Highlight contact networks in editing domain figure

---

## Computational Methods Note

**This analysis was only possible with AlphaFold3:**
- Accurate ligand positioning
- Reliable inter-atomic distances
- Multi-condition predictions (ancestral/modern, ±Zn)

Traditional crystallography would require 100+ structures across all conditions. AF3 provided comprehensive coverage in silico.
