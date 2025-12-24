# BioRender Schematic Instructions
## aaRS Evolution - Mechanism Diagrams

**Generated:** 2025-12-18
**Access:** https://app.biorender.com

---

## 🎨 Figure 1B: Domain Architecture Comparison

### Concept
Show side-by-side comparison of ProRS (has editing domain) vs ThrRS (no editing domain)

### Icons to Search & Use
1. **"Aminoacyl tRNA synthetase (cartoon)"** - Use for base enzyme shape
2. **"Protein domain"** or **"Generic enzyme"** - Use for modular domains
3. **"Enzyme (subunit, cartoon)"** - Purple oval shapes for domains

### Layout
```
ProRS (Modern):
┌────────────┬──────────┬────────────┐
│ Catalytic  │   INS    │  Editing   │
│  Domain    │ Domain   │  Domain    │
│  (500 aa)  │          │  (300 aa)  │
└────────────┴──────────┴────────────┘
     ↓             ↓           ↓
  [Purple]     [Teal]     [Orange]

ThrRS (Modern):
┌────────────┐
│ Catalytic  │
│  Domain    │
│  (663 aa)  │
└────────────┘
     ↓
  [Blue]
```

### Color Scheme
- **ProRS Catalytic**: Purple (#9b59b6)
- **ProRS INS**: Teal (#16a085)
- **ProRS Editing**: Orange (#f39c12)
- **ThrRS Catalytic**: Blue (#3498db)

### Labels
- Add text: "ProRS: 3 domains"
- Add text: "ThrRS: 1 domain"
- Arrow annotation: "Editing domain compensates for promiscuity"

### Key Message
ProRS evolved editing domain while ThrRS evolved zinc discrimination

---

## 🎨 Figure 2A: Double Sieve Mechanism

### Concept
Flow diagram showing two-stage quality control in ProRS

### Icons to Search & Use
1. **"Generic enzyme 1 (schematic)"** - Shows enzyme with binding pocket
2. **"Enzyme-ligand complex"** - For substrate binding
3. **"Protease-substrate binding"** - For sequential steps
4. Search **"filter"**, **"sieve"**, **"quality control"** for filtering icons

### Layout
```
      [Amino Acid Pool]
            ↓
    ┌───────────────┐
    │   SIEVE 1     │
    │  Catalytic    │  ← COARSE FILTER
    │   Domain      │     Lets errors through
    └───────────────┘
            ↓
      ┌─────┬─────┐
      │     │     │
    PRO   THR   ALA
    (OK)  (ERROR)(ERROR)
      │     │     │
      ↓     ↓     ↓
    ┌───────────────┐
    │   SIEVE 2     │
    │   Editing     │  ← FINE FILTER
    │   Domain      │     Catches errors
    └───────────────┘
            ↓
      ┌─────┬─────┐
      │     │     │
    PRO    [X]   [X]
  (PASS) (HYDROLYZED)
```

### Colors
- **Sieve 1 (Catalytic)**: Purple
- **Sieve 2 (Editing)**: Orange
- **PRO (correct)**: Green circles
- **THR, ALA (errors)**: Red circles
- **Hydrolyzed**: Gray X marks

### Arrows & Labels
- Downward arrows between stages
- Label: "Coarse discrimination" at Sieve 1
- Label: "Fine discrimination" at Sieve 2
- Add tRNA molecules (search "tRNA") at each stage

### Key Message
Two-stage quality control: catalytic site is coarse, editing site is fine

---

## 🎨 Figure 3A: Zinc Coordination Chemistry

### Concept
Show how Zn2+ coordinates with different amino acids

### Icons to Search & Use
1. Search **"zinc"**, **"metal ion"**, **"coordination"**
2. Search **"chemical bond"**, **"covalent bond"**
3. Use **"atom"** or **"molecule"** icons for Zn sphere

### Layout (Side-by-side comparison)

**Panel 1: THR Bidentate Coordination**
```
        O-H (hydroxyl)
         │
         │ ← Coordination bond 1
       [Zn²⁺]
         │ ← Coordination bond 2
         │
        NH₂ (amino)

Label: "Bidentate (2 bonds)"
```

**Panel 2: ILE Cannot Coordinate**
```
      CH₃-CH₂ (hydrophobic)
          │
          │  ← No coordination
        [Zn²⁺]  (gap)

          X

Label: "Hydrophobic - rejected"
```

**Panel 3: SER Bidentate (The Trap!)**
```
        O-H (hydroxyl)
         │
         │ ← Coordination bond 1
       [Zn²⁺]
         │ ← Coordination bond 2
         │
        NH₂ (amino)

Label: "Also bidentate (trapped!)"
```

### Colors
- **Zn ion**: Gray sphere
- **Coordination bonds**: Dashed lines (blue)
- **THR**: Green chemical structure
- **ILE**: Red chemical structure (with X)
- **SER**: Orange chemical structure (with warning symbol)

### Key Message
Zn discriminates by coordination geometry - THR and SER both coordinate (trap!)

---

## 🎨 Figure 4A: Zinc Filter Mechanism

### Concept
Show Zn as gatekeeper in active site pocket

### Icons to Search & Use
1. **"Generic enzyme"** with active site
2. Search **"gate"**, **"checkpoint"**, **"barrier"**
3. Search **"active site"**, **"binding pocket"**

### Layout
```
                 ENTRANCE
                     ↓
    ┌────────────────────────────────┐
    │    ACTIVE SITE POCKET          │
    │                                │
    │         [Zn²⁺] FILTER         │
    │            ⬇                   │
    │    ┌──────────────┐           │
    │    │ CATALYTIC    │           │
    │    │ RESIDUES     │           │
    │    └──────────────┘           │
    └────────────────────────────────┘

ACCEPTED ✓:              REJECTED ✗:
• THR (hydroxyl)         • ILE (methyl)
• SER (hydroxyl)         • VAL (methyl)
                         • LEU (methyl)
```

### Colors
- **Pocket**: Light purple background
- **Zn sphere**: Gray with glow
- **Accepted AAs**: Green with checkmarks
- **Rejected AAs**: Red with X marks

### Annotations
- Arrow pointing to Zn: "Steric filter"
- Label: "Requires hydroxyl for coordination"
- Label: "Hydrophobics cannot enter"

### Key Message
Zn acts as molecular gatekeeper based on chemical properties

---

## 🎨 Figure 5A: The Zinc Trap Concept

### Concept
Show both THR and SER coordinating Zn identically, requiring editing domain

### Icons to Search & Use
1. **"Enzyme-ligand complex"**
2. Search **"trap"**, **"escape route"**
3. Search **"error correction"**, **"proofreading"**

### Layout
```
    ZN FILTER (Catalytic Site)
    ┌─────────────────────────┐
    │     [Zn²⁺]              │
    │      │  │                │
    │    THR  SER              │
    │     ✓   ?                │  ← Both coordinate!
    └─────────────────────────┘
            │    │
            │    └──────┐
            │           ↓
            │    ┌──────────────┐
            │    │ SER TRAPPED! │
            ↓    └──────────────┘
    ┌─────────────────────────┐
    │  EDITING DOMAIN         │
    │  (Escape Route)         │
    │                         │
    │  THR → Hydrolyze ✓     │
    │  SER → Hydrolyze ✓     │
    └─────────────────────────┘
            ↓
        [Corrected]
```

### Colors
- **Zn filter**: Purple box
- **THR**: Green (with checkmark)
- **SER**: Orange (with question mark, then warning)
- **Editing domain**: Orange box (safety net)
- **Final state**: Green (corrected)

### Annotations
- Label at Zn site: "Both coordinate Zn identically!"
- Arrow to editing: "Mandatory correction"
- Label at editing: "Chemical backup required"

### Key Message
Zn filter fails against SER → editing domain is mandatory, not optional

---

## 📋 General BioRender Tips

### Workflow
1. Go to https://app.biorender.com
2. Create new project: "aaRS Evolution Schematics"
3. For each figure:
   - Start with blank canvas
   - Search for icons using terms above
   - Arrange according to layouts
   - Add text labels and annotations
   - Use arrow tools for flow
   - Color using hex codes provided

### Text Settings
- **Title font**: Bold, 18-20 pt
- **Labels**: Regular, 12-14 pt
- **Annotations**: Italic, 10-12 pt

### Export Settings
- **Format**: PNG (high resolution) + PDF (vector)
- **Resolution**: 300 DPI minimum
- **Size**: 3000 x 3000 px for square diagrams
- **Background**: White

### Color Consistency (Match Data Figures)
- **Cognate/Correct**: #2ecc71 (green)
- **Error/Trapped**: #f39c12 (orange)
- **Rejected**: #e74c3c (red)
- **Zn/Metal**: #95a5a6 (gray)
- **ProRS theme**: Purple/Orange
- **ThrRS theme**: Blue

---

## 📁 Save Locations

After creating in BioRender, export and save to:

```
/storage/kiran-stuff/aaRS/phase2/figures/biorender/

├── fig1b_domain_architecture.png
├── fig1b_domain_architecture.pdf
├── fig2a_double_sieve.png
├── fig2a_double_sieve.pdf
├── fig3a_zn_coordination.png
├── fig3a_zn_coordination.pdf
├── fig4a_zn_filter_mechanism.png
├── fig4a_zn_filter_mechanism.pdf
├── fig5a_zinc_trap_concept.png
└── fig5a_zinc_trap_concept.pdf
```

---

## ✅ Quality Checklist

Before exporting each schematic:

- [ ] All icons clearly visible
- [ ] Text is readable (not too small)
- [ ] Colors match overall figure scheme
- [ ] Arrows show clear direction of flow
- [ ] Key message is immediately obvious
- [ ] No clutter - keep it simple
- [ ] Labels are accurate (check spelling)
- [ ] Legend included if needed
- [ ] Background is white (for print)
- [ ] High resolution (300 DPI)

---

## 🎯 Priority Order

Based on manuscript flow:

1. **Figure 1B** - Domain architecture (context)
2. **Figure 2A** - Double sieve (ProRS solution)
3. **Figure 3A** - Zn coordination (chemical basis)
4. **Figure 4A** - Zn filter (ThrRS solution)
5. **Figure 5A** - Zinc trap (why editing still needed)

---

## 💡 Pro Tips

1. **Use templates**: BioRender has "protein structure" templates
2. **Layer management**: Keep related elements grouped
3. **Alignment tools**: Use grid and snap-to-align
4. **Color picker**: Save custom colors for consistency
5. **Duplicate elements**: Copy-paste for symmetry
6. **Preview**: Check at actual print size before exporting

---

**Estimated time per schematic:** 30-60 minutes
**Total estimated time:** 2-4 hours

**Note:** These schematics complement the data visualization panels already generated. They provide the mechanistic context that data alone cannot show.

---

**Created:** 2025-12-18
**For:** aaRS Evolution manuscript (Cell/NSMB)
