# *Hemipyrellia ligurriens* — figure & statistics source code

Analysis source code for the proteomic and lipidomic characterisation of
*Hemipyrellia ligurriens* (Diptera: Calliphoridae) larval **whole-body
extracts (WE)** and **excretory–secretory (ES) products**.

This repository contains the Python that regenerates the main-text
figures (Fig 1–4) and computes the inferential statistics reported in the
manuscript. Every figure and statistic is produced deterministically from
the processed data; no values are hand-edited.

## Contents

```
.
├── Supplementary.xlsx       MaxQuant proteinGroups output (raw input)
└── scripts/
    ├── make_fig1_proteome_hires.py        Fig 1 — total-proteome composite (A–D)
    ├── make_fig1_volcano_hires.py         Fig 1 panel C — differential-abundance volcano
    ├── make_fig2_scatter_hires.py         Fig 1 panel D — pairwise ES-vs-WE scatter
    ├── make_fig2_lowmw_combined_hires.py   Fig 2 — < 20 kDa overview + family composition
    ├── make_fig3_family_heatmap_hires.py   Fig 3 — family composition × compartment
    ├── make_fig4_fatty_acid_hires.py      Fig 3 — fatty-acid composition (GC-MS)
    ├── make_fig4_hacat_mtt_hires.py       Fig 4 — HaCaT MTT cytocompatibility
    └── stats_analysis.py                  Inferential statistics (all assays)
```

## Requirements

Python 3.9+ and:

```bash
pip install numpy pandas matplotlib scipy openpyxl matplotlib-venn adjustText Pillow
```

## Usage

Each script is self-contained and writes PLOS-ready 300-dpi TIFFs (LZW-compressed, <=2250x2625 px) to a `figures/`
directory (created automatically). Run any script directly:

```bash
python3 scripts/make_fig4_hacat_mtt_hires.py
python3 scripts/stats_analysis.py
```

Most scripts accept `--dpi` and `--out`; the volcano additionally accepts
`--no-heading`. Run a script with `-h` for details.

## What the statistics script computes

`scripts/stats_analysis.py` prints:

1. **Proteome compartment statistics** (from `Supplementary.xlsx`) —
   hypergeometric overlap test, two-proportion z-test for compartment
   exclusivity, Jaccard coefficient, and strict differential-abundance
   counts (|log2 FC| ≥ 1 and combined intensity ≥ 1e6).
2. **HaCaT MTT** — mean, SD and one-sample two-tailed t-test of
   %-viability versus the 100 % untreated control.
3. **Inhibition zones** — one-way ANOVA and 95 % CI across WE / ES /
   antibiotic control.

## Data and methodological notes

- **Protein identities are putative**, assigned by homology to the
  UniProt *Calliphoridae* database. Species names in UniProt descriptions
  reflect the closest-homolog match and do **not** indicate species of
  origin.
- The LC-MS/MS quantification is **single-replicate**; differential
  abundance is exploratory.
- **MTT** — the manuscript p-values were computed in GraphPad Prism 8 on
  absorbance data; the one-sample t-test on normalised %-viability in
  `stats_analysis.py` is a reproducible approximation and may differ
  slightly.
- **Inhibition zones** — only published mean ± SD (n = 3) are available,
  so triplicates are reconstructed from mean and SD before ANOVA; the
  ANOVA is indicative rather than a re-analysis of raw measurements.
