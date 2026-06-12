"""
Statistical analysis for the proteomic, lipidomic and functional-assay
data of *Hemipyrellia ligurriens* larval extracts (PLOS ONE, PONE-D-26-13228).

This script consolidates the inferential statistics underlying the
manuscript and the reviewer-response re-analysis. It is a faithful
extraction of the computations in the original `reviewer_response.py`
(PART G), with repo-relative paths and explicit methodological caveats.

It computes and prints:

  1. Proteome compartment statistics (from Supplementary.xlsx)
       - Hypergeometric test for WE/ES overlap significance
       - Two-proportion z-test for compartment-exclusivity asymmetry
       - Jaccard overlap coefficient
       - Strict differential-abundance counts (|log2FC| >= 1 and
         combined intensity >= 1e6), matching the volcano classification
  2. HaCaT MTT cytocompatibility (Fig 4)
       - mean, SD and one-sample two-tailed t-test of %-viability vs the
         100 % untreated control, per concentration, for WE and ES
  3. Agar-well-diffusion inhibition zones (manuscript Table 3)
       - one-way ANOVA across WE / ES / antibiotic control and 95 % CI

IMPORTANT CAVEATS (methodological honesty):
  * MTT — the manuscript p-values were computed in GraphPad Prism 8 on the
    absorbance data; the one-sample t-test on normalised %-viability
    reproduced here is an approximation and may differ slightly.
  * Inhibition zones — only published mean ± SD (n = 3) are available, so
    triplicates are reconstructed symmetrically from mean and SD before
    ANOVA; the ANOVA is therefore indicative, not a re-analysis of raw
    replicate measurements.

Usage:
    python3 scripts/stats_analysis.py
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom, f_oneway, ttest_1samp, norm

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "Supplementary.xlsx"


def rule(title: str) -> None:
    print("\n" + "=" * 70)
    print(title)
    print("=" * 70)


# ===================================================================
# 1. Proteome compartment statistics (from Supplementary.xlsx)
# ===================================================================
def proteome_stats() -> None:
    rule("1. PROTEOME COMPARTMENT STATISTICS  (Supplementary.xlsx)")
    df = pd.read_excel(SRC, sheet_name="proteinGroups", header=2)
    df = df.rename(columns={
        "Whole-body extracts -Intensity ": "WE",
        "Excretory-secretory (ES) products -Intensity ": "ES",
    })

    WE_n = int((df["WE"] > 0).sum())
    ES_n = int((df["ES"] > 0).sum())
    shared = int(((df["WE"] > 0) & (df["ES"] > 0)).sum())
    total = len(df)
    we_only = WE_n - shared
    es_only = ES_n - shared
    print(f"WE detected           : {WE_n:,}")
    print(f"ES detected           : {ES_n:,}")
    print(f"Shared (WE ∩ ES)      : {shared:,}")
    print(f"WE-exclusive          : {we_only:,}")
    print(f"ES-exclusive          : {es_only:,}")
    print(f"Union (universe)      : {total:,}")

    # --- Hypergeometric test for overlap significance ---
    # M = universe, n = WE detected, N = ES sample, k = observed overlap.
    overlap_p = float(hypergeom.sf(shared - 1, total, WE_n, ES_n))
    expected = ES_n * WE_n / total
    fold = shared / expected
    print("\n-- Hypergeometric overlap test --")
    print(f"  expected overlap    : {expected:.0f}")
    print(f"  observed overlap    : {shared:,}")
    print(f"  fold enrichment     : {fold:.2f}")
    print(f"  P(X >= k)           : {overlap_p:.3e}")
    print("  NB: P collapses because the union itself is the universe; the")
    print("      informative test is the exclusivity asymmetry below.")

    # --- Two-proportion z-test on exclusivity asymmetry ---
    p1 = we_only / WE_n
    p2 = es_only / ES_n
    pooled = (we_only + es_only) / (WE_n + ES_n)
    se = np.sqrt(pooled * (1 - pooled) * (1 / WE_n + 1 / ES_n))
    z = (p2 - p1) / se
    exclusive_p = 2 * (1 - norm.cdf(abs(z)))
    print("\n-- Two-proportion z-test (compartment-exclusivity asymmetry) --")
    print(f"  WE-exclusive fraction : {p1:.4f}  ({we_only:,}/{WE_n:,})")
    print(f"  ES-exclusive fraction : {p2:.4f}  ({es_only:,}/{ES_n:,})")
    print(f"  Z                     : {z:.2f}")
    print(f"  P                     : {exclusive_p:.3e}"
          + ("  (P < 1e-80)" if exclusive_p < 1e-80 else ""))

    # --- Jaccard ---
    jaccard = shared / (WE_n + ES_n - shared)
    print(f"\n  Jaccard overlap coefficient : {jaccard:.4f}")

    # --- Strict differential-abundance counts (volcano criteria) ---
    quant = df[(df["WE"] > 0) & (df["ES"] > 0)].copy()
    pc = min(df.loc[df["WE"] > 0, "WE"].min(),
             df.loc[df["ES"] > 0, "ES"].min()) / 10
    quant["log2FC"] = np.log2((quant["ES"] + pc) / (quant["WE"] + pc))
    quant["log10_sum"] = np.log10(quant["WE"] + quant["ES"] + pc)
    INT_THR = np.log10(1e6)
    up_es = int(((quant["log2FC"] >= 1) & (quant["log10_sum"] >= INT_THR)).sum())
    up_we = int(((quant["log2FC"] <= -1) & (quant["log10_sum"] >= INT_THR)).sum())
    print("\n-- Differential abundance (|log2FC| >= 1 and intensity >= 1e6) --")
    print(f"  co-quantified         : {len(quant):,}")
    print(f"  Up in ES              : {up_es:,}")
    print(f"  Up in WE              : {up_we:,}")


# ===================================================================
# 2. HaCaT MTT cytocompatibility (Fig 4)
# ===================================================================
# Normalised %-viability triplicates (24 h, n = 3); control = 100 %.
MTT = {
    "WE": {
        0.1: [103.034555, 99.769364, 105.249952],
        1.0: [110.629714, 113.165481, 114.109245],
        10:  [140.115769, 137.286181, 128.102683],
        100: [214.225574, 219.065924, 191.507431],
    },
    "ES": {
        0.1: [102.81719, 78.981019, 72.148860],
        1.0: [103.74255, 77.782218, 70.548219],
        10:  [96.40140, 75.704296, 76.270508],
        100: [86.98334, 63.516484, 62.224890],
    },
}


def sig(p):
    if p < 0.0001: return "****"
    if p < 0.001:  return "***"
    if p < 0.01:   return "**"
    if p < 0.05:   return "*"
    return "ns"


def mtt_stats() -> None:
    rule("2. HaCaT MTT CYTOCOMPATIBILITY  (one-sample t-test vs 100 % control)")
    print("Caveat: manuscript p-values were computed in GraphPad Prism 8 on")
    print("absorbance data; the values below are a reproducible approximation")
    print("from normalised %-viability and may differ slightly.\n")
    print(f"{'Cond':<5}{'Conc':>7}{'Mean':>9}{'SD':>8}{'t':>8}{'p':>11}{'sig':>6}")
    for cond, doses in MTT.items():
        for conc, reps in doses.items():
            reps = np.asarray(reps, float)
            m, s = reps.mean(), reps.std(ddof=1)
            t, p = ttest_1samp(reps, popmean=100.0)
            print(f"{cond:<5}{conc:>7}{m:>9.2f}{s:>8.2f}{t:>8.2f}"
                  f"{p:>11.4g}{sig(p):>6}")


# ===================================================================
# 3. Agar-well-diffusion inhibition zones (manuscript Table 3)
# ===================================================================
# Published mean (mm) and SD, n = 3. Raw triplicates were not reported;
# they are reconstructed symmetrically from mean and SD for an indicative
# one-way ANOVA across WE / ES / antibiotic control.
INHIBITION = {
    "Bacillus subtilis":          {"WE": (10.77, 0.75), "ES": (16.18, 0.61), "Ctrl": (35.09, 1.53)},
    "Staphylococcus epidermidis": {"WE": (0, 0),        "ES": (0, 0),        "Ctrl": (26.41, 1.00)},
    "Staphylococcus aureus":      {"WE": (0, 0),        "ES": (0, 0),        "Ctrl": (30.29, 2.14)},
    "Pseudomonas aeruginosa":     {"WE": (13.40, 3.44), "ES": (0, 0),        "Ctrl": (12.90, 1.22)},
    "Escherichia coli":           {"WE": (0, 0),        "ES": (0, 0),        "Ctrl": (28.41, 2.15)},
    "Proteus vulgaris":           {"WE": (0, 0),        "ES": (0, 0),        "Ctrl": (38.70, 3.81)},
    "Candida albicans":           {"WE": (0, 0),        "ES": (0, 0),        "Ctrl": (37.23, 2.01)},
}


def synth_triplicate(mean, sd):
    """Symmetric triplicate [m-d, m, m+d] whose sample SD equals `sd`."""
    if sd == 0:
        return [mean, mean, mean]
    d = sd * np.sqrt(3 / 2)
    return [mean - d, mean, mean + d]


def ci95(m, s):
    if s == 0:
        return (m, m)
    half = 1.96 * s / np.sqrt(3)
    return (m - half, m + half)


def inhibition_stats() -> None:
    rule("3. INHIBITION-ZONE STATISTICS  (one-way ANOVA + 95% CI, Table 3)")
    print("Caveat: triplicates reconstructed from published mean ± SD (n = 3);")
    print("ANOVA is indicative, not a re-analysis of raw measurements.\n")
    print(f"{'Organism':<28}{'WE':>14}{'ES':>14}{'Ctrl':>14}{'F':>9}{'p':>11}")
    for org, vals in INHIBITION.items():
        we = synth_triplicate(*vals["WE"])
        es = synth_triplicate(*vals["ES"])
        ct = synth_triplicate(*vals["Ctrl"])
        try:
            F, p = f_oneway(we, es, ct)
        except Exception:
            F, p = np.nan, np.nan
        we_ci = ci95(*vals["WE"])
        es_ci = ci95(*vals["ES"])
        we_s = f"{vals['WE'][0]:.1f}±{vals['WE'][1]:.2f}"
        es_s = f"{vals['ES'][0]:.1f}±{vals['ES'][1]:.2f}"
        ct_s = f"{vals['Ctrl'][0]:.1f}±{vals['Ctrl'][1]:.2f}"
        print(f"{org:<28}{we_s:>14}{es_s:>14}{ct_s:>14}{F:>9.2f}{p:>11.3g}")


def main() -> None:
    proteome_stats()
    mtt_stats()
    inhibition_stats()
    print()


if __name__ == "__main__":
    main()
