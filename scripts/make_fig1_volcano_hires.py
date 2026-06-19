"""
High-resolution re-render of Fig 1 (volcano-style differential-abundance
plot, ES vs WE).

Reproduces the volcano figure from `make_figures.py` exactly (same data,
classification, palette and label set) but renders at a configurable high
DPI (default 600) with enlarged axis-label and tick-label fonts for
figures destined for print where small text must remain legible.

Non-destructive: writes a separate output file (default
`figures/Fig1_volcano.tif`) and never touches the 320-dpi
`figures/Fig1_volcano.png` consumed by the HTML pipeline.

Usage:
    python3 scripts/make_fig1_volcano_hires.py [--dpi 600] [--out PATH]
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from adjustText import adjust_text

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "Supplementary.xlsx"
FIG_DIR = ROOT / "figures"


def build_dataframe() -> pd.DataFrame:
    df = pd.read_excel(SRC, sheet_name="proteinGroups", header=2)
    df = df.rename(columns={
        "Whole-body extracts -Intensity ": "WE",
        "Excretory-secretory (ES) products -Intensity ": "ES",
        "Sequence coverage (%)": "Coverage",
        "Molecular weight (kDa)": "MW",
        "Sequence length": "Length",
    })

    def parse_name(s):
        if not isinstance(s, str):
            return ("", "")
        gn = re.search(r"GN=(\S+)", s)
        desc = re.search(r"\|[A-Z0-9_]+\s+(.*?)\s+OS=", s)
        return (gn.group(1) if gn else "", desc.group(1) if desc else "")

    df["Gene"] = [p[0] for p in df["Protein name"].apply(parse_name)]
    df["Description"] = [p[1] for p in df["Protein name"].apply(parse_name)]

    # Pseudo-count to handle zeros (one tenth of smallest non-zero across both)
    nz_min = min(df.loc[df["WE"] > 0, "WE"].min(),
                 df.loc[df["ES"] > 0, "ES"].min())
    PC = nz_min / 10.0
    df["WE_p"] = df["WE"] + PC
    df["ES_p"] = df["ES"] + PC
    df["log2FC"] = np.log2(df["ES_p"] / df["WE_p"])
    df["log10_sum"] = np.log10(df["WE_p"] + df["ES_p"])

    FC_THR = 1.0
    INT_THR = np.log10(1e6)

    def classify(r):
        if r["WE"] == 0 and r["ES"] > 0:
            return "ES-exclusive"
        if r["ES"] == 0 and r["WE"] > 0:
            return "WE-exclusive"
        if r["log2FC"] >= FC_THR and r["log10_sum"] >= INT_THR:
            return "Up in ES"
        if r["log2FC"] <= -FC_THR and r["log10_sum"] >= INT_THR:
            return "Up in WE"
        return "Non-significant"

    df["Class"] = df.apply(classify, axis=1)
    return df


PALETTE = {
    "Up in ES":        "#D72660",
    "Up in WE":        "#1F77B4",
    "ES-exclusive":    "#F08C1E",
    "WE-exclusive":    "#6A3FA0",
    "Non-significant": "#B8BFC9",
}
FC_THR = 1.0
INT_THR = np.log10(1e6)

AMP_LABEL = {
    "P86471": "Lucifensin", "A0A7D5FFX3": "Lucifensin",
    "A0A0L0C905": "Defensin", "A0A0L0BRJ2": "Defensin", "P10891": "Phormicin",
    "C0HJX9": "Diptericin", "P18684": "Diptericin-D",
    "A0A0L0BNM1": "Lysozyme", "C0HLB7": "Lysozyme 2",
    "D9J148": "Lysozyme 2 (put.)", "D9J143": "Lysozyme 1B",
    "A0A0L0CIT1": "Lysozyme D",
    "A0A0L0C6T8": "Attacin-A", "A0A0L0BXA2": "Attacin C",
    "A0A0L0C799": "Attacin C", "A0A0L0C779": "Attacin C",
    "A0A7D5FH37": "Chymotrypsin", "Q25237": "Serine proteinase",
    "Q25231": "Serine proteinase", "Q9GSL6": "Ser. protease K16",
    "Q9GSN2": "Ser. protease K13", "Q9GSM7": "Ser. protease K16",
    "Q9GSM9": "Ser. protease K15", "Q9GSM8": "Ser. protease K3",
    "Q9GSL7": "Ser. protease K12", "Q9GSL9": "Ser. protease K14",
    "A0A0L0BYT7": "Ser/Thr kinase", "A0A0L0C1F0": "Ser/Thr kinase ATM",
    "A0A0L0BYH9": "Metallothionein-1", "A0A0L0BYA8": "Metallothionein-1",
    "A0A0U3AZ32": "Metallothionein-like",
    "D9J150": "Ferritin", "A0A0L0BQ22": "Cu-Zn SOD",
    "A0A0L0BQK7": "Insulin-like", "A0A0L0C4Y9": "Insulin-like 7",
}
AMP_PRIORITY = set(AMP_LABEL.keys())


def render(df: pd.DataFrame, dpi: int, out: Path,
           label_fs: float, tick_fs: float, heading: bool = True) -> None:
    # Enlarged axis-label and tick-label fonts (the rest of the styling
    # matches make_figures.py's global rcParams).
    mpl.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 11,
        "axes.titlesize": 14,
        "axes.titleweight": "bold",
        "axes.labelsize": label_fs,
        "axes.labelweight": "bold",
        "axes.linewidth": 1.2,
        "axes.edgecolor": "#222222",
        "xtick.labelsize": tick_fs,
        "ytick.labelsize": tick_fs,
        "xtick.major.width": 1.1,
        "ytick.major.width": 1.1,
        "xtick.major.size": 4.5,
        "ytick.major.size": 4.5,
        "legend.frameon": False,
        "savefig.bbox": "tight",
        "savefig.facecolor": "white",
        "axes.spines.top": False,
        "axes.spines.right": False,
    })

    fig, ax = plt.subplots(figsize=(7.5, 5.6))
    quant = df[(df["WE"] > 0) & (df["ES"] > 0)].copy()

    for cls in ["Non-significant", "Up in WE", "Up in ES"]:
        sub = quant[quant["Class"] == cls]
        ax.scatter(sub["log2FC"], sub["log10_sum"],
                   s=14 if cls == "Non-significant" else 24,
                   c=PALETTE[cls],
                   alpha=0.45 if cls == "Non-significant" else 0.85,
                   edgecolors="white", linewidths=0.3,
                   label=f"{cls} (n={len(sub):,})",
                   zorder=2 if cls == "Non-significant" else 3)

    ax.axvline(FC_THR, color="#444", lw=0.9, ls="--", alpha=0.6)
    ax.axvline(-FC_THR, color="#444", lw=0.9, ls="--", alpha=0.6)
    ax.axhline(INT_THR, color="#444", lw=0.9, ls="--", alpha=0.6)

    def _primary_acc(s):
        return str(s).split(";")[0]

    quant["_acc"] = quant["Protein IDs"].apply(_primary_acc)
    quant["_is_priority"] = quant["_acc"].isin(AMP_PRIORITY)

    sig = quant[quant["Class"].isin(["Up in ES", "Up in WE"])].copy()
    top_es = sig[sig["Class"] == "Up in ES"].nlargest(5, "log2FC")
    top_we = sig[sig["Class"] == "Up in WE"].nsmallest(5, "log2FC")
    top_fc = pd.concat([top_es, top_we])
    amp_pri = sig[sig["_is_priority"]].copy()
    top_fc_only = top_fc[~top_fc["_acc"].isin(amp_pri["_acc"])]

    texts = []
    for _, r in top_fc_only.iterrows():
        texts.append(ax.text(r["log2FC"], r["log10_sum"], r["_acc"],
                             fontsize=10, fontweight="normal",
                             color="#555", ha="center", va="center",
                             alpha=0.85))
    for _, r in amp_pri.iterrows():
        lab = AMP_LABEL.get(r["_acc"], r["_acc"])
        texts.append(ax.text(r["log2FC"], r["log10_sum"], lab,
                             fontsize=11, fontweight="bold",
                             color="#111", ha="center", va="center",
                             bbox=dict(boxstyle="round,pad=0.20",
                                       fc="#FFF6E0", ec="#C68C2D",
                                       lw=0.7, alpha=0.92)))
    adjust_text(texts, ax=ax,
                arrowprops=dict(arrowstyle="-", color="#777", lw=0.55, alpha=0.7),
                expand_points=(1.8, 2.2), expand_text=(1.3, 1.4),
                force_text=(0.6, 0.9), force_points=(0.5, 0.8))

    ax.set_xlabel("log$_2$ fold change (ES / WE)")
    ax.set_ylabel("log$_{10}$ combined intensity")
    if heading:
        ax.set_title(r"Differential abundance of $\mathit{Hemipyrellia\ ligurriens}$ proteins:"
                     "\nExcretory–Secretory (ES) vs Whole-body Extracts (WE)",
                     pad=14)

    ymin, ymax = ax.get_ylim()
    ax.set_ylim(ymin, ymax + 0.5)
    ax.text(0.02, 0.97, "← Enriched in WE", transform=ax.transAxes,
            ha="left", va="top", fontsize=12, fontweight="bold",
            color=PALETTE["Up in WE"],
            bbox=dict(boxstyle="round,pad=0.35", fc="white",
                      ec=PALETTE["Up in WE"], lw=1.2))
    ax.text(0.98, 0.97, "Enriched in ES →", transform=ax.transAxes,
            ha="right", va="top", fontsize=12, fontweight="bold",
            color=PALETTE["Up in ES"],
            bbox=dict(boxstyle="round,pad=0.35", fc="white",
                      ec=PALETTE["Up in ES"], lw=1.2))

    ex_es = (df["Class"] == "ES-exclusive").sum()
    ex_we = (df["Class"] == "WE-exclusive").sum()
    ax.text(0.01, 0.01,
            f"Co-quantified proteins shown = {len(quant):,}\n"
            f"Fraction-exclusive (not shown): ES-only n={ex_es:,}, WE-only n={ex_we:,}",
            transform=ax.transAxes, ha="left", va="bottom",
            fontsize=9.5, style="italic", color="#555",
            bbox=dict(boxstyle="round,pad=0.3", fc="#f6f7fa", ec="none"))

    ax.legend(loc="lower right", fontsize=11, frameon=True, framealpha=0.9,
              handletextpad=0.4, edgecolor="#ddd")
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=dpi, facecolor="white",
                pil_kwargs=({"compression": "tiff_lzw"}
                if str(out).lower().endswith((".tif", ".tiff")) else {}))
    plt.close(fig)

    from PIL import Image
    im = Image.open(out)
    print(f"wrote {out}")
    print(f"  dpi={im.info.get('dpi')}  px={im.size}  "
          f"print={im.size[0]/dpi:.2f}x{im.size[1]/dpi:.2f} in")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dpi", type=int, default=300)
    ap.add_argument("--out", type=Path,
                    default=FIG_DIR / "Fig2.tif")
    ap.add_argument("--label-fontsize", type=float, default=12,
                    help="axis x/y label font size (default 16; orig 12)")
    ap.add_argument("--tick-fontsize", type=float, default=11,
                    help="axis tick-number font size (default 14; orig 11)")
    ap.add_argument("--heading", action="store_true",
                    help="include the figure title (off by default; the "
                         "title belongs in the figure caption)")
    args = ap.parse_args()

    df = build_dataframe()
    render(df, args.dpi, args.out, args.label_fontsize, args.tick_fontsize,
           heading=args.heading)


if __name__ == "__main__":
    main()
