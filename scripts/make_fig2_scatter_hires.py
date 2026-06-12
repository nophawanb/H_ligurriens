"""
High-resolution re-render of Fig 2 (pairwise log10 intensity scatter,
ES vs WE) with density colouring.

Reproduces the scatter from `make_figures.py` exactly (same data, density
estimate, fold-change guides and Table-1 AMP highlights) but renders at a
configurable high DPI (default 600) and OMITS both the figure title
heading and the below-graph legend / n-count statistics (intended for
use as a composed panel / for journals that set titles and keys in the
caption rather than on the artwork).

Non-destructive: writes a separate output file (default
`figures/Fig2_scatter_600dpi.png`) and never touches the 320-dpi
`figures/Fig2_scatter.png` consumed by the HTML pipeline.

Usage:
    python3 scripts/make_fig2_scatter_hires.py [--dpi 600] [--out PATH]
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.stats import gaussian_kde

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "Supplementary.xlsx"
FIG_DIR = ROOT / "figures"

PALETTE = {
    "Up in ES":        "#D72660",
    "Up in WE":        "#1F77B4",
    "ES-exclusive":    "#F08C1E",
    "WE-exclusive":    "#6A3FA0",
    "Non-significant": "#B8BFC9",
}

AMP_PRIORITY = {
    "P86471", "A0A7D5FFX3", "A0A0L0C905", "A0A0L0BRJ2", "P10891",
    "C0HJX9", "P18684", "A0A0L0BNM1", "C0HLB7", "D9J148", "D9J143",
    "A0A0L0CIT1", "A0A0L0C6T8", "A0A0L0BXA2", "A0A0L0C799", "A0A0L0C779",
    "A0A7D5FH37", "Q25237", "Q25231", "Q9GSL6", "Q9GSN2", "Q9GSM7",
    "Q9GSM9", "Q9GSM8", "Q9GSL7", "Q9GSL9", "A0A0L0BYT7", "A0A0L0C1F0",
    "A0A0L0BYH9", "A0A0L0BYA8", "A0A0U3AZ32", "D9J150", "A0A0L0BQ22",
    "A0A0L0BQK7", "A0A0L0C4Y9",
}


def _primary_acc(s):
    return str(s).split(";")[0]


def build_dataframe() -> pd.DataFrame:
    df = pd.read_excel(SRC, sheet_name="proteinGroups", header=2)
    df = df.rename(columns={
        "Whole-body extracts -Intensity ": "WE",
        "Excretory-secretory (ES) products -Intensity ": "ES",
    })
    return df


def render(df: pd.DataFrame, dpi: int, out: Path,
           label_fs: float, tick_fs: float) -> None:
    mpl.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 11,
        "axes.labelsize": label_fs,
        "axes.labelweight": "bold",
        "axes.linewidth": 1.2,
        "axes.edgecolor": "#222222",
        "xtick.labelsize": tick_fs,
        "ytick.labelsize": tick_fs,
        "xtick.major.width": 1.1,
        "ytick.major.width": 1.1,
        "legend.frameon": False,
        "savefig.bbox": "tight",
        "savefig.facecolor": "white",
        "axes.spines.top": False,
        "axes.spines.right": False,
    })

    fig, ax = plt.subplots(figsize=(8, 7.6))
    both = df[(df["WE"] > 0) & (df["ES"] > 0)].copy()
    x = np.log10(both["WE"]); y = np.log10(both["ES"])
    xy = np.vstack([x, y])
    dens = gaussian_kde(xy)(xy)
    order = dens.argsort()
    sc = ax.scatter(x.values[order], y.values[order], c=dens[order],
                    cmap="viridis", s=12, edgecolors="none", alpha=0.85)
    cb = plt.colorbar(sc, ax=ax, fraction=0.04, pad=0.02)

    lo = min(x.min(), y.min()) - 0.2
    hi = max(x.max(), y.max()) + 0.2
    ax.plot([lo, hi], [lo, hi], color="#222", lw=1.0, ls="-", label="y = x")
    ax.plot([lo, hi], [lo + 1, hi + 1], color=PALETTE["Up in ES"],
            lw=0.9, ls="--", label="2× enriched in ES")
    ax.plot([lo, hi], [lo - 1, hi - 1], color=PALETTE["Up in WE"],
            lw=0.9, ls="--", label="2× enriched in WE")

    es_only = df[(df["WE"] == 0) & (df["ES"] > 0)]
    we_only = df[(df["ES"] == 0) & (df["WE"] > 0)]
    ax.scatter(np.full(len(es_only), lo + 0.15), np.log10(es_only["ES"]),
               marker="<", s=22, c=PALETTE["ES-exclusive"], alpha=0.7,
               label=f"ES-exclusive (n={len(es_only)})", clip_on=False)
    ax.scatter(np.log10(we_only["WE"]), np.full(len(we_only), lo + 0.15),
               marker="v", s=22, c=PALETTE["WE-exclusive"], alpha=0.7,
               label=f"WE-exclusive (n={len(we_only)})", clip_on=False)

    both["_acc"] = both["Protein IDs"].apply(_primary_acc)
    amp_hits = both[both["_acc"].isin(AMP_PRIORITY)]
    if len(amp_hits) > 0:
        ax.scatter(np.log10(amp_hits["WE"]), np.log10(amp_hits["ES"]),
                   s=64, facecolors="none", edgecolors="#C68C2D",
                   linewidths=1.4, zorder=5,
                   label=f"Table 1 AMP candidate (n={len(amp_hits)})")

    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    ax.set_xlabel(r"log$_{10}$ intensity — Whole-body Extracts (WE)")
    ax.set_ylabel(r"log$_{10}$ intensity — Excretory–Secretory (ES)")
    # NOTE: figure title heading intentionally omitted (per request).
    # NOTE: the below-graph legend (line-style key + n-count statistics)
    # is intentionally omitted (per request); these details belong in the
    # figure caption rather than on the artwork.
    ax.tick_params(axis="both", labelsize=tick_fs)
    cb.set_label("Point density", fontsize=label_fs)
    cb.ax.tick_params(labelsize=tick_fs)
    ax.grid(alpha=0.15, ls=":")
    fig.tight_layout()
    fig.savefig(out, dpi=dpi, facecolor="white")
    plt.close(fig)

    from PIL import Image
    im = Image.open(out)
    print(f"wrote {out}")
    print(f"  dpi={im.info.get('dpi')}  px={im.size}  "
          f"print={im.size[0]/dpi:.2f}x{im.size[1]/dpi:.2f} in")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dpi", type=int, default=600)
    ap.add_argument("--out", type=Path,
                    default=FIG_DIR / "Fig2_scatter_600dpi.png")
    ap.add_argument("--label-fontsize", type=float, default=15)
    ap.add_argument("--tick-fontsize", type=float, default=13)
    args = ap.parse_args()
    render(build_dataframe(), args.dpi, args.out,
           args.label_fontsize, args.tick_fontsize)


if __name__ == "__main__":
    main()
