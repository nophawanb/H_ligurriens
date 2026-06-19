"""
High-resolution re-render of Fig 4 (HaCaT MTT cytocompatibility,
*Hemipyrellia ligurriens* larval extracts).

Faithful copy of `scripts/make_fig4_hacat_mtt.py` (identical data,
statistics, styling and layout) with three deliberate differences,
matching the established hi-res criteria used for Fig 1 / Fig 2 / Fig 3:

  1. Renders at a configurable high DPI (default 600).
  2. OMITS the figure heading — both the master `suptitle` and the
     subtitle methods line beneath it (these belong in the caption).
     The WE/ES series legend and the significance-key footer are KEPT,
     as they are required to read the plot (not redundant overlays).
  3. Uses repo-relative paths (the original hard-codes another machine's
     absolute path) and writes a SEPARATE output file
     (`figures/Fig4_hacat_mtt.tif`), never touching the 320-dpi
     `figures/Fig4_hacat_mtt.png` consumed downstream.

Usage:
    python3 scripts/make_fig4_hacat_mtt_hires.py [--dpi 600] [--out PATH]
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent


# ===================== Data (sheet `HL`) =====================
# Each row: (concentration µg/mL, [rep1, rep2, rep3], reported_p)
WE = [
    (0.1, [103.034555, 99.769364, 105.249952], 0.16696),
    (1.0, [110.629714, 113.165481, 114.109245], 0.000262),
    (10,  [140.115769, 137.286181, 128.102683], 0.000633),
    (100, [214.225574, 219.065924, 191.507431], 0.000218),
]
ES = [
    (0.1, [102.81719,  78.981019, 72.148860], 0.17400),
    (1.0, [103.74255,  77.782218, 70.548219], 0.18809),
    (10,  [96.40140,   75.704296, 76.270508], 0.06479),
    (100, [86.98334,   63.516484, 62.224890], 0.02244),
]


def sig_marker(p):
    if p < 0.0001: return "****"
    if p < 0.001:  return "***"
    if p < 0.01:   return "**"
    if p < 0.05:   return "*"
    return "ns"


# ===================== Design tokens =====================
C_WE       = "#1F77B4"
C_ES       = "#D72660"
INK        = "#1F2D3D"
SUB_INK    = "#5D6675"
DIVIDER    = "#D8DCE3"
MUTED      = "#8B93A1"


def series_arrays(data):
    concs = np.array([c for c, _, _ in data], dtype=float)
    means = np.array([np.mean(reps) for _, reps, _ in data])
    sds   = np.array([np.std(reps, ddof=1) for _, reps, _ in data])
    reps  = [reps for _, reps, _ in data]
    ps    = [p for _, _, p in data]
    return concs, means, sds, reps, ps


def draw_series(ax, x, means, sds, reps, ps, color, marker,
                sig_va, sig_offset):
    ax.plot(x, means, "-", color=color, lw=1.8, zorder=3)
    ax.errorbar(x, means, yerr=sds, fmt="none", ecolor=color,
                elinewidth=1.2, capsize=4, capthick=1.2, zorder=4)
    ax.plot(x, means, marker, markersize=8.5,
            markerfacecolor=color, markeredgecolor="white",
            markeredgewidth=1.4, zorder=5)
    for xi, rs in zip(x, reps):
        ax.plot([xi] * len(rs), rs, "o", markersize=4.6,
                markerfacecolor="white", markeredgecolor=color,
                markeredgewidth=1.2, zorder=4)
    for xi, mi, si, pi in zip(x, means, sds, ps):
        if pi >= 0.05:
            label = "ns"; fs = 9.5; col = SUB_INK
        else:
            label = sig_marker(pi); fs = 12; col = color
        ax.text(xi, mi + sig_offset * (si + 6), label,
                ha="center", va=sig_va, fontsize=fs,
                fontweight="bold", color=col, zorder=6)


def render(dpi: int, out: Path) -> None:
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "axes.unicode_minus": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })

    fig = plt.figure(figsize=(7.5, 4.6), dpi=dpi)
    fig.patch.set_facecolor("white")

    # NOTE: master suptitle + subtitle methods line intentionally
    # omitted (per request); these belong in the figure caption.
    ax = fig.add_axes([0.09, 0.20, 0.86, 0.68])

    we_x, we_m, we_s, we_r, we_p = series_arrays(WE)
    es_x, es_m, es_s, es_r, es_p = series_arrays(ES)

    ax.axhline(100, color="#1F2D3D", lw=1.0, ls="--", zorder=1)
    ax.axhline(80,  color="#B85569", lw=0.9, ls=":",  zorder=1, alpha=0.85)

    draw_series(ax, we_x, we_m, we_s, we_r, we_p, C_WE, "o",
                sig_va="bottom", sig_offset=+1.0)
    draw_series(ax, es_x, es_m, es_s, es_r, es_p, C_ES, "s",
                sig_va="top",    sig_offset=-1.0)

    LABEL_X = 200
    ax.text(LABEL_X, 100, "100 % untreated baseline",
            ha="center", va="center", fontsize=9, color=SUB_INK,
            fontstyle="italic",
            bbox=dict(boxstyle="round,pad=0.30",
                      facecolor="white", edgecolor="none"), zorder=2)
    ax.text(LABEL_X, 80, "80 % cytotoxicity threshold",
            ha="center", va="center", fontsize=9, color="#B85569",
            fontstyle="italic",
            bbox=dict(boxstyle="round,pad=0.30",
                      facecolor="white", edgecolor="none"), zorder=2)

    ax.set_xscale("log")
    ax.set_xlim(0.06, 280)
    ax.set_xticks([0.1, 1, 10, 100])
    ax.set_xticklabels(["0.1", "1", "10", "100"], fontsize=11, color=INK)
    ax.set_xlabel(r"Extract concentration (µg mL$^{-1}$)  ·  log$_{10}$ scale",
                  fontsize=11.5, color=INK, labelpad=8)
    ax.set_ylabel("HaCaT viability (% of untreated control)",
                  fontsize=11.5, color=INK, labelpad=6)
    ax.set_ylim(40, 250)
    ax.set_yticks([50, 75, 100, 125, 150, 175, 200, 225, 250])
    ax.tick_params(axis="x", which="both", length=3,
                   color=MUTED, labelcolor=INK)
    ax.tick_params(axis="y", length=3, color=MUTED, labelcolor=SUB_INK)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.spines["left"].set_color(MUTED); ax.spines["left"].set_linewidth(0.8)
    ax.spines["bottom"].set_color(MUTED); ax.spines["bottom"].set_linewidth(0.8)
    ax.grid(axis="y", color=DIVIDER, alpha=0.4, lw=0.6, zorder=0)
    ax.set_axisbelow(True)

    # Legend — centred below the plot (kept: needed to read the plot).
    ax.plot([], [], "-o", color=C_WE, markersize=8,
            markerfacecolor=C_WE, markeredgecolor="white",
            markeredgewidth=1.2, lw=1.8, label="WE (whole-body extract)")
    ax.plot([], [], "-s", color=C_ES, markersize=8,
            markerfacecolor=C_ES, markeredgecolor="white",
            markeredgewidth=1.2, lw=1.8, label="ES (excretory–secretory)")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.16),
              ncol=2, frameon=False, fontsize=10.5,
              handlelength=2.4, handletextpad=0.6, columnspacing=2.5)

    # Significance-key footer (kept: explains on-plot markers). Only the
    # markers actually present in this figure are listed; ** (P<0.01) and
    # **** (P<0.0001) do not occur here and are omitted per request.
    ax.text(0.5, -0.27,
            "*  P < 0.05      ***  P < 0.001      ns: non-significant      "
            "Open circles = individual biological-replicate values (n = 3)",
            transform=ax.transAxes, ha="center", va="top",
            fontsize=9.2, color=SUB_INK, fontstyle="italic")

    out.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out, dpi=dpi, bbox_inches="tight",
                facecolor="white", edgecolor="none",
                pil_kwargs=({"compression": "tiff_lzw"}
                if str(out).lower().endswith((".tif", ".tiff")) else {}))
    plt.close()

    from PIL import Image
    im = Image.open(out)
    print(f"Wrote {out} ({out.stat().st_size/1024:.1f} KB)")
    print(f"  dpi={im.info.get('dpi')}  px={im.size}  "
          f"print={im.size[0]/dpi:.2f}x{im.size[1]/dpi:.2f} in")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dpi", type=int, default=300)
    ap.add_argument("--out", type=Path,
                    default=ROOT / "figures" / "Fig6.tif")
    args = ap.parse_args()
    render(args.dpi, args.out)


if __name__ == "__main__":
    main()
