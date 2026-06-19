"""
High-resolution re-render of Fig 4 (fatty-acid composition of
*Hemipyrellia ligurriens* larval whole-body extract, GC-MS).

Faithful copy of `scripts/make_fig4_fatty_acid.py` (identical data,
Table-2 sanity checks, layout and styling) with three deliberate
differences, matching the established hi-res criteria used for
Fig 1 / Fig 2 / Fig 3 / Fig 4 (HaCaT MTT):

  1. Renders at a configurable high DPI (default 600).
  2. OMITS the figure heading — the master `suptitle` (GC-MS / AOAC
     line). The panel titles "A." / "B." and the SFA/MUFA/PUFA class
     legend are KEPT (panel structure + data, not redundant overlays).
  3. Uses repo-relative paths (the original hard-codes another machine's
     absolute path) and writes a SEPARATE output file
     (`figures/Fig4_fatty_acid_600dpi.tif`), never touching the 320-dpi
     `figures/Fig4_fatty_acid.png` consumed downstream.

Usage:
    python3 scripts/make_fig4_fatty_acid_hires.py [--dpi 600] [--out PATH]
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

ROOT = Path(__file__).resolve().parent.parent

# ===================== Data (manuscript Table 2) =====================
# Tuple: (display name, % w/w of total extract, class)
FA_DATA = [
    # SFA
    ("Lauric (C12:0)",                        0.04, "SFA"),
    ("Myristic (C14:0)",                      0.20, "SFA"),
    ("Pentadecanoic (C15:0)",                 0.06, "SFA"),
    ("Palmitic (C16:0)",                      5.90, "SFA"),
    ("Heptadecanoic (C17:0)",                 0.10, "SFA"),
    ("Stearic (C18:0)",                       1.70, "SFA"),
    ("Arachidic (C20:0)",                     0.08, "SFA"),
    ("Heneicosanoic (C21:0)",                 0.27, "SFA"),
    ("Behenic (C22:0)",                       0.04, "SFA"),
    ("Lignoceric (C24:0)",                    0.02, "SFA"),
    # MUFA
    ("Myristoleic (C14:1)",                   0.07, "MUFA"),
    ("Palmitoleic (C16:1n7)",                 2.30, "MUFA"),
    ("cis-10-Heptadecenoic (C17:1n10)",       0.08, "MUFA"),
    ("cis-9-Oleic (C18:1n9c)",                6.26, "MUFA"),
    ("cis-11-Eicosenoic (C20:1n9)",           0.03, "MUFA"),
    ("Erucic (C22:1n9)",                      0.03, "MUFA"),
    ("Nervonic (C24:1n9)",                    0.03, "MUFA"),
    # PUFA
    ("cis-9,12-Linoleic (C18:2n6c)",          4.62, "PUFA"),
    ("gamma-Linolenic (C18:3n6)",             0.27, "PUFA"),
    ("alpha-Linolenic (C18:3n3)",             0.27, "PUFA"),
    ("cis-11,14-Eicosadienoic (C20:2)",       0.10, "PUFA"),
    ("cis-8,11,14-Eicosatrienoic (C20:3n6)",  0.09, "PUFA"),
    ("cis-11,14,17-Eicosatrienoic (C20:3n3)", 0.02, "PUFA"),
    ("Arachidonic (C20:4n6)",                 2.96, "PUFA"),
    ("Eicosapentaenoic (C20:5n3)",            0.49, "PUFA"),
    ("Docosahexaenoic (C22:6n3)",             0.17, "PUFA"),
]


# ===================== Design tokens =====================
C_SFA  = "#1F77B4"
C_MUFA = "#E58A1A"
C_PUFA = "#1A8A5D"
INK      = "#1F2D3D"
SUB_INK  = "#5D6675"
TITLE    = "#0F1E33"
MUTED    = "#8B93A1"
C_NONFA  = "#D8DCE3"
CLS_COL  = {"SFA": C_SFA, "MUFA": C_MUFA, "PUFA": C_PUFA}


def compute_totals():
    class_totals = {"SFA": 0.0, "MUFA": 0.0, "PUFA": 0.0}
    for _, pct, cls in FA_DATA:
        class_totals[cls] += pct
    expected = {"SFA": 8.41, "MUFA": 8.80, "PUFA": 8.99}
    for cls, exp in expected.items():
        if abs(class_totals[cls] - exp) > 0.02:
            raise AssertionError(
                f"Class total mismatch for {cls}: "
                f"{class_totals[cls]:.2f} vs Table 2 {exp:.2f}")
    total_fa = sum(class_totals.values())
    print(f"SFA = {class_totals['SFA']:.2f} %  "
          f"MUFA = {class_totals['MUFA']:.2f} %  "
          f"PUFA = {class_totals['PUFA']:.2f} %  "
          f"Total = {total_fa:.2f} %  (n = {len(FA_DATA)} fatty acids)")
    return class_totals, total_fa


def render(dpi: int, out: Path) -> None:
    class_totals, total_fa = compute_totals()

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans"],
        "axes.unicode_minus": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })

    fig = plt.figure(figsize=(13.0, 9.5), dpi=dpi)
    fig.patch.set_facecolor("white")
    gs = gridspec.GridSpec(
        nrows=1, ncols=2,
        width_ratios=[0.60, 1.00],
        wspace=0.30,
        left=0.04, right=0.97, top=0.92, bottom=0.06,
    )

    # --------- Panel A: SFA / MUFA / PUFA + Non-FA donut ---------
    non_fa  = 100.0 - total_fa
    sizes   = [class_totals["SFA"], class_totals["MUFA"],
               class_totals["PUFA"], non_fa]
    labels  = ["SFA", "MUFA", "PUFA", "Non-FA fraction"]
    colors  = [CLS_COL["SFA"], CLS_COL["MUFA"], CLS_COL["PUFA"], C_NONFA]
    ax = fig.add_subplot(gs[0, 0])
    ax.pie(sizes, colors=colors, startangle=90, counterclock=False,
           wedgeprops=dict(width=0.36, edgecolor="white", linewidth=2.2))
    ax.text(0, 0.18, "Total identified FA",
            ha="center", va="center", fontsize=10, color=SUB_INK)
    ax.text(0, 0.00, f"{total_fa:.2f} %",
            ha="center", va="center", fontsize=20, fontweight="bold",
            color=TITLE)
    ax.text(0, -0.18, "of extract  (w/w)",
            ha="center", va="center", fontsize=9.5, color=SUB_INK)

    ax.set_xlim(-1.4, 1.4); ax.set_ylim(-1.4, 1.4)
    ax.set_aspect("equal")
    ax.set_title("A", fontsize=14, fontweight="bold", color=TITLE,
                 pad=10, loc="left")

    band_r = 1.0 - 0.36 / 2
    cum = 0.0
    for cls, val in zip(labels, sizes):
        frac_mid = (cum + val / 2) / sum(sizes)
        theta = 90 - 360 * frac_mid
        rad = np.deg2rad(theta)
        tx, ty = band_r * np.cos(rad), band_r * np.sin(rad)
        is_nonfa = (cls == "Non-FA fraction")
        text_col = SUB_INK if is_nonfa else "white"
        ax.text(tx, ty + 0.05, cls,
                ha="center", va="center", fontsize=10.5, color=text_col,
                fontweight="bold",
                fontstyle="italic" if is_nonfa else "normal", zorder=4)
        ax.text(tx, ty - 0.07, f"{val:.2f} %",
                ha="center", va="center", fontsize=10, color=text_col,
                fontweight="bold" if not is_nonfa else "normal", zorder=4)
        cum += val

    ax.text(0.5, -0.04,
            "% w/w of total extract  ·  identified FA vs non-FA residual",
            transform=ax.transAxes, ha="center", va="top",
            fontsize=9.0, color=SUB_INK, clip_on=False)

    # --------- Panel B: ALL 26 fatty acids ---------
    ax = fig.add_subplot(gs[0, 1])
    order_idx = sorted(
        range(len(FA_DATA)),
        key=lambda i: (["SFA", "MUFA", "PUFA"].index(FA_DATA[i][2]),
                       -FA_DATA[i][1]),
    )
    names_o = [FA_DATA[i][0] for i in order_idx]
    vals_o  = [FA_DATA[i][1] for i in order_idx]
    class_o = [FA_DATA[i][2] for i in order_idx]
    y       = np.arange(len(order_idx))[::-1]
    colors_o = [CLS_COL[c] for c in class_o]

    ax.barh(y, vals_o, color=colors_o, edgecolor="white", linewidth=0.6,
            zorder=2)
    for yi, v in zip(y, vals_o):
        ax.text(v + max(vals_o) * 0.015, yi, f"{v:.2f}",
                ha="left", va="center", fontsize=8.8, color=INK)
    for i in range(1, len(order_idx)):
        if class_o[i] != class_o[i - 1]:
            ax.axhline(y[i] + 0.5, color="#C9CFDB", lw=0.7, alpha=0.95,
                       zorder=1)

    ax.set_yticks(y)
    ax.set_yticklabels(names_o, fontsize=9.0, color=INK)
    ax.set_xlim(0, max(vals_o) * 1.18)
    ax.set_xlabel("% w/w of total extract", fontsize=10.5, color=INK,
                  labelpad=6)
    ax.tick_params(axis="x", length=3, color=MUTED, labelcolor=SUB_INK)
    ax.tick_params(axis="y", length=0)
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    ax.spines["bottom"].set_color(MUTED)
    ax.spines["bottom"].set_linewidth(0.8)
    ax.set_title("B", fontsize=14, fontweight="bold", color=TITLE,
                 pad=10, loc="left")

    leg_handles = [
        plt.Rectangle((0, 0), 1, 1, facecolor=CLS_COL["SFA"],
                      edgecolor="white", lw=0.6),
        plt.Rectangle((0, 0), 1, 1, facecolor=CLS_COL["MUFA"],
                      edgecolor="white", lw=0.6),
        plt.Rectangle((0, 0), 1, 1, facecolor=CLS_COL["PUFA"],
                      edgecolor="white", lw=0.6),
    ]
    ax.legend(
        leg_handles,
        [f"SFA ({class_totals['SFA']:.2f} %)",
         f"MUFA ({class_totals['MUFA']:.2f} %)",
         f"PUFA ({class_totals['PUFA']:.2f} %)"],
        loc="lower right", frameon=False, fontsize=9.5,
        handlelength=1.0, handleheight=1.0, labelcolor=INK,
    )

    # --------- Master heading intentionally omitted (per request) -----
    # The GC-MS / AOAC 996.06 master suptitle is not drawn; it belongs in
    # the figure caption. bbox_inches="tight" trims the freed whitespace.

    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=dpi, bbox_inches="tight",
                facecolor="white", edgecolor="none",
                pil_kwargs=({"compression": "tiff_lzw"}
                if str(out).lower().endswith((".tif", ".tiff")) else {}))
    plt.close(fig)

    from PIL import Image
    im = Image.open(out)
    print(f"Wrote {out} ({out.stat().st_size/1024:.1f} KB)")
    print(f"  dpi={im.info.get('dpi')}  px={im.size}  "
          f"print={im.size[0]/dpi:.2f}x{im.size[1]/dpi:.2f} in")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dpi", type=int, default=600)
    ap.add_argument("--out", type=Path,
                    default=ROOT / "figures" / "Fig3_fatty_acid_600dpi.tif")
    args = ap.parse_args()
    render(args.dpi, args.out)


if __name__ == "__main__":
    main()
