"""
High-resolution re-render of Fig 1 (total-proteome overview composite,
4 panels A–D) of the *Hemipyrellia ligurriens* larval extracts.

Faithful copy of `scripts/make_fig1_proteome.py` (identical data,
computations, layout and styling) with deliberate differences matching
the established hi-res criteria used across the figure set:

  1. Renders at a configurable high DPI (default 600).
  2. Each panel keeps ONLY its letter label (A / B / C / D). The
     figure-level master heading (suptitle + methods subtitle) and the
     four per-panel descriptive titles ("Compartment partition …",
     "Annotation status …", "Volcano-style plot …", "Pairwise scatter …")
     are omitted; these belong in the figure caption. In-panel data
     annotations (Venn counts, union-n, bar %, legend strip) are KEPT.
  3. Panels C and D embed HEADINGLESS source renders so each panel is
     itself clean:
        C ← figures/_embed_volcano_noheading.tif
            (make_fig1_volcano_hires.py --no-heading)
        D ← figures/Fig2_scatter_600dpi.tif
            (make_fig2_scatter_hires.py — already headingless / no legend)
  4. Uses repo-relative paths (the original hard-codes another machine's
     absolute path) and writes a SEPARATE output file
     (`figures/Fig1_proteome_overview_600dpi.tif`), never touching the
     320-dpi original consumed downstream.

Usage:
    python3 scripts/make_fig1_proteome_hires.py [--dpi 600] [--out PATH]
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import openpyxl
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.patheffects as pe
from matplotlib.patches import Rectangle
from matplotlib_venn import venn2, venn2_circles

ROOT = Path(__file__).resolve().parent.parent
SRC  = ROOT / "Supplementary.xlsx"
# Embedded panel sources: headingless high-res renders.
PANEL_C_SRC = ROOT / "figures" / "_embed_volcano_noheading.tif"
PANEL_D_SRC = ROOT / "figures" / "Fig2_scatter_600dpi.tif"

# Fallbacks to the original 320-dpi renders if the headingless sources
# are not present (keeps the script runnable standalone).
PANEL_C_FALLBACK = ROOT / "figures" / "Fig1_volcano.png"
PANEL_D_FALLBACK = ROOT / "figures" / "Fig2_scatter.png"


# =================== 1. Data ===================
def is_uncharacterised(name):
    if not isinstance(name, str):
        return True
    n = name.lower()
    return "uncharacterized protein" in n or "uncharacterised protein" in n


def compute_counts():
    wb = openpyxl.load_workbook(SRC, read_only=True, data_only=True)
    ws = wb["proteinGroups"]
    we_total = es_total = 0
    we_char = we_unchar = es_char = es_unchar = 0
    both = we_only = es_only = 0
    for row in ws.iter_rows(min_row=4, values_only=True):
        if not row or row[0] is None:
            continue
        name = row[1]
        try:
            we_int = float(row[9] or 0) if row[9] not in (None, "") else 0.0
            es_int = float(row[10] or 0) if row[10] not in (None, "") else 0.0
        except (TypeError, ValueError):
            we_int = es_int = 0.0
        in_we = we_int > 0
        in_es = es_int > 0
        if not (in_we or in_es):
            continue
        unchar = is_uncharacterised(name)
        if in_we:
            we_total += 1
            if unchar: we_unchar += 1
            else:      we_char += 1
        if in_es:
            es_total += 1
            if unchar: es_unchar += 1
            else:      es_char += 1
        if in_we and in_es:
            both += 1
        elif in_we:
            we_only += 1
        elif in_es:
            es_only += 1
    return dict(we_total=we_total, es_total=es_total,
                we_char=we_char, we_unchar=we_unchar,
                es_char=es_char, es_unchar=es_unchar,
                both=both, we_only=we_only, es_only=es_only)


def pct(part, whole):
    return 100.0 * part / whole if whole else 0.0


# =================== 2. Design tokens ===================
C_WE       = "#1F77B4"
C_ES       = "#D72660"
C_SHARED   = "#6A3FA0"
C_CHAR     = "#2E5F7A"
C_UNCHAR   = "#C9B07A"
INK        = "#1F2D3D"
SUB_INK    = "#5D6675"
TITLE_NAVY = "#0F1E33"
DIVIDER    = "#D8DCE3"

SUBTITLE_FS      = 10
PANEL_LABEL_FS   = 18
VENN_COUNT_FS    = 14
VENN_PCT_FS      = 9
VENN_SETLABEL_FS = 11
BAR_VAL_FS       = 10
TICK_FS          = 10
AXIS_FS          = 10.5
LEGEND_FS        = 11


# =================== 3. Panel renderers ===================
def venn_partition(ax, we_only_n, both_n, es_only_n,
                   we_only_pct, both_pct, es_only_pct, total_label):
    v = venn2(subsets=(we_only_n, es_only_n, both_n),
              set_labels=("", ""), ax=ax)
    if v.get_patch_by_id("10"):
        p = v.get_patch_by_id("10")
        p.set_color(C_WE); p.set_alpha(0.78); p.set_edgecolor("none")
    if v.get_patch_by_id("01"):
        p = v.get_patch_by_id("01")
        p.set_color(C_ES); p.set_alpha(0.78); p.set_edgecolor("none")
    if v.get_patch_by_id("11"):
        p = v.get_patch_by_id("11")
        p.set_color(C_SHARED); p.set_alpha(0.55); p.set_edgecolor("none")

    centers = v.centers
    radii   = v.radii
    c0, c1 = centers[0], centers[1]
    cx_we = getattr(c0, "x", c0[0] if hasattr(c0, "__getitem__") else c0)
    cy_we = getattr(c0, "y", c0[1] if hasattr(c0, "__getitem__") else c0)
    cx_es = getattr(c1, "x", c1[0] if hasattr(c1, "__getitem__") else c1)
    cy_es = getattr(c1, "y", c1[1] if hasattr(c1, "__getitem__") else c1)
    r_we, r_es = float(radii[0]), float(radii[1])
    y_top = max(cy_we + r_we, cy_es + r_es)
    y_bot = min(cy_we - r_we, cy_es - r_es)
    x_lft = cx_we - r_we
    x_rgt = cx_es + r_es
    span_x = x_rgt - x_lft
    span_y = y_top - y_bot

    text_stroke = [pe.withStroke(linewidth=2.4, foreground="#1A2330")]
    pct_stroke  = [pe.withStroke(linewidth=2.0, foreground="white")]
    for sid in ("10", "01", "11"):
        lbl = v.get_label_by_id(sid)
        if lbl is None:
            continue
        x, y = lbl.get_position()
        lbl.set_visible(False)
        count_n = {"10": we_only_n, "01": es_only_n, "11": both_n}[sid]
        pct_v   = {"10": we_only_pct, "01": es_only_pct, "11": both_pct}[sid]
        t1 = ax.text(x, y + 0.045 * span_y, f"{count_n:,}",
                     ha="center", va="center",
                     fontsize=VENN_COUNT_FS, fontweight="bold",
                     color="white", zorder=5)
        t1.set_path_effects(text_stroke)
        t2 = ax.text(x, y - 0.055 * span_y,
                     f"{pct_v:.1f}%" + (" shared" if sid == "11" else ""),
                     ha="center", va="center",
                     fontsize=VENN_PCT_FS, color="#1A2330", zorder=5,
                     fontweight="bold")
        t2.set_path_effects(pct_stroke)

    venn2_circles(subsets=(we_only_n, es_only_n, both_n),
                  ax=ax, linewidth=1.0, color="#4A5260")

    pad_x = 0.40 * span_x
    pad_y_top = 0.22 * span_y
    pad_y_bot = 0.18 * span_y
    ax.set_xlim(x_lft - pad_x, x_rgt + pad_x)
    ax.set_ylim(y_bot - pad_y_bot, y_top + pad_y_top)

    # Panel descriptive title omitted (per request) — letter label only.
    # Data subtitle (union n) kept as an in-panel annotation.
    total_n = we_only_n + both_n + es_only_n
    ax.text(0.5, 0.94, f"union n = {total_n:,}   ·   {total_label}",
            ha="center", va="top",
            fontsize=SUBTITLE_FS, color=SUB_INK, style="italic",
            transform=ax.transAxes)

    ax.text(0.10, 0.05, "Whole-body\nExtracts (WE)",
            ha="center", va="bottom",
            fontsize=VENN_SETLABEL_FS, color=INK, fontweight="bold",
            linespacing=1.20, transform=ax.transAxes)
    ax.text(0.90, 0.05, "Excretory–\nSecretory (ES)",
            ha="center", va="bottom",
            fontsize=VENN_SETLABEL_FS, color=INK, fontweight="bold",
            linespacing=1.20, transform=ax.transAxes)

    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)


def composition_bars(ax, fractions, char_pcts, unchar_pcts, totals):
    y_pos = np.arange(len(fractions))[::-1]
    bar_h = 0.55
    for i, frac in enumerate(fractions):
        ax.barh(y_pos[i], char_pcts[i], height=bar_h,
                color=C_CHAR, edgecolor="white", linewidth=1.2, zorder=2)
        ax.barh(y_pos[i], unchar_pcts[i], left=char_pcts[i], height=bar_h,
                color=C_UNCHAR, edgecolor="white", linewidth=1.2, zorder=2)
        if char_pcts[i] > 10:
            ax.text(char_pcts[i] / 2, y_pos[i], f"{char_pcts[i]:.1f}%",
                    ha="center", va="center", fontsize=BAR_VAL_FS,
                    color="white", fontweight="bold", zorder=3)
        if unchar_pcts[i] > 8:
            ax.text(char_pcts[i] + unchar_pcts[i] / 2, y_pos[i],
                    f"{unchar_pcts[i]:.1f}%", ha="center", va="center",
                    fontsize=BAR_VAL_FS, color=INK, fontweight="bold",
                    zorder=3)
        ax.text(101, y_pos[i], f"n = {totals[i]:,}",
                ha="left", va="center", fontsize=BAR_VAL_FS, color=SUB_INK)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(fractions, fontsize=TICK_FS + 1,
                       color=INK, fontweight="bold")
    ax.set_xlim(0, 120)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xticklabels(["0%", "25%", "50%", "75%", "100%"], fontsize=TICK_FS)
    ax.set_xlabel("Share of detected proteins",
                  fontsize=AXIS_FS, color=SUB_INK, labelpad=4)
    ax.set_ylim(-0.6, len(fractions) - 0.4)
    # Panel descriptive title omitted (per request) — letter label only.
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.spines["left"].set_color(DIVIDER)
    ax.spines["bottom"].set_color(DIVIDER)
    ax.grid(axis="x", color=DIVIDER, alpha=0.6, lw=0.7, zorder=1)
    ax.set_axisbelow(True)


def embed_image(ax, src, fallback):
    use = src if src.exists() else fallback
    if not use.exists():
        ax.text(0.5, 0.5, f"missing:\n{src.name}",
                transform=ax.transAxes, ha="center", va="center",
                color="#C62828", fontsize=11)
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(False)
        return
    ax.imshow(mpimg.imread(use))
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    # Panel descriptive title omitted (per request) — letter label only.


# =================== 4. Build figure ===================
def render(dpi: int, out: Path) -> None:
    d = compute_counts()
    union_full = d["we_only"] + d["both"] + d["es_only"]
    WE_CHAR_PCT   = pct(d["we_char"], d["we_total"])
    WE_UNCHAR_PCT = pct(d["we_unchar"], d["we_total"])
    ES_CHAR_PCT   = pct(d["es_char"], d["es_total"])
    ES_UNCHAR_PCT = pct(d["es_unchar"], d["es_total"])
    WE_ONLY_PCT_F = pct(d["we_only"], union_full)
    BOTH_PCT_F    = pct(d["both"], union_full)
    ES_ONLY_PCT_F = pct(d["es_only"], union_full)

    print("=== Total proteome (computed from Supplementary.xlsx) ===")
    print(f"WE n = {d['we_total']:,} | ES n = {d['es_total']:,} | "
          f"union n = {union_full:,}")
    print(f"Venn:  WE-only {d['we_only']:,} | shared {d['both']:,} | "
          f"ES-only {d['es_only']:,}")

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "axes.unicode_minus": False,
        "axes.edgecolor": DIVIDER,
        "axes.linewidth": 1.0,
        "xtick.color": INK,
        "ytick.color": INK,
        "xtick.major.size": 3,
        "ytick.major.size": 3,
        "axes.titleweight": "bold",
    })

    fig = plt.figure(figsize=(13.0, 29.0), dpi=dpi)
    fig.patch.set_facecolor("white")
    # Master heading (suptitle + methods subtitle) omitted (per request).

    gs = fig.add_gridspec(
        nrows=4, ncols=2,
        width_ratios=[1.0, 1.0],
        height_ratios=[0.55, 0.06, 2.40, 2.40],
        wspace=0.22, hspace=0.22,
        left=0.06, right=0.97, top=0.965, bottom=0.02,
    )

    ax_a      = fig.add_subplot(gs[0, 0])
    ax_b      = fig.add_subplot(gs[0, 1])
    ax_legend = fig.add_subplot(gs[1, :])
    ax_c      = fig.add_subplot(gs[2, :])
    ax_d      = fig.add_subplot(gs[3, :])

    venn_partition(ax_a, d["we_only"], d["both"], d["es_only"],
                   WE_ONLY_PCT_F, BOTH_PCT_F, ES_ONLY_PCT_F, "full proteome")
    composition_bars(ax_b, fractions=["WE", "ES"],
                     char_pcts=[WE_CHAR_PCT, ES_CHAR_PCT],
                     unchar_pcts=[WE_UNCHAR_PCT, ES_UNCHAR_PCT],
                     totals=[d["we_total"], d["es_total"]])
    embed_image(ax_c, PANEL_C_SRC, PANEL_C_FALLBACK)
    embed_image(ax_d, PANEL_D_SRC, PANEL_D_FALLBACK)

    for ax, letter in [(ax_a, "A"), (ax_b, "B"), (ax_c, "C"), (ax_d, "D")]:
        ax.text(-0.10, 1.07, letter, transform=ax.transAxes,
                ha="left", va="center",
                fontsize=PANEL_LABEL_FS, fontweight="bold", color=TITLE_NAVY)

    ax_legend.axis("off")
    swatch_y, swatch_h, swatch_w = 0.55, 0.45, 0.018
    items = [
        (0.06, C_WE,     "WE only"),
        (0.18, C_SHARED, "Shared (WE ∩ ES)"),
        (0.36, C_ES,     "ES only"),
        (0.54, C_CHAR,   "Characterized"),
        (0.72, C_UNCHAR, "Uncharacterized"),
    ]
    for x, color, label in items:
        ax_legend.add_patch(Rectangle(
            (x, swatch_y - swatch_h / 2), swatch_w, swatch_h,
            facecolor=color, edgecolor="none",
            transform=ax_legend.transAxes, clip_on=False))
        ax_legend.text(x + swatch_w + 0.007, swatch_y, label,
                       transform=ax_legend.transAxes,
                       ha="left", va="center", fontsize=LEGEND_FS, color=INK)

    out.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out, dpi=dpi, bbox_inches="tight",
                facecolor="white", edgecolor="none",
                pil_kwargs=({"compression": "tiff_lzw"}
                if str(out).lower().endswith((".tif", ".tiff")) else {}))
    plt.close()

    from PIL import Image
    im = Image.open(out)
    print(f"Wrote {out} ({out.stat().st_size/1024:.1f} KB)")
    print(f"  dpi={im.info.get('dpi')}  px={im.size}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dpi", type=int, default=600)
    ap.add_argument("--out", type=Path,
                    default=ROOT / "figures" / "Fig1_total_proteome_characterization_600dpi.tif")
    args = ap.parse_args()
    render(args.dpi, args.out)


if __name__ == "__main__":
    main()
