"""
High-resolution re-render of Fig 2 (< 20 kDa subset combined overview +
family composition) of the *Hemipyrellia ligurriens* larval extracts.

Faithful copy of `scripts/make_fig2_lowmw_combined.py` (identical data,
computations, layout and styling) with deliberate differences matching
the established hi-res criteria used across the figure set:

  1. Renders at a configurable high DPI (default 600).
  2. Each panel keeps ONLY its letter label (A / B / C). The figure-level
     master heading (suptitle + subtitle) and the three per-panel
     descriptive titles ("Compartment partition …", "Annotation status …",
     "Family composition × compartment partition") are omitted; these
     belong in the figure caption. In-panel data annotations (Venn counts,
     union-n, bar %, KDE median, lollipop n, legend strip) are KEPT.
  3. Uses repo-relative paths (the original hard-codes another machine's
     absolute path) and writes a SEPARATE output file
     (`figures/Fig2_lowmw_combined.tif`), never touching the
     320-dpi original consumed downstream.

Usage:
    python3 scripts/make_fig2_lowmw_combined_hires.py [--dpi 600] [--out PATH]
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import openpyxl
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.patches import Rectangle
from matplotlib_venn import venn2, venn2_circles
from scipy.stats import gaussian_kde

ROOT = Path(__file__).resolve().parent.parent
SRC  = ROOT / "Supplementary.xlsx"

LMW_CUTOFF = 21.0


# =================== 1. Data ===================
def is_uncharacterised(name):
    if not isinstance(name, str):
        return True
    n = name.lower()
    return "uncharacterized protein" in n or "uncharacterised protein" in n


def primary_acc(s):
    return str(s).split(";")[0]


TABLE1_FAMILIES = {
    "P86471": "Lucifensin", "A0A7D5FFX3": "Lucifensin",
    "P10891": "Phormicin",
    "A0A0L0C905": "Defensin (other)", "A0A0L0BRJ2": "Defensin (other)",
    "C0HJX9": "Diptericin", "P18684": "Diptericin-D",
    "A0A0L0C6T8": "Attacin", "A0A0L0BXA2": "Attacin",
    "A0A0L0C799": "Attacin", "A0A0L0C779": "Attacin",
    "A0A0L0BNM1": "Lysozyme", "C0HLB7": "Lysozyme",
    "D9J148": "Lysozyme", "D9J143": "Lysozyme", "A0A0L0CIT1": "Lysozyme",
    "A0A7D5FH37": "Protease (serine)",
    "Q25237": "Protease (serine)", "Q25231": "Protease (serine)",
    "Q9GSL6": "Protease (serine)", "Q9GSN2": "Protease (serine)",
    "Q9GSM7": "Protease (serine)", "Q9GSM9": "Protease (serine)",
    "Q9GSM8": "Protease (serine)", "Q9GSL7": "Protease (serine)",
    "Q9GSL9": "Protease (serine)",
    "A0A0L0BYT7": "Ser/Thr kinase", "A0A0L0C1F0": "Ser/Thr kinase",
    "A0A0L0BYH9": "Metallothionein", "A0A0L0BYA8": "Metallothionein",
    "A0A0U3AZ32": "Metallothionein",
    "D9J150": "Ferritin", "A0A0L0BQ22": "Cu-Zn SOD",
    "A0A0L0BQK7": "Insulin-like", "A0A0L0C4Y9": "Insulin-like",
}


def keyword_family(name):
    if not isinstance(name, str):
        return "Uncharacterized"
    n = name.lower()
    if "uncharacterized protein" in n or "uncharacterised protein" in n:
        return "Uncharacterized"
    if "lysozyme" in n:
        return "Lysozyme"
    if any(k in n for k in ["serine protease", "trypsin", "chymotrypsin",
                            "metalloprotease", "metallopeptidase",
                            "peptidase", "cathepsin"]):
        return "Protease (other)"
    if any(k in n for k in ["protease inhibitor", "serpin", "kazal"]):
        return "Protease inhibitor"
    if any(k in n for k in ["heat shock", "chaperonin", "chaperone"]):
        return "Heat-shock / chaperone"
    if any(k in n for k in ["ribosomal", "translation", "elongation factor"]):
        return "Ribosomal / translation"
    if any(k in n for k in ["actin", "tubulin", "myosin", "tropomyosin",
                            "troponin"]):
        return "Cytoskeleton / structural"
    if any(k in n for k in ["nadh", "atp synthase", "cytochrome",
                            "succinate dehydrogenase", "oxidase"]):
        return "Energy / electron transport"
    if any(k in n for k in ["transport", "carrier", "binding protein"]):
        return "Transport / binding"
    if any(k in n for k in ["histone", "chitin", "cuticle"]):
        return "Cuticle / chromatin"
    return "Other characterized"


ROW_ORDER = [
    ("Lucifensin",                  "AMP / immune-effector"),
    ("Phormicin",                   None),
    ("Defensin (other)",            None),
    ("Diptericin",                  None),
    ("Diptericin-D",                None),
    ("Attacin",                     None),
    ("Lysozyme",                    None),
    ("Protease (serine)",           "Proteolytic"),
    ("Ser/Thr kinase",              None),
    ("Protease (other)",            None),
    ("Protease inhibitor",          None),
    ("Metallothionein",             "Antioxidant / metal-binding"),
    ("Ferritin",                    None),
    ("Cu-Zn SOD",                   None),
    ("Insulin-like",                "Signalling"),
    ("Heat-shock / chaperone",      "Other characterized"),
    ("Ribosomal / translation",     None),
    ("Cytoskeleton / structural",   None),
    ("Energy / electron transport", None),
    ("Transport / binding",         None),
    ("Cuticle / chromatin",         None),
    ("Other characterized",         None),
    ("Uncharacterized",             "Uncharacterized"),
]
ROW_LABELS = [r[0] for r in ROW_ORDER]
COMP_LABELS = ["WE-only", "Shared", "ES-only"]


def pct(part, whole):
    return 100.0 * part / whole if whole else 0.0


def safe_log2_es_over_we(we_n, es_n):
    if we_n == 0 and es_n == 0:
        return np.nan
    return float(np.log2((es_n + 0.5) / (we_n + 0.5)))


def compute():
    wb = openpyxl.load_workbook(SRC, read_only=True, data_only=True)
    ws = wb["proteinGroups"]
    lmw_we = lmw_es = 0
    lmw_we_char = lmw_we_unchar = lmw_es_char = lmw_es_unchar = 0
    lmw_both = lmw_we_only = lmw_es_only = 0
    matrix = {r: {c: 0 for c in COMP_LABELS} for r in ROW_LABELS}

    for row in ws.iter_rows(min_row=4, values_only=True):
        if not row or row[0] is None:
            continue
        acc_field, name = row[0], row[1]
        try:
            we_int = float(row[9] or 0) if row[9] not in (None, "") else 0.0
            es_int = float(row[10] or 0) if row[10] not in (None, "") else 0.0
        except (TypeError, ValueError):
            we_int = es_int = 0.0
        try:
            mw = float(row[5]) if row[5] not in (None, "") else None
        except (TypeError, ValueError):
            mw = None
        if mw is None or mw >= LMW_CUTOFF:
            continue
        in_we = we_int > 0
        in_es = es_int > 0
        if not (in_we or in_es):
            continue
        unchar = is_uncharacterised(name)
        if in_we:
            lmw_we += 1
            if unchar: lmw_we_unchar += 1
            else:      lmw_we_char += 1
        if in_es:
            lmw_es += 1
            if unchar: lmw_es_unchar += 1
            else:      lmw_es_char += 1
        if in_we and in_es:
            lmw_both += 1
        elif in_we:
            lmw_we_only += 1
        elif in_es:
            lmw_es_only += 1
        acc = primary_acc(acc_field)
        fam = TABLE1_FAMILIES.get(acc) or keyword_family(name)
        if fam not in matrix:
            fam = "Other characterized"
        if we_int > 0 and es_int > 0:
            matrix[fam]["Shared"] += 1
        elif we_int > 0:
            matrix[fam]["WE-only"] += 1
        else:
            matrix[fam]["ES-only"] += 1

    return dict(
        lmw_we=lmw_we, lmw_es=lmw_es,
        lmw_we_char=lmw_we_char, lmw_we_unchar=lmw_we_unchar,
        lmw_es_char=lmw_es_char, lmw_es_unchar=lmw_es_unchar,
        lmw_both=lmw_both, lmw_we_only=lmw_we_only, lmw_es_only=lmw_es_only,
        matrix=matrix,
    )


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
MUTED      = "#8B93A1"
HAIR       = "#D8DCE3"


# =================== 3. Panel renderers ===================
def venn_panel(ax, we_only_n, both_n, es_only_n,
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
        if lbl is None: continue
        x, y = lbl.get_position()
        lbl.set_visible(False)
        count_n = {"10": we_only_n, "01": es_only_n, "11": both_n}[sid]
        pct_v   = {"10": we_only_pct, "01": es_only_pct, "11": both_pct}[sid]
        t1 = ax.text(x, y + 0.075 * span_y, f"{count_n:,}",
                     ha="center", va="center", fontsize=12,
                     fontweight="bold", color="white", zorder=5)
        t1.set_path_effects(text_stroke)
        t2 = ax.text(x, y - 0.085 * span_y,
                     f"{pct_v:.1f}%",
                     ha="center", va="center", fontsize=8,
                     color="#1A2330", zorder=5, fontweight="bold")
        t2.set_path_effects(pct_stroke)

    venn2_circles(subsets=(we_only_n, es_only_n, both_n),
                  ax=ax, linewidth=1.0, color="#4A5260")

    pad_x = 0.40 * span_x
    pad_y_top = 0.22 * span_y
    pad_y_bot = 0.18 * span_y
    ax.set_xlim(x_lft - pad_x, x_rgt + pad_x)
    ax.set_ylim(y_bot - pad_y_bot, y_top + pad_y_top)

    total_n = we_only_n + both_n + es_only_n
    ax.text(0.5, 0.99, f"union n = {total_n:,}   ·   {total_label}",
            ha="center", va="top", fontsize=10, color=SUB_INK,
            style="italic", transform=ax.transAxes)
    ax.text(0.10, 0.04, "Whole-body\nExtracts (WE)",
            ha="center", va="bottom", fontsize=10.5, color=INK,
            fontweight="bold", linespacing=1.20, transform=ax.transAxes)
    ax.text(0.90, 0.04, "Excretory–\nSecretory (ES)",
            ha="center", va="bottom", fontsize=10.5, color=INK,
            fontweight="bold", linespacing=1.20, transform=ax.transAxes)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)


def annot_bars(ax, fractions, char_pcts, unchar_pcts, totals):
    y_pos = np.arange(len(fractions))[::-1]
    bar_h = 0.55
    for i in range(len(fractions)):
        ax.barh(y_pos[i], char_pcts[i], height=bar_h,
                color=C_CHAR, edgecolor="white", linewidth=1.2, zorder=2)
        ax.barh(y_pos[i], unchar_pcts[i], left=char_pcts[i], height=bar_h,
                color=C_UNCHAR, edgecolor="white", linewidth=1.2, zorder=2)
        if char_pcts[i] > 10:
            ax.text(char_pcts[i] / 2, y_pos[i], f"{char_pcts[i]:.1f}%",
                    ha="center", va="center", fontsize=10,
                    color="white", fontweight="bold", zorder=3)
        if unchar_pcts[i] > 8:
            ax.text(char_pcts[i] + unchar_pcts[i] / 2, y_pos[i],
                    f"{unchar_pcts[i]:.1f}%", ha="center", va="center",
                    fontsize=10, color=INK, fontweight="bold", zorder=3)
        ax.text(101, y_pos[i], f"n = {totals[i]:,}",
                ha="left", va="center", fontsize=10, color=SUB_INK)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(fractions, fontsize=11, color=INK, fontweight="bold")
    ax.set_xlim(0, 120)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xticklabels(["0%", "25%", "50%", "75%", "100%"], fontsize=10)
    ax.set_xlabel("Share of detected proteins",
                  fontsize=10.5, color=SUB_INK, labelpad=4)
    ax.set_ylim(-0.6, len(fractions) - 0.4)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.spines["left"].set_color(DIVIDER)
    ax.spines["bottom"].set_color(DIVIDER)
    ax.grid(axis="x", color=DIVIDER, alpha=0.6, lw=0.7, zorder=1)
    ax.set_axisbelow(True)


# =================== 4. Build figure ===================
def render(dpi: int, out: Path) -> None:
    d = compute()
    matrix = d["matrix"]
    union_lmw = d["lmw_we_only"] + d["lmw_both"] + d["lmw_es_only"]
    LMW_WE_CHAR_PCT   = pct(d["lmw_we_char"], d["lmw_we"])
    LMW_WE_UNCHAR_PCT = pct(d["lmw_we_unchar"], d["lmw_we"])
    LMW_ES_CHAR_PCT   = pct(d["lmw_es_char"], d["lmw_es"])
    LMW_ES_UNCHAR_PCT = pct(d["lmw_es_unchar"], d["lmw_es"])
    WE_ONLY_PCT_L = pct(d["lmw_we_only"], union_lmw)
    BOTH_PCT_L    = pct(d["lmw_both"], union_lmw)
    ES_ONLY_PCT_L = pct(d["lmw_es_only"], union_lmw)

    M = np.array([[matrix[fam][c] for c in COMP_LABELS] for fam in ROW_LABELS],
                 dtype=float)
    totals = M.sum(axis=1)
    with np.errstate(invalid="ignore"):
        fam_pct = np.where(totals[:, None] > 0,
                           100.0 * M / totals[:, None], 0.0)
    log2_ratio = np.array([
        safe_log2_es_over_we(M[i, 0], M[i, 2]) for i in range(M.shape[0])
    ])

    print("=== < 20 kDa subset (overview) ===")
    print(f"WE n = {d['lmw_we']:,} | ES n = {d['lmw_es']:,} | "
          f"union n = {union_lmw:,}")
    print(f"Venn: WE-only {d['lmw_we_only']:,} | shared {d['lmw_both']:,} | "
          f"ES-only {d['lmw_es_only']:,}")
    print(f"Family-row total: {int(M.sum())}")

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "axes.unicode_minus": False,
        "axes.edgecolor": DIVIDER,
        "axes.linewidth": 1.0,
        "xtick.color": INK,
        "ytick.color": INK,
        "axes.titleweight": "bold",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })

    fig = plt.figure(figsize=(6.9, 8.75), dpi=dpi)
    fig.patch.set_facecolor("white")
    # Master heading (suptitle + subtitle) omitted (per request).

    gs_master = fig.add_gridspec(
        nrows=3, ncols=1,
        height_ratios=[0.55, 1.40, 0.07],
        hspace=0.18,
        left=0.06, right=0.97, top=0.95, bottom=0.04,
    )
    gs_overview = gs_master[0].subgridspec(
        nrows=1, ncols=2, width_ratios=[1.0, 1.0], wspace=0.30,
    )
    ax_a = fig.add_subplot(gs_overview[0, 0])
    ax_b = fig.add_subplot(gs_overview[0, 1])

    gs_main = gs_master[1].subgridspec(
        nrows=2, ncols=2,
        width_ratios=[1.0, 0.26],
        height_ratios=[0.14, 1.0],
        wspace=0.025, hspace=0.06,
    )
    ax_kde   = fig.add_subplot(gs_main[0, 0])
    ax_main  = fig.add_subplot(gs_main[1, 0])
    ax_right = fig.add_subplot(gs_main[1, 1])
    ax_legend = fig.add_subplot(gs_master[2])

    # --------- Panel A: Venn (descriptive title omitted) ---------
    venn_panel(ax_a, d["lmw_we_only"], d["lmw_both"], d["lmw_es_only"],
               WE_ONLY_PCT_L, BOTH_PCT_L, ES_ONLY_PCT_L, "MW < 21 kDa subset")

    # --------- Panel B: Annotation bars (descriptive title omitted) ---------
    annot_bars(ax_b, ["WE", "ES"],
               [LMW_WE_CHAR_PCT, LMW_ES_CHAR_PCT],
               [LMW_WE_UNCHAR_PCT, LMW_ES_UNCHAR_PCT],
               [d["lmw_we"], d["lmw_es"]])

    # --------- Panel C: family composition + KDE + lollipop ---------
    n_rows = len(ROW_LABELS)
    y_pos  = np.arange(n_rows)[::-1]
    bar_h  = 0.74
    GROUP_BG_TINT = {
        "AMP / immune-effector":       "#FDE9E9",
        "Proteolytic":                 "#EAF0FA",
        "Antioxidant / metal-binding": "#EAF4EA",
        "Signalling":                  "#FAF1E1",
        "Other characterized":         "#F4F5F8",
        "Uncharacterized":             "#EEEEEE",
    }
    row_group_of = []
    current = None
    for _, g in ROW_ORDER:
        if g is not None:
            current = g
        row_group_of.append(current)

    ax = ax_main
    for i, g in enumerate(row_group_of):
        ax.axhspan(y_pos[i] - 0.5, y_pos[i] + 0.5, xmin=0, xmax=1,
                   facecolor=GROUP_BG_TINT.get(g, "white"),
                   edgecolor="none", alpha=0.45, zorder=0)
    ax.barh(y_pos, fam_pct[:, 0], left=0, height=bar_h,
            color=C_WE, edgecolor="white", linewidth=1.0, zorder=2)
    ax.barh(y_pos, fam_pct[:, 1], left=fam_pct[:, 0], height=bar_h,
            color=C_SHARED, edgecolor="white", linewidth=1.0, zorder=2)
    ax.barh(y_pos, fam_pct[:, 2], left=fam_pct[:, 0] + fam_pct[:, 1],
            height=bar_h, color=C_ES, edgecolor="white", linewidth=1.0,
            zorder=2)
    MIN_LBL = 7.0
    for i in range(n_rows):
        if M[i, 0] > 0 and fam_pct[i, 0] >= MIN_LBL:
            ax.text(fam_pct[i, 0] / 2, y_pos[i], f"{int(M[i, 0])}",
                    ha="center", va="center", fontsize=8.5,
                    color="white", fontweight="bold", zorder=3)
        if M[i, 1] > 0 and fam_pct[i, 1] >= MIN_LBL:
            ax.text(fam_pct[i, 0] + fam_pct[i, 1] / 2, y_pos[i],
                    f"{int(M[i, 1])}", ha="center", va="center",
                    fontsize=8.5, color="white", fontweight="bold", zorder=3)
        if M[i, 2] > 0 and fam_pct[i, 2] >= MIN_LBL:
            ax.text(fam_pct[i, 0] + fam_pct[i, 1] + fam_pct[i, 2] / 2,
                    y_pos[i], f"{int(M[i, 2])}", ha="center", va="center",
                    fontsize=8.5, color="white", fontweight="bold", zorder=3)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(ROW_LABELS, fontsize=10, color=INK)
    ax.set_xlim(0, 100)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xticklabels(["0", "25", "50", "75", "100"], fontsize=9.5,
                       color=SUB_INK)
    ax.set_xlabel("% within family (compartment partition)",
                  fontsize=10.5, color=INK, labelpad=6)
    ax.set_ylim(-0.6, n_rows - 0.4)
    ax.tick_params(axis="x", length=3, color=MUTED, labelcolor=SUB_INK)
    ax.tick_params(axis="y", length=0)
    ax.axvline(50, color=HAIR, lw=0.6, ls=":", zorder=1)
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    ax.spines["bottom"].set_color(MUTED)
    ax.spines["bottom"].set_linewidth(0.8)
    for idx, (_, group) in enumerate(ROW_ORDER):
        if group is None or idx == 0:
            continue
        ax.axhline(y_pos[idx] + 0.5, color="#C9CFDB", lw=0.6, alpha=0.9,
                   zorder=1)

    # Right: log2(ES/WE) lollipop
    ax = ax_right
    ax.set_xlim(-3.6, 3.6); ax.set_ylim(-0.6, n_rows - 0.4)
    ax.axvline(0,  color=MUTED, lw=0.8, zorder=1)
    ax.axvline(-1, color=HAIR, lw=0.5, ls=":", zorder=1)
    ax.axvline( 1, color=HAIR, lw=0.5, ls=":", zorder=1)
    for i in range(n_rows):
        r = log2_ratio[i]
        if np.isnan(r): continue
        r_clipped = max(-3.5, min(3.5, r))
        colour = C_ES if r > 0 else C_WE if r < 0 else MUTED
        ax.hlines(y_pos[i], 0, r_clipped, color=colour, lw=1.4, zorder=2)
        ax.plot(r_clipped, y_pos[i], "o", markersize=5.5, color=colour,
                markeredgecolor="white", markeredgewidth=0.8, zorder=3)
    for i in range(n_rows):
        ax.text(3.9, y_pos[i], f"n={int(totals[i]):,}",
                ha="left", va="center", fontsize=9, color=SUB_INK)
    ax.set_xticks([-2, 0, 2])
    ax.set_xticklabels(["−2", "0", "+2"], fontsize=9, color=SUB_INK)
    ax.set_xlabel(r"log$_2$(ES/WE)", fontsize=9.5, color=SUB_INK, labelpad=4)
    ax.set_yticks([])
    ax.tick_params(axis="x", length=3, color=MUTED, labelcolor=SUB_INK)
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    ax.spines["bottom"].set_color(MUTED)
    ax.spines["bottom"].set_linewidth(0.8)

    # Top: KDE of %ES-only across families
    ax = ax_kde
    es_pct_array = fam_pct[:, 2][totals > 0]
    y_max = 1.0
    if len(es_pct_array) >= 3:
        kde = gaussian_kde(es_pct_array, bw_method=0.35)
        xs = np.linspace(0, 100, 400)
        ys = kde(xs)
        y_max = float(ys.max()) * 2.4
        ax.fill_between(xs, 0, ys, color=C_ES, alpha=0.18,
                        edgecolor="none", zorder=1)
        ax.plot(xs, ys, color=C_ES, lw=1.4, zorder=2)
        median_es = float(np.median(es_pct_array))
        ax.axvline(median_es, color=C_ES, lw=0.9, ls="--", alpha=0.7, zorder=3)
        ax.text(0.5, y_max * 0.92,
                f"median = {median_es:.1f} %  (n = {len(es_pct_array)} families)",
                fontsize=9, color=C_ES, fontweight="bold", ha="left", va="top")
    ax.axvline(50, color=HAIR, lw=0.6, ls=":", zorder=1)
    ax.text(51.5, y_max * 0.92, "50 %  (equal partition)",
            ha="left", va="top", fontsize=8.5, color=SUB_INK)
    ax.set_xlim(0, 100); ax.set_xticks([])
    ax.set_ylim(0, y_max); ax.set_yticks([])
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    ax.spines["bottom"].set_color(HAIR)
    ax.spines["bottom"].set_linewidth(0.5)

    # Panel letters only (descriptive panel-C title omitted).
    ax_a.text(-0.10, 1.04, "A", transform=ax_a.transAxes,
              ha="left", va="bottom", fontsize=18, fontweight="bold",
              color=TITLE_NAVY)
    ax_b.text(-0.10, 1.04, "B", transform=ax_b.transAxes,
              ha="left", va="bottom", fontsize=18, fontweight="bold",
              color=TITLE_NAVY)
    ax_kde.text(-0.10, 1.55, "C", transform=ax_kde.transAxes,
                ha="left", va="bottom", fontsize=18, fontweight="bold",
                color=TITLE_NAVY)

    # Bottom legend
    ax_legend.axis("off")
    items_main = [
        (0.01, C_WE,     "WE only"),
        (0.15, C_SHARED, "Shared (WE ∩ ES)"),
        (0.39, C_ES,     "ES only"),
        (0.53, C_CHAR,   "Characterized"),
        (0.73, C_UNCHAR, "Uncharacterized"),
    ]
    for x, color, label in items_main:
        ax_legend.add_patch(Rectangle(
            (x, 0.30), 0.018, 0.42, facecolor=color, edgecolor="none",
            transform=ax_legend.transAxes, clip_on=False))
        ax_legend.text(x + 0.025, 0.50, label,
                       transform=ax_legend.transAxes, ha="left", va="center",
                       fontsize=9, color=INK)

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
    ap.add_argument("--dpi", type=int, default=300)
    ap.add_argument("--out", type=Path,
                    default=ROOT / "figures" / "Fig4.tif")
    args = ap.parse_args()
    render(args.dpi, args.out)


if __name__ == "__main__":
    main()
