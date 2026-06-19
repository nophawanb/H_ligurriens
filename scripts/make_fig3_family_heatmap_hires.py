"""
High-resolution re-render of Fig 3 (family composition × compartment
partition, < 20 kDa subset of the *Hemipyrellia ligurriens* larval
proteome).

Faithful copy of `scripts/make_fig3_family_stackedbars.py` (identical
data, taxonomy, layout and styling) with three deliberate differences,
matching the established hi-res criteria used for Fig 1 / Fig 2:

  1. Renders at a configurable high DPI (default 600).
  2. OMITS the figure heading — both the master `suptitle` and the
     subtitle statistics line beneath it (these belong in the caption).
  3. Uses repo-relative paths (the original hard-codes another machine's
     absolute path) and writes a SEPARATE output file
     (`figures/Fig3_family_heatmap.tif`), never touching the
     320-dpi `figures/Fig3_family_heatmap.png` consumed downstream.

Usage:
    python3 scripts/make_fig3_family_heatmap_hires.py [--dpi 600] [--out PATH]
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import openpyxl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle
from scipy.stats import gaussian_kde

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "Supplementary.xlsx"

LMW_CUTOFF = 21.0

# ===================== Family taxonomy =====================
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


def primary_acc(s):
    return str(s).split(";")[0]


def compute_matrix():
    row_group_of = []
    current = None
    for _, g in ROW_ORDER:
        if g is not None:
            current = g
        row_group_of.append(current)

    wb = openpyxl.load_workbook(SRC, read_only=True, data_only=True)
    ws = wb["proteinGroups"]
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
        if we_int <= 0 and es_int <= 0:
            continue
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

    M = np.array([[matrix[fam][c] for c in COMP_LABELS] for fam in ROW_LABELS],
                 dtype=float)
    return M, row_group_of


def safe_log2_es_over_we(we_n, es_n):
    if we_n == 0 and es_n == 0:
        return np.nan
    a = (es_n + 0.5)
    b = (we_n + 0.5)
    return float(np.log2(a / b))


# ===================== Design tokens =====================
C_WE      = "#1F77B4"
C_SHARED  = "#6A3FA0"
C_ES      = "#D72660"
INK       = "#1F2D3D"
SUB_INK   = "#5D6675"
MUTED     = "#8B93A1"
TITLE     = "#0F1E33"
HAIR      = "#D8DCE3"

GROUP_COLOR = {
    "AMP / immune-effector":       "#B85569",
    "Proteolytic":                 "#5C6F9C",
    "Antioxidant / metal-binding": "#5C8A5C",
    "Signalling":                  "#B58E47",
    "Other characterized":         "#7F8896",
    "Uncharacterized":             "#A3A3A3",
}

GROUP_BG_TINT = {
    "AMP / immune-effector":       "#FDE9E9",
    "Proteolytic":                 "#EAF0FA",
    "Antioxidant / metal-binding": "#EAF4EA",
    "Signalling":                  "#FAF1E1",
    "Other characterized":         "#F4F5F8",
    "Uncharacterized":             "#EEEEEE",
}


def render(dpi: int, out: Path) -> None:
    M, row_group_of = compute_matrix()
    totals = M.sum(axis=1)
    with np.errstate(invalid="ignore"):
        pct = np.where(totals[:, None] > 0, 100.0 * M / totals[:, None], 0.0)
    log2_ratio = np.array([
        safe_log2_es_over_we(M[i, 0], M[i, 2]) for i in range(M.shape[0])
    ])
    print(f"Grand total in < 20 kDa subset: {int(M.sum())}")

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans"],
        "axes.unicode_minus": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })

    fig = plt.figure(figsize=(7.5, 8.5), dpi=dpi)
    fig.patch.set_facecolor("white")
    fig.set_size_inches(7.5, 8.5)
    gs = gridspec.GridSpec(
        nrows=3, ncols=2,
        width_ratios=[1.0, 0.26],
        height_ratios=[0.16, 1.0, 0.08],
        wspace=0.025, hspace=0.10,
        left=0.27, right=0.95, top=0.92, bottom=0.05,
    )

    ax_top   = fig.add_subplot(gs[0, 0])
    ax_main  = fig.add_subplot(gs[1, 0])
    ax_right = fig.add_subplot(gs[1, 1])
    ax_leg   = fig.add_subplot(gs[2, 0])

    n_rows = len(ROW_LABELS)
    y_pos  = np.arange(n_rows)[::-1]
    bar_h  = 0.74

    # --------- Main: 100 % stacked horizontal bars ---------
    ax = ax_main
    for i, g in enumerate(row_group_of):
        ax.axhspan(y_pos[i] - 0.5, y_pos[i] + 0.5,
                   xmin=0, xmax=1,
                   facecolor=GROUP_BG_TINT.get(g, "white"),
                   edgecolor="none", alpha=0.45, zorder=0)

    ax.barh(y_pos, pct[:, 0], left=0,                  height=bar_h,
            color=C_WE,     edgecolor="white", linewidth=1.0, zorder=2)
    ax.barh(y_pos, pct[:, 1], left=pct[:, 0],          height=bar_h,
            color=C_SHARED, edgecolor="white", linewidth=1.0, zorder=2)
    ax.barh(y_pos, pct[:, 2], left=pct[:, 0]+pct[:, 1], height=bar_h,
            color=C_ES,     edgecolor="white", linewidth=1.0, zorder=2)

    MIN_LBL = 7.0
    for i in range(n_rows):
        if M[i, 0] > 0 and pct[i, 0] >= MIN_LBL:
            ax.text(pct[i, 0]/2, y_pos[i], f"{int(M[i,0])}",
                    ha="center", va="center",
                    fontsize=8.5, color="white", fontweight="bold", zorder=3)
        if M[i, 1] > 0 and pct[i, 1] >= MIN_LBL:
            ax.text(pct[i, 0] + pct[i, 1]/2, y_pos[i], f"{int(M[i,1])}",
                    ha="center", va="center",
                    fontsize=8.5, color="white", fontweight="bold", zorder=3)
        if M[i, 2] > 0 and pct[i, 2] >= MIN_LBL:
            ax.text(pct[i, 0] + pct[i, 1] + pct[i, 2]/2, y_pos[i],
                    f"{int(M[i,2])}",
                    ha="center", va="center",
                    fontsize=8.5, color="white", fontweight="bold", zorder=3)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(ROW_LABELS, fontsize=10, color=INK)
    ax.set_xlim(0, 100)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xticklabels(["0", "25", "50", "75", "100"],
                       fontsize=9.5, color=SUB_INK)
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
        ax.axhline(y_pos[idx] + 0.5, color="#C9CFDB", lw=0.6,
                   alpha=0.9, zorder=1)

    # --------- Right: total n + log2(ES/WE) lollipop ---------
    ax = ax_right
    ax.set_xlim(-3.6, 3.6); ax.set_ylim(-0.6, n_rows - 0.4)
    ax.axvline(0,  color=MUTED, lw=0.8, zorder=1)
    ax.axvline(-1, color=HAIR, lw=0.5, ls=":", zorder=1)
    ax.axvline( 1, color=HAIR, lw=0.5, ls=":", zorder=1)

    for i in range(n_rows):
        r = log2_ratio[i]
        if np.isnan(r):
            continue
        r_clipped = max(-3.5, min(3.5, r))
        colour = C_ES if r > 0 else C_WE if r < 0 else MUTED
        ax.hlines(y_pos[i], 0, r_clipped, color=colour, lw=1.4, zorder=2)
        ax.plot(r_clipped, y_pos[i], "o",
                markersize=5.5, color=colour,
                markeredgecolor="white", markeredgewidth=0.8, zorder=3)

    for i in range(n_rows):
        total = int(totals[i])
        ax.text(3.9, y_pos[i], f"n={total:,}",
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

    # --------- Top: KDE of %ES-only across families ---------
    ax = ax_top
    es_pct_array = pct[:, 2][totals > 0]
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
                fontsize=9, color=C_ES, fontweight="bold",
                ha="left", va="top")
    ax.axvline(50, color=HAIR, lw=0.6, ls=":", zorder=1)
    ax.text(51.5, y_max * 0.92, "50 %  (equal partition)",
            ha="left", va="top", fontsize=8.5, color=SUB_INK)

    ax.set_xlim(0, 100); ax.set_xticks([])
    ax.set_ylim(0, y_max)
    ax.set_yticks([])
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    ax.spines["bottom"].set_color(HAIR)
    ax.spines["bottom"].set_linewidth(0.5)

    # --------- Bottom: compartment legend ---------
    ax = ax_leg
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.axis("off")
    items = [("WE-only", C_WE), ("Shared (WE ∩ ES)", C_SHARED),
             ("ES-only", C_ES)]
    for k, (lbl, col) in enumerate(items):
        x = 0.05 + k * 0.30
        ax.add_patch(Rectangle((x, 0.40), 0.025, 0.30,
                               facecolor=col, edgecolor="none",
                               transform=ax.transAxes))
        ax.text(x + 0.030, 0.55, lbl,
                transform=ax.transAxes, ha="left", va="center",
                fontsize=10, color=INK, fontweight="bold")

    # --------- Heading intentionally omitted (per request) ---------
    # The master suptitle and the subtitle statistics line
    # ("n = ... proteins · MW < 21 kDa · per-row partition ...") are not
    # drawn; these details belong in the figure caption. bbox_inches=
    # "tight" trims the freed top whitespace.

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
    ap.add_argument("--dpi", type=int, default=300)
    ap.add_argument("--out", type=Path,
                    default=ROOT / "figures" / "Fig3_family_heatmap.tif")
    args = ap.parse_args()
    render(args.dpi, args.out)


if __name__ == "__main__":
    main()
