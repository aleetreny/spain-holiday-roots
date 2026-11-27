"""Shared typography, plotting helpers and read-only inputs for the journal figures."""

from functools import lru_cache
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
import pandas as pd

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
MAIN = Path("/Users/alejandrotreny/Proyectos/spain-holiday-roots")
JOURNAL = MAIN.with_name("spain-holiday-roots-journal")
ANALYSIS = JOURNAL / "journal/analysis"
HEAVY = JOURNAL / "_local/journal_data"
A = "A_calendar"
A2 = "A2_confirmation2025"
B = "B_cohorts"
C = "C_municipal"
D = "D_villages"

INK, MUTED, NEUTRAL = "#222222", "#666666", "#444444"
INLAND, COAST = "#c0623b", "#2b6f9e"
EXODUS, RECENT = "#7b3fa0", "#0f8f7a"
SENDER, MAP_ZERO, RECEIVER = "#7b3fa0", "#efefef", "#1b7837"
ORDINARY, BAND, GRID = "#bdbdbd", "#e6e6e6", "#ececec"
MARK = 4.0         # point-estimate marker size (pt)
SMALL = 3.5        # background/observation marker size (pt)

# Type scale (pt) used by every figure: ticks, annotations and legends; axis labels;
# panel titles; panel letters.
TEXT, LABEL, TITLE, LETTER = 7, 7.5, 8.5, 10
MM = 1 / 25.4
LETTER_GAP = 4.2   # mm from a panel letter to its title
TITLE_RISE = 3.2   # mm from the axes top to the letter/title baseline

FONTS = Path("/System/Library/Fonts/Supplemental")
for name in ("Arial.ttf", "Arial Bold.ttf", "Arial Italic.ttf"):
    font_manager.fontManager.addfont(FONTS / name)
mpl.rcParams.update({
    "font.family": "Arial", "font.size": TEXT, "pdf.fonttype": 42,
    "axes.titlesize": TITLE, "axes.labelsize": LABEL, "axes.linewidth": 0.5,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "axes.grid.axis": "y", "axes.axisbelow": True,
    "grid.color": GRID, "grid.linewidth": 0.5,
    "xtick.labelsize": TEXT, "ytick.labelsize": TEXT,
    "xtick.major.pad": 2.5, "ytick.major.pad": 2.5, "axes.labelpad": 3,
    "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "xtick.color": INK, "ytick.color": INK, "axes.labelcolor": INK,
    "axes.edgecolor": NEUTRAL, "text.color": INK,
    "lines.linewidth": 0.85, "lines.markersize": 3.8,
    "legend.fontsize": TEXT, "legend.frameon": False, "legend.handletextpad": .35,
    "legend.columnspacing": 1.4, "legend.handlelength": 1.4, "legend.borderaxespad": 0,
    "legend.borderpad": 0,
    "mathtext.fontset": "custom", "mathtext.rm": "Arial", "mathtext.it": "Arial:italic",
    "mathtext.bf": "Arial:bold",
    "savefig.dpi": 300, "figure.facecolor": "white", "savefig.facecolor": "white",
})


def read(arm, name, **kwargs):
    return pd.read_csv(ANALYSIS / arm / "out" / f"{name}.csv", **kwargs)


def across_years(name):
    return pd.concat([read(A, name), read(A2, name)], ignore_index=True)


def figure(height, width=170):
    return plt.figure(figsize=(width * MM, height * MM))


def place(fig, left, top, width, height, **kwargs):
    """Add axes from millimetres measured from the figure's top-left corner."""
    W, H = fig.get_figwidth() / MM, fig.get_figheight() / MM
    return fig.add_axes([left / W, 1 - (top + height) / H, width / W, height / H], **kwargs)


def panel(ax, letter, title=None, *, at, rise=TITLE_RISE):
    """Panel letter at `at` mm from the figure's left edge and the title beside it.

    Both sit on one baseline `rise` mm above the axes top, so every panel letter in a
    column shares its x position and every panel in a row shares its baseline.
    """
    fig = ax.figure
    x = at / (fig.get_figwidth() / MM)
    kwargs = dict(xycoords=(fig.transFigure, ax.transAxes), textcoords="offset points",
                  va="baseline", ha="left", annotation_clip=False)
    ax.annotate(letter, (x, 1), xytext=(0, rise * 72 * MM),
                fontsize=LETTER, fontweight="bold", **kwargs)
    if title:
        ax.annotate(title, (x, 1), xytext=(LETTER_GAP * 72 * MM, rise * 72 * MM),
                    fontsize=TITLE, **kwargs)


def legend(ax, handles, *, loc="upper left", anchor=(0, 1), ncol=1, **kwargs):
    return ax.legend(handles=handles, loc=loc, bbox_to_anchor=anchor, ncol=ncol, **kwargs)


def reference(ax, value=0, *, vertical=False):
    line = ax.axvline if vertical else ax.axhline
    line(value, color=NEUTRAL, linestyle=(0, (3, 3)), lw=0.55, zorder=1)


def interval(ax, position, point, lo, hi, *, horizontal=True, color=NEUTRAL,
             marker="o", filled=True, size=MARK):
    if horizontal:
        ax.hlines(position, lo, hi, color=color, lw=0.85, zorder=3)
        x, y = point, position
    else:
        ax.vlines(position, lo, hi, color=color, lw=0.85, zorder=3)
        x, y = position, point
    ax.plot(x, y, linestyle="none", marker=marker, ms=size, mfc=color if filled else "white",
            mec=color, mew=0.75, zorder=4)


def forest(ax, labels, limits, ticks):
    ax.set_yticks(range(len(labels)), labels)
    ax.set_ylim(len(labels) - 0.45, -0.55)
    ax.set_xlim(*limits)
    ax.set_xticks(ticks)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    reference(ax, vertical=True)


def save(fig, number, *tables):
    OUT.mkdir(parents=True, exist_ok=True)
    name = f"fig{number}"
    pd.concat(tables, ignore_index=True, sort=False).to_csv(OUT / f"{name}_data.csv", index=False)
    for suffix in ("pdf", "png"):
        fig.savefig(OUT / f"{name}.{suffix}")
    plt.close(fig)


@lru_cache
def nodes():
    return pd.read_parquet(MAIN / "results/b01_networks/nodes.parquet")


@lru_cache
def map_geometry():
    """Use the coursework EPSG:3042 zone geometry and mainland node universe."""
    import geopandas as gpd

    geo = gpd.read_file(MAIN / "clean_data/v2/zones/municipios_mitma.gpkg")[["id", "geometry"]]
    geo = geo.loc[geo.id.isin(nodes().id)].copy()
    assert geo.crs.to_epsg() == 3042 and len(geo) == len(nodes())
    geo["prov"] = geo.id.str[:2]
    provinces = geo.dissolve("prov")
    # Shared-edge simplification preserves an uninterrupted province basemap.
    if hasattr(provinces.geometry, "simplify_coverage") and provinces.geometry.is_valid_coverage():
        provinces.geometry = provinces.geometry.simplify_coverage(500)
    else:
        provinces.geometry = provinces.geometry.simplify(300)
    geo.geometry = geo.geometry.simplify(350)
    return geo, provinces


def map_axes(ax):
    ax.set(xlim=(-30_000, 1_150_000), ylim=(3_960_000, 4_880_000), aspect="equal")
    ax.set_axis_off()


def decimals(n):
    """Tick formatter with a fixed number of decimals and a true minus sign."""
    from matplotlib.ticker import FuncFormatter
    return FuncFormatter(lambda value, _: f"{value:.{n}f}".replace("-", "−"))
