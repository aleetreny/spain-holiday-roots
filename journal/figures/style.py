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

FONTS = Path("/System/Library/Fonts/Supplemental")
for name in ("Arial.ttf", "Arial Bold.ttf"):
    font_manager.fontManager.addfont(FONTS / name)
mpl.rcParams.update({
    "font.family": "Arial", "font.size": 8, "pdf.fonttype": 42,
    "axes.titlesize": 8.5, "axes.labelsize": 8, "axes.linewidth": 0.5,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "axes.grid.axis": "y", "axes.axisbelow": True,
    "grid.color": GRID, "grid.linewidth": 0.5,
    "xtick.labelsize": 7, "ytick.labelsize": 7,
    "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "xtick.color": INK, "ytick.color": INK, "axes.labelcolor": INK,
    "axes.edgecolor": NEUTRAL, "text.color": INK,
    "lines.linewidth": 0.85, "lines.markersize": 3.8,
    "legend.fontsize": 7, "legend.frameon": False,
    "savefig.dpi": 300, "figure.facecolor": "white", "savefig.facecolor": "white",
})


def read(arm, name, **kwargs):
    return pd.read_csv(ANALYSIS / arm / "out" / f"{name}.csv", **kwargs)


def across_years(name):
    return pd.concat([read(A, name), read(A2, name)], ignore_index=True)


def figure(height, width=170):
    return plt.figure(figsize=(width / 25.4, height / 25.4))


def panel(ax, letter, title=None, x=-0.045, y=1.03):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=10, fontweight="bold",
            va="bottom", ha="left")
    if title:
        ax.set_title(title, loc="left", pad=6)


def reference(ax, value=0, *, vertical=False):
    line = ax.axvline if vertical else ax.axhline
    line(value, color=NEUTRAL, linestyle=(0, (3, 3)), lw=0.55, zorder=1)


def interval(ax, position, point, lo, hi, *, horizontal=True, color=NEUTRAL,
             marker="o", filled=True, size=4.2):
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
