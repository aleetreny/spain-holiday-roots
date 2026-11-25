"""F1 drafts: four candidate hero figures (three maps + BOSAM).

Drafts only. Shared visual system: serif type (CMU Serif, as LaTeX), sea blue = coast,
terracotta = inland, warm grey = base/ordinary, dark ink for text.
Run: .venv/bin/python -I pipeline/f01_hero_drafts.py
"""
from pathlib import Path
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "reports" / "figures_drafts"
OUT.mkdir(parents=True, exist_ok=True)
NET = ROOT / "results" / "b01_networks"

BLUE, TERRA, GREY, INK, MUTED = "#1f6e8c", "#c0623b", "#b5b1aa", "#1d1d1b", "#6b6862"
LAND, BORDER = "#f1eee8", "#ffffff"
plt.rcParams.update({"font.family": "STIX Two Text", "font.size": 9, "axes.titlesize": 10,
                     "text.color": INK, "axes.edgecolor": INK, "savefig.dpi": 300})
EXCL = {"07", "35", "38", "51", "52"}
NIGHTS = {"feb": "2024-02-14", "aug": "2024-08-15", "xmas": "2024-12-25"}
LABEL = {"feb": "An ordinary night\nWednesday 14 February 2024",
         "aug": "A holiday night\nThursday 15 August 2024",
         "xmas": "Christmas Eve\n24–25 December 2024"}
MIN_KM, MIN_W = 100, 100

nodes = pd.read_parquet(NET / "nodes.parquet").set_index("idx").sort_index()
X, Y, COAST = nodes.x.values, nodes.y.values, nodes.coastal_true.values.astype(bool)

def provinces():
    g = gpd.read_file(ROOT / "clean_data/v2/zones/municipios_mitma.gpkg")
    g["prov"] = g["id"].str[:2]
    g = g[~g.prov.isin(EXCL)]
    p = g.dissolve("prov").reset_index()
    p["geometry"] = p.geometry.simplify(800)
    return p

PROV = provinces()

def base(ax):
    PROV.plot(ax=ax, color=LAND, edgecolor=BORDER, linewidth=0.5)
    PROV.dissolve().boundary.plot(ax=ax, color="#cfcac1", linewidth=0.5)
    ax.set_axis_off(); ax.set_aspect("equal")
    ax.set_xlim(-30000, 1_150_000); ax.set_ylim(3_960_000, 4_880_000)

def edges(night):
    e = pd.read_parquet(NET / f"edges_{NIGHTS[night]}__full.parquet")
    e = e[e.i != e.j].copy()
    e["km"] = np.hypot(X[e.i] - X[e.j], Y[e.i] - Y[e.j]) / 1000
    e["coast"] = COAST[e.j]
    return e

def arcs(e, bend=0.18, n=24):
    """Quadratic Bezier arcs, bending to the right of travel (direction is readable)."""
    x0, y0, x1, y1 = X[e.i], Y[e.i], X[e.j], Y[e.j]
    mx, my = (x0 + x1) / 2, (y0 + y1) / 2
    dx, dy = x1 - x0, y1 - y0
    cx, cy = mx + bend * dy, my - bend * dx
    t = np.linspace(0, 1, n)[None, :]
    px = (1 - t) ** 2 * x0[:, None] + 2 * (1 - t) * t * cx[:, None] + t ** 2 * x1[:, None]
    py = (1 - t) ** 2 * y0[:, None] + 2 * (1 - t) * t * cy[:, None] + t ** 2 * y1[:, None]
    return np.stack([px, py], axis=-1)

def draw_flows(ax, e, wref, lw_max=1.6):
    e = e.sort_values("w")
    s = np.sqrt(e.w.values / wref).clip(0, 1)
    cols = np.where(e.coast.values, BLUE, TERRA)
    rgba = matplotlib.colors.to_rgba_array(cols)
    rgba[:, 3] = 0.06 + 0.70 * s
    ax.add_collection(LineCollection(arcs(e), colors=rgba, linewidths=0.10 + lw_max * s,
                                     capstyle="round", zorder=3))

def legend(ax, loc="lower right"):
    h = [Line2D([], [], color=BLUE, lw=2, label="to the coast"),
         Line2D([], [], color=TERRA, lw=2, label="inland")]
    ax.legend(handles=h, loc=loc, frameon=False, fontsize=8.5, handlelength=1.8)

def stats(e):
    far = e[e.km >= MIN_KM]
    return far.w.sum(), far[far.coast].w.sum() / far.w.sum()

E = {k: edges(k) for k in NIGHTS}

# ---------- Draft A: two Spains side by side ----------
def draft_a():
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.6))
    sel = {k: E[k][(E[k].km >= MIN_KM) & (E[k].w >= MIN_W)] for k in ("feb", "aug")}
    wref = sel["aug"].w.quantile(0.995)
    for ax, k in zip(axs, ("feb", "aug")):
        base(ax); draw_flows(ax, sel[k], wref)
        tot, cs = stats(E[k])
        ax.set_title(LABEL[k], loc="left", fontsize=9.5, linespacing=1.3)
        ax.text(0.0, -0.02, f"{tot/1e6:.1f} million people slept at least {MIN_KM} km from home\n"
                f"{cs:.0%} of them at the coast, {1-cs:.0%} inland",
                transform=ax.transAxes, va="top", fontsize=8, color=MUTED)
    legend(axs[1])
    fig.savefig(OUT / "F1_A_two_spains.png", bbox_inches="tight"); plt.close(fig)

# ---------- Draft B: August fans from the five largest cities ----------
def draft_b():
    big = nodes.sort_values("population", ascending=False).head(5)
    e = E["aug"]; e = e[(e.km >= MIN_KM) & e.i.isin(big.index) & (e.w >= 20)]
    fig, ax = plt.subplots(figsize=(6.4, 5.0)); base(ax)
    draw_flows(ax, e, e.w.quantile(0.995), lw_max=2.2)
    ax.scatter(X[big.index], Y[big.index], s=26, color=INK, zorder=5, edgecolor="white", lw=1.2)
    for i, r in big.iterrows():
        ax.text(X[i] + 14000, Y[i] + 12000, r["name"].split(" agreg")[0], fontsize=8.5, zorder=6,
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.7, pad=0.6))
    cs = e[e.coast].w.sum() / e.w.sum()
    ax.set_title("15 August 2024: where the residents of the five largest cities slept\n"
                 f"(at least {MIN_KM} km from home; {cs:.0%} at the coast, {1-cs:.0%} inland)",
                 loc="left", fontsize=9.5)
    legend(ax)
    fig.savefig(OUT / "F1_B_city_fans.png", bbox_inches="tight"); plt.close(fig)

# ---------- Draft C: one coastal town, three nights ----------
def draft_c():
    tgt = nodes.index[nodes.id == "29094"][0]
    fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.6))
    vmax = max(E[k][E[k].j == tgt].w.max() for k in NIGHTS)
    for ax, k in zip(axs, ("feb", "aug", "xmas")):
        base(ax)
        v = E[k][(E[k].j == tgt)]
        same = nodes.prov.values[v.i] == nodes.prov[tgt]
        ax.scatter(X[v.i], Y[v.i], s=60 * np.sqrt(v.w / vmax) * 6, c=np.where(same, GREY, TERRA),
                   alpha=0.75, lw=0, zorder=3)
        ax.scatter([X[tgt]], [Y[tgt]], marker="*", s=70, color=BLUE, zorder=5, edgecolor="white", lw=0.6)
        ax.set_title(LABEL[k].replace("\n", "\n"), loc="left", fontsize=8.5)
        ax.text(0, -0.02, f"{v.w.sum():,.0f} visitors\n{same @ v.w.values / v.w.sum():.0%} from its own province",
                transform=ax.transAxes, va="top", fontsize=7.5, color=MUTED)
    h = [Line2D([], [], marker="*", ls="", color=BLUE, ms=9, label="the coastal town"),
         Line2D([], [], marker="o", ls="", color=GREY, ms=6, label="visitors from its own province"),
         Line2D([], [], marker="o", ls="", color=TERRA, ms=6, label="visitors from elsewhere")]
    fig.legend(handles=h, loc="lower center", ncol=3, frameon=False, fontsize=8, bbox_to_anchor=(0.5, -0.12))
    fig.savefig(OUT / "F1_C_one_town.png", bbox_inches="tight"); plt.close(fig)

# ---------- Draft D: BOSAM, February vs August ----------
def draft_d():
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.8))
    tags = {"feb": "ordinary_2024-02-14", "aug": "summer_2024-08-15"}
    for ax, k in zip(axs, ("feb", "aug")):
        ed = pd.read_parquet(ROOT / f"results/b02_bosam_edges_{tags[k]}__w10.parquet")
        order = pd.read_csv(ROOT / f"results/b02_bosam_order_{tags[k]}__w10.csv")
        coast_pos = COAST[order.sort_values("position").idx.values]
        touch = coast_pos[ed.row.values] | coast_pos[ed["col"].values]
        for m, c, z in ((~touch, TERRA, 2), (touch, BLUE, 3)):
            ax.scatter(ed["col"].values[m], ed.row.values[m], s=0.02, color=c, lw=0, rasterized=True, zorder=z)
        N = len(order)
        ax.set_xlim(0, N); ax.set_ylim(N, 0); ax.set_aspect("equal")
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values(): s.set_linewidth(0.4); s.set_color(MUTED)
        # strip: which positions are coastal zones
        ax.vlines(np.where(coast_pos)[0], N + 25, N + 95, color=BLUE, lw=0.25, clip_on=False)
        kmean = 2 * len(ed) / 2 / N
        ax.set_title(LABEL[k], loc="left", fontsize=9.5)
        ax.text(0, -0.09, f"{len(ed)//2:,} links; nodes sorted by degree (low to high)",
                transform=ax.transAxes, va="top", fontsize=7.5, color=MUTED)
    h = [Line2D([], [], marker="s", ls="", color=BLUE, ms=6, label="link touches a coastal zone"),
         Line2D([], [], marker="s", ls="", color=TERRA, ms=6, label="link between inland zones")]
    fig.legend(handles=h, loc="lower center", ncol=2, frameon=False, fontsize=8, bbox_to_anchor=(0.5, -0.06))
    fig.savefig(OUT / "F1_D_bosam.png", bbox_inches="tight"); plt.close(fig)

if __name__ == "__main__":
    for f in (draft_a, draft_b, draft_c, draft_d):
        f(); print("done", f.__name__)
