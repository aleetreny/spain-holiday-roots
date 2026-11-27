"""Map the measurement geography, the model dyad and source coverage."""

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from matplotlib import cm
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.patches import FancyArrowPatch
from matplotlib.ticker import PercentFormatter

from style import (
    A, C, HEAVY, MAP_ZERO, MUTED, NEUTRAL, ORDINARY, RECEIVER, SENDER,
    across_years, figure, map_axes, map_geometry, nodes, panel, read, save,
)

HOLIDAYS = {"EASTER": ("D", "Good Friday"), "ALLSAINTS": ("s", "All Saints"),
            "XMAS": ("^", "Christmas Eve")}


def map_data(manifest, provinces):
    """Population-weighted province balances and the 12 largest excess dyads."""
    base = HEAVY / A / "panel"
    daily = pq.read_table(base / "panel_zones2435_daily.parquet",
                          filters=[("date", "=", "2024-08-15")]).to_pandas()
    daily = daily.merge(nodes()[["id", "population"]], left_on="zone_id", right_on="id",
                        validate="one_to_one")
    balance = daily.assign(prov=daily.id.str[:2]).groupby("prov")[["net_flow", "population"]].sum()
    balance["net_per_resident"] = balance.net_flow / balance.population

    summer = manifest.query("year == 2024 and night_class == 'AUG'").file_date.tolist()
    ordinary = manifest.query("year == 2024 and night_class == 'O'").file_date.tolist()
    dates = summer + ordinary
    city_flows = pq.read_table(base / "panel_cities137_prov.parquet",
                               columns=["date", "city", "dprov", "flow_all"],
                               filters=[("date", "in", dates)]).to_pandas()
    assert set(city_flows.date) == set(dates) and city_flows.flow_all.notna().all()
    totals = city_flows.assign(summer=city_flows.date.isin(summer)).pivot_table(
        index=["city", "dprov"], columns="summer", values="flow_all", aggfunc="sum", fill_value=0,
    )
    corridors = pd.DataFrame({"august_mean": totals[True] / len(summer),
                              "ordinary_mean": totals[False] / len(ordinary)})
    corridors["excess_persons"] = corridors.august_mean - corridors.ordinary_mean
    corridors = corridors.nlargest(12, "excess_persons").reset_index()
    corridors["n_august"], corridors["n_ordinary"] = len(summer), len(ordinary)
    corridors = corridors.merge(nodes()[["id", "name", "x", "y"]],
                                left_on="city", right_on="id", validate="many_to_one")
    centres = provinces.geometry.representative_point()
    corridors["dest_x"] = corridors.dprov.map(centres.x)
    corridors["dest_y"] = corridors.dprov.map(centres.y)
    assert corridors[["dest_x", "dest_y"]].notna().all().all()
    return balance, corridors


def draw_map(fig, ax, balance, corridors, provinces):
    panel(ax, "a", x=-.02, y=1.015)
    ax.text(.032, 1.026, "Overnight balance and excess August flows",
            fontsize=8.5, transform=ax.transAxes, va="bottom")
    cmap = LinearSegmentedColormap.from_list("balance", [SENDER, MAP_ZERO, RECEIVER])
    magnitude = np.ceil(balance.net_per_resident.abs().max() * 10) / 10
    norm = TwoSlopeNorm(vmin=-magnitude, vcenter=0, vmax=magnitude)
    provinces.join(balance).plot(ax=ax, column="net_per_resident", cmap=cmap, norm=norm,
                                edgecolor=ORDINARY, linewidth=0.35)
    map_axes(ax)
    max_flow = corridors.excess_persons.max()
    for row in corridors.sort_values("excess_persons").itertuples():
        arrow = FancyArrowPatch((row.x, row.y), (row.dest_x, row.dest_y),
                                connectionstyle="arc3,rad=0.18", arrowstyle="-|>",
                                mutation_scale=5, linewidth=1.65 * row.excess_persons / max_flow,
                                color=NEUTRAL, alpha=0.78, shrinkA=1, shrinkB=1, zorder=3)
        ax.add_patch(arrow)

    offsets = {"Madrid": (-5, 6), "Barcelona": (-9, 7), "València": (5, -6),
               "Valencia": (5, -6), "Sevilla": (-3, -8), "Zaragoza": (-5, 6),
               "Bilbao": (-4, 6)}
    labels = (corridors.groupby("name").excess_persons.sum().nlargest(3).index)
    for row in corridors.drop_duplicates("name").query("name in @labels").itertuples():
        dx, dy = offsets.get(row.name, (4, 4))
        ax.plot(row.x, row.y, "o", color=NEUTRAL, ms=3.6, zorder=4)
        ax.annotate(row.name, (row.x, row.y), xytext=(dx, dy), textcoords="offset points",
                    ha="right" if dx < 0 else "left", va="center", fontsize=7, zorder=5,
                    bbox={"fc": "white", "ec": "none", "alpha": .87, "pad": .65})
    ax.plot([4_000, 204_000], [4_023_000] * 2, color=NEUTRAL, lw=0.85)
    ax.text(104_000, 3_997_000, "200 km", ha="center", va="top", fontsize=7)

    cax = fig.add_axes([.106, .135, .36, .018])
    colorbar = fig.colorbar(
        cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax,
        orientation="horizontal", ticks=[-magnitude, 0, magnitude], format=PercentFormatter(1),
    )
    colorbar.solids.set_rasterized(False)
    cax.set_xlabel("Net overnight balance per resident")
    fig.text(.285, .034, "Arrows show the 12 largest excess flows", ha="center", fontsize=7)


def draw_dyad(ax):
    panel(ax, "b", "One city, two destination provinces", x=-.08)
    ax.set(xlim=(0, 1), ylim=(0, 1))
    ax.set_axis_off()
    box = {"boxstyle": "round,pad=0.34,rounding_size=0.18", "fc": "white", "ec": NEUTRAL, "lw": .55}
    ax.text(.045, .63, "City", ha="center", va="center", bbox=box, zorder=4)
    for y, province in ((.77, "Province p"), (.48, "Province q")):
        ax.text(.98, y, province, ha="right", va="center", bbox=box, zorder=4)
        for rad, ls, colour in ((.085, "--", MUTED), (-.085, "-", NEUTRAL)):
            ax.add_patch(FancyArrowPatch((.105, .63), (.73, y), connectionstyle=f"arc3,rad={rad}",
                                        arrowstyle="-|>", linestyle=ls, color=colour, lw=.85,
                                        mutation_scale=6, shrinkA=0, shrinkB=0))
    ax.plot([.05, .16], [.985, .985], "--", color=MUTED, lw=.85)
    ax.text(.19, .985, "City residents born there", fontsize=7, va="center")
    ax.plot([.05, .16], [.87, .87], color=NEUTRAL, lw=.85)
    ax.text(.19, .87, "City residents sleeping there", fontsize=7, va="center")
    ax.text(.01, .26, "Flow ~ birth share", weight="bold", fontsize=8, va="top")
    ax.text(.01, .16, "Compare cities of the same province against\n"
            "the same destination, with distance and income.", fontsize=7, va="top")


def draw_coverage(ax, manifest, months):
    panel(ax, "c", "Coverage by source", x=-.08)
    ax.set(xlim=(0, 366), ylim=(-.95, 3.62))
    ax.set_yticks([3, 2, 1, 0], ["2022", "2023", "2024", "2025"])
    starts = pd.date_range("2024-01-01", periods=12, freq="MS").day_of_year - 1
    ax.set_xticks(starts[::2], ["Jan", "Mar", "May", "Jul", "Sep", "Nov"])
    ax.spines[["left", "bottom"]].set_visible(False)
    ax.tick_params(length=0, pad=3)
    ax.grid(False)
    ax.fill_betweenx([-.12, .13], 182, 366, facecolor="#f5f5f5",
                     edgecolor=ORDINARY, linewidth=.35, hatch="////", zorder=0)
    for year, nights in manifest.loc[manifest.file_present].groupby("year"):
        y = 2025 - year
        for klass, subset in nights.groupby("coverage_class"):
            if klass in HOLIDAYS:
                ax.plot(subset.day, [y] * len(subset), linestyle="none", marker=HOLIDAYS[klass][0],
                        ms=3.7, color=NEUTRAL, zorder=5)
            else:
                col = NEUTRAL if klass == "AUG" else ORDINARY if klass == "O" else "#dedede"
                width = .42 if klass in ("AUG", "O") else .25
                ax.vlines(subset.day, y - .105, y + .105, color=col, linewidth=width)
        available = months.loc[months.year.eq(year)]
        ax.hlines(np.full(len(available), y - .27), available.start, available.end - 3, color=MUTED, lw=.8)
    ax.text(362, .28, "New method", ha="right", fontsize=7, color=MUTED)
    ax.text(0, -.75, "Daily ticks   MITMA     Monthly bars   INE", fontsize=7, color=MUTED)


def main():
    _, provinces = map_geometry()
    manifest = across_years("night_manifest")
    manifest["coverage_class"] = manifest.night_class.str.replace("2025", "", regex=False)
    # Project month/day onto a leap year so every row aligns with the month ticks.
    manifest["day"] = pd.to_datetime("2024-" + manifest.file_date.str[5:]).dt.day_of_year - 1
    months = read(C, "od_inventory").assign(year=lambda d: d.month.str[:4].astype(int))
    months = months.query("2022 <= year <= 2025").copy()
    months["start"] = pd.to_datetime("2024-" + months.month.str[-2:] + "-01").dt.day_of_year - 1
    months["end"] = months.start + pd.to_datetime(months.month).dt.days_in_month
    balance, corridors = map_data(manifest, provinces)
    fig = figure(120)
    draw_map(fig, fig.add_axes([.028, .185, .52, .695]), balance, corridors, provinces)
    draw_dyad(fig.add_axes([.604, .555, .366, .34]))
    draw_coverage(fig.add_axes([.604, .196, .366, .272]), manifest, months)
    x = .604
    for mark, label, width, col in (
        ("|", "Ordinary", .125, ORDINARY), ("|", "August window", .20, NEUTRAL),
        ("D", "Good Friday", .142, NEUTRAL), ("s", "All Saints", .118, NEUTRAL),
        ("^", "Christmas Eve", .16, NEUTRAL),
    ):
        y = .118 if mark == "|" else .069
        if mark == "D":
            x = .604
        fig.add_artist(Line2D([x, x], [y, y], transform=fig.transFigure, linestyle="none",
                              marker=mark, color=col, ms=3.7, mew=.65))
        fig.text(x + .017, y, label, fontsize=7, va="center")
        x += width
    save(fig, 1, balance.reset_index().assign(panel="a", series="province balance"),
         corridors.assign(panel="a", series="excess corridor"),
         manifest.loc[manifest.file_present, ["file_date", "year", "day", "is_valid", "coverage_class", "method"]]
         .assign(panel="c", series="MITMA night"),
         months.assign(panel="c", series="INE month"))


if __name__ == "__main__":
    main()
