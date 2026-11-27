"""Map the measurement geography, the model dyad and source coverage."""

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from matplotlib import cm, patheffects
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.patches import FancyArrowPatch
from matplotlib.ticker import PercentFormatter

from style import (
    A, C, HEAVY, INK, LABEL, MAP_ZERO, MUTED, NEUTRAL, ORDINARY, RECEIVER, SENDER, SMALL, TEXT,
    across_years, figure, map_axes, map_geometry, nodes, panel, place, read, save,
)

TICKS = {"O": "#a8a8a8", "AUG": INK, "other": "#e4e4e4"}
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


HALO = [patheffects.withStroke(linewidth=2.2, foreground="white")]
LEFT, RIGHT = 2.0, 108.0      # mm: panel-letter columns


def arrow_width(excess, max_flow):
    return 2.4 * excess / max_flow


def draw_map(fig, ax, balance, corridors, provinces):
    panel(ax, "a", "Overnight balance and excess August flows", at=LEFT)
    cmap = LinearSegmentedColormap.from_list("balance", [SENDER, MAP_ZERO, RECEIVER])
    magnitude = np.ceil(balance.net_per_resident.abs().max() * 10) / 10
    norm = TwoSlopeNorm(vmin=-magnitude, vcenter=0, vmax=magnitude)
    provinces.join(balance).plot(ax=ax, column="net_per_resident", cmap=cmap, norm=norm,
                                edgecolor="white", linewidth=0.45)
    map_axes(ax)
    max_flow = corridors.excess_persons.max()
    for row in corridors.sort_values("excess_persons").itertuples():
        width = arrow_width(row.excess_persons, max_flow)
        ax.add_patch(FancyArrowPatch(
            (row.x, row.y), (row.dest_x, row.dest_y), connectionstyle="arc3,rad=0.18",
            arrowstyle="-|>,head_length=0.55,head_width=0.32", mutation_scale=3 + 1.5 * width,
            linewidth=width, color=INK, shrinkA=2, shrinkB=1, zorder=3, joinstyle="miter",
            path_effects=[patheffects.withStroke(linewidth=width + 1.3, foreground="white"),
                          patheffects.Normal()],
        ))
    offsets = {"Madrid": (6, 5, "left"), "Barcelona": (-5, 6, "right"), "Sevilla": (-5, -6, "right"),
               "Valencia": (6, 0, "left"), "València": (6, 0, "left")}
    labels = corridors.groupby("name").excess_persons.sum().nlargest(3).index
    for row in corridors.drop_duplicates("name").query("name in @labels").itertuples():
        dx, dy, ha = offsets.get(row.name, (5, 5, "left"))
        ax.plot(row.x, row.y, "o", color=INK, ms=4.2, mec="white", mew=.7, zorder=6)
        ax.annotate(row.name, (row.x, row.y), xytext=(dx, dy), textcoords="offset points",
                    ha=ha, va="center", fontsize=TEXT, zorder=7, path_effects=HALO)
    ax.plot([20_000, 220_000], [4_000_000] * 2, color=NEUTRAL, lw=0.85, solid_capstyle="butt")
    ax.text(120_000, 4_012_000, "200 km", ha="center", va="bottom", fontsize=TEXT)
    return cmap, norm, magnitude, max_flow


def map_keys(fig, cmap, norm, magnitude, max_flow):
    """Colour bar and arrow-width key on one row below the map."""
    cax = place(fig, 8, 93, 42, 2.2)
    colorbar = fig.colorbar(
        cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation="horizontal",
        ticks=np.linspace(-magnitude, magnitude, 5), format=PercentFormatter(1, decimals=0),
    )
    colorbar.solids.set_rasterized(False)
    colorbar.outline.set_linewidth(.4)
    cax.tick_params(length=2, width=.5)
    cax.set_title("Net overnight balance per resident", fontsize=TEXT, pad=3, loc="left")
    key = place(fig, 60, 93, 40, 6.5)
    key.set(xlim=(0, 1), ylim=(0, 1))
    key.set_axis_off()
    key.set_title("Excess people per night", fontsize=TEXT, pad=3, loc="left")
    for x, value in ((0, 40_000), (.5, 80_000)):
        width = arrow_width(value, max_flow)
        key.add_patch(FancyArrowPatch(
            (x, .72), (x + .2, .72), arrowstyle="-|>,head_length=0.55,head_width=0.32",
            mutation_scale=3 + 1.5 * width, linewidth=width, color=INK, shrinkA=0, shrinkB=0,
        ))
        key.text(x + .24, .72, f"{value:,}", va="center", fontsize=TEXT)


def draw_dyad(ax):
    panel(ax, "b", "One city, two destination provinces", at=RIGHT)
    ax.set(xlim=(0, 1), ylim=(0, 1))
    ax.set_axis_off()
    box = {"boxstyle": "round,pad=0.38,rounding_size=0.15", "fc": "white", "ec": NEUTRAL, "lw": .6}
    origin, ends = (.08, .76), ((.70, .97), (.70, .55))
    ax.text(*origin, "City", ha="center", va="center", bbox=box, zorder=4, fontsize=LABEL)
    for (x, y), name in zip(ends, ("p", "q")):
        ax.text(x + .015, y, f"Province $\\it{{{name}}}$", ha="left", va="center", bbox=box,
                zorder=4, fontsize=LABEL)
        for rad, ls, colour in ((.10, (0, (3, 2)), MUTED), (-.10, "-", INK)):
            ax.add_patch(FancyArrowPatch((origin[0] + .07, origin[1]), (x - .01, y),
                                         connectionstyle=f"arc3,rad={rad}",
                                         arrowstyle="-|>,head_length=0.5,head_width=0.25",
                                         linestyle=ls, color=colour, lw=.85,
                                         mutation_scale=7, shrinkA=0, shrinkB=0))
    for y, ls, colour, text in ((.33, (0, (3, 2)), MUTED, "City residents born there (census)"),
                                (.21, "-", INK, "City residents sleeping there (MITMA)")):
        ax.plot([0, .09], [y, y], ls=ls, color=colour, lw=.85)
        ax.text(.12, y, text, fontsize=TEXT, va="center")
    ax.text(0, .04, "Flow ~ birth share", weight="bold", fontsize=LABEL, va="baseline")
    ax.text(0, -.03, "Compare cities of the same province against\n"
            "the same destination, with distance and income.", fontsize=TEXT, va="top")


def draw_coverage(fig, ax, manifest, months):
    panel(ax, "c", "Coverage by source", at=RIGHT)
    ax.set(xlim=(0, 366), ylim=(3.45, -.5))
    ax.set_yticks([0, 1, 2, 3], ["2022", "2023", "2024", "2025"])
    starts = pd.date_range("2024-01-01", periods=12, freq="MS").day_of_year - 1
    mids = starts + pd.date_range("2024-01-01", periods=12, freq="MS").days_in_month / 2
    ax.set_xticks(mids, list("JFMAMJJASOND"))
    ax.set_xticks(list(starts) + [366], minor=True)
    ax.tick_params(axis="x", which="major", length=0, pad=2)
    ax.tick_params(axis="x", which="minor", length=2, width=.5)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.grid(False)
    ax.fill_between([182, 366], 3 - .2, 3 + .3, facecolor="none", edgecolor=ORDINARY,
                    linewidth=0, hatch="//////", zorder=0)
    for year, nights in manifest.loc[manifest.file_present].groupby("year"):
        y = year - 2022
        for klass, subset in nights.groupby("coverage_class"):
            if klass in HOLIDAYS:
                ax.plot(subset.day, [y - .02] * len(subset), linestyle="none", marker=HOLIDAYS[klass][0],
                        ms=SMALL, color=INK, mec="white", mew=.4, zorder=5)
            else:
                col = TICKS[klass] if klass in TICKS else TICKS["other"]
                ax.vlines(subset.day, y - .2, y + .16, color=col, linewidth=.45)
        available = months.loc[months.year.eq(year)]
        ax.hlines(np.full(len(available), y + .3), available.start + 1.5, available.end - 1.5,
                  color=MUTED, lw=1.1)
    handles = [
        Line2D([], [], ls="", marker="|", ms=6, mew=.9, color=TICKS["O"], label="Ordinary night"),
        Line2D([], [], ls="", marker="|", ms=6, mew=.9, color=TICKS["AUG"], label="August window"),
        Line2D([], [], ls="", marker="|", ms=6, mew=.9, color=TICKS["other"], label="Other MITMA night"),
        Line2D([], [], color=MUTED, lw=1.1, label="INE month"),
    ] + [Line2D([], [], ls="", marker=mark, ms=SMALL, color=INK, label=label)
         for mark, label in HOLIDAYS.values()] + [
        Patch(facecolor="white", edgecolor=ORDINARY, hatch="//////", lw=0, label="New MITMA method"),
    ]
    order = [0, 1, 2, 3, 4, 5, 6, 7]
    ax.legend(handles=[handles[i] for i in order], loc="upper left", bbox_to_anchor=(-.1, -.2),
              ncol=2, handlelength=1.1, columnspacing=1.0, labelspacing=.5)


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
    fig = figure(106)
    keys = draw_map(fig, place(fig, 4, 10, 98, 76.4), balance, corridors, provinces)
    map_keys(fig, *keys)
    draw_dyad(place(fig, 110, 10, 57, 34))
    draw_coverage(fig, place(fig, 110, 59, 57, 26), manifest, months)
    save(fig, 1, balance.reset_index().assign(panel="a", series="province balance"),
         corridors.assign(panel="a", series="excess corridor"),
         manifest.loc[manifest.file_present, ["file_date", "year", "day", "is_valid", "coverage_class", "method"]]
         .assign(panel="c", series="MITMA night"),
         months.assign(panel="c", series="INE month"))


if __name__ == "__main__":
    main()
