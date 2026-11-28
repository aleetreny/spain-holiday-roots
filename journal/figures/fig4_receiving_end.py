"""Net receiver geography and arrivals/departures in the most emptied fifth."""

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from matplotlib import patheffects
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import PercentFormatter

from style import (
    A, D, GRID, HEAVY, INK, LABEL, MAP_ZERO, MARK, MUTED, NEUTRAL, RECEIVER, SENDER, SMALL, TEXT,
    decimals, figure, map_axes, map_geometry, panel, place, read, save,
)

NIGHTS = {"feb14": ("2024-02-14", "13–14 February"), "aug15": ("2024-08-15", "14–15 August"),
          "dec25": ("2024-12-25", "Christmas Eve")}
CLASSES = ["Ordinary, Nov–Mar", "Ordinary, October", "Ordinary, Apr–Jun, Sep", "Fri–Sun",
           "August window", "Other August", "Christmas Eve", "Other holidays", "Other nights"]
LEFT, RIGHT = 2.0, 106.0   # mm: panel-letter columns
HALO = [patheffects.withStroke(linewidth=2.2, foreground="white")]


def draw_maps(fig, zone_results):
    geo, provinces = map_geometry()
    q5 = zone_results.query("loss_q == 5")
    rows = []
    for i, (stub, (date, label)) in enumerate(NIGHTS.items()):
        ax = place(fig, 4 + i * 56.5, 18, 53, 41.3)
        metric = f"imp_{stub}"
        data = geo.merge(q5[["id", metric]], on="id", how="left", validate="one_to_one")
        data["receives"] = data[metric].map({1.: "Net receiver", 0.: "Net sender"})
        for mask, colour in [(data.receives.isna(), MAP_ZERO),
                             (data.receives.eq("Net sender"), SENDER),
                             (data.receives.eq("Net receiver"), RECEIVER)]:
            data.loc[mask].plot(ax=ax, color=colour, linewidth=0)
        provinces.boundary.plot(ax=ax, color="white", linewidth=.4)
        map_axes(ax)
        share = q5[metric].mean()
        ax.annotate(label, (.5, 1), xycoords="axes fraction", xytext=(0, 10.5),
                    textcoords="offset points", ha="center", va="baseline", fontsize=LABEL)
        ax.annotate(f"{share:.1%} net receivers", (.5, 1), xycoords="axes fraction",
                    xytext=(0, 2.5), textcoords="offset points", ha="center", va="baseline",
                    fontsize=TEXT, color=MUTED)
        if i == 0:
            panel(ax, "a", "Net receivers in the most emptied fifth", at=LEFT, rise=10.5)
        rows.append(data[["id", "receives"]].assign(date=date, share_net_receiver=share,
                                                   panel="a", series="zone map"))
    ax.legend(handles=[Patch(color=c, label=t) for c, t in
                       [(RECEIVER, "Net receiver"), (SENDER, "Net sender"), (MAP_ZERO, "Other zones")]],
              loc="upper center", bbox_to_anchor=(.5, 0), bbox_transform=fig.transFigure, ncol=3,
              handlelength=1.1, handleheight=.9, columnspacing=1.8)
    ax.get_legend().set_bbox_to_anchor((.5, 1 - 61 / fig.get_figheight() / 25.4), fig.transFigure)
    return pd.concat(rows, ignore_index=True)


def calendar_data(zones):
    data = read(A, "descriptive_q5_villages_daily").copy()
    single = zones.query("loss_q == 5 and n_municipalities == 1")[["id"]]
    daily = pq.read_table(
        HEAVY / A / "panel/panel_zones2435_daily.parquet",
        columns=["zone_id", "date", "visitors_in", "residents_away"],
        filters=[("zone_id", "in", single.id.tolist()), ("date", "in", data.date.tolist())],
    ).to_pandas()
    daily["is_receiver"] = daily.visitors_in.gt(daily.residents_away)
    single_stats = daily.groupby("date").is_receiver.agg(["mean", "size"])
    assert single_stats["size"].eq(len(single)).all() and len(single_stats) == len(data)
    data = data.merge(single_stats, left_on="date", right_index=True, validate="one_to_one").rename(
        columns={"mean": "single_mun_share", "size": "single_mun_n"}
    )
    bucket = data.desc_bucket
    holiday = bucket.str.startswith(("ALLSAINTS", "EASTER", "EPIPHANY", "NEW_YEAR", "OTHER_holiday",
                                     "PUENTE", "XMAS_NEWYEAR"))
    month = pd.to_datetime(data.evening_date).dt.month
    ordinary = bucket.str.startswith("O_ordinary")
    data["class_index"] = np.select(
        [ordinary & month.isin([11, 12, 1, 2, 3]), ordinary & month.eq(10), ordinary,
         bucket.str.startswith("WEEKEND"), bucket.str.startswith("AUG_peak"), bucket.str.startswith("AUG"),
         bucket.eq("XMAS_eve_1225"), holiday],
        [0, 1, 2, 3, 4, 5, 6, 7], default=8,
    )
    data["night_category"] = data.class_index.map(dict(enumerate(CLASSES))).str.replace("\n", " ")
    return data


def draw_calendar(ax, data, boot):
    panel(ax, "b", "Net receiver share across nights", at=LEFT)
    rng = np.random.default_rng(475)
    data["jitter"] = rng.uniform(-.07, .07, len(data))
    for prop, shift, face, edge, alpha in (
        ("q5_share_net_receiver", -.2, NEUTRAL, NEUTRAL, .25),
        ("single_mun_share", .2, "none", MUTED, .45),
    ):
        ax.scatter(data[prop], data.class_index + shift + data.jitter, s=SMALL ** 2,
                   facecolor=face, edgecolor=edge, linewidth=.5, alpha=alpha)
    baselines = boot.query(
        "mode == 'level' and cls == 'Q5' and sample == 'all' and metric == 'share_net_importer' "
        "and period in ['2022-02-16', '2023-02-15', '2024-02-14', '2024-10-mean12']"
    ).copy()
    assert len(baselines) == 4
    ax.set_yticks(range(len(CLASSES)), CLASSES)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.grid(axis="y", visible=False)
    ax.grid(axis="x", visible=True)
    ax.set(ylim=(len(CLASSES) - .5, -.6), xlim=(0, 1.02), xticks=[0, .25, .5, .75, 1])
    ax.xaxis.set_major_formatter(PercentFormatter(1, decimals=0))
    ax.set_xlabel("Zones receiving more people than they send")
    ax.legend(handles=[
        Line2D([], [], ls="", marker="o", ms=SMALL, color=NEUTRAL, alpha=.6, label="All 475 zones"),
        Line2D([], [], ls="", marker="o", ms=SMALL, mfc="none", mec=MUTED, label="84 single-municipality zones"),
    ], loc="upper center", bbox_to_anchor=(.5, -.25), ncol=2)
    return baselines.assign(panel="b", series="selected ordinary baseline")


def draw_flows(ax, flows):
    panel(ax, "c", "Arrivals and departures per resident", at=RIGHT)
    width = .36
    for row in flows.itertuples():
        j = [x[0] for x in NIGHTS.values()].index(row.period)
        arriving = row.metric == "arrivals_per_res"
        x = j + (-width / 2 - .02 if arriving else width / 2 + .02)
        colour = RECEIVER if arriving else SENDER
        ax.bar(x, row.point, width=width, color=colour, edgecolor="white", linewidth=.6, zorder=2)
        ax.errorbar(x, row.point, yerr=[[row.point - row.ci_lo], [row.ci_hi - row.point]], fmt="none",
                    ecolor=INK, elinewidth=.7, capsize=0, zorder=3)
        ax.annotate(f"{row.point:.2f}", (x, row.ci_hi), xytext=(0, 2.5), textcoords="offset points",
                    ha="center", va="bottom", fontsize=TEXT)
    ax.legend(handles=[Patch(color=RECEIVER, label="Arrivals"),
                       Patch(color=SENDER, label="Departures")],
              loc="upper right", bbox_to_anchor=(1, 1.02), frameon=False, handlelength=1.1, handleheight=.9)
    ax.set_xticks(range(3), ["13–14 Feb", "14–15 Aug", "Christmas Eve"])
    ax.tick_params(axis="x", length=0, pad=4)
    ax.set(xlim=(-.55, 2.55), ylim=(0, .85), yticks=[0, .2, .4, .6, .8])
    ax.yaxis.set_major_formatter(decimals(1))
    ax.set_ylabel("People per resident")


def draw_inland(ax, controlled):
    """Adjusted inland changes as an aligned table under panel c."""
    ax.set_axis_off()
    ax.set(xlim=(0, 1), ylim=(0, 1))
    lookup = controlled.set_index(["period", "outcome"]).coef
    ax.text(0, 1, "Inland zones, adjusted for covariates", fontsize=TEXT, va="top", weight="bold")
    columns = (.62, .95)
    ax.text(0, .74, "Change vs 13–14 Feb", fontsize=TEXT, va="top", color=MUTED)
    for x, head in zip(columns, ("Arrivals", "Departures")):
        ax.text(x, .74, head, ha="right", va="top", fontsize=TEXT, color=MUTED)
    ax.plot([0, .95], [.5, .5], color=GRID, lw=.6)
    for k, (period, day) in enumerate([("Aug15_minus_Feb14", "14–15 Aug"),
                                       ("Dec25_minus_Feb14", "Christmas Eve")]):
        y = .42 - k * .25
        ax.text(0, y, day, va="top", fontsize=TEXT)
        for x, outcome in zip(columns, ("arrivals_per_res_diff", "departures_per_res_diff")):
            value = lookup.loc[period, outcome]
            ax.text(x, y, f"{value:+.3f}".replace("-", "−"), ha="right", va="top", fontsize=TEXT)


def main():
    zones = read(D, "d_zones_main_nights", dtype={"id": str, "prov": str})
    flows = read(D, "fig_receiving_end_tidy").query(
        "panel == 'C_decomp_quintile_level' and sample == 'all' and group == 'Q5' "
        "and stat == 'mean' and metric != 'net_per_res'"
    )
    dates = [value[0] for value in NIGHTS.values()]
    flows = flows.loc[flows.period.isin(dates)]
    controlled = read(D, "d_decomp_controlled_inland").query(
        "sample == 'inland' and period in ['Aug15_minus_Feb14', 'Dec25_minus_Feb14'] "
        "and outcome in ['arrivals_per_res_diff', 'departures_per_res_diff']"
    )
    data = calendar_data(zones)
    fig = figure(140)
    maps = draw_maps(fig, zones)
    baselines = draw_calendar(place(fig, 30, 78, 66, 44), data, read(D, "d_quintiles_main_boot"))
    draw_flows(place(fig, 118, 78, 50, 26), flows)
    draw_inland(place(fig, 110.2, 113, 58, 17), controlled)
    save(fig, 4, maps, data.assign(panel="b", series="nightly net receiver shares"), baselines,
         flows.assign(panel="c", series="zone means"), controlled.assign(panel="c", series="inland adjusted"))


if __name__ == "__main__":
    main()
