"""Net receiver geography and arrivals/departures in the most emptied fifth."""

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from matplotlib.patches import Rectangle
from matplotlib.ticker import PercentFormatter

from style import (
    A, D, HEAVY, MAP_ZERO, MUTED, NEUTRAL, ORDINARY, RECEIVER, SENDER,
    figure, map_axes, map_geometry, panel, read, save,
)

NIGHTS = {"feb14": ("2024-02-14", "13–14 February"), "aug15": ("2024-08-15", "14–15 August"),
          "dec25": ("2024-12-25", "Christmas Eve")}
CLASSES = ["Ordinary", "Fri–Sun", "August\nwindow", "Other\nAugust",
           "Christmas\nEve", "Other\nholidays", "Other\nnights"]


def draw_maps(fig, zone_results):
    geo, provinces = map_geometry()
    q5 = zone_results.query("loss_q == 5")
    rows = []
    for i, (stub, (date, label)) in enumerate(NIGHTS.items()):
        ax = fig.add_axes([.035 + i * .315, .535, .301, .318])
        metric = f"imp_{stub}"
        data = geo.merge(q5[["id", metric]], on="id", how="left", validate="one_to_one")
        data["receives"] = data[metric].map({1.: "Net receiver", 0.: "Net sender"})
        for mask, colour in [(data.receives.isna(), MAP_ZERO),
                             (data.receives.eq("Net sender"), SENDER),
                             (data.receives.eq("Net receiver"), RECEIVER)]:
            data.loc[mask].plot(ax=ax, color=colour)
        provinces.boundary.plot(ax=ax, color=ORDINARY, linewidth=.34)
        map_axes(ax)
        share = q5[metric].mean()
        ax.set_title(f"{label}   {share:.1%}", pad=2)
        rows.append(data[["id", "receives"]].assign(date=date, share_net_receiver=share,
                                                   panel="a", series="zone map"))
    fig.text(.017, .938, "a", fontsize=10, weight="bold")
    fig.text(.045, .941, "Net receivers in the most emptied fifth", fontsize=8.5)
    for x, colour, text in [(.19, RECEIVER, "Net receiver"), (.43, SENDER, "Net sender"),
                            (.65, MAP_ZERO, "Other zones")]:
        fig.patches.append(Rectangle(
            (x, .517), .019, .016, transform=fig.transFigure, color=colour, linewidth=0))
        fig.text(x + .025, .525, text, fontsize=7, va="center")
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
    data["class_index"] = np.select(
        [bucket.str.startswith("O_ordinary"), bucket.str.startswith("WEEKEND"),
         bucket.str.startswith("AUG_peak"), bucket.str.startswith("AUG"),
         bucket.eq("XMAS_eve_1225"), holiday],
        [0, 1, 2, 3, 4, 5], default=6,
    )
    data["night_category"] = data.class_index.map(dict(enumerate(CLASSES))).str.replace("\n", " ")
    return data


def draw_calendar(ax, data, boot):
    panel(ax, "b", "Net receiver share across nights", x=-.17)
    rng = np.random.default_rng(475)
    data["jitter"] = rng.uniform(-.125, .125, len(data))
    for prop, shift, face, edge, alpha in (
        ("q5_share_net_receiver", -.16, NEUTRAL, NEUTRAL, .25),
        ("single_mun_share", .16, "none", MUTED, .32),
    ):
        ax.scatter(data.class_index + shift + data.jitter, data[prop], s=3.5 ** 2,
                   facecolor=face, edgecolor=edge, linewidth=.5, alpha=alpha)
    baselines = boot.query(
        "mode == 'level' and cls == 'Q5' and sample == 'all' and metric == 'share_net_importer' "
        "and period in ['2022-02-16', '2023-02-15', '2024-02-14', '2024-10-mean12']"
    ).copy()
    lo, hi = baselines["mean"].min(), baselines["mean"].max()
    assert len(baselines) == 4
    ax.vlines(-.49, lo, hi, color=NEUTRAL, lw=.8)
    ax.hlines([lo, hi], -.55, -.44, color=NEUTRAL, lw=.8)
    for y, align in [(lo, "top"), (hi, "bottom")]:
        ax.text(-.49, y, f"{y:.0%}", ha="center", va=align, fontsize=7, color=NEUTRAL,
                bbox={"fc": "white", "ec": "none", "pad": .25})
    ax.set_xticks(range(len(CLASSES)), CLASSES)
    ax.tick_params(axis="x", length=0)
    ax.set(xlim=(-.81, 6.55), ylim=(0, 1.07), yticks=[0, .25, .5, .75, 1])
    ax.yaxis.set_major_formatter(PercentFormatter(1, decimals=0))
    ax.text(.97, -.27, "●  All 475 zones    ○  Single-municipality zones",
            transform=ax.transAxes, ha="right", fontsize=7, va="top")
    ax.text(.97, -.38, "Bracket: selected ordinary baselines", color=MUTED,
            transform=ax.transAxes, ha="right", fontsize=7, va="top")
    return baselines.assign(panel="b", series="selected ordinary baseline")


def draw_flows(ax, flows, controlled):
    panel(ax, "c", "Arrivals and departures per resident", x=-.17)
    for row in flows.itertuples():
        j = [x[0] for x in NIGHTS.values()].index(row.period)
        arriving = row.metric == "arrivals_per_res"
        x = j + (-.12 if arriving else .12)
        ax.errorbar(x, row.point, yerr=[[row.point - row.ci_lo], [row.ci_hi - row.point]],
                    fmt="o" if arriving else "s", color=NEUTRAL, ms=3.5, mew=.7,
                    elinewidth=.85, capsize=2.7, capthick=.65, zorder=3)
    ax.set_xticks(range(3), ["13–14 Feb", "14–15 Aug", "Christmas Eve"])
    ax.set(xlim=(-.45, 2.42), ylim=(0, .82), yticks=[0, .2, .4, .6, .8])
    ax.text(.045, .92, "●  Arrivals\n■  Departures", transform=ax.transAxes, va="top", fontsize=7)
    lookup = controlled.set_index(["period", "outcome"]).coef
    for j, (period, day) in enumerate([("Aug15_minus_Feb14", "August"), ("Dec25_minus_Feb14", "Christmas Eve")]):
        arrivals = lookup.loc[period, "arrivals_per_res_diff"]
        departures = lookup.loc[period, "departures_per_res_diff"]
        ax.text(.0, -.38 - j * .14, f"{day}   {arrivals:+.3f} arrivals, {departures:+.3f} departures",
                transform=ax.transAxes, fontsize=7)
    ax.text(.0, -.24, "Inland change vs February, adjusted for covariates",
            transform=ax.transAxes, fontsize=7, weight="bold")


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
    fig = figure(120)
    maps = draw_maps(fig, zones)
    baselines = draw_calendar(fig.add_axes([.083, .174, .435, .259]), data,
                              read(D, "d_quintiles_main_boot"))
    draw_flows(fig.add_axes([.599, .174, .369, .259]), flows, controlled)
    save(fig, 4, maps, data.assign(panel="b", series="nightly net receiver shares"), baselines,
         flows.assign(panel="c", series="zone means"), controlled.assign(panel="c", series="inland adjusted"))


if __name__ == "__main__":
    main()
