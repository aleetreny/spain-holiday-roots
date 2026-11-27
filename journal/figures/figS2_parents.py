"""Parental ties and the publication threshold in INE municipal travel data."""

from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import pandas as pd

from style import (
    BAND, C, INK, MARK, NEUTRAL, decimals, figure, panel, place, read, reference, save,
)

LEFT, AX_LEFT, AX_WIDTH = 2.0, 20.0, 147.0   # mm
ENVELOPE = "#c8c8c8"


def title_legend(ax, handles):
    ax.legend(handles=handles, loc="lower right", bbox_to_anchor=(1, 1), ncol=len(handles),
              borderaxespad=.55)


def main():
    calendar = read(C, "calendar", dtype={"baseline": str}).query("treatment == 'publication'").copy()
    calendar["date"] = pd.to_datetime(calendar.month)
    calendar["baseline_month"] = calendar.baseline.str[-2:]
    for field in ("estimate", "lo95", "hi95"):
        calendar[f"{field}_pp"] = calendar[field] * 100
    estimates = read(C, "estimates", dtype={"reference": str}).query(
        "sample == 'all_years' and treatment in ['interval_lower', 'interval_upper']"
    )
    fig = figure(113)
    for i, (month, label) in enumerate([("02", "February"), ("10", "October")]):
        ax = place(fig, AX_LEFT, 10 + 36 * i, AX_WIDTH, 25)
        points = calendar.loc[calendar.baseline_month.eq(month)].set_index("date").sort_index()
        points = points.reindex(pd.date_range("2019-07", "2026-09", freq="MS"))
        ax.fill_between(points.index, points.lo95_pp, points.hi95_pp, color=BAND, lw=0)
        ax.plot(points.index, points.estimate_pp, color=NEUTRAL)
        aug = points.dropna(subset=["estimate_pp"]).loc[lambda d: d.index.month == 8]
        ax.plot(aug.index, aug.estimate_pp, "o", color=INK, ms=MARK)
        reference(ax)
        panel(ax, "ab"[i], f"Monthly publication probability vs {label}", at=LEFT)
        if i == 0:
            title_legend(ax, [Line2D([], [], color=NEUTRAL, label="Monthly change"),
                              Line2D([], [], ls="", marker="o", ms=MARK, color=INK, label="August"),
                              Patch(color=BAND, label="95% CI")])
        ax.set_ylim(-6, 8.6)
        ax.set_yticks([-4, 0, 4, 8])
        ax.yaxis.set_major_formatter(decimals(0))
        ax.set_ylabel("Change (pp)")
        years = pd.to_datetime([f"{year}-01" for year in range(2020, 2027)])
        ax.set_xticks(years, years.year)
        ax.set_xlim(pd.Timestamp("2019-06"), pd.Timestamp("2026-10"))
    ax = place(fig, AX_LEFT, 86, AX_WIDTH, 17)
    panel(ax, "c", "August trip-count bounds under suppression", at=LEFT)
    reference(ax, vertical=True)
    for y, month in enumerate(["02", "10"]):
        rows = estimates.loc[estimates.reference.eq(month)].set_index("treatment")
        lower, upper = rows.loc["interval_lower"], rows.loc["interval_upper"]
        ax.plot([lower.lo95, upper.hi95], [y, y], color=ENVELOPE, lw=4, solid_capstyle="butt")
        ax.hlines(y, lower.estimate, upper.estimate, color=INK, lw=.85)
        ax.plot([lower.estimate, upper.estimate], [y] * 2, "o", color=INK, ms=MARK)
    title_legend(ax, [Line2D([], [], ls="-", marker="o", ms=MARK, lw=.85, color=INK,
                             label="Estimated bounds"),
                      Line2D([], [], color=ENVELOPE, lw=4, solid_capstyle="butt",
                             label="Outer 95% envelope")])
    ax.set(yticks=[0, 1], yticklabels=["vs February", "vs October"], xlim=(-105, 110),
           ylim=(1.55, -.57), xlabel="Change in trips per city–village pair")
    ax.xaxis.set_major_formatter(decimals(0))
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.grid(axis="y", visible=False)
    save(fig, "S2", calendar.assign(panel=calendar.baseline_month.map({"02": "a", "10": "b"})),
         estimates.assign(panel="c"))


if __name__ == "__main__":
    main()
