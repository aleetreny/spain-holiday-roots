"""Parental ties and the publication threshold in INE municipal travel data."""

import pandas as pd

from style import BAND, C, MUTED, NEUTRAL, figure, panel, read, reference, save


def main():
    calendar = read(C, "calendar", dtype={"baseline": str}).query("treatment == 'publication'").copy()
    calendar["date"] = pd.to_datetime(calendar.month)
    calendar["baseline_month"] = calendar.baseline.str[-2:]
    for field in ("estimate", "lo95", "hi95"):
        calendar[f"{field}_pp"] = calendar[field] * 100
    estimates = read(C, "estimates", dtype={"reference": str}).query(
        "sample == 'all_years' and treatment in ['interval_lower', 'interval_upper']"
    )
    fig = figure(109)
    for i, (month, label) in enumerate([("02", "February"), ("10", "October")]):
        ax = fig.add_axes([.105, .705 - .274 * i, .854, .184])
        points = calendar.loc[calendar.baseline_month.eq(month)].set_index("date").sort_index()
        points = points.reindex(pd.date_range("2019-07", "2026-09", freq="MS"))
        ax.fill_between(points.index, points.lo95_pp, points.hi95_pp, color=BAND, lw=0)
        ax.plot(points.index, points.estimate_pp, color=NEUTRAL)
        aug = points.dropna(subset=["estimate_pp"]).loc[lambda d: d.index.month == 8]
        ax.plot(aug.index, aug.estimate_pp, "o", color=NEUTRAL, ms=3.6)
        reference(ax)
        panel(ax, "ab"[i], f"Monthly publication probability vs {label}", x=-.093)
        ax.set_ylim(-6, 8.6)
        ax.set_yticks([-4, 0, 4])
        ax.set_ylabel("Change (pp)")
        years = pd.to_datetime([f"{year}-01" for year in range(2020, 2027)])
        ax.set_xticks(years, years.year)
        ax.set_xlim(pd.Timestamp("2019-06"), pd.Timestamp("2026-10"))
    ax = fig.add_axes([.105, .136, .854, .177])
    panel(ax, "c", "August trip-count bounds under suppression", x=-.093)
    reference(ax, vertical=True)
    for y, month in enumerate(["02", "10"]):
        rows = estimates.loc[estimates.reference.eq(month)].set_index("treatment")
        lower, upper = rows.loc["interval_lower"], rows.loc["interval_upper"]
        ax.hlines(y, lower.lo95, upper.hi95, color="#aaaaaa", lw=4)
        ax.hlines(y, lower.estimate, upper.estimate, color=NEUTRAL, lw=1)
        ax.plot([lower.estimate, upper.estimate], [y] * 2, "o", color=NEUTRAL, ms=3.6)
    ax.set(yticks=[0, 1], yticklabels=["February", "October"], xlim=(-105, 110), ylim=(1.55, -.57),
           xlabel="Change in trips per city–village pair")
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.text(.02, .063, "Outer 95% envelope", fontsize=7, color=MUTED, transform=ax.transAxes)
    save(fig, "S2", calendar.assign(panel=calendar.baseline_month.map({"02": "a", "10": "b"})),
         estimates.assign(panel="c"))


if __name__ == "__main__":
    main()
