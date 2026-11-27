"""Nightly roots estimates, pre-registered August test and cohort calendar."""

import numpy as np
import pandas as pd

from fig1_design import HOLIDAYS
from style import (
    A, A2, BAND, EXODUS, MUTED, NEUTRAL, ORDINARY, RECENT,
    across_years, figure, interval, panel, read, reference, save,
)

BRIDGE = pd.Timestamp("2025-07-01")


def nightly(name):
    data = across_years(name).query("sample == 'all'").copy()
    data["date"] = pd.to_datetime(data.date)
    data["klass"] = data.night_class.str.replace("2025", "", regex=False)
    data["calendar_day"] = pd.to_datetime("2024-" + data.date.dt.strftime("%m-%d")).dt.day_of_year - 1
    return data.sort_values("date")


def pooled_calendar(fig, data):
    means = []
    labelled = set()
    for i, (year, rows) in enumerate(data.groupby("year")):
        ax = fig.add_axes([.122, .807 - i * .094, .846, .084])
        plot = rows.set_index("calendar_day").reindex(range(366))
        x = plot.index
        ax.fill_between(x, plot.ci_lo, plot.ci_hi, color=BAND, alpha=.72, lw=0)
        changing = plot.klass.ne("O")
        ax.plot(x, plot.beta.where(changing), color=NEUTRAL, lw=.78)
        ordinary = rows.query("klass == 'O'")
        ax.plot(ordinary.calendar_day, ordinary.beta, linestyle="none", marker="o",
                ms=3.5, mfc="white" if year == 2025 else ORDINARY, mec=ORDINARY, mew=.55)
        avg = ordinary.beta.mean()
        means.append({"year": year, "ordinary_mean": avg, "n_ordinary": len(ordinary)})
        reference(ax, avg)
        aug = rows.query("klass == 'AUG'")
        ax.plot(aug.calendar_day, aug.beta, linestyle="none", marker="o", ms=3.6,
                mec=NEUTRAL, mew=.7, mfc="white" if year == 2025 else NEUTRAL)
        for klass, (mark, label) in HOLIDAYS.items():
            holidays = rows.loc[rows.klass.eq(klass)]
            ax.plot(holidays.calendar_day, holidays.beta, marker=mark, ms=4.0,
                    linestyle="none", mec=NEUTRAL, mew=.7,
                    mfc="white" if year == 2025 else NEUTRAL, zorder=5)
            if len(holidays) and klass not in labelled:
                anchor = holidays.iloc[0]
                offsets = {"EASTER": (5, 6), "ALLSAINTS": (-4, 7), "XMAS": (-6, 7)}
                dx, dy = offsets[klass]
                ax.annotate(label, (anchor.calendar_day, anchor.beta),
                            xytext=(dx, dy), textcoords="offset points",
                            ha="right" if dx < 0 else "left", fontsize=7, zorder=6,
                            bbox={"fc": "white", "ec": "none", "alpha": .9, "pad": .6})
                labelled.add(klass)
        ax.set(xlim=(0, 365), ylim=(.27, .96), yticks=[.4, .7])
        ax.text(-.045, .53, str(year), transform=ax.transAxes, ha="right", va="center", fontsize=7)
        ax.tick_params(axis="y", pad=2)
        # Place the year outside the tick labels.
        ax.yaxis.set_ticks_position("right")
        ax.tick_params(axis="y", labelright=False, right=False)
        ax.text(1.008, (.4 - .27) / (.96 - .27), "0.4", transform=ax.transAxes, va="center", fontsize=7)
        ax.text(1.008, (.7 - .27) / (.96 - .27), "0.7", transform=ax.transAxes, va="center", fontsize=7)
        ax.spines["left"].set_visible(False)
        if year == 2025:
            ax.axvline(182, color=NEUTRAL, linestyle=":", lw=.75)
            ax.text(185, .91, "1 July: new method", va="top", fontsize=7)
            starts = pd.date_range("2024-01-01", periods=12, freq="MS")
            ax.set_xticks(starts.day_of_year - 1, starts.strftime("%b"))
        else:
            ax.set_xticks([])
        if i == 0:
            panel(ax, "a", x=-.125, y=1.13)
            ax.text(0, 1.16, "Nightly roots elasticity", transform=ax.transAxes,
                    fontsize=8.5, va="bottom")
    return pd.DataFrame(means).assign(panel="a", series="ordinary reference")


def august_test(ax):
    panel(ax, "b", "August premium and ordinary-night placebos", x=-.125)
    confirmed = read(A, "preregistered_tests").query("hypothesis == 'P1' and split == 'confirmation_O'")
    holdout = read(A2, "preregistered_tests").query("hypothesis == 'P1'")
    tests = pd.concat([confirmed, holdout], ignore_index=True)
    tests["label"] = ["2022", "2023", "2024", "2022–24", "2025"]
    tests["lo"] = tests.ci_wcr_lo.fillna(tests.ci_crv_lo)
    tests["hi"] = tests.ci_wcr_hi.fillna(tests.ci_crv_hi)
    tests["interval"] = np.where(tests.ci_wcr_lo.isna(), "CRV1", "wild-cluster bootstrap")
    positions = {"2022": 0, "2023": 1, "2024": 2, "2022–24": 3.15, "2025": 4.3}
    tests["x"] = tests.label.map(positions)
    placebos = across_years("placebo_windows").query(
        "window_type == 'placebo' and split in ['confirmation_O', 'confirmation_ii']"
    ).copy()
    placebos["x"] = placebos.year.astype(str).map(positions)
    rng = np.random.default_rng(20261124)
    placebos["jittered_x"] = placebos.x + rng.uniform(-.15, .15, len(placebos))
    ax.plot(placebos.jittered_x, placebos.delta_P1, linestyle="none", marker="o",
            color=ORDINARY, alpha=.61, ms=3.5, mew=0)
    for row in tests.itertuples():
        interval(ax, row.x, row.est, row.lo, row.hi, horizontal=False, filled=row.label != "2025")
        annual = placebos.loc[placebos.year.astype(str).eq(row.label), "delta_P1"]
        assert annual.empty or row.est > annual.max()
    reference(ax)
    ax.set(xlim=(-.42, 4.72), ylim=(min(-.09, placebos.delta_P1.min() - .01), .25),
           yticks=[-.05, .05, .15, .25])
    ax.set_ylabel("Change in roots elasticity\nAugust vs ordinary nights", fontsize=7.2)
    ax.yaxis.set_label_coords(-.085, .46)
    ax.set_xticks(list(positions.values()), ["2022", "2023", "2024", "Pooled 2022–24", "2025 hold-out"])
    ax.text(2.76, -.047, "Grey dots: ordinary-night placebos", fontsize=7, color=MUTED)
    return tests.assign(panel="b", series="August contrast"), placebos.assign(panel="b", series="placebo")


def cohort_calendar(ax, data):
    panel(ax, "c", "Roots elasticity by arrival cohort", x=-.125)
    spans = []
    for year in range(2022, 2026):
        spans.extend([(f"{year}-08-10", f"{year}-08-19", "August window"),
                      (f"{year - 1}-12-25", f"{year}-01-07", "Christmas to Epiphany")])
    spans = pd.DataFrame(spans, columns=["start", "end", "series"])
    for row in spans.itertuples():
        ax.axvspan(pd.Timestamp(row.start), pd.Timestamp(row.end), facecolor=BAND, alpha=.55, lw=0)
    groups = {"pre1981": (EXODUS, "Arrived before 1981"),
              "recent_2016_2023": (RECENT, "Arrived 2016–2023")}
    for i, (term, (colour, label)) in enumerate(groups.items()):
        rows = data.loc[data.term.eq(term)]
        train = rows.loc[rows.year.lt(2025)]
        holdout = rows.loc[rows.year.eq(2025)]
        ax.plot(train.date, train.beta, linestyle="none", marker="o",
                color=colour, alpha=.14, ms=3.5, mew=0)
        ax.plot(holdout.date, holdout.beta, linestyle="none", marker="o", ms=3.5,
                mfc="none", mec=colour, mew=.5, alpha=.32)
        for _, segment in rows.groupby("segment"):
            ax.plot(segment.date, segment.rolling_median, color=colour, lw=.85)
        x = .44 + i * .26
        ax.plot([x, x + .035], [1.14] * 2, color=colour, transform=ax.transAxes, clip_on=False)
        ax.text(x + .049, 1.14, label, transform=ax.transAxes, va="center", fontsize=7)
    ax.axvline(BRIDGE, color=NEUTRAL, linestyle=":", lw=.75)
    ax.set_xlim(pd.Timestamp("2022-01-01"), pd.Timestamp("2026-01-01"))
    ax.set_ylim(-.03, .33)
    ax.set_yticks([0, .1, .2, .3])
    ax.set_xticks(pd.to_datetime([f"{y}-01-01" for y in range(2022, 2026)]),
                  ["2022", "2023", "2024", "2025"])
    return spans.assign(panel="c")


def main():
    pooled = nightly("nightly_pooled")
    cohorts = nightly("nightly_cohorts").query("model == 'cohort_joint'").copy()
    # Split at gaps and the method change so smoothed lines stop at coverage boundaries.
    for _, rows in cohorts.groupby("term"):
        gap = rows.date.diff().gt(pd.Timedelta(days=1)) | rows.date.eq(BRIDGE)
        cohorts.loc[rows.index, "segment"] = gap.cumsum()
        cohorts.loc[rows.index, "rolling_median"] = rows.beta.groupby(gap.cumsum()).transform(
            lambda series: series.rolling(7, center=True, min_periods=7).median()
        )

    fig = figure(150)
    refs = pooled_calendar(fig, pooled)
    contrasts, placebos = august_test(fig.add_axes([.122, .258, .846, .179]))
    ax = fig.add_axes([.122, .062, .846, .125])
    spans = cohort_calendar(ax, cohorts)
    save(fig, 2, pooled.assign(panel="a", series="nightly pooled"), refs, contrasts, placebos,
         cohorts.assign(panel="c", series="nightly cohort"), spans)


if __name__ == "__main__":
    main()
