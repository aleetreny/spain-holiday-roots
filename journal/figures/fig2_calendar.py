"""Nightly roots estimates, pre-registered August test and cohort calendar."""

import numpy as np
import pandas as pd
from matplotlib import patheffects
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from fig1_design import HOLIDAYS
from style import (
    A, A2, BAND, EXODUS, GRID, INK, LABEL, MARK, NEUTRAL, ORDINARY, RECENT, SMALL, TEXT,
    across_years, decimals, figure, interval, panel, place, read, reference, save,
)

BRIDGE = pd.Timestamp("2025-07-01")
HALO = [patheffects.withStroke(linewidth=2.2, foreground="white")]


def nightly(name):
    data = across_years(name).query("sample == 'all'").copy()
    data["date"] = pd.to_datetime(data.date)
    data["klass"] = data.night_class.str.replace("2025", "", regex=False)
    data["calendar_day"] = pd.to_datetime("2024-" + data.date.dt.strftime("%m-%d")).dt.day_of_year - 1
    return data.sort_values("date")


LEFT, AX_LEFT, AX_WIDTH = 2.0, 20.0, 138.0   # mm
ROW_TOP, ROW_HEIGHT, ROW_GAP = 11.0, 16.0, 3.2


def title_legend(ax, handles, ncol):
    """Key on the panel-title baseline, right-aligned with the axes."""
    ax.legend(handles=handles, loc="lower right", bbox_to_anchor=(1, 1), ncol=ncol,
              borderaxespad=.55, handlelength=1.2, columnspacing=1.1)


def pooled_calendar(fig, data):
    means = []
    labelled = set()
    starts = pd.date_range("2024-01-01", periods=12, freq="MS")
    for i, (year, rows) in enumerate(data.groupby("year")):
        ax = place(fig, AX_LEFT, ROW_TOP + i * (ROW_HEIGHT + ROW_GAP), AX_WIDTH, ROW_HEIGHT)
        plot = rows.set_index("calendar_day").reindex(range(366))
        x = plot.index
        ax.fill_between(x, plot.ci_lo, plot.ci_hi, color=BAND, alpha=.8, lw=0)
        changing = plot.klass.ne("O")
        ax.plot(x, plot.beta.where(changing), color=NEUTRAL, lw=.75)
        ordinary = rows.query("klass == 'O'")
        ax.plot(ordinary.calendar_day, ordinary.beta, linestyle="none", marker="o",
                ms=SMALL, mfc="white" if year == 2025 else ORDINARY, mec=ORDINARY, mew=.55)
        avg = ordinary.beta.mean()
        means.append({"year": year, "ordinary_mean": avg, "n_ordinary": len(ordinary)})
        reference(ax, avg)
        aug = rows.query("klass == 'AUG'")
        ax.plot(aug.calendar_day, aug.beta, linestyle="none", marker="o", ms=MARK,
                mec=INK, mew=.7, mfc="white" if year == 2025 else INK, zorder=5)
        for klass, (mark, label) in HOLIDAYS.items():
            holidays = rows.loc[rows.klass.eq(klass)]
            ax.plot(holidays.calendar_day, holidays.beta, marker=mark, ms=MARK,
                    linestyle="none", mec=INK, mew=.7,
                    mfc="white" if year == 2025 else INK, zorder=5)
            if len(holidays) and klass not in labelled:
                anchor = holidays.iloc[0]
                ax.annotate(label, (anchor.calendar_day, anchor.beta), xytext=(0, 6),
                            textcoords="offset points", ha="center", va="bottom",
                            fontsize=TEXT, zorder=6, path_effects=HALO)
                labelled.add(klass)
        ax.set(xlim=(0, 366), ylim=(.27, .9), yticks=[.4, .6, .8])
        ax.yaxis.set_major_formatter(decimals(1))
        ax.text(1.012, .5, str(year), transform=ax.transAxes, ha="left", va="center",
                fontsize=LABEL, weight="bold")
        ax.set_xticks(starts.day_of_year - 1, starts.strftime("%b") if year == 2025 else [])
        ax.tick_params(axis="x", length=2 if year == 2025 else 0)
        if year == 2025:
            ax.axvline(182, color=NEUTRAL, linestyle=":", lw=.75)
            ax.text(185, .86, "1 July 2025: new MITMA method", va="top", fontsize=TEXT,
                    path_effects=HALO)
        if i == 0:
            panel(ax, "a", "Nightly roots elasticity", at=LEFT)
            title_legend(ax, [
                Line2D([], [], ls="", marker="o", ms=SMALL, mfc=ORDINARY, mec=ORDINARY, label="Ordinary night"),
                Line2D([], [], ls="", marker="o", ms=MARK, color=INK, label="August window"),
                Line2D([], [], color=NEUTRAL, lw=.75, label="Other nights"),
                Patch(color=BAND, label="95% CI"),
            ], ncol=4)
        if i == 1:
            first = ax
    fig.text(7 / 170, (first.get_position().y0 + first.get_position().y1) / 2 - .5 * (ROW_HEIGHT + ROW_GAP) / fig.get_figheight() / 25.4,
             "Roots elasticity", rotation=90, ha="center", va="center", fontsize=LABEL)
    return pd.DataFrame(means).assign(panel="a", series="ordinary reference")


def august_test(ax):
    panel(ax, "b", "August premium and ordinary-night placebos", at=LEFT)
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
            color=ORDINARY, alpha=.7, ms=SMALL, mew=0)
    for row in tests.itertuples():
        interval(ax, row.x, row.est, row.lo, row.hi, horizontal=False, filled=row.label != "2025",
                 color=INK)
        annual = placebos.loc[placebos.year.astype(str).eq(row.label), "delta_P1"]
        assert annual.empty or row.est > annual.max()
    reference(ax)
    ax.axvline(2.575, color=GRID, lw=.6, zorder=0)
    ax.set(xlim=(-.42, 4.72), ylim=(min(-.09, placebos.delta_P1.min() - .01), .25),
           yticks=[-.05, 0, .05, .10, .15, .20, .25])
    ax.yaxis.set_major_formatter(decimals(2))
    ax.set_ylabel("Change in roots elasticity")
    ax.yaxis.set_label_coords(-13 / AX_WIDTH, .5)
    ax.set_xticks(list(positions.values()), ["2022", "2023", "2024", "Pooled 2022–24", "2025 hold-out"])
    ax.tick_params(axis="x", length=0, pad=4)
    title_legend(ax, [
        Line2D([], [], ls="", marker="o", ms=MARK, color=INK, label="August vs ordinary, 2022–24"),
        Line2D([], [], ls="", marker="o", ms=MARK, mfc="white", mec=INK, label="2025 hold-out"),
        Line2D([], [], ls="", marker="o", ms=SMALL, color=ORDINARY, label="Ordinary-night placebo"),
    ], ncol=3)
    return tests.assign(panel="b", series="August contrast"), placebos.assign(panel="b", series="placebo")


def cohort_calendar(ax, data):
    panel(ax, "c", "Roots elasticity by arrival cohort", at=LEFT)
    spans = []
    for year in range(2022, 2026):
        spans.extend([(f"{year}-08-10", f"{year}-08-19", "August window"),
                      (f"{year - 1}-12-25", f"{year}-01-07", "Christmas to Epiphany")])
    spans = pd.DataFrame(spans, columns=["start", "end", "series"])
    for row in spans.itertuples():
        ax.axvspan(pd.Timestamp(row.start), pd.Timestamp(row.end), facecolor=BAND, lw=0, zorder=0)
    groups = {"pre1981": (EXODUS, "Arrived before 1981"),
              "recent_2016_2023": (RECENT, "Arrived 2016–2023")}
    for term, (colour, label) in groups.items():
        rows = data.loc[data.term.eq(term)]
        train = rows.loc[rows.year.lt(2025)]
        holdout = rows.loc[rows.year.eq(2025)]
        ax.plot(train.date, train.beta, linestyle="none", marker="o",
                color=colour, alpha=.16, ms=SMALL, mew=0)
        ax.plot(holdout.date, holdout.beta, linestyle="none", marker="o", ms=SMALL,
                mfc="none", mec=colour, mew=.5, alpha=.35)
        for _, segment in rows.groupby("segment"):
            ax.plot(segment.date, segment.rolling_median, color=colour, lw=1.0)
    title_legend(ax, [Line2D([], [], color=c, lw=1.0, label=l) for c, l in groups.values()]
                 + [Patch(color=BAND, label="August window; Christmas–Epiphany")], ncol=3)
    ax.axvline(BRIDGE, color=NEUTRAL, linestyle=":", lw=.75)
    ax.text(BRIDGE + pd.Timedelta(days=6), .315, "New method", va="top", fontsize=TEXT,
            path_effects=HALO)
    ax.set_xlim(pd.Timestamp("2022-01-01"), pd.Timestamp("2026-01-01"))
    ax.set_ylim(-.01, .33)
    ax.set_yticks([0, .1, .2, .3])
    ax.yaxis.set_major_formatter(decimals(1))
    ax.set_ylabel("Cohort term")
    ax.yaxis.set_label_coords(-13 / AX_WIDTH, .5)
    ax.set_xticks(pd.to_datetime([f"{y}-01-01" for y in range(2022, 2027)]),
                  ["2022", "2023", "2024", "2025", "2026"])
    ax.set_xticks(pd.to_datetime([f"{y}-07-01" for y in range(2022, 2026)]), minor=True)
    ax.tick_params(axis="x", which="minor", length=1.5, width=.5)
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

    fig = figure(182)
    refs = pooled_calendar(fig, pooled)
    contrasts, placebos = august_test(place(fig, AX_LEFT, 101, AX_WIDTH, 34))
    spans = cohort_calendar(place(fig, AX_LEFT, 150, AX_WIDTH, 24), cohorts)
    save(fig, 2, pooled.assign(panel="a", series="nightly pooled"), refs, contrasts, placebos,
         cohorts.assign(panel="c", series="nightly cohort"), spans)


if __name__ == "__main__":
    main()
