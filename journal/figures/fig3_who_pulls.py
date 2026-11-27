"""August roots contrasts by cohort, destination, and age band."""

import subprocess
import sys

from matplotlib.lines import Line2D
import pandas as pd

from style import (
    A, A2, B, COAST, EXODUS, HERE, INLAND, NEUTRAL, OUT, RECENT,
    figure, forest, interval, panel, read, save,
)

GROUPS = {
    "P2": ("Arrived before 1981", EXODUS),
    "S5": ("Arrived 2016–2023", RECENT),
    "S4": ("Birthplace ≤ 2,000", EXODUS),
    "urban": ("Birthplace > 10,000", RECENT),
    "S6": ("Pre-1981 minus recent", NEUTRAL),
}
BASELINES = {"Aug15_vs_Feb14": ("o", -.16, "13–14 February"),
             "Aug15_vs_OctMean12": ("s", .16, "October weekdays")}


def cohort_data():
    first = read(A, "preregistered_tests").query(
        "split == 'confirmation_O' and year_scope == 'equal_year_avg'"
    ).assign(period="2022–24")
    second = read(A2, "preregistered_tests").assign(period="2025")
    urban_file = OUT / "fig3_urban_contrast.csv"
    if not urban_file.exists():
        subprocess.run([sys.executable, "-I", "-B", str(HERE / "urban_contrast.py")], check=True)
    urban = pd.read_csv(urban_file, dtype={"period": str}).assign(hypothesis="urban")
    preregistered = pd.concat([first, second], ignore_index=True).rename(
        columns={"ci_wcr_lo": "ci_lo", "ci_wcr_hi": "ci_hi"}
    ).assign(interval_type="wild-cluster bootstrap")
    urban = urban.rename(columns={"ci_crv_lo": "ci_lo", "ci_crv_hi": "ci_hi"}).assign(
        interval_type="CRV1"
    )
    fields = ["hypothesis", "term", "period", "est", "ci_lo", "ci_hi", "G", "n_obs", "interval_type"]
    return pd.concat([preregistered[fields], urban[fields]], ignore_index=True).query(
        "hypothesis in @GROUPS"
    )


def draw_cohorts(ax, data):
    panel(ax, "a", "Arrival cohort and birthplace size", x=-.61)
    forest(ax, [pair[0] for pair in GROUPS.values()], (-.087, .166), [-.05, 0, .05, .10, .15])
    for row in data.itertuples():
        i = list(GROUPS).index(row.hypothesis)
        filled = row.period != "2025"
        interval(ax, i + (-.16 if filled else .16), row.est, row.ci_lo, row.ci_hi,
                 color=GROUPS[row.hypothesis][1], filled=filled)
    ax.set_xlabel("August vs ordinary nights")
    ax.legend(handles=[
        Line2D([], [], ls="", marker="o", color=NEUTRAL, label="2022–24 confirmation"),
        Line2D([], [], ls="", marker="o", mfc="white", color=NEUTRAL, label="2025 hold-out"),
    ], loc="upper center", bbox_to_anchor=(.32, -1.165), ncol=2, columnspacing=1, handletextpad=.25)


def draw_destinations(ax, data):
    panel(ax, "b", "By destination type in 2024", x=-.41)
    kinds = {"all": ("All", NEUTRAL), "inland": ("Inland", INLAND),
             "coastal": ("Coastal", COAST), "common_19_provs": ("Inland − coastal", NEUTRAL)}
    group_names = ["pre1981", "born_village_le2k"]
    rows = [(group, kind) for group in group_names for kind in kinds]
    positions = [0, 1, 2, 3, 4.8, 5.8, 6.8, 7.8]
    forest(ax, [], (-.087, .241), [-.05, 0, .05, .10, .15, .20])
    ax.set_yticks(positions, [kinds[kind][0] for _, kind in rows])
    ax.set_ylim(8.36, -1.10)
    for pos, text in [(-.76, "Arrived before 1981"), (4.03, "Birthplace ≤ 2,000")]:
        ax.text(-.092, pos, text, ha="right", va="center", weight="bold", fontsize=7)
    for row in data.itertuples():
        y = positions[rows.index((row.group, row.dest_type))]
        marker, shift, _ = BASELINES[row.baseline]
        interval(ax, y + shift, row.point, row.ci_lo, row.ci_hi,
                 color=kinds[row.dest_type][1], marker=marker)
    ax.set_xlabel("August vs February or October")
    ax.legend(handles=[
        Line2D([], [], ls="", marker=shape, color=NEUTRAL, label=label)
        for shape, _, label in BASELINES.values()
    ], loc="upper center", bbox_to_anchor=(.31, -.17), ncol=2, columnspacing=1, handletextpad=.25)


def draw_age(ax, data):
    panel(ax, "c", "Within-age arrival gap", x=-.61)
    forest(ax, ["Aged 43–64", "Aged 65+"], (-.087, .166), [-.05, 0, .05, .10, .15])
    for row in data.itertuples():
        marker, shift, _ = BASELINES[row.contrast]
        position = 0 if row.term.startswith("within_43_64") else 1
        interval(ax, position + shift, row.delta, row.ci_wcr_lo, row.ci_wcr_hi, marker=marker)
    ax.set_xlabel("August vs February or October")


def main():
    cohorts = cohort_data()
    tidy = read(B, "fig_who_pulls_tidy")
    destinations = tidy.loc[
        tidy.group.isin(["pre1981", "born_village_le2k"]) & tidy.baseline.isin(BASELINES)
        & (tidy.dest_type.isin(["all", "inland", "coastal"]) & tidy.panel.eq("A_main_contrasts")
           | tidy.dest_type.eq("common_19_provs") & tidy.panel.eq("B_inland_x_coast"))
    ].copy()
    age = read(B, "b4_age_matched").query(
        "spec == 'A2_age_matched_binary_43_64_and_65plus' and sample == 'all' and is_within_band_diff"
    )
    assert (len(cohorts), len(destinations), len(age)) == (10, 16, 4)
    fig = figure(110)
    draw_cohorts(fig.add_axes([.196, .533, .282, .357]), cohorts)
    draw_destinations(fig.add_axes([.646, .23, .324, .66]), destinations)
    draw_age(fig.add_axes([.196, .23, .282, .148]), age)
    fig.text(.5, .042, "Change in roots elasticity", ha="center", fontsize=8, weight="bold")
    save(fig, 3, cohorts.assign(panel="a"), destinations.assign(panel="b"), age.assign(panel="c"))


if __name__ == "__main__":
    main()
