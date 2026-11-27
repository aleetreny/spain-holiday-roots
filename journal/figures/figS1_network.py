"""Coursework network diagnostics at the ≥10-person link threshold."""

from matplotlib.lines import Line2D
import pandas as pd

from style import (
    COAST, INK, INLAND, LABEL, MAIN, MARK, MUTED, NEUTRAL, ORDINARY, SMALL, TEXT,
    decimals, figure, nodes, panel, place, save,
)

DATES = {"2024-02-14": "13–14 February", "2024-08-15": "14–15 August",
         "2024-11-01": "All Saints"}
LEFT, RIGHT = 2.0, 88.0   # mm: panel-letter columns
RANGE = "#dddddd"


def zone_groups():
    """Destination-type roles shared with the main figures: inland, coast, large cities."""
    data = nodes()[["idx", "population", "coastal_true"]].copy()
    big = data.population.ge(100_000)
    data["group"] = "inland"
    data.loc[~big & data.coastal_true, "group"] = "coast"
    data.loc[big, "group"] = "city"
    return data.set_index("idx").group


GROUPS = {"inland": (INLAND, .28, SMALL, "Inland zone"),
          "coast": (COAST, .55, SMALL, "Coastal zone"),
          "city": (INK, .95, MARK, "City over 100,000")}


def degree(ax, points, title, letter, at, y_label=False):
    rho = points.kin_w10.corr(points.kout_w10, method="spearman")
    panel(ax, letter, title, at=at)
    ax.plot([1, 3000], [1, 3000], ls=(0, (3, 3)), lw=.55, color=NEUTRAL, zorder=1)
    groups = points.idx.map(zone_groups())
    for group, (colour, alpha, size, _) in GROUPS.items():
        subset = points.loc[groups.eq(group)]
        ax.plot(subset.kout_w10 + 1, subset.kin_w10 + 1, linestyle="none", marker="o",
                ms=size, mew=.3 if group == "city" else 0, mec="white", mfc=colour, alpha=alpha,
                zorder=2 + list(GROUPS).index(group))
    ax.set(xscale="log", yscale="log", xlim=(.8, 3200), ylim=(.8, 3200),
           xlabel="Out-degree + 1")
    ax.set_ylabel("In-degree + 1" if y_label else "")
    ax.text(.97, .04, f"Spearman ρ = {rho:.2f}", ha="right", va="bottom", transform=ax.transAxes,
            fontsize=TEXT)
    return points.assign(panel=letter, spearman_rho=rho, zone_group=groups)


def metric(ax, data, title, letter, at, real, null, label, low=None, high=None):
    panel(ax, letter, title, at=at)
    for y, row in enumerate(data.itertuples()):
        if low:
            ax.plot([getattr(row, low), getattr(row, high)], [y, y], color=RANGE, lw=5.5,
                    solid_capstyle="butt", zorder=1)
        ax.plot(getattr(row, null), y, marker="|", ms=8, mew=1.1, color=MUTED, ls="none", zorder=2)
        ax.plot(getattr(row, real), y, "o", ms=MARK, color=INK, mec="white", mew=.5, zorder=3)
    ax.set_yticks(range(len(data)), list(DATES.values()))
    ax.set_ylim(len(data) - .45, -.55)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.grid(axis="y", visible=False)
    ax.grid(axis="x", visible=True)
    ax.xaxis.set_major_formatter(decimals(2))
    ax.set_xlabel(label)
    columns = ["date", real, null] + ([low, high] if low else [])
    return data[columns].assign(panel=letter)


def main():
    result = MAIN / "results"
    degree_data = pd.read_parquet(result / "p_node_metrics.parquet", columns=[
        "date", "idx", "kin_w10", "kout_w10",
    ])
    structure = pd.read_csv(result / "b02_structural_constraints.csv").query(
        "backbone == 'w10'"
    ).set_index("date").loc[list(DATES)].reset_index()
    fig = figure(124)
    emitted = []
    for i, (date, title) in enumerate(list(DATES.items())[:2]):
        ax = place(fig, 20 + i * 86, 10, 60, 54)
        emitted.append(degree(ax, degree_data.loc[degree_data.date.eq(date)], title, "ab"[i],
                              (LEFT, RIGHT)[i], i == 0))
        if i == 0:
            ax.legend(handles=[Line2D([], [], ls="", marker="o", ms=MARK, mfc=c, mec="white", mew=.3,
                                      label=l) for c, _, _, l in GROUPS.values()],
                      loc="upper left", bbox_to_anchor=(.03, .97))
    for i, (title, real, null, label, low, high) in enumerate([
        ("Degree assortativity", "r_real", "maxrand_conn_r_mean", "Assortativity r",
         "r_min_maxdisassortative", "r_max_maxassortative"),
        ("Second-order mixing", "rcpres_R_avg_real", "rcpres_R_avg_mean", "Mean second-order mixing",
         None, None),
    ]):
        ax = place(fig, 40 + i * 86, 85, 40, 17)
        emitted.append(metric(ax, structure, title, "cd"[i], (LEFT, RIGHT)[i], real, null, label,
                              low, high))
    fig.legend(handles=[
        Line2D([], [], ls="", marker="o", ms=MARK, color=INK, label="Observed"),
        Line2D([], [], ls="", marker="|", ms=8, mew=1.1, color=MUTED, label="Rewired mean (50 runs)"),
        Line2D([], [], color=RANGE, lw=5.5, solid_capstyle="butt", label="Range reachable by degree-preserving rewiring (c)"),
    ], loc="upper center", bbox_to_anchor=(.5, 1 - 113.5 / 124), ncol=3, handlelength=2)
    share = structure.maxrand_conn_r_mean.iloc[:2].diff().iloc[-1] / structure.r_real.iloc[:2].diff().iloc[-1]
    fig.text(.5, 1 - 119 / 124, f"Rewiring reproduces {share:.0%} of the February-to-August "
             "assortativity shift.", ha="center", va="top", fontsize=TEXT, color=MUTED)
    save(fig, "S1", *emitted, pd.DataFrame([{"panel": "c", "share_shift_reproduced": share}]))


if __name__ == "__main__":
    main()
