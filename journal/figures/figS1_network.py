"""Coursework network diagnostics at the ≥10-person link threshold."""

import pandas as pd

from style import MAIN, MUTED, NEUTRAL, ORDINARY, figure, panel, save

DATES = {"2024-02-14": "13–14 February", "2024-08-15": "14–15 August",
         "2024-11-01": "All Saints"}


def degree(ax, points, title, letter, y_label=False):
    rho = points.kin_w10.corr(points.kout_w10, method="spearman")
    panel(ax, letter, title, x=-.21)
    ax.plot([1, 3000], [1, 3000], ls="--", lw=.6, color=ORDINARY)
    ax.plot(points.kout_w10 + 1, points.kin_w10 + 1, linestyle="none", marker="o",
            ms=3.5, mew=0, mfc=NEUTRAL, alpha=.18)
    ax.set(xscale="log", yscale="log", xlim=(.8, 3200), ylim=(.8, 3200),
           xlabel="Out-degree + 1")
    if y_label:
        ax.set_ylabel("In-degree + 1")
    ax.text(.035, .94, f"Spearman ρ = {rho:.2f}", va="top", transform=ax.transAxes, fontsize=7.5,
            bbox={"fc": "white", "ec": "none", "pad": 1})
    return points.assign(panel=letter, spearman_rho=rho)


def metric(ax, data, title, letter, real, null):
    panel(ax, letter, title, x=-.21)
    for field, marker, line, fc in [(real, "o", "-", NEUTRAL), (null, "s", "--", "white")]:
        ax.plot(range(len(data)), data[field], marker=marker, linestyle=line,
                color=NEUTRAL, mfc=fc, mew=.7, ms=4)
    ax.set_xticks(range(len(data)), ["February", "August", "All Saints"])
    return data[["date", real, null]].assign(panel=letter)


def main():
    result = MAIN / "results"
    degree_data = pd.read_parquet(result / "p_node_metrics.parquet", columns=[
        "date", "idx", "kin_w10", "kout_w10",
    ])
    structure = pd.read_csv(result / "b02_structural_constraints.csv").query(
        "backbone == 'w10'"
    ).set_index("date").loc[list(DATES)].reset_index()
    fig = figure(105)
    emitted = []
    for i, (date, title) in enumerate(list(DATES.items())[:2]):
        emitted.append(degree(fig.add_axes([.111 + i * .479, .55, .37, .34]),
                              degree_data.loc[degree_data.date.eq(date)], title, "ab"[i], i == 0))
    for i, (title, real, null) in enumerate([
        ("Degree assortativity", "r_real", "maxrand_conn_r_mean"),
        ("Second-order mixing", "rcpres_R_avg_real", "rcpres_R_avg_mean"),
    ]):
        ax = fig.add_axes([.111 + i * .479, .147, .37, .25])
        emitted.append(metric(ax, structure, title, "cd"[i], real, null))
        ax.text(.34, .86, "●  Observed    □  Rewired",
                transform=ax.transAxes, va="bottom", color=MUTED, fontsize=7)
    share = structure.maxrand_conn_r_mean.iloc[:2].diff().iloc[-1] / structure.r_real.iloc[:2].diff().iloc[-1]
    fig.text(.111, .052, f"Rewiring reproduces {share:.0%} of the February-to-August assortativity shift.", fontsize=7)
    save(fig, "S1", *emitted, pd.DataFrame([{"panel": "c", "share_shift_reproduced": share}]))


if __name__ == "__main__":
    main()
