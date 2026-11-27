"""Compare zero-cell rules in predicted-flow units rather than transformed slopes."""

from matplotlib.lines import Line2D

from style import (
    B, EXODUS, MARK, NEUTRAL, RECENT, decimals, figure, forest, interval, panel, place, read, save,
)

RULES = {
    "haldane_0.5_main": "Add 0.5 to count",
    "plus_1.0_log1p": "Add 1 to count",
    "plus_0.1_pop_plus1": "Add 0.1 to count",
    "log_share_plus_1e4": "Add 0.0001 to share",
    "log_share_plus_1e5": "Add 0.00001 to share",
    "asinh_1000_share": "asinh(1,000 × share)",
    "R2_112_large_cities_haldane0.5": "112 largest cities, +0.5",
}
TERMS = {
    "pre1981": ("Arrived before 1981", EXODUS),
    "recent_2016_2023": ("Arrived 2016–2023", RECENT),
    "born_village_le2k": ("Birthplace ≤ 2,000", EXODUS),
    "born_urban_gt10k": ("Birthplace > 10,000", RECENT),
}
BASELINES = {"Aug15_vs_Feb14": ("o", -.15, "13–14 February"),
             "Aug15_vs_OctMean12": ("s", .15, "October weekdays")}


def main():
    data = read(B, "b5_zero_cells_common_scale").query("sample == 'all' and term in @TERMS and not model.str.contains('complete')").copy()
    fig = figure(116)
    estimate = "pct_flow_1sd_rawshare_fw"
    extrema = min(data[f"{estimate}_ci_lo"]), max(data[f"{estimate}_ci_hi"])
    for i, (term, (title, colour)) in enumerate(TERMS.items()):
        x, y = i % 2, i // 2
        ax = place(fig, 42 + x * 84, 10 + y * 47, 40, 34)
        panel(ax, "abcd"[i], title, at=2 + x * 84)
        forest(ax, list(RULES.values()), (extrema[0] - 1, extrema[1] + 1), [-10, 0, 10, 20])
        ax.xaxis.set_major_formatter(decimals(0))
        if y == 1:
            ax.set_xlabel("August change in predicted flow (%)\nper one-SD rise in raw birth share")
        rows = data.loc[data.term.eq(term)].set_index(["rule", "contrast"])
        for j, rule in enumerate(RULES):
            for contrast, (marker, shift, _) in BASELINES.items():
                row = rows.loc[(rule, contrast)]
                interval(ax, j + shift, row[estimate], row[f"{estimate}_ci_lo"], row[f"{estimate}_ci_hi"],
                         color=colour, marker=marker)
    fig.legend(handles=[Line2D([], [], ls="", marker=mark, ms=MARK, color=NEUTRAL, label=label)
                        for mark, _, label in BASELINES.values()],
               loc="upper center", bbox_to_anchor=(.5, 1 - 110 / 116), ncol=2)
    save(fig, "S5", data.assign(panel=data.term.map(dict(zip(TERMS, "abcd")))))


if __name__ == "__main__":
    main()
