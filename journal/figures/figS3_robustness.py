"""Sensitivity of the pooled August and pre-1981 contrasts."""

import numpy as np
import pandas as pd

from style import (
    A, INK, MARK, ORDINARY, SMALL, decimals, figure, forest, interval, panel, place, read, save,
)

LABELS = {
    "main": "Main specification",
    "1_year_matched_income_cohort": "Year-matched income",
    "1_year_matched_roots_income": "Year-matched census and income",
    "2_road_travel_time": "Road travel time",
    "3_drop_BCN_only": "Exclude Barcelona",
    "3_drop_MAD_and_BCN": "Exclude Madrid and Barcelona",
    "3_drop_MAD_only": "Exclude Madrid",
    "5_coastal_zones_19": "Coastal zones",
    "5_inland_provinces_only_28": "Inland provinces",
    "5_inland_zones_all47": "Inland zones",
    "6_zero_112_large_cities": "112 largest cities",
    "6_zero_asinh_1000_share": "asinh(1,000 × share)",
    "6_zero_plus_0.1": "Add 0.1 to count",
    "6_zero_plus_1.0": "Add 1 to count",
    "6_zero_share_plus_1e4": "Add 0.0001 to share",
    "6_zero_share_plus_1e5": "Add 0.00001 to share",
    "7_seasonal_ref_Jun_Sep_only": "June/September reference",
    "8_equal_weight_by_night": "Equal weight per night",
    "8_single_pooled_delta_across_yrs": "Single pooled contrast",
    "9_drop_3_artefact_dest_zones": "Exclude three anomalous zones",
}


def main():
    main_tests = read(A, "preregistered_tests").query(
        "split == 'confirmation_O' and year_scope == 'equal_year_avg' and hypothesis in ['P1','P2']"
    ).assign(check="main")
    checks = read(A, "robustness").query("hypothesis in ['P1','P2']")
    data = pd.concat([main_tests, checks], ignore_index=True)
    lodo = read(A, "robustness_lodo_dest", dtype={"dropped_dest_prov": str})
    fig = figure(117)
    limits, ticks = (-.155, .3), [-.1, 0, .1, .2, .3]
    pitch, top = 4.6, 10
    counts = {hyp: sum(key in set(data.loc[data.hypothesis.eq(hyp), "check"]) for key in LABELS)
              for hyp in ("P1", "P2")}
    strip_top = top + max(counts.values()) * pitch + 9
    for i, (hyp, title) in enumerate([("P1", "P1 · Pooled roots"), ("P2", "P2 · Arrived before 1981")]):
        left = 44 + i * 84
        rows = data.loc[data.hypothesis.eq(hyp)].set_index("check")
        keys = [key for key in LABELS if key in rows.index]
        ax = place(fig, left, top, 38, len(keys) * pitch)
        panel(ax, "ab"[i], title, at=2 + i * 84)
        forest(ax, [LABELS[key] for key in keys], limits, ticks)
        for y, row in enumerate(rows.loc[keys].itertuples()):
            interval(ax, y, row.est, row.ci_wcr_lo, row.ci_wcr_hi, marker="D" if y == 0 else "o",
                     color=INK)
        ax.xaxis.set_major_formatter(decimals(1))
        strip = place(fig, left, strip_top, 38, 8)
        forest(strip, ["Leave one destination\nprovince out (47 fits)"], limits, ticks)
        rng = np.random.default_rng(20261124)
        lodo[f"{hyp}_strip_y"] = rng.uniform(-.27, .27, len(lodo))
        strip.plot(lodo[f"{hyp}_est"], lodo[f"{hyp}_strip_y"], linestyle="none", marker="o",
                   color=ORDINARY, ms=SMALL, alpha=.7, mew=0)
        strip.plot(rows.loc["main", "est"], 0, "o", color=INK, ms=MARK)
        strip.set_ylim(.5, -.5)
        strip.xaxis.set_major_formatter(decimals(1))
        strip.set_xlabel("Change in roots elasticity\nAugust vs ordinary nights")
    save(fig, "S3", data.assign(panel=data.hypothesis.map({"P1": "a", "P2": "b"})),
         lodo.assign(panel="a,b", series="destination leave-one-out"))


if __name__ == "__main__":
    main()
