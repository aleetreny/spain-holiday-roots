"""Outcome-blind within-corridor matched placebo exposure."""
import numpy as np
import pandas as pd


def matched_placebos(frame, h1, h2):
    features = np.column_stack([np.log1p(frame["pop"]), frame["log_km"], frame["housing_share"]])
    tolerances = np.array([0.25, 0.10, 0.10])
    scaled = features / tolerances
    usable = np.isfinite(scaled).all(axis=1)
    groups = frame.assign(first=h1 > 0).groupby(["corridor", "coastal", "first"]).indices
    paired = []
    for dyads in groups.values():
        treated = dyads[(h2[dyads] > 0) & usable[dyads]]
        controls = dyads[(h2[dyads] == 0) & usable[dyads]]
        used = np.zeros(len(controls), dtype=bool)
        possibilities = []
        for t in treated:
            diff = scaled[controls] - scaled[t]
            candidates = np.flatnonzero((np.abs(diff) <= 1).all(axis=1))
            costs = (diff[candidates] ** 2).sum(axis=1)
            possibilities.append((len(candidates), t, candidates[np.argsort(costs, kind="stable")]))
        for _, t, candidates in sorted(possibilities, key=lambda row: (row[0], row[1])):
            available = candidates[~used[candidates]]
            if len(available):
                c = available[0]
                used[c] = True
                paired.append((t, controls[c]))
    matched = np.array(paired)
    original, placebos = matched.T
    diff = abs(features[original] - features[placebos])
    balance = pd.DataFrame({
        "covariate": ["log1p_population", "log_distance", "nonprimary_housing_share"],
        "caliper": tolerances, "mean_abs_difference": diff.mean(axis=0),
        "p95_abs_difference": np.quantile(diff, 0.95, axis=0),
        "matched": len(original), "eligible_ties": int((h2 > 0).sum()),
    })
    return original, placebos, balance
