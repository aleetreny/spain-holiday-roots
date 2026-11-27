"""Within-dyad roots and monthly travel, on extensive and intensive margins."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import splu
from scipy.stats import t

from pilot_data import OUT, ine_counts, municipal_dyads

REGRESSORS = ["parents", "first", "log_km", "log_km_sq"]


class TwoWayEffects:
    def __init__(self, groups):
        n = len(groups)
        columns = groups.copy()
        columns[:, 1] += columns[:, 0].max() + 1
        full = sparse.csr_matrix(
            (np.ones(2 * n), (np.repeat(np.arange(n), 2), columns.ravel()))
        )
        normal = full.T @ full
        _, labels = connected_components(normal, directed=False)
        # Omit one intercept per connected component to make the system full rank.
        keep = np.ones(full.shape[1], dtype=bool)
        keep[np.unique(labels, return_index=True)[1]] = False
        self.design = full[:, keep].tocsr()
        self.degrees = self.design.shape[1]
        self.solve = splu((self.design.T @ self.design).tocsc())

    def residualize(self, values):
        coefficients = self.solve.solve(self.design.T @ values)
        return values - self.design @ coefficients


class WithinFit:
    def __init__(self, base, regressors, effects, cluster, cluster_universe):
        groups = np.column_stack([pd.factorize(base[name])[0] for name in effects])
        self.absorb = TwoWayEffects(groups)
        self.x = self.absorb.residualize(base[regressors].to_numpy(dtype=float))
        self.cluster = pd.Categorical(base[cluster], categories=cluster_universe).codes
        self.n_clusters = len(np.unique(self.cluster))
        self.cluster_universe = cluster_universe
        self.inverse = np.linalg.inv(self.x.T @ self.x)

    def estimate(self, outcome):
        y = self.absorb.residualize(outcome.astype(float).reshape(-1, 1)).ravel()
        beta = self.inverse @ (self.x.T @ y)
        scores = (self.x @ self.inverse[:, 0]) * (y - self.x @ beta)
        influences = np.bincount(self.cluster, weights=scores, minlength=len(self.cluster_universe))
        clusters = self.n_clusters
        n = len(y)
        df_model = self.absorb.degrees + len(beta)
        se = np.sqrt((clusters / (clusters - 1)) * ((n - 1) / (n - df_model)) * (influences @ influences))
        radius = t.ppf(0.975, clusters - 1) * se
        return {"estimate": beta[0], "se": se, "lo95": beta[0] - radius, "hi95": beta[0] + radius,
                "n": n, "destination_provinces": clusters}


def main():
    base = municipal_dyads()
    flows = ine_counts(base)
    clusters = sorted(base["dprov"].unique())
    effects = ["city_province", "dest"]
    # Use a static origin × village frame, independent of the flow outcome.
    full_fit = WithinFit(base, REGRESSORS, effects, "dprov", clusters)
    sample_stats = [{
        "cities": base["orig"].nunique(), "villages": base["dest"].nunique(), "dyads": len(base),
        "parents_positive": int(base["parents"].sum()), "first_positive": int(base["first"].sum()),
        "parents_two_households": int(base["g2_households"].ge(2).sum()),
        "months": len(flows),
    }]
    rows = []
    for month, outcome in flows.items():
        year, month_number = month.split("-")
        references = ["2019-10"] if year == "2019" else [f"{year}-02"]
        if month_number == "08" and year != "2019":
            references.append(f"{year}-10")
        for reference in references:
            if reference not in flows or reference == month:
                continue
            baseline = flows[reference]
            published = np.isfinite(outcome)
            published_baseline = np.isfinite(baseline)
            both = published & published_baseline
            detection = published.astype(float) - published_baseline
            intensive = np.log(outcome[both]) - np.log(baseline[both])
            rows.append({"month": month, "reference": reference, "outcome": "published",
                         **full_fit.estimate(detection)})
            rows.append({"month": month, "reference": reference, "outcome": "log_count_published_both",
                         **WithinFit(base.loc[both], REGRESSORS, effects, "dprov", clusters).estimate(intensive)})
    # Use an equal-weight mean of annual differences and keep the same dyads across years.
    for reference_number in ["02", "10"]:
        years = [str(y) for y in range(2022, 2027) if f"{y}-{reference_number}" in flows]
        target = np.column_stack([flows[f"{y}-08"] for y in years])
        reference = np.column_stack([flows[f"{y}-{reference_number}"] for y in years])
        published_diff = np.isfinite(target).astype(float) - np.isfinite(reference)
        intensive_diff = np.log(target) - np.log(reference)
        fixed = np.isfinite(intensive_diff).all(axis=1)
        fixed_fit = WithinFit(base.loc[fixed], REGRESSORS, effects, "dprov", clusters)
        for outcome, fit, change in [
            ("published", full_fit, published_diff),
            ("log_count_common_pairs", fixed_fit, intensive_diff[fixed]),
        ]:
            year_centred = np.arange(len(years)) - (len(years) - 1) / 2
            for kind, y in [
                ("mean", change.mean(axis=1)),
                ("linear_trend_per_year", change @ year_centred / (year_centred @ year_centred)),
            ]:
                rows.append({"month": f"{years[0]}-{years[-1]} Aug {kind}",
                             "reference": reference_number, "outcome": outcome,
                             **fit.estimate(y)})
    # Check whether the baseline-dependent results survive larger thresholds and better
    # sampled ties. Thresholds are applied to observable counts only.
    robustness = []
    target = np.column_stack([flows[f"{y}-08"] for y in range(2022, 2027)])
    reference = np.column_stack([flows[f"{y}-02"] for y in range(2022, 2027)])
    for cutoff in [30, 60, 100]:
        appears = (target >= cutoff).astype(float) - (reference >= cutoff)
        robustness.append({"specification": f"publication_at_{cutoff}_or_more",
                           **full_fit.estimate(appears.mean(axis=1))})
        paired = ((target >= cutoff) & (reference >= cutoff)).all(axis=1)
        fit = WithinFit(base.loc[paired], REGRESSORS, effects, "dprov", clusters)
        contrast = (np.log(target[paired]) - np.log(reference[paired])).mean(axis=1)
        robustness.append({"specification": f"common_pairs_at_{cutoff}_or_more",
                           **fit.estimate(contrast)})

    base["two_parent_households"] = base["g2_households"].ge(2).astype(float)
    fit = WithinFit(base, ["two_parent_households", *REGRESSORS], effects, "dprov", clusters)
    detection = np.isfinite(target).astype(float) - np.isfinite(reference)
    robustness.append({"specification": "additional_2plus_households_vs_one",
                       **fit.estimate(detection.mean(axis=1))})

    OUT.mkdir(exist_ok=True)
    pd.DataFrame(robustness).to_csv(OUT / "municipal_robustness.csv", index=False, float_format="%.8g")
    pd.DataFrame(sample_stats).to_csv(OUT / "municipal_sample.csv", index=False)
    pd.DataFrame(rows).to_csv(OUT / "municipal_calendar.csv", index=False, float_format="%.8g")


if __name__ == "__main__":
    main()
