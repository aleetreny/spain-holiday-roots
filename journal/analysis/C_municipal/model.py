"""Differenced fixed effects and interval identification for monthly trip counts."""
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import splu
from scipy.stats import t


class FixedEffects:
    """Project on origin-city × destination-province and destination effects."""

    def __init__(self, frame, balanced):
        corridor = pd.factorize(frame.corridor)[0]
        village = pd.factorize(frame.dest)[0]
        province = pd.factorize(frame.dprov)[0]
        self.groups = [corridor, village, province]
        self.group_sizes = [np.bincount(group) for group in self.groups]
        self.rank = len(set(corridor)) + len(set(village)) - len(set(province))
        self.balanced = balanced
        if balanced:
            # Each destination province is a rectangular city × village block.
            counts = frame.groupby("dprov").agg(dyads=("orig", "size"), cities=("orig", "nunique"),
                                                villages=("dest", "nunique"))
            assert (counts.dyads == counts.cities * counts.villages).all()
            return

        indices = np.column_stack([corridor, village + self.group_sizes[0].size])
        design = sparse.csr_matrix((np.ones(indices.size),
                                   (np.repeat(np.arange(len(frame)), 2), indices.ravel())))
        _, components = connected_components(design.T @ design, directed=False)
        keep = np.ones(design.shape[1], dtype=bool)
        keep[np.unique(components, return_index=True)[1]] = False
        self.design = design[:, keep]
        self.rank = self.design.shape[1]
        self.solve = splu((self.design.T @ self.design).tocsc()).solve

    def residual(self, values):
        if not self.balanced:
            return values - self.design @ self.solve(self.design.T @ values)
        means = [np.bincount(group, weights=values) / sizes
                 for group, sizes in zip(self.groups, self.group_sizes)]
        return values - means[0][self.groups[0]] - means[1][self.groups[1]] + means[2][self.groups[2]]


class Design:
    """Frisch–Waugh–Lovell regression and two-way province-cluster inference."""

    def __init__(self, frame, balanced=True):
        self.n = len(frame)
        self.effects = FixedEffects(frame, balanced)
        d = frame.log_km.to_numpy()
        controls = np.column_stack([d, d ** 2, d * frame.coastal.to_numpy()])
        controls = np.column_stack([self.effects.residual(column) for column in controls.T])
        u, singular, _ = np.linalg.svd(controls, full_matrices=False)
        self.controls = u[:, singular > 1e-8]
        self.clusters = [pd.factorize(x)[0] for x in [frame.oprov, frame.dprov,
                                                      frame.oprov + frame.dprov]]
        self.g = [int(group.max()) + 1 for group in self.clusters]
        assert min(self.g) > 1
        self.critical = t.ppf(0.975, min(self.g[:2]) - 1)

    def residual(self, y):
        within = self.effects.residual(y)
        return within - self.controls @ (self.controls.T @ within)

    def fit(self, regressors):
        return Fit(self, regressors)


class Fit:
    def __init__(self, design, regressors):
        self.design = design
        self.x = np.column_stack([design.residual(column) for column in regressors.T])
        rank = self.x.shape[1]
        self.focus_identified = True
        if design.effects.balanced:
            self.inverse = np.linalg.inv(self.x.T @ self.x)
            self.weights = self.x @ self.inverse
        else:
            # In the selected panel some resamples have a redundant first-generation
            # control. Discard only numerically null directions and check estimability.
            u, s, v = np.linalg.svd(self.x, full_matrices=False)
            keep = s > max(1e-10, s.max() * 1e-10)
            rank = keep.sum()
            self.focus_identified = bool(np.all(np.abs(v[~keep, 0]) < 1e-8))
            self.weights = (u[:, keep] / s[keep]) @ v[keep]
        self.df_model = design.effects.rank + design.controls.shape[1] + rank

    def coefficients(self, y_residual):
        return self.weights.T @ y_residual

    def summary(self, y_residual, focus=0):
        beta = self.coefficients(y_residual)
        residual = y_residual - self.x @ beta
        weights = (self.weights[:, focus[0]] - self.weights[:, focus[1]]
                   if isinstance(focus, tuple) else self.weights[:, focus])
        scores = weights * residual
        components = []
        for group, g in zip(self.design.clusters, self.design.g):
            sums = np.bincount(group, weights=scores, minlength=g)
            components.append((sums @ sums) * g / (g - 1))
        variance = components[0] + components[1] - components[2]
        df_scale = (self.design.n - 1) / (self.design.n - self.df_model)
        return dict(estimate=weights @ y_residual, variance_cluster=max(0, variance) * df_scale,
                    variance_cluster_unclipped=variance * df_scale,
                    variance_origin=components[0] * df_scale,
                    variance_destination=components[1] * df_scale,
                    variance_intersection=components[2] * df_scale,
                    dyads=self.design.n, origin_provinces=self.design.g[0],
                    destination_provinces=self.design.g[1], critical=self.design.critical)

    def interval_endpoints(self, low, high):
        """Sharp box bounds on the OLS slope with unrestricted cell completions."""
        positive = self.weights[:, 0] >= 0
        return (np.where(positive, low, high), np.where(positive, high, low))


def finish(summary, refits):
    """Add survey bootstrap variance to conditional cluster variance."""
    estimate = summary["estimate"]
    survey = np.asarray(refits, dtype=float)
    valid = survey[np.isfinite(survey)]
    var_cluster = summary["variance_cluster"]
    conditional_radius = summary["critical"] * np.sqrt(var_cluster)
    if len(valid) < 0.95 * len(survey):
        return dict(**summary, variance_survey=np.nan, se=np.nan, lo95=np.nan, hi95=np.nan,
                    survey_draws=len(survey), identified_draws=len(valid),
                    survey_bootstrap_mean=np.nan, lo95_conditional=estimate - conditional_radius,
                    hi95_conditional=estimate + conditional_radius)
    survey = valid
    var_survey = survey.var(ddof=1)
    se = np.sqrt(var_cluster + var_survey)
    radius = summary["critical"] * se
    return dict(**summary, variance_survey=var_survey, se=se,
                lo95=estimate - radius, hi95=estimate + radius,
                survey_draws=len(refits), identified_draws=len(survey),
                survey_bootstrap_mean=survey.mean(), lo95_conditional=estimate - conditional_radius,
                hi95_conditional=estimate + conditional_radius)
