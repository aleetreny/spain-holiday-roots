"""Compare August prediction with and without survey birthplace features."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd
from scipy.special import xlogy
from sklearn.ensemble import HistGradientBoostingClassifier, HistGradientBoostingRegressor
from sklearn.metrics import average_precision_score

from pilot_data import OUT, ine_counts, municipal_dyads

SEED = 20261127
GEOGRAPHY = ["log_km", "dest_pop", "orig_pop", "income", "coast_km", "historical_loss"]
HISTORY = ["february", "previous_august", "previous_october"]
ROOTS = ["log_g1", "log_g2", "log_parent_households"]


def features(base, flows, year):
    table = base.copy()
    table["february"] = np.log1p(flows[f"{year}-02"])
    table["previous_august"] = np.log1p(flows[f"{year-1}-08"])
    table["previous_october"] = np.log1p(flows[f"{year-1}-10"])
    table["count"] = flows[f"{year}-08"]
    return table


def loss_difference(base_loss, roots_loss, groups):
    group = pd.factorize(groups)[0]
    n = np.bincount(group)
    total_base = np.bincount(group, weights=base_loss)
    total_roots = np.bincount(group, weights=roots_loss)
    rng = np.random.default_rng(SEED)
    draws = rng.integers(len(n), size=(1999, len(n)))
    improvements = 100 * (1 - total_roots[draws].sum(axis=1) / total_base[draws].sum(axis=1))
    return {
        "base_loss": total_base.sum() / n.sum(), "roots_loss": total_roots.sum() / n.sum(),
        "reduction_pct": 100 * (1 - total_roots.sum() / total_base.sum()),
        "lo95_pct": np.quantile(improvements, 0.025), "hi95_pct": np.quantile(improvements, 0.975),
        "dyads": len(group), "destination_provinces": len(n),
    }


def fit(train, tests, cols):
    config = dict(max_iter=120, max_leaf_nodes=15, min_samples_leaf=100, learning_rate=0.1,
                  early_stopping=False, random_state=SEED)
    published = train["count"].notna()
    clf = HistGradientBoostingClassifier(**config).fit(train[cols], published)
    reg = HistGradientBoostingRegressor(loss="poisson", **config).fit(
        train.loc[published, cols], train.loc[published, "count"],
    )
    losses = {}
    for year, table in tests.items():
        observable = table["count"].notna()
        p = np.clip(clf.predict_proba(table[cols])[:, 1], 1e-9, 1 - 1e-9)
        loss_pub = np.where(observable, -np.log(p), -np.log1p(-p))
        y = table.loc[observable, "count"].to_numpy()
        mu = reg.predict(table.loc[observable, cols])
        loss_count = 2 * (xlogy(y, y / mu) - y + mu)
        losses[year] = (loss_pub, loss_count, average_precision_score(observable, p))
    return losses


def main():
    base = municipal_dyads()
    flows = ine_counts(base)
    base["log_g1"] = np.log1p(base["g1"])
    base["log_g2"] = np.log1p(base["g2"])
    base["log_parent_households"] = np.log1p(base["g2_households"])
    table_cols = GEOGRAPHY + ROOTS + ["dprov"]
    train = pd.concat([features(base[table_cols], flows, y) for y in [2022, 2023, 2024]], ignore_index=True)
    tests = {y: features(base[table_cols], flows, y) for y in [2025, 2026]}
    rows = []
    for label, cols in [("geography", GEOGRAPHY), ("geography_and_history", GEOGRAPHY + HISTORY)]:
        base_fit = fit(train, tests, cols)
        roots_fit = fit(train, tests, cols + ROOTS)
        for year, table in tests.items():
            base_pub, base_count, base_ap = base_fit[year]
            roots_pub, roots_count, roots_ap = roots_fit[year]
            pub = table["count"].notna()
            for metric, loss0, loss1, clusters in [
                ("publication_log_loss", base_pub, roots_pub, table["dprov"]),
                ("published_count_poisson_deviance", base_count, roots_count, table.loc[pub, "dprov"]),
            ]:
                rows.append({"year": year, "baseline": label, "metric": metric,
                             **loss_difference(loss0, loss1, clusters),
                             "base_ap": base_ap if metric == "publication_log_loss" else np.nan,
                             "roots_ap": roots_ap if metric == "publication_log_loss" else np.nan})
    pd.DataFrame(rows).to_csv(OUT / "august_prediction.csv", index=False, float_format="%.8g")


if __name__ == "__main__":
    main()
