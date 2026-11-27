"""Run Analysis C with the existing isolated project interpreter."""
import argparse
from dataclasses import dataclass
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd

from falsification import matched_placebos
from inputs import MONTHS, OUT, PRIVATE, Survey, THRESHOLD, fixed_frame, n1, panel
from margins import audit
from model import Design, finish


@dataclass
class Outcome:
    sample: str
    reference: str
    years: str
    treatment: str
    residual: np.ndarray
    low: np.ndarray | None = None
    high: np.ndarray | None = None


def regressors(g1, g2, h1, h2):
    # One basis for presence and a separate one for dose.
    control = [h1 > 0, np.log1p(10_000 * g1)]
    return {
        "primary": np.column_stack([h2 > 0, *control]),
        "share": np.column_stack([h2 > 0, np.log1p(10_000 * g2), *control]),
        "households": np.column_stack([h2 > 0, h2 >= 2, *control]),
    }


def annual_pairs(reference):
    return [(m, m[:4] + "-" + reference) for m in MONTHS
            if m.endswith("-08") and m[:4] + "-" + reference in MONTHS]


def year_sets(reference):
    pairs = annual_pairs(reference)
    yield "all_years", pairs
    yield "post_2021", [p for p in pairs if p[0][:4] >= "2022"]
    yield "2020_2021", [p for p in pairs if p[0][:4] in {"2020", "2021"}]


def create_outcomes(design, flows, selected, selected_design):
    outcomes = []
    selected_outcomes = []
    high_value = THRESHOLD - 1
    for reference in ["02", "10"]:
        for sample, pairs in year_sets(reference):
            labels = " ".join(a[:4] for a, _ in pairs)
            raw = {name: np.zeros(design.n) for name in ["publication", "fill_0", "fill_29",
                                                         "interval_low", "interval_high"]}
            logged = np.zeros(selected.sum())
            for a, b in pairs:
                av, bv = flows[a], flows[b]
                ap, bp = av >= 0, bv >= 0
                a0, b0 = np.maximum(av, 0), np.maximum(bv, 0)
                a29, b29 = np.where(ap, av, high_value), np.where(bp, bv, high_value)
                raw["publication"] += ap.astype(float) - bp
                raw["fill_0"] += a0 - b0
                raw["fill_29"] += a29 - b29
                raw["interval_low"] += a0 - b29
                raw["interval_high"] += a29 - b0
                logged += np.log(av[selected]) - np.log(bv[selected])
            raw = {k: v / len(pairs) for k, v in raw.items()}
            for t in ["publication", "fill_0", "fill_29"]:
                outcomes.append(Outcome(sample, reference, labels, t, design.residual(raw[t])))
            for t in ["interval_lower", "interval_upper"]:
                outcomes.append(Outcome(sample, reference, labels, t, np.empty(0),
                                        raw["interval_low"], raw["interval_high"]))
            selected_outcomes.append(Outcome(sample, reference, labels, "selected_log_count",
                                             selected_design.residual(logged / len(pairs))))
    return outcomes, selected_outcomes


def fit_outcomes(fit, outcomes, full=True):
    if not fit.focus_identified:
        if full:
            raise ValueError("Parental tie coefficient is unidentified on the observed frame")
        return [np.nan] * len(outcomes)
    summaries = []
    for o in outcomes:
        y = o.residual
        if o.low is not None:
            low, high = fit.interval_endpoints(o.low, o.high)
            y = low if o.treatment == "interval_lower" else high
            if not full:
                summaries.append(fit.weights[:, 0] @ y)
                continue
            y = fit.design.residual(y)
        summaries.append(fit.summary(y) if full else fit.coefficients(y)[0])
    return summaries


def calendar_outcomes(design, flows):
    matrix = np.empty((design.n, 3 * len(MONTHS)), order="F")
    for i, month in enumerate(MONTHS):
        v = flows[month]
        matrix[:, i] = design.residual((v >= 0).astype(float))
        matrix[:, i + len(MONTHS)] = design.residual(np.maximum(v, 0))
        matrix[:, i + 2 * len(MONTHS)] = design.residual(np.where(v >= 0, v, THRESHOLD - 1))
    return matrix


def export_calendar(fit, calendar, draws):
    rows = []
    for k, treatment in enumerate(["publication", "fill_0", "fill_29"]):
        column = {m: k * len(MONTHS) + i for i, m in enumerate(MONTHS)}
        for m in MONTHS:
            for reference in ["02", "10"]:
                b = m[:4] + "-" + reference
                if b not in column or m == b:
                    continue
                j, l = column[m], column[b]
                s = fit.summary(calendar[:, j] - calendar[:, l])
                rows.append(dict(month=m, baseline=b, treatment=treatment,
                                 pandemic_year=m[:4] in {"2020", "2021"},
                                 **finish(s, draws[:, j] - draws[:, l])))
    pd.DataFrame(rows).to_csv(OUT / "calendar.csv", index=False)


def export_estimates(file, outcomes, summaries, boot):
    rows = [dict(sample=o.sample, reference=o.reference, years=o.years,
                 treatment=o.treatment, **finish(s, b))
            for o, s, b in zip(outcomes, summaries, boot.T)]
    pd.DataFrame(rows).to_csv(OUT / file, index=False)


def model_check(fit, selected_fit, design, selected_design, frame, x):
    # Compare the fast balanced projection against an independent explicit FE solve.
    rebuilt = Design(frame, balanced=False)
    probe = np.random.default_rng(100).normal(size=len(frame))
    np.testing.assert_allclose(design.residual(probe), rebuilt.residual(probe), atol=1e-9)
    assert np.max(np.abs(design.controls.T @ fit.x)) < 1e-9
    assert np.max(np.abs(selected_design.controls.T @ selected_fit.x)) < 1e-9
    lo, hi = -np.abs(probe), np.abs(probe)
    lower, upper = fit.interval_endpoints(lo, hi)
    assert fit.weights[:, 0] @ lower <= 0 <= fit.weights[:, 0] @ upper
    assert np.max(np.abs(fit.weights.T @ x - np.eye(x.shape[1]))) < 1e-9


def main(args):
    OUT.mkdir(exist_ok=True)
    PRIVATE.mkdir(parents=True, exist_ok=True)
    frame, codes = fixed_frame()
    survey = Survey(frame, codes)
    observed = survey.exposures(np.ones(len(survey.total), dtype=int))
    g1, g2, h1, h2 = observed
    flows = panel(frame)
    if not args.skip_margin_audit:
        audit(frame, codes, flows, h2)
    selected = np.logical_and.reduce([v >= 0 for v in flows.values()])
    pd.DataFrame([dict(
        cities=frame.orig.nunique(), villages=frame.dest.nunique(), dyads=len(frame),
        origin_provinces=frame.oprov.nunique(), destination_provinces=frame.dprov.nunique(),
        months=len(flows), households=len(survey.total),
        city_born_households=int((survey.local > 0).sum()),
        parent_ties=int((h2 > 0).sum()), parent_ties_2plus_households=int((h2 >= 2).sum()),
        first_gen_ties=int((h1 > 0).sum()), first_gen_ties_2plus_households=int((h1 >= 2).sum()),
        coastal_villages=frame.loc[frame.coastal, "dest"].nunique(),
        coastal_cities=int(pd.read_parquet(PRIVATE / "municipal_geometry.parquet")
                           .loc[frame.orig.unique(), "coastal"].sum()),
        missing_housing_villages=frame.loc[frame.housing_share.isna(), "dest"].nunique(),
        distances_present=int(frame.km.notna().sum()), min_km=frame.km.min(), max_km=frame.km.max(),
        always_published=int(selected.sum()), always_published_tied=int((selected & (h2 > 0)).sum()),
    )]).to_csv(OUT / "sample.csv", index=False)
    private_exp = frame[["orig", "dest"]].assign(g1=g1, g2=g2, g1_households=h1, g2_households=h2)
    private_exp.to_parquet(PRIVATE / "exposures.parquet", index=False)

    print("Fit main and selected designs", flush=True)
    design = Design(frame)
    selected_design = Design(frame.loc[selected], balanced=False)
    x = regressors(*observed)
    base = design.fit(x["primary"])
    selected_fit = selected_design.fit(x["primary"][selected])
    if args.check:
        model_check(base, selected_fit, design, selected_design, frame, x["primary"])
    outcomes, selected_outcomes = create_outcomes(design, flows, selected, selected_design)
    summaries = fit_outcomes(base, outcomes) + fit_outcomes(selected_fit, selected_outcomes)

    print("Prepare calendar and dose diagnostics", flush=True)
    calendar = calendar_outcomes(design, flows)
    support = np.union1d(survey.g1_rows, survey.g2_rows)
    sparse_calendar = calendar[support]
    dose_outcomes = [o for o in outcomes if o.sample == "all_years"
                     and o.treatment in {"publication", "fill_0", "fill_29"}]
    dose_models = ["share", "households"]
    dose_fits = [design.fit(x[name]) for name in dose_models]
    dose_summaries = [[fit.summary(o.residual, focus=1) for o in dose_outcomes] for fit in dose_fits]
    original, placebo, balance = matched_placebos(frame, h1, h2)
    balance.to_csv(OUT / "match_balance.csv", index=False)
    # Include original tie and the matched placebo together. Also separate the matched true ties
    # from those outside common support, to compare like with like.
    matched_tie = np.zeros(len(frame))
    matched_tie[original] = 1
    placebo_indicator = np.zeros(len(frame))
    placebo_indicator[placebo] = 1

    def matched_design(primary):
        retained_match = matched_tie * primary[:, 0]
        unmatched = primary[:, 0] - retained_match
        return design.fit(np.column_stack([retained_match, placebo_indicator, unmatched, primary[:, 1:]]))

    match = matched_design(x["primary"])
    match_summaries = [[match.summary(o.residual, focus=j) for o in dose_outcomes]
                       for j in [0, 1, (0, 1)]]
    rng = np.random.default_rng(n1.C.SEED)
    boot = np.empty((args.draws, len(summaries)))
    cal_boot = np.empty((args.draws, calendar.shape[1]))
    dose_boot = np.empty((len(dose_models), args.draws, len(dose_outcomes)))
    match_boot = np.empty((len(match_summaries), args.draws, len(dose_outcomes)))
    dose_y = np.column_stack([o.residual for o in dose_outcomes])
    for draw in range(args.draws):
        bx = regressors(*survey.exposures(survey.draw(rng)))
        bfit = design.fit(bx["primary"])
        sfit = selected_design.fit(bx["primary"][selected])
        boot[draw] = fit_outcomes(bfit, outcomes, full=False) + fit_outcomes(
            sfit, selected_outcomes, full=False)
        # Residualised outcomes are orthogonal to nuisance controls. Their cross-product
        # with the exposure can be evaluated only where any household reported a tie.
        cal_boot[draw] = (bfit.inverse[0] @ bx["primary"][support].T) @ sparse_calendar
        for k, name in enumerate(dose_models):
            dose_boot[k, draw] = design.fit(bx[name]).coefficients(dose_y)[1]
        match_b = matched_design(bx["primary"])
        paired = match_b.coefficients(dose_y)[:2]
        match_boot[:2, draw, :] = paired
        match_boot[2, draw, :] = paired[0] - paired[1]
        if (draw + 1) % 50 == 0:
            print(f"Refitted {draw + 1}/{args.draws} household resamples", flush=True)

    export_estimates("estimates.csv", outcomes + selected_outcomes, summaries, boot)
    export_calendar(base, calendar, cal_boot)
    for file, terms, estimates, draws in [
        ("dose.csv", dose_models, dose_summaries, dose_boot),
        ("falsification.csv", ["matched_true", "matched_placebo", "true_minus_placebo"],
         match_summaries, match_boot),
    ]:
        rows = []
        for term, term_estimates, term_draws in zip(terms, estimates, draws):
            for o, summary, bs in zip(dose_outcomes, term_estimates, term_draws.T):
                rows.append(dict(term=term, reference=o.reference, years=o.years,
                                 treatment=o.treatment, **finish(summary, bs)))
        pd.DataFrame(rows).to_csv(OUT / file, index=False)
    print("Analysis C complete", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--draws", type=int, default=399)
    parser.add_argument("--skip-margin-audit", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.draws < 2:
        parser.error("--draws must be at least two")
    main(args)
