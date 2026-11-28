"""Compare the zone-excluded rerun against the original Analysis A outputs."""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_noart import DEFAULT_A, OUT

KEYS = ["split", "year_scope", "hypothesis"]
SAMPLE = [
    "family", "years_included", "holiday_class", "n_holiday_nights",
    "n_ordinary_nights", "model", "term", "G", "n_blocks",
]
METRICS = [
    "est", "ci_wcr_lo", "ci_wcr_hi", "ci_crv_lo", "ci_crv_hi",
    "p_holm", "p_holm_hybrid_S",
]


def paired_rows(before, after, keys, sample):
    before = before.set_index(keys)
    after = after.set_index(keys)
    if not before.index.is_unique or not after.index.is_unique:
        raise ValueError(f"Published or rerun rows are duplicated for {keys}")
    if set(before.index) != set(after.index):
        raise ValueError(f"Published and rerun rows differ for {keys}")
    after = after.reindex(before.index)
    if not before[sample].equals(after[sample]):
        raise ValueError(f"Published and rerun samples differ for {keys}")
    return before, after


def compare_tests(before, after):
    import numpy as np

    before, after = paired_rows(before, after, KEYS, SAMPLE)
    result = before[SAMPLE].copy()
    for metric in METRICS:
        result[f"{metric}_a"] = before[metric]
        result[f"{metric}_noart"] = after[metric]
        result[f"abs_diff_{metric}"] = (after[metric] - before[metric]).abs()

    # A's tests are two-sided. A Holm rejection supports the stated hypothesis
    # only if the estimate also has the pre-registered sign (negative for S5).
    direction = np.where(before.index.get_level_values("hypothesis") == "S5", -1, 1)
    for p_col, label in (("p_holm", "decision"), ("p_holm_hybrid_S", "hybrid_s_decision")):
        for suffix, frame in (("a", before), ("noart", after)):
            result[f"{label}_{suffix}"] = np.select(
                [
                    frame[p_col].isna(),
                    frame[p_col].le(0.05) & ((frame["est"] * direction) > 0),
                    frame[p_col].le(0.05),
                ],
                ["", "supported", "opposite_sign"],
                default="not_detected",
            )
        result[f"{label}_changed"] = result[f"{label}_a"].ne(result[f"{label}_noart"])

    return result.reset_index()


def main():
    import pandas as pd

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-a", type=Path, default=DEFAULT_A)
    published = parser.parse_args().analysis_a / "out"

    original = pd.read_csv(published / "preregistered_tests.csv")
    rerun = pd.read_csv(OUT / "preregistered_tests_noart.csv")
    comparison = compare_tests(original, rerun)
    comparison.to_csv(OUT / "comparison_noart.csv", index=False)

    before, after = paired_rows(
        pd.read_csv(published / "placebo_summary.csv"),
        pd.read_csv(OUT / "placebo_summary_noart.csv"),
        ["split", "year", "hypothesis"], ["n_placebo_windows"],
    )
    print("August ranks changed:", before["aug_rank_of_total"].ne(after["aug_rank_of_total"]).sum())
    windows_a, windows_noart = paired_rows(
        pd.read_csv(published / "placebo_windows.csv"),
        pd.read_csv(OUT / "placebo_windows_noart.csv"),
        ["split", "year", "window_type", "window_idx"],
        ["start_date", "end_date", "n_win", "n_ref"],
    )
    print("Placebo and August windows matched:", len(windows_noart))
    columns = [f"delta_P{i}" for i in range(1, 4)]
    print("Largest absolute window shift by test:")
    print((windows_a[columns] - windows_noart[columns]).abs().max().to_string())


if __name__ == "__main__":
    main()
