# Analysis A with the three artefact zones excluded

This rerun implements the flow exclusion in §1 of the pre-registration for 2022–2024 (Deviation 2). I excluded `31905_AM` (Beriáin), `19103_AM` (Checa) and `40149_AM` (Olombrada). These are all inland destinations. I checked both the original city lists and the observed panel origins. The three zones occur in neither the 137-city pooled panel nor the 675-city cohort panel.

Both original panels contain `flow_inland_noart`. I independently reaggregated the 315 nights used by the models directly from the daily files, excluding both endpoints before grouping by city and destination province. The largest absolute difference from the saved inland or coastal values was less than 10⁻¹² persons in either panel. Applying A's coverage and rolling-volume rules to the 1,085 available files with the zones removed yielded the same exclusions. I therefore used A's night manifest, including the original split assignments. Following A2's flow-only sensitivity convention, I kept the province-centroid distance and other fixed model covariates as in A.

The wrapper loads both panels through A's `use_noart=True` option and calls A's original stacked PPML and placebo routines. It uses equal weight across years, the same 9,999 Webb draws, the same two-sided restricted score tests clustered jointly by origin province across all stacked nights, and the same Holm families. The pooled and cohort models have 22 and 44 origin-province clusters, respectively. The reported intervals for the equal-year averages are 95% WCR intervals obtained by test inversion.

For the August hypotheses, the `confirmation_O` split compares previously seen August nights with new ordinary nights. This is confirmation against new reference nights. For Christmas Eve and Good Friday, A also reports a strict confirmation split and a mixed secondary Holm family using those strict rows.

## Results

**No Holm decision changed.** This holds for the confirmation, discovery and all-nights splits, the additional splits that A reported, and the alternative Holm family using strict confirmation for S1 and S2. On the confirmation split, P2, S1 and S3–S6 remain supported in their specified directions. P1, P3 and S2 remain undetected at 0.05.

| Test | Published estimate [95% WCR CI] | Excluded-zone estimate [95% WCR CI] | Absolute estimate change | Holm p, published → excluded | Decision in both runs |
|---|---:|---:|---:|---:|---|
| P1 | +0.106379 [+0.010216, +0.196998] | +0.106258 [+0.009956, +0.197010] | 0.000120 | 0.0664 → 0.0668 | Not detected |
| P2 | +0.070750 [+0.041447, +0.092359] | +0.070669 [+0.041349, +0.092200] | 0.000081 | 0.0003 → 0.0003 | Supported |
| P3 | +0.122136 [-0.037296, +0.279604] | +0.122136 [-0.037296, +0.279604] | 0.000000 | 0.1116 → 0.1116 | Not detected |
| S1 | +0.088230 [+0.049328, +0.176708] | +0.088192 [+0.049273, +0.176633] | 0.000039 | 0.0006 → 0.0006 | Supported |
| S2 | +0.036610 [-0.040537, +0.109626] | +0.036610 [-0.040663, +0.109652] | 0.000001 | 0.2493 → 0.2489 | Not detected |
| S3 | +0.073122 [+0.022432, +0.134884] | +0.073233 [+0.022325, +0.135046] | 0.000111 | 0.0165 → 0.0171 | Supported |
| S4 | +0.073142 [+0.055170, +0.094618] | +0.073094 [+0.055064, +0.094628] | 0.000047 | 0.0006 → 0.0006 | Supported |
| S5 | -0.044017 [-0.069728, -0.011301] | -0.043998 [-0.069717, -0.011272] | 0.000019 | 0.0165 → 0.0171 | Supported |
| S6 | +0.114767 [+0.071190, +0.139198] | +0.114667 [+0.071113, +0.139084] | 0.000100 | 0.0006 → 0.0006 | Supported |

The largest absolute change in a confirmation estimate was **0.000120** for P1 (+0.106379 to +0.106258). Its rerun estimate matches the artefact-destination check in A's `robustness.csv` exactly. P3 uses city–province pairs in coastal provinces, which excludes all three affected provinces from its support in both runs.

Across the three main splits, the largest absolute change among equal-year-average estimates was **0.000350** for discovery S1. Including the year-specific rows, the largest was **0.000420** for all-nights P1 in 2024. The largest change to a 95% WCR interval endpoint was **0.000391**, for the lower endpoint of discovery P1. The largest absolute change in a Holm p-value was **0.0008**, for all-nights P1 (0.0844 to 0.0852). The paired CSV includes absolute differences for both interval endpoints and for the Holm p-values.

## Placebo windows

All 444 placebo windows and six August windows match A's dates and sample sizes. August ranks first for P1–P3 in each confirmation year (1 of 77, 74 and 40) and in each all-nights year (1 of 89, 85 and 85). All 18 ranks agree with A. The largest absolute shift in a window contrast was 0.000459 for P1 in the 2024 all-nights placebo diagnostic.

## Reproduction

Analysis A was read from commit `3d16db437f6ccedf4262fa42cec433ba80ef0a5b` of the journal checkout. The wrapper reads A's manifest and the saved 137-city and 675-city parquets under the journal checkout's `_local/journal_data/A_calendar/panel/`. It writes its outputs in this directory. Run both commands from the repository root.

```sh
.venv/bin/python -I journal/analysis/A_calendar_noart/run_noart.py
.venv/bin/python -I journal/analysis/A_calendar_noart/compare_noart.py
```

[Tests](out/preregistered_tests_noart.csv) · [Paired comparison](out/comparison_noart.csv) · [Placebo summary](out/placebo_summary_noart.csv) · [Individual placebo windows](out/placebo_windows_noart.csv)
