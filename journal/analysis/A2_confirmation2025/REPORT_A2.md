# Confirmation II | The August roots premium in 2025

27 November 2026. The binding rules are in the [Confirmation II addendum](../PREREGISTRATION.md#addendum-confirmation-ii-frozen-27-nov-2026-before-any-2025-mitma-file-is-downloaded-or-opened), frozen in commit `8506c36`. The comparison is between unseen August nights in 2025 and ordinary nights in the same year.

## Data and reproduction

I reused the conversion, calendar, quality rules, panels, fixed-effect PPML and restricted score bootstrap from Analysis A (commit `3d16db4`). I ran its original reproduction suite before opening the 2025 data. The 49 recorded checks matched Analysis A, with a maximum absolute error against the coursework of 8.22 × 10⁻¹⁵. For the night of 13–14 February 2024, the pooled elasticity was 0.532359 and the pre-1981 term was 0.095693. [Reproduction checks](out/reproduction_check.csv).

I verified the recorded SHA-256 and size of each of the 202 supplied files before conversion. All 196 standard files passed the coverage and volume rules. The lowest residence-zone coverage was 100% of the median. The largest deviation from the centred rolling population median was 1.18%. The six recalculations were reserved for the paired bridge. The two April blackout files were absent from the ministry inventory and were flagged as incidents. Other unavailable dates are marked as missing. [Night manifest](out/night_manifest.csv), [Analysis A inventory audit](../../research/design_first_principles.md), [MITMA RSS](https://movilidad-opendata.mitma.es/RSS.xml), accessed 27 November 2026.

The evening preceding the file date sets calendar eligibility. The public-holiday screen covers national holidays and the five specified regions. The 2025 ordinary sample includes 61 old-method nights and one new-method night, 15–16 October. Andalusia's observed National Day holiday rules out the evenings of 13 and 14 October. The August sample covers the nights of 9–10, 13–14, 14–15, 15–16 and 16–17 August. The Easter file covers the night of 17–18 April, beginning on Holy Thursday. [Holiday table](out/boe_holidays_2025.csv), [Contrast dates](out/contrast_nights.csv), [BOE 2025 labour calendar](https://www.boe.es/diario_boe/txt.php?id=BOE-A-2024-21316), accessed 27 November 2026.

## Step 1 | Method bridge

On 1 July 2025, the ministry changed how it filters static devices from the mobile-phone sample. It recalculated six earlier nights with the new procedure and cautions that the adjustment varies across origin–destination pairs. [MITMA methodology note, §§2–6](https://movilidad-opendata.mitma.es/LEEME_CambioMetodolog%C3%ADa_IndicadoresBigdata_posteriores_a_20250701.pdf), accessed 27 November 2026.

| Term | Mean absolute difference | Upper tolerance | Verdict |
| --- | ---: | ---: | --- |
| Pooled roots elasticity | 0.013817 | 0.035 | Pass |
| Pre-1981 arrivals | 0.005424 | 0.024 | Pass |

Both differences are below their thresholds. The preregistered rule therefore uses all 62 eligible ordinary nights. These six pairs assess comparability for the two specified coefficients. They cover no August nights and provide no direct calibration of the inland–coastal contrast or the other cohort terms. [Bridge summary](out/method_bridge_summary.csv).

| Night ending | Pooled old → new | Pre-1981 old → new | All person-nights, old → new (millions) | Away person-nights, old → new (millions) | Away pairs changed |
| --- | --- | --- | --- | --- | --- |
| 2024-10-15 | 0.48072 → 0.48335 | 0.09230 → 0.10579 | 44.371 → 44.041 | 2.628 → 2.743 | 99.974% |
| 2024-10-18 | 0.47332 → 0.47503 | 0.09497 → 0.09771 | 44.362 → 44.027 | 2.783 → 2.955 | 99.974% |
| 2024-10-19 | 0.51740 → 0.51287 | 0.10363 → 0.10692 | 44.353 → 43.996 | 3.955 → 4.182 | 99.971% |
| 2025-01-13 | 0.59114 → 0.61289 | 0.11506 → 0.10866 | 44.405 → 44.077 | 1.933 → 2.044 | 99.973% |
| 2025-01-19 | 0.54511 → 0.57415 | 0.11174 → 0.10920 | 44.387 → 44.033 | 3.670 → 3.841 | 99.980% |
| 2025-02-08 | 0.54282 → 0.56606 | 0.10628 → 0.10220 | 44.381 → 44.035 | 3.309 → 3.513 | 99.978% |

Pairs are directed flows between distinct mainland zones with positive weight under either method. For this diagnostic I counted a missing weight as zero and a pair as changed when the absolute weight difference exceeded 10⁻⁸ persons. Across the six bridges, total mainland person-nights changed by −0.80% to −0.74%, while person-nights away from the residence zone changed by +4.37% to +6.18%. A small movement in the two coefficients can therefore coexist with substantial changes in the flow weights. [Paired coefficients and flow totals](out/method_bridge_dates.csv).

## Step 2 | Stacked tests

I used Analysis A's night-specific controls and fixed effects, and its common-support inland–coastal specification for P3. The outcomes use the published flow weights, including weights below 10. The intervals invert the two-sided wild cluster restricted score test with 9,999 Webb draws. The pooled and cohort models cluster jointly across nights by origin province, using 22 and 44 clusters respectively. Holm adjustment covers three primary tests and the four available secondary tests. The interval for each coefficient is at the individual 95% level. The main pooled roots measure comes from the 2024 annual census. The arrival and birthplace covariates come from the 2024 full-variable microdata and retain Analysis A's cohort definitions. [Tests](out/preregistered_tests.csv), [Kline and Santos (2012)](https://doi.org/10.1515/2156-6674.1006), [INE table 68545](https://www.ine.es/jaxiT3/Tabla.htm?t=68545), [Annual Census and microdata](https://www.ine.es/dyngs/INEbase/es/operacion.htm?c=Estadistica_C&cid=1254736176992&menu=resultados&idp=1254735572981), accessed 27 November 2026.

| Test | Holiday nights | Ordinary nights | Estimate | 95% WCR CI | Raw p | Holm p |
| --- | --- | --- | --- | --- | --- | --- |
| P1 Pooled August premium | 5 | 62 | +0.134 | [+0.037, +0.233] | 0.0150 | 0.0300 |
| P2 Pre-1981 arrivals | 5 | 62 | +0.057 | [+0.022, +0.076] | 0.0033 | 0.0099 |
| P3 Inland minus coast | 5 | 62 | +0.129 | [−0.055, +0.319] | 0.1364 | 0.1364 |
| S2 Easter pooled premium | 1 | 62 | +0.076 | [−0.013, +0.157] | 0.0765 | 0.1530 |
| S4 Born in municipalities of ≤2,000 | 5 | 62 | +0.056 | [+0.031, +0.078] | 0.0010 | 0.0040 |
| S5 Recent arrivals, 2016–2023 | 5 | 62 | −0.033 | [−0.060, +0.006] | 0.0842 | 0.1530 |
| S6 Pre-1981 minus recent arrivals | 5 | 62 | +0.089 | [+0.034, +0.117] | 0.0010 | 0.0040 |

P1, P2, S4, S6 pass Holm at 0.05 in the prespecified direction. [Full test output](out/preregistered_tests.csv).

### Placebo windows and census sensitivity

The August ranks against all 54 rolling nine-night ordinary windows are 1/55 for P1, 1/55 for P2, 1/55 for P3. Each placebo window is compared with the remaining ordinary nights. Rank 1 denotes the largest contrast. The ranks are diagnostic. P3 still has a wide interval and misses the Holm threshold. August has five observed nights and the placebo windows have nine, as specified. [Window estimates](out/placebo_windows.csv), [ranks](out/placebo_summary.csv).

| Test | 2025 census estimate | 95% WCR CI | Holm p |
| --- | --- | --- | --- |
| P1 | +0.131 | [+0.032, +0.233] | 0.0332 |
| P3 | +0.132 | [−0.055, +0.323] | 0.1352 |
| S2 | +0.073 | [−0.016, +0.156] | 0.1684 |

[Census-year sensitivity](out/roots_year_sensitivity.csv).

The [pooled nightly series](out/nightly_pooled.csv) and the [cohort nightly series](out/nightly_cohorts.csv) cover each valid standard 2025 file, with the ministry method and per-night cluster intervals. The paired recalculations are in the [bridge table](out/method_bridge_dates.csv).

## Combined decision

The addendum requires the same sign in both periods and a Holm-adjusted result in the predicted direction in at least one period. The abstract can report P1 (pooled August premium), P2 (pre-1981 arrivals), S4 (born in municipalities of ≤2,000), S5 (recent arrivals, 2016–2023), S6 (pre-1981 minus recent arrivals). [Combined decisions](out/combined_decision.csv), [Analysis A tests](../A_calendar/out/preregistered_tests.csv).

| Test | 2022–2024 estimate | Holm p | 2025 estimate | Holm p | Same sign | Abstract |
| --- | --- | --- | --- | --- | --- | --- |
| P1 | +0.106 | 0.0664 | +0.134 | 0.0300 | Yes | Yes |
| P2 | +0.071 | 0.0003 | +0.057 | 0.0099 | Yes | Yes |
| P3 | +0.122 | 0.1116 | +0.129 | 0.1364 | Yes | No |
| S2 | +0.037 | 0.2493 | +0.076 | 0.1530 | Yes | Outside the August rule |
| S4 | +0.073 | 0.0006 | +0.056 | 0.0040 | Yes | Yes |
| S5 | −0.044 | 0.0165 | −0.033 | 0.1530 | Yes | Yes |
| S6 | +0.115 | 0.0006 | +0.089 | 0.0040 | Yes | Yes |

The 2022–2024 results use previously analysed August nights and unseen ordinary nights (Deviation 1). The 2025 contrast uses unseen nights on both sides. P1 qualifies through the 2025 result. S5 qualifies through the 2022–2024 result and has the predicted sign in 2025. The estimates describe associations between a city's birth-province composition and its overnight flows to other provinces. The census and aggregate mobility records do not identify the travellers or the destination villages of their families.

## Implementation decisions and limitations

- The main 2025 fits use Analysis A's existing artefact-zone-excluded flow column to implement §1. Analysis A's published headline estimates retain the three artefact zones as destinations. I left those published estimates intact for the combined rule and checked the effect of using A's original flow convention in 2025. The largest absolute change across the seven estimates was 0.00047. All seven Holm verdicts were the same. [Original A panel sensitivity](out/artefact_sensitivity.csv).
- I applied A's 28-observation centred rolling median, with a minimum of 14 available observations. Missing dates were skipped in the window. In the sparse post-change sample, the reference window can span gaps in calendar time. Both the 28- and 29-observation checks gave the same exclusions.
- I used the six recalculations as fixed paired bridge inputs. Their sparse spacing does not support a standalone quality median for the recalculated series. The three 2024 originals passed A's quality screen and the three 2025 originals passed the 2025 screen.
- The ministry incident page returned HTTP 403 on recheck. I used the two blackout exclusions documented in the pre-analysis MITMA inventory. The file coverage and volume screens ran on all supplied standard files. [Ministry incident page](https://www.transportes.gob.es/ministerio/proyectos-singulares/estudios-de-movilidad-con-big-data/opendata-movilidad), attempted 27 November 2026.
- The addendum changes the pooled exposure table year in the sensitivity. Income and the 2024 microdata cohort exposures stay fixed. The sensitivity reapplies Holm within each full family, with the fixed cohort p-values carried over.
- S1 and S3 have no available 2025 holiday files. The combined abstract rule applies to the August claims. S2 is reported as a secondary test.

## Reproduction

From this worktree, use the existing environment and the local data snapshots specified in the brief. The wrapper imports the copy of A in the worktree's ignored `_local/A_reference/` when present and otherwise uses the journal checkout. The reference version is `3d16db4`. The raw inputs, converted parquets and panels are kept in the journal checkout's private data area. The scripts read the supplied MITMA files without downloading them.

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
PYTHON=/Users/alejandrotreny/Proyectos/spain-holiday-roots/.venv/bin/python
A2=journal/analysis/A2_confirmation2025
"$PYTHON" -I -B "$A2/a2_data.py"
"$PYTHON" -I -B "$A2/a2_models.py" tests --workers 7
"$PYTHON" -I -B "$A2/a2_models.py" nightly --workers 5
"$PYTHON" -I -B "$A2/a2_models.py" placebos --workers 3
```

To rerun the original Analysis A panel treatment, add `--include-artefact-zones --output-dir <private-output-directory>` to the tests command.
