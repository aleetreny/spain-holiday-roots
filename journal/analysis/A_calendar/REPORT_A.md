# Analysis A: the roots calendar — report

Date: 27 November 2026. Binding design: `journal/analysis/PREREGISTRATION.md` and `design.yml` (journal worktree, frozen 27 Nov 2026). Every number below is read from a CSV in `out/` (named in brackets). Heavy outputs (night parquet dataset, panels) live in `_local/journal_data/A_calendar/` of the journal checkout and are not in git.

## 1. What was done

| Step | Script | Output | Wall time (approx., from file timestamps) |
|---|---|---|---|
| 1. BOE labour calendar 2022–2024 + night rules | `a01_calendar.py` | `out/boe_holidays_2022_2024.csv` | < 1 min (done 15:34) |
| 2. Night dataset (all 2022–2024 files) + validation | `a02_build_dataset.py` | `pernoct_mun_2022_2024/` (date-partitioned parquet), `out/dataset_validation.csv` | ≈ 15 min (15:35–15:50) |
| 3. Night manifest + panels | `a03_manifest_and_panel.py` | `out/night_manifest.csv`, `panel/` (137-city, 675-city, 2,435-zone totals) | ≈ 1 min |
| 4. Nightly series + reproduction | `a04_nightly_series.py` (solver `a05_solver.py`) | `out/reproduction_check.csv`, `out/nightly_pooled.csv`, `out/nightly_cohorts.csv` | ≈ 8 min (15:50–15:58) |
| 6. Descriptive extras | `a06_descriptive.py` | `out/descriptive_*.csv` | ≈ 1 min |
| 5a. Pre-registered tests P1–P3, S1–S6 | `a05_preregistered_tests.py` | `out/preregistered_tests.csv` | ≈ 50 min (12 workers, 16:02–16:52) |
| 5b. Placebo windows | same | `out/placebo_windows.csv`, `out/placebo_summary.csv` | ≈ 23 min (16:52–17:15) |
| 5c. Robustness (§7) | same, `--skip-main --skip-placebo --workers 8` | `out/robustness.csv`, `out/robustness_lodo_dest.csv` | ≈ 63 min (17:51–18:54) |

Run everything with `/Users/alejandrotreny/Proyectos/spain-holiday-roots/.venv/bin/python -I journal/analysis/A_calendar/<script>.py` (set `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1` for a05).

## 2. Data checks

**Dataset.** 1,085 MITMA daily files 2022–2024 converted with the coursework processing (`b00_common` zoning, same artefact-zone removal). On all 306 discovery dates the per-(r, s) sums equal the coursework parquets: no missing or extra pairs, max |diff| = 9.3e-10 persons (float rounding) [`dataset_validation.csv`]. No processing difference.

**Reproduction.** All 43 coursework reference values are reproduced to max |diff| = 8.2e-15 [`reproduction_check.csv`]: pooled β 0.5324 (14–15 Feb 2024 night file 2024-02-14) and 0.6287 (2024-08-15), with identical SEs; cohort_joint pre-1981 0.0957 / 0.1793 and the stacked Aug − Feb contrast 0.0836; recent arrivals −0.0716.

**Nights excluded (17 calendar dates)** [`night_manifest.csv`, column `exclusion_reasons`]:
- 11 dates have no file on the ministry server: 2023-10-26, 27, 30, 31, 11-01, 02, 03 (Orange outage, ministry incident list), 2024-04-18, 2024-11-10 (also incident list), 2024-11-09, 2024-11-26.
- 6 files present but excluded: 2024-04-04, 2024-04-19, 2024-11-11 (rule a, ministry incident list); 2024-10-29, 30, 31 (rule b, DANA window).
- Rules c (residence-zone coverage < 99 %) and d (total deviates > 10 % from 28-night rolling median) excluded nothing: minimum coverage 1.00, maximum deviation 0.57 %.
- Result: 1,079 valid nights; 280 ordinary (O), 27 AUG, 3 XMAS, 3 EASTER, 2 ALLSAINTS (no 2023 file), 764 other.

**Calendar sources.** BOE-A-2021-17113 (2022), BOE-A-2022-16755 (2023), BOE-A-2023-22014 (2024), each with the `txt.php` and ELI URLs in `boe_holidays_2022_2024.csv`. Incident list: MITMA open-data page (`transportes.gob.es/.../opendata-movilidad`), accessed 27 Nov 2026.

## 3. An important caveat about the "confirmation set"

All 27 AUG nights (files 10–18 August 2022, 2023, 2024), the 2024 Christmas, Easter and All Saints nights, and the 2022/2024 All Saints nights were already in the coursework cache, so they are **discovery** nights [`night_manifest.csv`, crosstab of `night_class` × `split`]. The only *new* holiday nights are Christmas Eve and Good Friday of 2022 and 2023. Therefore:

- For P1–P3 and S3–S6, "confirmation" means **discovery holiday nights compared with never-opened ordinary nights** (split `confirmation_O`: 212 new O nights). The reference side is out of sample; the August side is not.
- For S1 and S2 a strict version is also possible (new holiday nights vs new O nights, 2022–2023): split `strict_confirmation`.
- The split `unseen_aug_vs_conf_O` (August nights other than 15 Aug, 2022–2023) is reported but should not be sold as out of sample: those nights were in the cache and the coursework used them as sister nights.

This is a design limitation of the pre-registration, not a choice made after seeing results. It is the strongest point a referee could raise against "confirmation". It can only be fixed with 2025 August nights.

## 4. Pre-registered tests on the confirmation set

Stacked PPML, every FE and covariate interacted with night, equal-weight average of year-specific contrasts, wild cluster restricted score bootstrap by origin province (22 clusters pooled, 44 cohort), 9,999 Webb draws, 95 % CI by test inversion, Holm within family [`preregistered_tests.csv`, `split == confirmation_O`, `year_scope == equal_year_avg`].

| Test | Holiday vs O nights | Estimate | 95 % CI | p (WCR) | Holm p | Verdict |
|---|---|---|---|---|---|---|
| **P1** August premium, pooled | 27 vs 212 | +0.106 | [+0.010, +0.197] | 0.033 | **0.066** | **Not confirmed** at Holm 0.05 |
| **P2** pre-1981 cohort × AUG | 27 vs 212 | +0.071 | [+0.041, +0.092] | 0.0001 | **0.0003** | **Confirmed** |
| **P3** roots × AUG × inland | 27 vs 212 | +0.122 | [−0.037, +0.280] | 0.112 | **0.112** | **Not confirmed** |
| S1 Christmas Eve | 3 vs 212 | +0.088 | [+0.049, +0.177] | 0.0001 | 0.0006 | holds |
| S2 Good Friday | 3 vs 212 | +0.037 | [−0.041, +0.110] | 0.249 | 0.249 | not detected |
| S3 All Saints | 2 vs 131 | +0.073 | [+0.022, +0.135] | 0.0055 | 0.0165 | holds |
| S4 village-born (≤ 2,000) × AUG | 27 vs 212 | +0.073 | [+0.055, +0.095] | 0.0001 | 0.0006 | holds |
| S5 recent arrivals × AUG | 27 vs 212 | −0.044 | [−0.070, −0.011] | 0.0079 | 0.0165 | holds (decrease) |
| S6 pre-1981 − recent × AUG | 27 vs 212 | +0.115 | [+0.071, +0.139] | 0.0001 | 0.0006 | holds |

Strict confirmation for S1/S2 (new holiday nights 2022–2023 only): S1 +0.087 [+0.042, +0.176], p 0.0003; S2 +0.025 [−0.076, +0.103], p 0.49. With a mixed S family (strict S1/S2 + confirmation_O S3–S6, column `p_holm_hybrid_S`) the Holm p-values are S1 0.0012, S2 0.49, S3 0.0165, S4 0.0006, S5 0.0165, S6 0.0006: same conclusions.

**Year by year (P1)**: +0.100 (2022), +0.106 (2023), +0.113 (2024). Positive every year, each with a CRV1 SE ≈ 0.04. P2 by year: +0.072, +0.074, +0.065.

### Placebo windows (diagnostic, not a test)

For every block of 9 consecutive ordinary-eligible nights in the same year, the same contrast [`placebo_summary.csv`]. On the confirmation set the AUG contrast is larger than **every** placebo window, for P1, P2 and P3, in all three years: rank 1 of 77 (2022), 1 of 74 (2023), 1 of 40 (2024). Largest placebo for P1: +0.083, +0.057, +0.109; for P2: +0.025, +0.018, +0.008. Same on all nights (rank 1 of 89, 85, 85).

### Discovery and all nights (reported next to the confirmation set)

| Test | Discovery: est [CI], Holm p | All nights: est [CI], Holm p |
|---|---|---|
| P1 | +0.143 [+0.048, +0.242], 0.020 | +0.108 [+0.005, +0.201], 0.084 |
| P2 | +0.077 [+0.043, +0.102], 0.0003 | +0.074 [+0.043, +0.097], 0.0003 |
| P3 | +0.172 [−0.021, +0.352], 0.078 | +0.139 [−0.035, +0.312], 0.097 |
| S1 | +0.087 [+0.052, +0.183], 0.0008 | +0.090 [+0.054, +0.178], 0.0008 |
| S2 | +0.055 [−0.037, +0.141], 0.18 | +0.039 [−0.048, +0.114], 0.27 |
| S3 | +0.126 [+0.076, +0.176], 0.0015 | +0.076 [+0.028, +0.138], 0.0057 |
| S4 | +0.086 [+0.064, +0.108], 0.0006 | +0.079 [+0.060, +0.101], 0.0006 |
| S5 | −0.041 [−0.068, −0.006], 0.034 | −0.045 [−0.073, −0.011], 0.020 |
| S6 | +0.118 [+0.065, +0.150], 0.0006 | +0.120 [+0.072, +0.146], 0.0006 |

## 5. Honest verdict

- **P2 is confirmed** and is the strongest result: the pre-1981 (exodus) cohort's roots elasticity rises in August by about +0.07, with a tight interval, on 44 clusters, every year, and above every placebo window. The cohort contrast (S6, +0.115) and the village-born term (S4, +0.073) agree; recent arrivals move the other way (S5, −0.044).
- **P1 fails the pre-registered decision rule.** The pooled August premium is +0.106 with a 95 % interval above zero and raw p = 0.033, but the Holm-adjusted p is 0.066 (> 0.05). By the rule in §3 it does **not** enter the abstract. It is positive in each year and ranks first among placebo windows, so it is "suggestive, not confirmed". The pooled model has only 22 clusters, which makes it wide; the discovery estimate (+0.143) was larger, consistent with some winner's-curse shrinkage.
- **P3 fails.** The inland − coast difference is positive (+0.122) but its interval includes zero. The coursework's "inland only, coast not detected" stays a descriptive pattern.
- Christmas Eve (S1) and All Saints (S3) show a roots premium; Good Friday (S2) does not.
- Caveat of §3: the August side of every confirmation test was seen before. P2 is "confirmed against new ordinary nights", not against new August nights.

## 6. Robustness (§7)

All on the confirmation split (discovery AUG vs new O nights), equal-year average, WCR bootstrap unless marked [`robustness.csv`, `robustness_lodo_dest.csv`]. Runtime: 47 min for the 24 checks, then 16 min for leave-one-out (8 workers).

| Check | P1 est [95 % CI], p | P2 est [95 % CI], p | P3 est [95 % CI], p |
|---|---|---|---|
| Main (for reference) | +0.106 [+0.010, +0.197], 0.033 | +0.071 [+0.041, +0.092], 0.0001 | +0.122 [−0.037, +0.280], 0.11 |
| Year-matched roots and income | +0.107 [+0.013, +0.191], 0.032 | +0.071 [+0.041, +0.092], 0.0003 | — |
| Road travel time instead of km | +0.099 [−0.007, +0.202], 0.069 | — | +0.152 [+0.022, +0.347], 0.031 |
| Drop Madrid origins | +0.119 [+0.031, +0.183], 0.014 | +0.071 [+0.036, +0.093], 0.0001 | — |
| Drop Barcelona origins | +0.071 [+0.017, +0.185], 0.018 | +0.064 [+0.031, +0.089], 0.0007 | — |
| Drop Madrid and Barcelona | +0.096 [+0.035, +0.171], 0.009 | +0.061 [+0.024, +0.091], 0.001 | +0.001 [−0.126, +0.244], 0.99 |
| Inland-province destinations only (28) | +0.154 [+0.039, +0.245], 0.020 | +0.147 [+0.091, +0.184], 0.0002 | — |
| Inland zones only (47 provinces) | +0.138 [+0.024, +0.223], 0.021 | +0.131 [+0.093, +0.157], 0.0001 | — |
| Coastal zones only (19 provinces) | +0.011 [−0.135, +0.134], 0.87 | −0.005 [−0.028, +0.014], 0.57 | — |
| P3 on full support (47 vs 19) | — | — | +0.127 [−0.020, +0.275], 0.074 |
| Seasonal reference (June/Sept O only) | +0.088 [+0.014, +0.168], 0.022 | +0.060 [+0.036, +0.080], 0.0001 | +0.079 [−0.049, +0.233], 0.18 |
| Equal weight by night | +0.110 [−0.007, +0.200], 0.060 | +0.071 [+0.038, +0.093], 0.0009 | +0.132 [−0.058, +0.289], 0.14 |
| Single δ across years | +0.103 [+0.007, +0.193], 0.041 | +0.071 [+0.041, +0.093], 0.0003 | +0.128 [−0.030, +0.293], 0.097 |
| Drop 3 artefact destination zones | +0.106 [+0.011, +0.196], 0.034 | — | — |
| Leave one destination province out (47) | range +0.094 (drop 18) to +0.136 (drop 28) | range +0.066 (drop 18) to +0.076 (drop 43) | — |

**Zero-cell rules (cohort model; Table A2 rules)** — P2 / S4 / S5 / S6 estimates:
- +1.0: +0.094 / +0.083 / −0.063 / +0.158
- +0.1: +0.032 / +0.039 / −0.021 / +0.053
- share + 1e-4: +0.051 / +0.076 / −0.031 / +0.082
- share + 1e-5: +0.014 (CI −0.003 to +0.025, p 0.093) / +0.025 / −0.012 (p 0.070) / +0.025
- asinh(1000 × share): +0.122 / +0.164 / −0.079 / +0.201
- Haldane 0.5 on the 112 large cities: +0.077 / +0.073 / −0.077 / +0.154

For S4–S6 under the zero-cell rules the WCR interval could not be inverted, so the CI shown in the CSV is CRV1 (column `ci_wcr_lo` empty); the WCR p is still reported.

**Reading.** P2 is positive with an interval above zero in every check except two: coastal destinations only (−0.005, not detected) and the smallest zero-cell offset (share + 1e-5, +0.014, p = 0.09). It moves by little under leave-one-out. Its *size* depends on the zero-cell rule (+0.014 to +0.122), as in the coursework Table A2: report the sign, not the size, as robust. P1 stays between +0.07 and +0.15, but its interval touches zero with road time (p 0.069) and with equal weighting by night (p 0.060). The August premium sits entirely in inland destinations (P1 inland +0.14 to +0.15; coastal +0.01), which supports P3's direction, but the formal P3 interaction is only above zero with road time (+0.152, p 0.031) and vanishes when Madrid and Barcelona origins are dropped (+0.001).

## 7. Descriptive extras (no tests)

[`descriptive_calendar_summary.csv`, `descriptive_top_nights.csv`, `nightly_pooled.csv`]

- **Mean pooled β by night type:** ordinary Mon–Wed 0.521 (280 nights); non-holiday Thursday 0.539; weekends 0.570 (Fri), 0.570 (Sat), 0.592 (Sun); July 0.597; early August (1–8) 0.634; peak (10–18) 0.635; late August (18–31) 0.645; Christmas Eve 0.617; All Saints 0.598; New Year's Eve 0.569; Good Friday 0.565; Epiphany eve 0.648.
- **Highest β nights** are the end of August and the first Sunday evenings of September (e.g. file 2022-08-31: 0.702; 2023-09-04: 0.698; 2022-09-01: 0.695), and Sunday evenings in mid-January (2024-01-22: 0.685). So the pooled elasticity peaks at the *end* of the holidays, not at 15 August; the peak flows (largest totals) are 10–18 August.
- **Inland vs coastal destinations:** in August the inland β (≈ 0.70–0.71) is far above coastal (≈ 0.46–0.48); on ordinary nights 0.55 vs 0.46.
- **Pre-1981 term** peaks in the 10–18 August window (mean 0.178) against 0.104 on ordinary nights; the recent-arrival term is lowest in August (0.119–0.132) and highest at Christmas/New Year and Epiphany (≈ 0.20–0.22): recent arrivals go home at Christmas, the exodus generation in August.
- Village receiving end, nightly: `descriptive_q5_villages_daily.csv`.

## 8. Deviations from the pre-registration

1. **Confirmation set composition** (§3 above): the pre-registration's confirmation set contains no new AUG nights; tests labelled `confirmation_O` use discovery AUG nights against new O nights. Not a change of method, but must be stated in the paper.
2. **Extra splits** `strict_confirmation` and `unseen_aug_vs_conf_O` and the mixed-family Holm column `p_holm_hybrid_S` are additions; the pre-registered column is `p_holm`.
3. **Recent cohort label (not a deviation).** The coursework label `recent_2016_2024` in `v1_c2_main_contrasts.csv` is only a name: `pipeline/v1_verify_n1.py:565` counts arrival years 2016–2023, as pre-registered.
4. **Run interruption.** The first full run of `a05` stopped during the robustness stage (session change); robustness was re-run alone with 8 workers. Main tests and placebos were not re-run and are unaffected (same seeds).
