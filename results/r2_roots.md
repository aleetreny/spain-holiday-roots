# R2: does the holiday network align with the roots network more than the ordinary network does?

*Agent R2 (wave 2), 2026-11-24. Code: `pipeline/r2_roots.py` (one run, 633 s, log `results/r2_run.log`). Every number below is read from the `results/r2_*.csv` file named next to it. Decisions: `[R2]` lines in `.planning/2026-11-24-consolidation/decisions.md`.*

**Dates.** File date d is the night ending on d (config). "Aug 15" is the night of 14 to 15 August, "Nov 1" the night of 31 Oct to 1 Nov, "Dec 25" Christmas Eve, "Mar 29" the night of Holy Thursday to Good Friday.

## Answer in one paragraph
Residents of a city sleep more in the provinces where more of their fellow residents were born. This holds on every night, ordinary or holiday, and it holds when two cities of the same province are compared against the same destination, beyond distance and income. The roots elasticity is 0.53 on an ordinary February night and 0.63 on Aug 15. Holidays strengthen the link only modestly: about +0.05 to +0.10 against February, and +0.08 to +0.14 against October. Almost all of the summer increase comes from inland destinations: inland zones +0.16 against Feb 14, coastal zones −0.05. With 22 province clusters, the summer increase over an ordinary February night is not statistically secure (p = 0.18). It is secure against the October baseline (p = 0.006) and for inland provinces (p = 0.010). Christmas Eve is the only holiday whose increase passes against both baselines. Across destination provinces (pattern #7's index), August is not more roots-aligned than a weekday, because beaches take a large share of summer nights whatever a city's roots are. The honest headline is: **roots shape where Spaniards sleep every night; holidays, the family ones and inland summer, add a modest extra pull; beaches dominate the summer total.**

## 1. Design
- **Roots.** INE Censo anual, 1 Jan 2024, residents of each city by province of birth (table 68545, 151 cities), and the 47 × 47 province matrix (table 68523). Files: `r2_roots_city_prov.csv`, `r2_roots_prov_prov.csv`. One roots year (2024) is used for every night.
- **City set.** The 137 mainland census cities (`r2_city_set.csv`). Each is a single-municipality MITMA zone. The neighbouring-city design is identified only by the **112 cities in the 22 provinces with at least two census cities**: Madrid 23, Barcelona 19, Alicante, Cádiz and Málaga 8 each, and so on. Single-city provinces are absorbed by their pair effects.
- **Outcome.** For each single night (never a window mean), persons resident in city c sleeping in mainland zones of province p ≠ prov(c). The data are the b01 full edge lists. Nights that b01 did not save were built with `b00.build_edges`: 03-29, 11-03, 12-25 and October.
- **PPML (main specification).** The roots regressor is log(share of residents born in p). Controls are log km, (log km)², and 46 city-income × destination-province slopes (ADRH 2024, net income per person). Fixed effects: city, and origin province × destination province. The solver is our own; β equals pyfixest to 1e-6 on all five main nights (`r2_solver_check.csv`).
- **Inference.**
  - Clusters are origin provinces (G = 22).
  - Every p-value and CI marked "WCR" comes from the wild cluster restricted score bootstrap (Kline & Santos 2012), with Webb weights and 9,999 draws. CIs are obtained by test inversion.
  - CRV1 t-intervals are also in the CSVs.
  - No normal-approximation p-values are used.
- **Difference test (new).** A stacked two-night PPML with night-specific fixed effects and covariates, a common roots slope and a roots × holiday term δ. H0: δ = 0, tested by WCR. The unrestricted δ equals β_holiday − β_ordinary exactly (asserted).
- **QAP.**
  - Null: doubly-constrained gravity, fitted by IPF on rows, columns and distance-bin totals. Margins are matched to < 1e-9 (`r2_gravity_convergence.csv`).
  - City design: residuals are demeaned within origin province × destination province. Cities are permuted within their province.
  - Province design: simultaneous row and column permutation.
  - Partial QAP is **Dekker, Krackhardt & Snijders (2007) double semi-partialling** (it permutes the residual of roots on the ordinary-night residual). It is not Freedman–Lane.
  - **Difference QAP:** r(Z_holiday − Z_ordinary, Z_roots). Each test uses 2,000 permutations; p values are one-sided, with floor 1/2001.
- **Baselines.**
  - The 13 February Tue/Wed/Thu nights (Feb 14 is the main night).
  - The mean of those 13 nights.
  - The 15 October Tue/Wed/Thu nights and their mean (pattern #12).
  - Placebo: October mean − February mean = **−0.053**, WCR p = 0.12, CI [−0.127, 0.013] (`r2_diff_tests.csv`, baseline `placebo`). The ordinary baseline itself moves by about half the size of the holiday effects. This is why both baselines are reported.

## 2. Headline numbers (main specification, all destinations)
Single nights (`r2_ppml_nights.csv`, `r2_ppml_summary.csv`):

| Night | β roots | WCR 95 % CI | Sister median [range] |
|---|---:|---|---|
| Ordinary, Feb 14 | 0.532 | [0.382, 0.694] | 0.534 [0.515, 0.600], 13 nights |
| Ordinary, Feb mean (13) | 0.544 | [0.389, 0.716] | — |
| Ordinary, Oct mean (15) | 0.492 | [0.380, 0.621] | Oct nights 0.487 [0.453, 0.542] |
| Summer, Aug 15 | 0.629 | [0.520, 0.769] | 0.621 [0.596, 0.649], 9 nights |
| All Saints, Nov 1 | 0.584 | [0.488, 0.724] | 0.541 [0.539, 0.584], 3 nights; Nov 1 is the maximum |
| Christmas Eve, Dec 25 | 0.598 | [0.487, 0.787] | single night (appendix) |
| Good Friday, Mar 29 | 0.566 | [0.485, 0.761] | single night (appendix) |
| 2022: Feb 16 / Aug 15 / Nov 1 | 0.413 / 0.643 / 0.612 | [0.329, 0.541] / [0.526, 0.798] / [0.504, 0.765] | — |
| 2023: Feb 15 / Aug 15 | 0.538 / 0.649 | [0.377, 0.723] / [0.551, 0.804] | — |

Every night: WCR p(β = 0) ≤ 0.002. An elasticity of 0.53 means that a city with twice the share of natives of province p sends 2^0.53 = 1.44 times as many nights to p, compared with a neighbouring city of the same province with the same distance and income profile. In August the factor is 2^0.63 = 1.55.

**Difference holiday − ordinary** (`r2_diff_tests.csv`, `r2_diff_summary.csv`):

| Contrast | δ | WCR p | WCR 95 % CI | Over sister pairs: median [range]; pairs with δ > 0; pairs with WCR p < 0.05 |
|---|---:|---:|---|---|
| Aug 15 − Feb 14 | **+0.096** | 0.183 | [−0.031, 0.242] | 9 × 13 pairs: 0.080 [−0.005, 0.134]; 116/117; 1/117 |
| Aug 15 − Oct mean | +0.137 | 0.006 | [0.048, 0.238] | 15 Oct nights: 0.141 [0.086, 0.175]; 15/15; 14/15 |
| Nov 1 − Feb 14 | +0.052 | 0.221 | [−0.025, 0.160] | 3 × 13: 0.010 [−0.062, 0.070]; 27/39; 0/39 |
| Nov 1 − Oct mean | +0.092 | 0.001 | [0.034, 0.179] | 0.097 [0.042, 0.131]; 15/15; 14/15 |
| Dec 25 − Feb 14 | +0.065 | **0.019** | [0.013, 0.174] | 13 Feb nights: 0.064 [−0.003, 0.083]; 12/13; 8/13 |
| Dec 25 − Oct mean | +0.106 | <0.001 | [0.067, 0.210] | 0.110 [0.055, 0.144]; 15/15; 15/15 |
| Mar 29 − Feb 14 | +0.034 | 0.463 | [−0.053, 0.139] | 0.033 [−0.034, 0.052]; 11/13; 0/13 |
| Mar 29 − Oct mean | +0.075 | 0.006 | [0.026, 0.154] | 0.079 [0.024, 0.113]; 15/15; 13/15 |
| 2022 Aug 15 − Feb 16 | +0.230 | 0.002 | [0.112, 0.330] | — |
| 2022 Nov 1 − Feb 16 | +0.199 | 0.001 | [0.087, 0.322] | — |
| 2023 Aug 15 − Feb 15 | +0.111 | 0.141 | [−0.035, 0.232] | — |

- **Sign.** The summer increase has the right sign on 116 of 117 February pairs. It is significant on only 1 of 117, so it is consistent but small relative to province-level heterogeneity.
- **2022.** The large 2022 contrasts come from a low February 2022 baseline (0.413). Pattern #12 shows that Feb 2022 away volume was 1.17 M against 2.11 M in 2024. Treat 2022 as fragile.
- **Barcelona.** The summer difference leans on Barcelona (`r2_diff_jackknife.csv`). Dropping its 19 cities moves Aug 15 − Feb 14 from 0.096 to 0.037 for all destinations, and from 0.155 to 0.086 for inland destinations. No other single province moves it below 0.088 (all destinations).
- **Variants** (Aug 15 − Feb 14). With road time instead of km: +0.089 (p 0.33). Without income controls: +0.160 (p 0.16), against Feb 14 β = 0.492. About 40 % of the raw summer increase goes with income: poorer cities have more natives of other provinces.

## 3. Beach confound: inland-only, coastal-only, inland provinces
Same model, with the outcome restricted to particular destination zones (`r2_ppml_summary.csv`, `r2_diff_summary.csv`). As audit m7 warns, this changes the outcome population, so it decomposes the effect rather than testing the confound.

| Destinations | Feb 14 | Aug 15 | Nov 1 | Dec 25 | Mar 29 | Aug − Feb 14 (WCR p) | Aug − Oct mean (p) | Nov 1 − Oct mean (p) | Dec 25 − Feb 14 (p) |
|---|---:|---:|---:|---:|---:|---|---|---|---|
| All | 0.532 | 0.629 | 0.584 | 0.598 | 0.566 | +0.096 (0.18) | +0.137 (0.006) | +0.092 (0.001) | +0.065 (0.019) |
| Inland zones (coastal_true = 0) | 0.549 | 0.704 | 0.597 | 0.611 | 0.626 | +0.155 (0.086), CI [−0.013, 0.292] | +0.184 (0.001) | +0.076 (0.008) | +0.062 (0.034) |
| Coastal zones only | 0.501 | 0.451 | 0.491 | 0.550 | 0.346 | −0.050 (0.57), CI [−0.229, 0.153] | +0.006 (0.93) | +0.046 (0.26) | +0.050 (0.29) |
| 28 inland provinces | 0.492 | 0.719 | 0.577 | 0.602 | 0.633 | **+0.227 (0.010)**, CI [0.057, 0.345] | +0.208 (0.006) | +0.066 (0.014) | +0.110 (0.005) |

- **Sister pairs, inland zones.** Aug − Feb has median +0.138 [0.040, 0.197], positive on 117/117 pairs, and WCR p < 0.05 on 29 of 117 (25 %). Against October, 15/15 pairs are significant.
- **Sister pairs, coastal zones.** Aug − Feb has median −0.087, negative on 117/117 pairs, with none significant.
- **All Saints (Nov 1 − Feb 14) by destination.** Inland +0.048 (p 0.37), coastal −0.010, inland provinces +0.085 (p 0.10).

The summer roots premium is an inland phenomenon. At the coast, roots matter no more in August than on a February weekday.

## 4. QAP (`r2_qap.csv`, `r2_qap_summary.csv`)

| Statistic | Feb 14 | Aug 15 | Nov 1 | Dec 25 | Mar 29 |
|---|---:|---:|---:|---:|---:|
| City design, r(Z_night, Z_roots), within-province permutation | 0.228 | 0.483 | 0.354 | 0.449 | 0.424 |
| City design, DSP partial r given Feb 14 | — | 0.455 | 0.321 | 0.417 | 0.384 |
| City design, **difference r(Z_h − Z_Feb14, Z_roots)** | — | **0.105** (p < 0.001) | 0.053 (p 0.001) | 0.096 (p < 0.001) | 0.080 (p < 0.001) |
| Province 47 × 47, r(Z_night, Z_roots) | 0.539 | 0.540 | 0.560 | 0.690 | 0.542 |
| Province, DSP partial r given Feb 14 | — | 0.352 | 0.338 | 0.522 | 0.328 |
| Province, **difference r** | — | 0.061 (p 0.029) | 0.010 (p 0.37) | **0.161** (p < 0.001) | 0.010 (p 0.38) |

- **Sister robustness, city design.** Every holiday sister × February sister difference is positive and significant: summer 117/117, maximum p 0.001; All Saints 39/39, maximum p 0.033. The October pairs are also 15/15 for every holiday.
- **Sister robustness, province design.** Summer passes on 77 of 117 pairs (66 %, maximum p 0.33). All Saints and Good Friday mostly fail against February. Christmas passes 13/13.
- **Two tests, two nulls.** The city-design QAP rejects for every holiday, while the PPML difference test often does not. The QAP null ("cities of a province are exchangeable") forbids province-level heterogeneity in the roots slope. The cluster bootstrap allows it. The QAP is also unweighted (every city × province cell counts the same), while PPML weights cells by volume, which makes Barcelona heavy. I recommend the cluster bootstrap as the confirmatory test and the QAP as the network-science description.

## 5. Reconciliation with the pattern-mining roots index (pattern #7)
The index of `p_patterns.md` is reproduced exactly: max |diff| 1e-15 over 47 nights against `p_daily_metrics.roots_corr_partial`. File: `r2_index_bridge.csv`.

**What each statistic identifies.**
- **The index** is a *correlation* across all 137 cities × other provinces between log(1+flow) and log(1+born). It removes city means and log distance, and gives every cell equal weight. Its variation is between destination provinces: does a city send more nights to the provinces its residents come from, *compared with other provinces at the same distance*? A correlation falls whenever other causes add variance to flows. In August, beach provinces receive many nights regardless of roots.
- **The PPML elasticity** is a *slope* identified *within* origin province × destination province. It compares two cities of the same province sending to the same destination. Anything about the destination (beaches, size, prices) and about the corridor cancels. It asks: does a city with more natives of p send more nights to p *than its neighbour does*? It is weighted by volume.

**Bridge** (median over sister nights, all destinations; `r2_index_bridge.csv`):

| Statistic | Ordinary Feb | Oct | Summer | All Saints | Christmas (Dec 25) | Good Friday |
|---|---:|---:|---:|---:|---:|---:|
| Index: corr, between provinces (= pattern #7) | 0.590 | 0.630 | 0.603 | 0.681 | 0.768 | 0.668 |
| OLS slope, between provinces | 0.699 | 0.682 | 0.645 | 0.692 | 0.783 | 0.675 |
| corr, within origin prov × dest prov | 0.205 | 0.240 | 0.459 | 0.299 | 0.386 | 0.382 |
| OLS slope, within | 0.428 | 0.408 | 0.610 | 0.474 | 0.561 | 0.551 |
| PPML β (main spec) | 0.534 | 0.487 | 0.621 | 0.541 | 0.598 | 0.566 |

**Difference tests of the index.** These use a cluster (origin province) pairs bootstrap with 1,999 draws (`r2_index_diff_tests.csv`, `r2_index_diff_summary.csv`).
- Between-province index:
  - Aug 15 − Feb 14 = +0.009, CI [−0.039, 0.055], p 0.78. Over 117 sister pairs, the median is +0.013 and none is significant.
  - Against the Feb mean it is −0.048 (p 0.03), and against the Oct mean −0.082 (p < 0.001).
  - Inland destinations only: +0.116 (p < 0.001) against Feb 14, but −0.014 (p 0.53) against the Oct mean.
  - Coastal destinations only: −0.149 (p < 0.001).
  - Christmas: +0.167 against Feb 14 and +0.076 against the Oct mean (both p < 0.001).
  - All Saints: +0.086 (p < 0.001) against Feb 14 and −0.005 against the Oct mean.
- Within-province index:
  - Aug 15 − Feb 14 = +0.270, CI [0.156, 0.355]. All 117 sister pairs are significant.
  - Against the Oct mean: +0.104 (p 0.011).

So pattern #7 and R2 agree once each is read for what it measures:
- **Between provinces**, August is not more roots-aligned than a weekday. Beaches supply a large, roots-independent share of summer nights: the coastal-only index falls by 0.15. Christmas, and to a smaller degree All Saints and Good Friday against February, are more aligned.
- **Within provinces** (neighbouring cities, same destination), August shows the largest roots premium of any night. The premium is concentrated in inland destinations.
- **Formal test.** With 22 province clusters and volume weighting, the August premium over an ordinary February night is not statistically secure (p 0.18). It is secure against the October baseline, for inland provinces, and in the city-design QAP. Christmas Eve is the one holiday significant in every test.

## 6. Reproduced vs changed numbers relative to legacy

| Legacy claim (source) | Legacy value | Reproduced here | Clean value (this run) | Why it changed |
|---|---|---|---|---|
| City × province PPML, mid-Aug, c + pairFE (p72) | 0.721 | 0.721 (L0); 0.722 with correct mean (L1); 0.720 on 137 mainland cities (L2) | Aug 15 legacy spec 0.731; **main spec 0.629** [0.520, 0.769] | Log share + distance + income × destination; clusters 22 provinces instead of 151 cities (SE 0.023 → 0.052) |
| Same, Feb weekday (p72) | 0.633 | 0.633 (L0); 0.644 (L1); 0.651 (L2) | Feb 14 legacy spec 0.650; **main 0.532** [0.382, 0.694] | As above |
| Income-robust c + pairFE + p[linc], mid-Aug (Check 1A) | 0.717 (0.032) | — | 0.629 (0.052 CRV1) | Adds distance, 137 mainland cities, mainland destinations, single night |
| "Holiday jump 1.8–2×" (Check 4, zone level) | 0.29 → 0.52 | not reproduced (no code, audit M5) | Aug − Feb 14 = +0.096 (×1.18), p 0.18 | The difference was never tested |
| Inland zones, pairFE, mid-Aug / All Saints (Check 6) | 0.798 / 0.734 | — | 0.704 / 0.597 | Specification as above |
| Coastal only, pairFE, mid-Aug (Check 6) | 0.534 | — | 0.451 | — |
| "Removing coast strengthens the effect" (Check 6) | monotone | yes in level | Inland Aug − Feb +0.155 (p 0.086); inland provinces +0.227 (p 0.010) | Now a tested difference |
| Multi-year, mid-Aug 2022 / 2023 (Check 2) | 0.733 / 0.733 | — | Aug 15: 0.643 / 0.649 | Single night, main spec |
| Feb weekday 2022 (Check 2) | 0.589 | — | Feb 16 2022: 0.413 | Low 2022 baseline (pattern #12) |
| City neighbour QAP, Feb / mid-Aug (p93) | 0.330 / 0.559 | — | Feb 14 0.228 / Aug 15 0.483 | Single nights; roots permuted, gravity on 137 cities |
| City partial QAP, mid-Aug (README "0.49") | 0.493 | — | DSP 0.455 (Aug 15 given Feb 14) | Named correctly (DSP, not Freedman–Lane) |
| Province spatial-null r, Feb / mid-Aug (p93) | 0.573 / 0.545 | — | 0.539 / 0.540 | Single nights |
| Province partial QAP, Aug excess (p93) | 0.085, p 0.02 | — | Difference QAP 0.061, p 0.029 (sisters up to p 0.33) | Direct test of the difference |
| Echo, inland: b98 \| recent in + out, Feb / mid-Aug (Check 5) | −0.042 / +0.295 | — | Feb 14 −0.001 (0.031) / Aug 15 +0.346 (0.045) | See appendix |
| Inland deviance, birthplace vs SCI, mid-Aug (Check 5B) | 0.651 vs 0.428 | — | 0.630 vs 0.416 | — |
| Small-sample p-values (p82, "p = 1e-35") | — | — | Dropped; WCR bootstrap only | Audit M6 |

## 7. Appendix: old vs recent migration and Meta SCI
Province level, 47 × 47, PPML with origin province, destination province and distance-bin fixed effects, CRV1 by origin province (47 clusters). File: `r2_echo_sci_nights.csv`, summaries in `r2_echo_sci_summary.csv`. These regressors exist only between provinces, so in the city design they would be absorbed by the pair effects. The numbers below are for the main nights, inland destination provinces.

| Night | Deviance explained: birth 2024 | Deviance explained: SCI | Old (1988–95) vs recent (2021–24), joint | Birth 1998 given recent in + out |
|---|---:|---:|---|---|
| Feb 14 | 0.356 | 0.481 | old 0.100 (0.040), recent 0.744 (0.053) | −0.001 (0.031) |
| Aug 15 | **0.630** | 0.416 | old **0.532** (0.058), recent 0.349 (0.069) | **0.346** (0.045) |
| Nov 1 | 0.469 | 0.493 | old 0.258 (0.042), recent 0.582 (0.060) | 0.129 (0.034) |
| Dec 25 | 0.542 | 0.655 | old 0.231 (0.039), recent 0.693 (0.052) | 0.109 (0.032) |

Legacy Check 5 holds on single nights for inland destinations. Old ties and the 1998 birthplace stock matter in August, while recent migration and SCI dominate ordinary nights and Christmas. With all destinations included, summer flows are explained less well by every relation (birth 0.298, SCI 0.245).

## 8. Most defensible headline sentence
> *On every night, people sleep more in the provinces where their neighbours were born: comparing cities of the same province, a doubled share of natives of a province means about 1.4 times as many nights there on an ordinary February night and 1.55 times on 15 August, beyond distance and income. Holidays add a modest extra pull (Christmas Eve +0.07, p = 0.02; summer +0.10, p = 0.18 against February, p = 0.006 against October), concentrated in inland destinations (+0.23 in the inland provinces, p = 0.01) and absent at the coast. Across provinces, beaches take so many summer nights that August is no more roots-aligned than a weekday.*

Short Spanish version for the report: «Las raíces pesan todas las noches; las vacaciones, sobre todo las familiares y el verano en el interior, las refuerzan un poco; en verano, la playa domina el total.»

## Files
`r2_city_set.csv`, `r2_roots_city_prov.csv`, `r2_roots_prov_prov.csv`, `r2_solver_check.csv`, `r2_ppml_nights.csv`, `r2_ppml_summary.csv`, `r2_diff_tests.csv`, `r2_diff_summary.csv`, `r2_diff_jackknife.csv`, `r2_gravity_convergence.csv`, `r2_qap.csv`, `r2_qap_summary.csv`, `r2_legacy_reproduction.csv`, `r2_echo_sci_nights.csv`, `r2_echo_sci_summary.csv`, `r2_index_bridge.csv`, `r2_index_diff_tests.csv`, `r2_index_diff_summary.csv`, `r2_run.log`.
