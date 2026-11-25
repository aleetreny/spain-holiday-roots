# R1 (node level): do emptied villages gain a place in the holiday network?

*Agent R1, 2026-11-24. Code: [pipeline/r1_villages.py](/Users/alejandrotreny/Proyectos/spain-holiday-roots/pipeline/r1_villages.py), run with `.venv/bin/python -I pipeline/r1_villages.py` (log `results/r1_diag/run.log`). Every number below comes from the `results/r1_*.csv` file named next to it. Decisions are in `decisions.md`, in the lines tagged [R1].*

## Short answer
In network terms, no. Once size, zone aggregation, distance to a city and the coast are held fixed, emptied zones do not move towards the core on holidays. Their undirected coreness falls behind comparable zones, and degree-preserving rewiring reproduces that drop exactly. What does survive is a change in their role:
- Emptied zones flip from net exporters to net importers of people sleeping away from home.
- They close the visitor deficit they show on ordinary nights.
- They draw in-partners from more provinces than a gravity null predicts.

In summer the flip comes mostly from residents staying home, not from more visitors per resident. Visitors per resident beyond what the controls predict appear on family holidays: clearly at Christmas, weakly on All Saints. They also appear in the INE municipal data, in small municipalities.

## Definitions (fixed before estimation)
- **Loss** = `loss_peak_2021` of the node table = 1 − pop₂₀₂₁ / max(pop₁₉₀₀ … pop₁₉₆₀). The census columns are summed over the zone's INE municipalities (b00).
  - 56 of 2,435 zones have no loss and are excluded (54 of them single-municipality, mostly post-1960 municipalities).
  - The 3 `ARTEFACT_ZONES` are excluded too.
  - Samples: **all** n = 2,376; **single-municipality** n = 1,473; 47 provinces in each (`r1_sample.csv`).
  - The regressor is the percentile rank of loss, so a coefficient is the change from the least to the most emptied zone.
- **Node metrics, per night:**
  - s_in: visitors from other zones on the full matrix.
  - k_in and k_out: on the w ≥ 10 backbone; equal-edges as a check.
  - Coreness: b02, undirected projection.
  - In-coreness: igraph `mode="in"` on the directed w10 backbone.
  - n_prov: origin provinces other than the zone's own with ≥ 10 persons.
  - eff_prov: inverse Simpson index over origin provinces.
  - share_far: share of visitors from more than 100 km away.
  - away_share = s_out / pop.
  - net_per_res = (s_in − s_out) / pop; net_bal = (s_in − s_out) / (s_in + s_out); net importer = 1[s_in > s_out].
- **Controls:** log pop, log n_municipalities, log area, `coastal_true`, log(1 + km to the nearest zone with ≥ 100k people; 58 such zones), log(1 + km to the nearest coastal zone), and province fixed effects.
- **Inference:** CRV1 by province with G = 47. 95% CIs are wild-cluster percentile-t (Webb weights, B = 9,999), and p-values are restricted wild-cluster (p_WCR). My OLS matches pyfixest to 1e-8 in the coefficient (`r1_diag/check_fe_ols_vs_pyfixest.csv`).
- **Nights:**
  - Main pairs: 15 Aug − 14 Feb and 1 Nov − 14 Feb. The file date is the night *ending* on that date, so these are holiday-eve nights.
  - Sister ranges use every (holiday sister × ordinary sister) pair: summer 9 × 13 = 117, All Saints (1, 2, 3 Nov) 3 × 13 = 39.
  - Years: 2022 and 2023, each paired with that year's February Wednesday.
  - Placebo: each ordinary sister against 14 Feb (12 pairs).
  - Appendix: Christmas Eve night (25 Dec) and Good Friday (29 Mar), each against 14 Feb.
- **Change versus level.** Emptied zones have a low ordinary-night baseline, so a change score can show a "gain" that is only catch-up (the C2 problem, mirrored). I therefore also report the level on each night and an ANCOVA (holiday level given the ordinary level): `r1_models_levels_ancova_summary.csv`.

## Survive / not survive

Notation: β [95% wild-cluster CI], p_WCR. "Sisters" = how many sister pairs share the sign of the main estimate, and their range. Years = 2022 / 2023 coefficients. From `r1_models_summary.csv` unless another file is named.

| # | Claim | All zones (n = 2,376) | Single-municipality (n = 1,473) | Sisters, years | Verdict |
|---|---|---|---|---|---|
| 1 | Undirected k-core promotion (legacy "×4.5") | Aug Δcore **−9.4** [−13.6, −5.3]; Nov −2.4 [−4.4, −0.3] | −15.3 [−19.8, −10.7] | 117/117 negative (−12.3…−7.1); −11.9 / −11.0 | **Does not survive. Reversed.** |
| 1b | …beyond the degree-preserving rewired null | Aug excess −0.5 [−1.7, 0.8], p 0.44; rewired null alone predicts −8.9 | −0.3 [−1.6, 1.0], p 0.69 | main nights only | Coreness here is a function of degree (`r1_rewired_coreness_*`) |
| 2 | In-coreness gain (core built from visits only) | Aug +3.0 [0.8, 5.3], p 0.006; Nov +1.0 [−1.4, 3.3] | Aug +1.3 [−1.1, 3.8] | 117/117 positive (1.8…4.6); 3.0 / 3.2 | Survives in all zones only. **Not robust to aggregation.** |
| 2b | …beyond gravity (NB draws) | Aug +2.1 [0.9, 3.3], p 0.003; Nov +1.4 [0.4, 2.4] | Aug +2.7 [1.3, 4.2]; Nov +2.5 [1.1, 4.0] | main nights only (100 draws) | Survives the gravity null (`r1_gravity_summary.csv`) |
| 3 | Proportional visitor gain, log(s_in,h / s_in,o) | Aug **+0.46** [0.20, 0.73] (×1.58); Nov +0.54 [0.33, 0.75] | Aug +0.14 [−0.15, 0.43]; Nov +0.32 [0.09, 0.55] | Aug 117/117 (0.25…0.57); 0.47 / 0.53 | Survives in all zones and on All Saints. **Fails in single-municipality zones in summer.** |
| 3b | …but it is catch-up: level on each night | Feb −0.51 [−0.71, −0.31]; Aug −0.05 [−0.40, 0.30]; Nov +0.02 [−0.22, 0.26] | Feb −0.41; Aug −0.27 [−0.69, 0.14] | all sisters agree in sign on Feb | No excess level on holidays (`r1_models_levels_ancova_summary.csv`) |
| 3c | ANCOVA s_in (holiday given Feb) | Aug +0.28 [0.03, 0.53]; Nov +0.32 [0.13, 0.50] | Aug +0.02 [−0.27, 0.33]; Nov +0.16 [−0.05, 0.36] | Aug 117/117 positive | Weak; disappears in single-municipality zones |
| 4 | Proportional in-degree gain, log((k_in,h+1)/(k_in,o+1)) | Aug **+0.49** [0.27, 0.71]; Nov +0.38 [0.19, 0.57] | Aug +0.33 [0.05, 0.61]; Nov +0.29 [0.01, 0.55] | 117/117 (0.31…0.58); 0.40 / 0.53; equal-edges 0.48 | Survives. Absolute Δk_in is **−36** [−65, −5], so this is proportional only. |
| 4b | …beyond gravity (analytic NB) | Aug +0.11 [−0.08, 0.30]; Nov −0.06 [−0.21, 0.09] | Aug **+0.29** [0.11, 0.48], p 0.008; Nov +0.09 [−0.05, 0.23] | Aug single: 117/117 positive (0.05…0.36); 0.13 / 0.16 | In all zones gravity explains it. In single-municipality zones a summer excess survives. |
| 5 | Origin diversity beyond gravity, log eff_prov | Aug +0.19 [−0.02, 0.41], p 0.08 | Aug **+0.26** [0.02, 0.49], p 0.017 | 117/117 positive (0.20…0.39); 0.21 / 0.30 | Survives in summer (gravity predicts concentration, −0.25). Not on All Saints. |
| 6 | Net importer flip (level, probability) | Aug **+0.39** [0.23, 0.55]; Nov +0.33 [0.12, 0.54]; Feb −0.12 [−0.29, 0.04] | Aug +0.44 [0.21, 0.67]; Nov +0.31 [0.10, 0.53] | 9/9 and 3/3 sisters; 0.36 / 0.37 | **Survives** (`r1_models_levels_ancova_summary.csv`) |
| 6b | Net balance (in−out)/(in+out), change | Aug +0.42 [0.23, 0.61]; Nov +0.38 [0.25, 0.51] | Aug +0.24 [0.02, 0.47]; Nov +0.24 [0.09, 0.39] | 117/117 and 39/39; 0.39 / 0.44 | **Survives** |
| 6c | …but per resident, net (in−out)/pop | Aug level 0.00 [−0.20, 0.20]; Δ +0.03 [−0.16, 0.22] | Aug Δ −0.12 [−0.33, 0.09] | Aug sisters 78/117 positive | **In summer the flip is not more visitors per resident.** It comes from residents staying home: away_share level −0.72 [−0.86, −0.59] in log, about 51% fewer residents away. |
| 6d | Net per resident on family holidays | Nov Δ +0.07 [0.03, 0.12]; Christmas Δ **+0.19** [0.14, 0.24]; Good Friday +0.10 [−0.04, 0.24] | Nov +0.03 [−0.01, 0.08]; Christmas +0.15 [0.10, 0.19] | Nov 39/39; Christmas is a single night | Survives at Christmas in both samples, on All Saints in all zones only |
| 7 | Holiday-specific or "more people away"? (volume curve) | Gradient of log share of away nights: Aug −7.0 sd (sisters −7.7…−6.8); Nov +0.35 sd; Christmas −1.5 sd | Aug −5.2 sd; Nov +0.12 sd | all August nights extrapolated (100%) | August does **not** add a village effect. It sits below the curve: the beach takes the extra volume (`r1_volume_curve_summary.csv`). |
| 8 | INE 53001, municipalities, log(Aug/Feb tourists) | **+0.55** [0.43, 0.68], n = 5,831; ANCOVA +0.84 [0.69, 0.98] | Municipalities that are single-municipality zones: +0.41 [−0.04, 0.87]; ANCOVA +0.56 [0.05, 1.07] | 2022 0.46, 2023 0.67, 2025 0.27 | **Survives at municipal level.** It is concentrated in small municipalities: those inside '_AM' zones +0.68 [0.52, 0.83]; pop < 1,000 +0.68. |
| 8b | INE net (53001 arrivals − 53000 departures) per resident, Aug − Feb | **+2.85** [2.01, 3.73], n = 4,672 | single-municipality zones −0.00 [−0.53, 0.52]; in '_AM' zones +2.60 [1.66, 3.55] | 2022 3.26, 2023 3.57 | Survives, in small municipalities only (`r1_ine53001_models.csv`) |

Placebo ranges (12 ordinary-vs-ordinary pairs) are small for every outcome: for example log s_in [−0.03, 0.09], k_in [−0.04, 0.08], Δcore [−1.0, 0.5], net_bal [−0.04, 0.07].

## The steps in detail
1. **Coreness (task 3, second null).** For each main night, mean rewired coreness is an almost exact function of degree: R² on degree 0.9997 (Feb 14), 0.9999 (Aug 15) and 0.9998 (Nov 1). Real coreness is 0.96–0.99 (`r1_rewired_coreness_summary.csv`). So in these networks, coreness under rewiring is a function of degree, and a coreness result is a degree result.
   - The undirected projection mixes visits with residents' departures. Emptied zones send far fewer residents away in August: log(1 + k_out) level −0.79 [−1.03, −0.55] below comparable zones (`r1_models_levels_ancova_summary.csv`).
   - That is why their coreness lags even though their in-degree catches up.
   - Unadjusted quintile medians of Δcore (`r1_quintiles.csv`) are U-shaped: Q1 +42, Q3 +23, Q5 +34. Q5 has a median of 5 municipalities per zone and is 82% '_AM'. In single-municipality zones the medians fall monotonically: Q1 +45, Q5 +17.
   - The legacy "×4.5" came from aggregation and the low February baseline.
2. **Gravity null (task 3).** The fit is 3-way IPF: origin totals, destination totals, and 21 fixed distance bins up to 1,200 km. It was fitted separately on each of the 35 nights and converged in 32–45 iterations with relative margin error < 1e-8 (`r1_gravity_fit.csv`).
   - It preserves s_in and s_out exactly, so it cannot test visitor volume or net balance. It reproduces them by construction.
   - Binary metrics come from negative-binomial draws. Dispersion α is estimated by moments, 0.25–0.46.
   - The null under-predicts backbone edges: 25.4k against 35.8k on Feb 14, and 105.8k against 116.0k on Aug 15. No α closes the gap (the expected edge count is flat in α), so I compare only cross-sectional slopes.
   - Gravity predicts most of the all-zone proportional in-degree gain: expected +0.37 against observed +0.49. In single-municipality zones it predicts almost none (+0.04), so an excess of +0.29 remains.
3. **Net balance (pattern #6).** I reproduce the descriptive numbers (`r1_net_balance_by_quintile.csv`):
   - Emptied inland zones (loss ≥ 0.5, n = 455) that are net importers: **99.1%** on Aug 15, 94.3% on Nov 1, 22.6% on Feb 14.
   - In single-municipality zones (n = 76): 94.7%, 84.2% and 38.2%.
   - The gradient survives all controls (rows 6 and 6b).
   - The mechanism differs by holiday. In August, comparable zones receive the same visitors per resident while emptied zones' residents stay home (6c). At Christmas, and weakly on All Saints, visitors per resident do rise (6d).
4. **Volume curve (pattern #3).** For each of the 306 nights, the loss gradient comes from the same controlled regression. Inputs are `p_zone_daily.parquet`, asserted equal to my s_in and s_out on every overlapping night. The gradients are fitted on log(away persons) over 118 non-holiday 2024 nights (R² 0.83 for the log-share gradient).
   - Family holidays sit on or above the curve. On net per resident: All Saints +3.4 sd and Christmas +14.5 sd.
   - August sits far below it: −7.0 sd on log share and −9.0 sd on net per resident.
   - Every August night is an extrapolation beyond the largest ordinary volume, so the sd values are indicative, but the sign is unambiguous. The emptied-village share of away nights is not an August effect.
5. **Years (pattern #12).** Year comparisons use within-year ratios or shares. A log ratio of s_in equals a log ratio of shares up to a constant that the fixed effects absorb. Net balance is scale-free. For binary metrics, equal-edges gives the same signs (k_in 0.38 / 0.56).
6. **INE (task 4).** At true municipality level the gradient is larger than in MITMA zones, and it sits where MITMA cannot see. The emptied-quintile median municipality has 68 residents (`r1_ine53001_quintiles.csv`). MITMA's smallest zone has 1,504. Limitations:
   - INE counts tourists per month (domestic, other province or distance threshold), not nights.
   - Cells < 30 are suppressed. The per-resident results are bounded [0, 29]: 2.72–2.81 for Aug − Feb.
   - Coast and distance controls come from the parent MITMA zone.
   - There is no municipal area control.

## Not done or open
- No directed rewiring null for in-coreness, so row 2 is tested only against gravity.
- Pattern #7 (roots index) belongs to R2.
- Christmas and Good Friday are single nights with no sister ranges (appendix).

## Most defensible sentence
> *On holiday nights Spain's emptied villages switch from net exporters to net importers of people sleeping away from home (99% of emptied inland zones on 15 August, against 23% on an ordinary February night). This holds with size, zone aggregation, distance to cities, coast and province fixed, and in single-municipality zones. It is a change of role, not a rise to the core: in summer their visitors only catch up with similar places, and their residents stay home. More visitors per resident than similar places, beyond what those factors predict, appear only on family holidays (Christmas, All Saints) and in small municipalities.*
