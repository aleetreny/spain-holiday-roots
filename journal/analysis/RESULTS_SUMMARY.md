# Results summary (analysis phase closed, 27 November 2026)

One page for writing and figures. Every number comes from the report named in brackets; those reports trace each number to a CSV.
Reports: [A](A_calendar/REPORT_A.md) (2022–2024 roots calendar), [A2](A2_confirmation2025/REPORT_A2.md) (2025 hold-out), [B and D](REPORT_BD.md) (cohorts, villages), [C](C_municipal/REPORT_C.md) (parents' villages). Rules: [pre-registration](PREREGISTRATION.md) (two deviations logged).

## Data
- MITMA overnight stays: 1,079 valid nights 2022–2024 (17 dates missing or excluded by pre-registered rules) + 196 nights of 2025 (5 unseen August nights and 62 unseen ordinary nights in the test). Coursework numbers reproduced exactly (43 values, max difference 8e-15) [A, A2].
- Method change of 1 July 2025: on 6 dates computed both ways, the roots elasticity moves by 0.014 on average (tolerance 0.035) and the pre-1981 term by 0.005 (tolerance 0.024): **bridge passes** [A2].

## Pre-registered tests

August holiday nights (files 10–18 Aug) against ordinary weekday nights of the same year. "2022–24" = new ordinary nights against August nights already seen (Deviation 1). "2025" = unseen on both sides. Estimates with 95 % wild-cluster bootstrap intervals; Holm-adjusted p.

| | 2022–24 | 2025 hold-out | Abstract? |
|---|---|---|---|
| P1 pooled August premium | +0.106 [0.010, 0.197], p = 0.066 | **+0.134 [0.037, 0.233], p = 0.030** | **Yes** (via 2025, same sign) |
| P2 pre-1981 arrivals × August | **+0.071 [0.041, 0.092], p = 0.0003** | **+0.057 [0.022, 0.076], p = 0.010** | **Yes** (both) |
| P3 inland minus coast | +0.122 [−0.037, 0.280], p = 0.11 | +0.129 [−0.055, 0.319], p = 0.14 | No |
| S1 Christmas Eve | **+0.088 [0.049, 0.177], p = 0.0006** | no 2025 night | Secondary result |
| S2 Good Friday | +0.037 [−0.041, 0.110], p = 0.25 | +0.076 [−0.013, 0.157], p = 0.15 | No |
| S3 All Saints | **+0.073 [0.022, 0.135], p = 0.017** | no 2025 night | Secondary result |
| S4 born in a village ≤ 2,000 × August | **+0.073 [0.055, 0.095], p = 0.0006** | **+0.056 [0.031, 0.078], p = 0.004** | **Yes** |
| S5 arrived 2016–2023 × August | **−0.044 [−0.070, −0.011], p = 0.017** | −0.033 [−0.060, 0.006], p = 0.15 | Yes (2022–24; same sign 2025) |
| S6 pre-1981 minus recent | **+0.115 [0.071, 0.139], p = 0.0006** | **+0.089 [0.034, 0.117], p = 0.004** | **Yes** |

Placebo windows: in every year 2022–2025 the August contrast is above every block of 9 consecutive ordinary nights (P1–P3) [A, A2]. Census 2025 instead of 2024: P1 2025 = +0.131, p = 0.033 [A2].

**Reading.** On holiday, cities send relatively more overnight visitors to the provinces where more of their residents were born, beyond distance, income and the destination pair. The premium is carried by the exodus generation (arrived before 1981) and by people born in villages; recent arrivals move the other way in August. The premium shows on family holidays without beach (Christmas Eve, All Saints) but not detectably at Easter.

## Robustness and descriptive facts (not tests)
- P2 holds in every robustness check except coastal destinations only (−0.005) and the smallest zero-cell offset; its **size** depends on the zero-cell rule (+0.014 to +0.122): report the sign as robust, not the size [A].
- P1 is between +0.07 and +0.15 across checks; intervals touch zero with road time (p 0.069) and equal night weighting (p 0.060) [A].
- The August premium sits in **inland destinations** (+0.14 to +0.15) and not on the coast (+0.01), but the formal inland × coast test (P3) is not significant and vanishes without Madrid and Barcelona origins [A]. Say "concentrated inland", not "larger inland".
- Calendar: mean pooled elasticity 0.521 on ordinary Mon–Wed nights, 0.635 on 10–18 August, 0.645 in late August, 0.617 on Christmas Eve, 0.648 on Epiphany eve. The highest nights are the end of August and early-September Sunday evenings; peak volumes are 10–18 August [A]. The pre-1981 term peaks in August (0.178 vs 0.104 ordinary); the recent-arrival term peaks at Christmas, New Year and Epiphany (≈ 0.20–0.22) [A].

## Cohorts on the coursework nights (B)
- Against 14 Feb 2024 and against clean October 2024 weekdays: pre-1981 +0.084 / +0.087; inland +0.155 / +0.154; village-born +0.086 / +0.093. Direct inland × coast for pre-1981 on common 19-province support +0.114 / +0.128. Holds within age bands, without Madrid and Barcelona, in all province jackknives and in 117/117 sister pairs.
- Recent arrivals −0.071 vs February but −0.029 (n.s.) vs October: February has an unusually high recent-arrival slope. The pre-1981 minus recent gap is robust (+0.155 / +0.116).

## Receiving end (D)
- Most emptied fifth of zones (n = 475): net receivers 97–99 % on every August night of 2022–2024 (93/93 nights) and 99.6 % on Christmas Eve 2024. On ordinary nights the share ranges 24 % (Feb 2024) to 56 % (Feb 2022) — report the range.
- 14–15 Aug vs 13–14 Feb 2024, most emptied fifth: arrivals per resident +0.644, departures per resident −0.013, so almost all of the net gain is arrivals. Inland, with controls: arrivals +0.144 and fewer departures −0.076 (about two thirds arrivals), similar to Christmas Eve.
- Second-home attenuation (48 % August vs 2 % Christmas Eve) is not detected in single-municipality zones (+0.083, p = 0.18): SI only, labelled attenuation.

## Parents' villages (C) — failed the gate, SI only
- A survey-reported parental tie raises the August − February probability that INE publishes the city-to-village flow by +3.98 pp [2.24, 5.72] (vs October +3.08 pp [1.52, 4.63]).
- The change in trip counts is not identified under the suppression rule (interval envelope −65 to +93 trips). One sentence in the main text, details in the SI.

## Not in the main text
Network-structure diagnostics (SI), BOSAM, k-core, territories, prediction gain (2.1 % over geography, 0.4–0.6 % over lagged flows), SCI comparison, festivals.

## Abstract-level claims (decision rule satisfied)
1. Roots elasticity rises on August holiday nights (2025 hold-out +0.134).
2. The rise is driven by pre-1981 arrivals and village-born residents (both periods), and recent arrivals move the other way (S5, S6).
3. Christmas Eve and All Saints also show a premium (secondary, 2022–2024).
4. The most emptied zones become net receivers on 97–99 % of August nights (descriptive).
