# Analyses B and D: Who Pulls (Arrival Cohort & Birthplace Size) and the Receiving End (Emptied Villages)

**Date:** 27 November 2026  
**Scope:** Discovery cache only (`data/derived/pernoct_mun_2024_6m.parquet` and `data/derived/pernoct_mun_2022_2023.parquet`, 306 cached nights across 2022–2024). Zero files from the held-out confirmation set (`_local/external_journal/mitma_pernoct/`) were opened or inspected.  
**Code & outputs (this worktree):**
- Analysis B scripts & CSVs: [`journal/analysis/B_cohorts/`](B_cohorts/) (`b_common.py`, `run_b_cohorts.py`, 13 CSVs in `out/`)
- Analysis D scripts & CSVs: [`journal/analysis/D_villages/`](D_villages/) (`run_d_villages.py`, `d_inland_decomp.py`, 9 CSVs in `out/`)
- Heavy intermediates (not in git): `_local/journal_data/BD/` (`b2_sister_pairs.parquet`, `d_zone_daily_panel.parquet`, run logs)

All dates follow the binding night convention: a file dated $d$ is the night ending on the morning of $d$ (evening $d-1$). Throughout, estimates on city–province flows describe aggregate province-level ties rather than individual returns to a specific village; non-significant effects are described as *not detected*; non-primary housing regressions are treated strictly as *slope attenuation* rather than mediation or proof of where visitors sleep; and no causal claims are made.

---

## 0. Executive Summary and Verdicts

| Check | Core question | Headline result (95% CI) | Change vs coursework (`V1` / `R1`) | Verdict |
|---|---|---|---|---|
| **B.1** Audit recode | Do V1 estimates change when dropping special `VARANOP='7777'` records (11,350) and 3 artefact destination zones? | Exact reproduction `< 1e-15`. Combined recode shifts `cohort_joint` by `<= 0.0003` and `born_village_le2k` by `-0.0007` (`all`) / `-0.0016` (`inland`). Dropping DANA 29–31 Oct weakens `recent` vs Oct from `p=0.060` to `p=0.081`. | Tiny numeric changes; confirms V1 `cohort_joint` already excluded `7777`, whereas `birthsize_joint` and `E2a` had absorbed `7777`. | **Confirmed** |
| **B.2** Baselines, 2022–2023 & stacked | Do `pre-1981`, `village-born`, `recent`, and `pre-1981 − recent` hold against both `Feb14` and `OctMean12` and across 2022–2024? | `pre1981`: `+0.0836 [+0.0459, +0.1175]` (Feb14), `+0.0870 [+0.0551, +0.1149]` (Oct12). `village_le2k`: `+0.0857 [+0.0614, +0.1084]` (Feb14), `+0.0927 [+0.0705, +0.1180]` (Oct12). `pre1981 − recent` (S6): `+0.1550` / `+0.1161` (both `p=0.0001`). `recent` alone: `-0.0714` (`p=0.0077`) vs Feb14, **not detected** vs Oct12 (`-0.0291`, `p=0.081`) or in 2022 (`-0.0223`, `p=0.322`). | Hardens `pre1981`, `village_le2k`, and `S6` across 100% of sister pairs (`117/117`, `108/108`) and 3-year stacked PPML; downgrades standalone `recent` to baseline/year-dependent (confirmed by Feb-vs-Oct placebo). | **Confirmed (`pre1981`, `village`, `S6`); Baseline-dependent (`recent`)** |
| **B.3** Direct `inland × coast` | In one 4-block stacked model on common 19-province support, is the August lift larger inland than on the coast? | `pre1981` (`common_19_provs`): **`+0.1144 [+0.0631, +0.1630]`** (`p=0.0002`, Feb14) and **`+0.1276 [+0.0846, +0.1650]`** (`p=0.0001`, Oct12). `born_village_le2k` (`common_19_provs`): **`+0.0496 [-0.0048, +0.1079]`** (`p=0.075`, Feb14) and **`+0.0345 [-0.0058, +0.0843]`** (`p=0.092`, Oct12), vs `+0.1036` / `+0.0879` (`p <= 0.0002`) on full 47 inland vs 19 coast. | **Major nuance:** `pre1981` inland dominance is rock-solid even inside the 19 coastal provinces; `village_le2k` inland dominance is largely driven by the 28 landlocked provinces (`+0.164` / `+0.169`) while coastal zones also see a modest village-born lift (`+0.041`, `p=0.0086` vs Oct12). | **Confirmed (`pre1981`); Nuanced (`village_le2k`)** |
| **B.4** Age-matched contrasts | Is the cohort or birthplace effect merely an age effect? | Within `65+`: `pre1981 − post1981` = `+0.0947 [+0.0484, +0.1430]` (`p=0.0003`, Feb14) and `+0.0565 [+0.0285, +0.0851]` (`p=0.0004`, Oct12); inland `+0.1187` / `+0.0914`. Within `43–64`: `pre1981 − recent` = `+0.0395 [+0.0119, +0.0710]` (`p=0.0086`) and `+0.0339 [+0.0128, +0.0565]` (`p=0.0050`). Within both `43–64` and `65+`, `village − urban` = `+0.058..+0.084` (`p <= 0.0016`). | Audit drop of 3,245 `7777` records aged 65+ slightly *strengthens* the `65+` within-band gap (`+0.0926 -> +0.0965`). Adds direct within-band WCR intervals and confirms pattern in `43–64` and `43–59`/`60+`. | **Confirmed** |
| **B.5** Zero-cell rules (common scale) | Are effect sizes comparable across the 7 Table A2 zero-cell rules on a common predicted-flow contrast? | For `+1 within-FE SD` of raw share (`Aug15 vs Feb14`, `all`): `pre1981` is `+4.1%` to `+10.3%` across 6 of 7 rules (`+7.9% [+4.2%, +11.3%]` main) and `+8.2%` to `+14.9%` inland; `village_le2k` is `+8.8%` to `+17.4%` (`+13.6%` main). `log(share+1e-5)` shrinks sizes 3–4× (`+2.0%`, `p=0.091` `all`, but `+5.2%`, `p=0.0003` inland; `p=0.035` vs Oct12). | Common flow scale removes most of the >10× raw coefficient spread across `asinh`, `+1.0`, `+0.5`, `+0.1`, `+1e-4`, and `112 cities`, while flagging `1e-5` as an extreme pseudo-zero leverage outlier. 3-term partitions confirm a clean monotonic gradient. | **Confirmed (with `1e-5` sensitivity note)** |
| **B.6** Jackknives & drop Madrid/Barcelona | Does any province or Madrid/Barcelona drive the cohort results? | 44/44 origin provinces and 47/47 destination provinces keep `p <= 0.0007` (`pre1981`) and `p <= 0.0001` (`village_le2k`) on both baselines. Dropping **both** Madrid & Barcelona (`G=42`) leaves `pre1981`, `village_le2k`, `S6`, and `common_19_provs` `pre1981` inland×coast (`+0.1010` / `+0.1137`, `p <= 0.0005`) intact. | Adds full 44-origin and 47-destination jackknives plus joint Madrid+Barcelona exclusion with WCR CIs. | **Confirmed** |
| **D.1** Receiving end across 306 nights | Do the most emptied zones (Q5) turn into net receivers across all cached August nights and baselines? | All **93/93** cached August nights (2022–2024) have Q5 receiver share `96.6%–99.0%` (`86.9%–97.6%` in `single_mun`). `2024-08-15`: `98.95% [97.98%, 99.78%]`; `2024-12-25`: `99.58% [98.91%, 100.00%]`. | **Baseline nuance:** Q5 ordinary receiver share is `24.0%` on 14 Feb 2024 (`25.5%` on 15 Feb 2023) but `47.4%` on `2024-10-mean12` and `55.6%` on 16 Feb 2022. August reversal is massive from every baseline (`+51.6` to `+74.9 pp`). | **Confirmed** |
| **D.2** Arrivals vs departures decomposition | Does August work via fewer departures or more arrivals (compared with Christmas Eve)? | **Raw Q5** (`Aug15 − Feb14`): arrivals `+0.644 [+0.607, +0.683]`, departures `-0.013 [-0.017, -0.009]` (**98%** of Q5 net gain from arrivals; `Q5 − Q1` = **76.4%** arrivals, **23.6%** fewer departures). **Controlled `all`** (`n=2,376`): arrivals `-0.040` (n.s.), departures `-0.069` (`p=0.0001`). **Controlled `inland`** (`n=2,053`): arrivals **`+0.144 [+0.058, +0.224]`** (65.5%), departures **`-0.076 [-0.090, -0.060]`** (34.5%), net **`+0.221 [+0.130, +0.306]`** — matching Christmas Eve (`+0.200`: `+0.120` arrivals, `-0.080` departures). | **Major refinement vs coursework R1:** the coursework's all-zone null on August linear `arrivals_per_res` and `net_per_res` was caused by extreme-tail coastal resorts in Q1–Q2 swamping linear slopes. Among inland zones, August operates via **both** higher arrivals and fewer departures, just like Christmas Eve. | **Refined (major clarification)** |
| **D.3** Housing slope attenuation | Does non-primary housing (`share_np`) attenuate the emptying slope more in August than on Christmas Eve? | `all` (`n=2,376`, `log_s_in_ratio`): Aug attenuation `0.221` (48.1%) vs Xmas `0.012` (1.8%); paired **Aug − Xmas diff = `+0.209`** (zone boot `[+0.120, +0.292]`, prov boot `[+0.061, +0.361]`, `p <= 0.005`). `inland` (`n=2,053`): raw slopes identical (`0.713` vs `0.713`), diff **`+0.271`** (`[+0.198, +0.343]` zone / `[+0.161, +0.398]` prov). `single_mun` (`n=1,473`): diff **`+0.083`** (`[-0.039, +0.184]` zone, `p=0.179`, **not detected**). | Confirmed on full sample, inland sample, both baselines, and `log1p_np_per_res` (`+0.180`), with paired bootstrap CIs. **Caveat:** does not reach significance in `single_mun` alone because the largest non-primary housing shares sit in multi-municipality `_AM` zones. | **Confirmed (`all`, `inland`); Not detected (`single_mun`)** |

---

## 1. Analysis B — Who Pulls: Arrival Cohort and Birthplace Size

### 1.1 Sample, Specification, and Estimation Mechanics
- **Dyad frame:** 675 mainland municipalities $> 10,000$ residents that form their own MITMA zone, across $G = 44$ origin provinces with $\ge 2$ such cities (`pipeline/v1_verify_n1.py` imported read-only), paired with the 46 other mainland destination provinces ($675 \times 46 = 31,050$ directed city–province dyads per night; $12,372$ coastal-province dyads when restricted to the 19 coastal provinces).
- **Destination flow construction:** full MITMA *pernoctaciones* weights (`w >= 0`, no $w \ge 10$ threshold), excluding self-province flows and excluding the three binding artefact destination zones (`31905_AM`, `19103_AM`, `40149_AM`; pre-registration §1). Coastal vs inland destinations follow `coastal_true` from `nodes.parquet` (never the buffer flag).
- **PPML & inference:** exact block-diagonal stacked HDFE PPML (`FastStackedMultiBlock` in `b_common.py`) interacting city FE, origin-province × destination-province FE, log-distance + (log-distance)², and 46 income × destination-province slopes (ADRH 2024) with every block (night and/or destination stratum). Verified against `pyfixest.fepois` on both 2-night differences and 4-block `inland × coast` models (`|diff| < 1e-5`; `run_b_all.log`). All $p$-values and 95% confidence intervals use Kline–Santos (2012) wild cluster restricted (WCR) score bootstrap with 6-point Webb weights and 9,999 draws, clustered by origin province ($G = 44$, or $43$/$42$ in leave-out checks) jointly across all stacked blocks.
- **Baselines:** two binding ordinary references for 2024: (i) **13–14 February 2024** (`2024-02-14`, coursework main ordinary night) and (ii) **clean October 2024 weekday mean** (`2024-10-mean12`, mean of the 12 Mon/Tue/Wed-evening weekday files in October 2024 after dropping the 3 DANA dates `2024-10-29`, `2024-10-30`, `2024-10-31`). We also report `2024-02-mean13` (13 February weekday nights), `2024-10-mean10_noVal` (further dropping Valencia Day evenings 8–9 Oct, files `2024-10-09` and `2024-10-10`), `2024-08-mean9` (mean of 10–18 August), and full multi-night stacked PPML across 2022–2024.

### 1.2 B.1 — Audit Recode vs Coursework (`b1_v1_audit_comparison.csv`)
In the 2024 Census 10% microdata (`md_CensoAnual_muestra_2024.parquet`), there are **592,961** interprovincial mainland-born resident records across the 675 cities. Exactly **11,350** of them (`1.91%`) carry the special arrival-year code `VARANOP = '7777'` (and zero other non-dated codes appear), leaving **581,611** dated arrivals (`0001` pre-1981: 243,858; `0002`–`0005` 1981–2015: 234,441; `2016`–`2023` recent: 103,312). Among the 11,350 `7777` records:
- **1,266** (`1.03%` of raw village-born) have `CMUN_NAC = '991'` ($\le 2,000$),
- **8,856** (`2.39%` of raw urban-born) have `CMUN_NAC` $> 10,000$, and
- **3,245** (`4.42%` of V1's raw `post1981_ge65`) are aged 65+.

Per `PREREGISTRATION.md` §1, only explicit dated arrival codes enter any group; `7777` is dropped from both `cohort_joint` and `birthsize_joint` (and from all age-matched cells in B.4). Table 1 decomposes the exact shift against `results/v1_c2_main_contrasts.csv` (source: `out/b1_v1_audit_comparison.csv`).

| Variant | Sample | `pre1981` ($\Delta$ vs V1) | `recent_2016_2023` ($\Delta$ vs V1) | `born_village_le2k` ($\Delta$ vs V1) | `born_urban_gt10k` ($\Delta$ vs V1) |
|---|---|---:|---:|---:|---:|
| `v1_exact_repro` (Feb14) | `all` | `+0.08365` (`4e-16`) | `-0.07156` (`-1e-17`) | `+0.08646` (`-3e-16`) | `-0.06156` (`3e-16`) |
| `drop_7777_only` (Feb14) | `all` | `+0.08365` (`4e-16`) | `-0.07156` (`-1e-17`) | `+0.08580` (`-0.00066`) | `-0.06124` (`+0.00032`) |
| `drop_artefact_dest_only` (Feb14) | `all` | `+0.08357` (`-0.00008`) | `-0.07142` (`+0.00014`) | `+0.08639` (`-0.00007`) | `-0.06126` (`+0.00030`) |
| **`prereg_binding_both` (Feb14)** | `all` | **`+0.08357`** (`-0.00008`, `p=0.0002`) | **`-0.07142`** (`+0.00014`, `p=0.0077`) | **`+0.08574`** (`-0.00072`, `p=0.0001`) | **`-0.06091`** (`+0.00065`, `p=0.0226`) |
| **`prereg_binding_both` (Feb14)** | `inland` | **`+0.15445`** (`+0.00009`, `p=0.0001`) | **`-0.07346`** (`+0.00028`, `p=0.0096`) | **`+0.12190`** (`-0.00156`, `p=0.0001`) | **`-0.04120`** (`+0.00212`, `p=0.1787`) |
| **`prereg_binding_oct12_noDANA`** | `all` | **`+0.08696`** (`+0.00133`, `p=0.0001`) | **`-0.02911`** (`+0.00312`, `p=0.0810`) | **`+0.09266`** (`-0.00033`, `p=0.0001`) | **`-0.01126`** (`+0.00617`, `p=0.5901`) |
| **`prereg_binding_oct12_noDANA`** | `inland` | **`+0.15401`** (`+0.00107`, `p=0.0001`) | **`-0.03197`** (`+0.00225`, `p=0.0760`) | **`+0.12872`** (`-0.00190`, `p=0.0001`) | **`+0.00788`** (`+0.00675`, `p=0.6985`) |

*Verdict (B.1):* **Confirmed.** Coursework `cohort_joint` was completely immune to `7777` (exact zero shift because V1 selected `VARANOP = '0001'` and `2016..2023` explicitly). Dropping `7777` from `birthsize_joint` and excluding the 3 artefact destination zones moves point estimates by at most `0.0021`. Excluding the 3 DANA nights from the October reference leaves `pre1981` and `born_village_le2k` slightly sharper (`p = 0.0001`) while moving `recent_2016_2023` from `p = 0.060` (`-0.0322`) to `p = 0.081` (`-0.0291`).

### 1.3 B.2 — Both Baselines, Multi-Night Stacked PPML, Sister Grids, and 2022–2023 (`b2_main_contrasts.csv`, `b2_stacked_multinight.csv`, `b2_sister_summary.csv`, `b2_key_single_nights.csv`)

Single-night baseline levels (`b2_key_single_nights.csv`, `all`) explain why `Feb14` and `OctMean12` agree on `pre1981` and `born_village_le2k` but disagree on `recent_2016_2023` and `born_urban_gt10k`:
- `pre1981`: `0.0955` (14 Feb) and `0.0921` (Oct12) $\to$ `0.1791` (15 Aug); inland `0.1027` / `0.1031` $\to$ `0.2571`.
- `born_village_le2k`: `0.0358` (14 Feb) and `0.0288` (Oct12) $\to$ `0.1215` (15 Aug); inland `0.0466` / `0.0398` $\to$ `0.1685`.
- `recent_2016_2023`: `0.1961` (14 Feb) vs `0.1538` (Oct12) $\to$ `0.1247` (15 Aug). Ordinary February has an elevated `recent` (and `born_urban_gt10k`: `0.2806` vs `0.2310` in Oct12) slope, directly visible in the **ordinary-night placebo** (`Placebo_OctMean12_vs_FebMean13`, `b2_main_contrasts.csv`): `pre1981` `-0.0041` (`p=0.657`) and `born_village_le2k` `-0.0075` (`p=0.234`) show zero drift between October and February, whereas `recent_2016_2023` shifts `-0.0395 [-0.0588, -0.0207]` (`p=0.0002`) and `born_urban_gt10k` shifts `-0.0462 [-0.0759, -0.0188]` (`p=0.0008`).

| Contrast / Stacked Spec | Sample | `pre1981` (95% WCR CI; $p$) | `recent_2016_2023` (95% WCR CI; $p$) | `pre1981 − recent` (S6) (95% WCR CI; $p$) | `born_village_le2k` (95% WCR CI; $p$) |
|---|---|---|---|---|---|
| **2024 Aug15 vs Feb14** | `all` | `+0.0836 [+0.0459, +0.1175]` (`0.0002`) | `-0.0714 [-0.1138, -0.0201]` (`0.0077`) | `+0.1550 [+0.0933, +0.2019]` (`0.0001`) | `+0.0857 [+0.0614, +0.1084]` (`0.0001`) |
| **2024 Aug15 vs OctMean12** | `all` | `+0.0870 [+0.0551, +0.1149]` (`0.0001`) | `-0.0291 [-0.0547, +0.0038]` (`0.0810`) | `+0.1161 [+0.0725, +0.1470]` (`0.0001`) | `+0.0927 [+0.0705, +0.1180]` (`0.0001`) |
| 2024 Aug15 vs OctMean10_noVal | `all` | `+0.0878 [+0.0557, +0.1156]` (`0.0001`) | `-0.0300 [-0.0567, +0.0042]` (`0.0840`) | `+0.1178 [+0.0727, +0.1498]` (`0.0001`) | `+0.0942 [+0.0719, +0.1193]` (`0.0001`) |
| **2024 Aug15 vs Feb14** | `inland` | `+0.1545 [+0.1108, +0.1906]` (`0.0001`) | `-0.0735 [-0.1220, -0.0193]` (`0.0096`) | `+0.2279 [+0.1565, +0.2821]` (`0.0001`) | `+0.1219 [+0.0932, +0.1486]` (`0.0001`) |
| **2024 Aug15 vs OctMean12** | `inland` | `+0.1540 [+0.1148, +0.1875]` (`0.0001`) | `-0.0320 [-0.0589, +0.0039]` (`0.0760`) | `+0.1860 [+0.1339, +0.2213]` (`0.0001`) | `+0.1287 [+0.1022, +0.1576]` (`0.0001`) |
| 2024 Aug15 vs Feb14 | `coastal` | `-0.0173 [-0.0591, +0.0284]` (`0.4259`) | `-0.0459 [-0.0927, +0.0172]` (`0.1335`) | `+0.0287 [-0.0340, +0.0785]` (`0.3154`) | `+0.0183 [-0.0188, +0.0548]` (`0.3132`) |
| 2024 Aug15 vs OctMean12 | `coastal` | `+0.0030 [-0.0240, +0.0324]` (`0.8246`) | `-0.0184 [-0.0469, +0.0145]` (`0.2397`) | `+0.0214 [-0.0132, +0.0538]` (`0.1962`) | **`+0.0408 [+0.0107, +0.0744]`** (`0.0086`) |
| 2024 Aug15 vs Feb14 | `inlandprov` (28) | `+0.1851 [+0.1224, +0.2362]` (`0.0001`) | `-0.0559 [-0.1042, +0.0070]` (`0.0789`) | `+0.2411 [+0.1421, +0.3133]` (`0.0001`) | `+0.1637 [+0.1168, +0.2006]` (`0.0001`) |
| 2024 Aug15 vs OctMean12 | `inlandprov` (28) | `+0.1695 [+0.1123, +0.2182]` (`0.0001`) | `-0.0250 [-0.0521, +0.0115]` (`0.1698`) | `+0.1945 [+0.1218, +0.2472]` (`0.0001`) | `+0.1693 [+0.1262, +0.2040]` (`0.0001`) |
| **Stacked 2024 Aug9 vs Feb13** | `all` | `+0.0792 [+0.0427, +0.1104]` (`0.0001`) | `-0.0723 [-0.1099, -0.0292]` (`0.0008`) | `+0.1514 [+0.0957, +0.1907]` (`0.0001`) | `+0.0807 [+0.0580, +0.1034]` (`0.0001`) |
| **Stacked 2024 Aug9 vs Oct12** | `all` | `+0.0822 [+0.0504, +0.1082]` (`0.0001`) | `-0.0340 [-0.0580, -0.0034]` (`0.0292`) | `+0.1162 [+0.0729, +0.1446]` (`0.0001`) | `+0.0871 [+0.0655, +0.1108]` (`0.0001`) |
| Stacked 2024 Aug9 vs Oct12 | `inland` | `+0.1456 [+0.1071, +0.1760]` (`0.0001`) | `-0.0358 [-0.0615, -0.0005]` (`0.0471`) | `+0.1813 [+0.1290, +0.2149]` (`0.0001`) | `+0.1229 [+0.0968, +0.1502]` (`0.0001`) |
| **Stacked 2022 Aug9 vs Feb12** | `all` | `+0.0696 [+0.0344, +0.0955]` (`0.0006`) | `-0.0210 [-0.0484, +0.0121]` (`0.2722`) | `+0.0906 [+0.0352, +0.1251]` (`0.0017`) | `+0.0796 [+0.0565, +0.1062]` (`0.0001`) |
| **Stacked 2023 Aug9 vs Feb12** | `all` | `+0.0759 [+0.0400, +0.1032]` (`0.0003`) | `-0.0489 [-0.0827, -0.0107]` (`0.0082`) | `+0.1247 [+0.0716, +0.1608]` (`0.0001`) | `+0.0859 [+0.0587, +0.1060]` (`0.0001`) |
| 2024 Mar29 (Good Fri) vs Feb14 | `all` | `+0.0566 [+0.0263, +0.0839]` (`0.0009`) | `-0.0706 [-0.1080, -0.0225]` (`0.0029`) | `+0.1271 [+0.0676, +0.1711]` (`0.0002`) | `+0.0671 [+0.0471, +0.0846]` (`0.0001`) |
| 2024 Nov01 (All Saints) vs Feb14 | `all` | `+0.0375 [+0.0047, +0.0678]` (`0.0277`) | `-0.0322 [-0.0641, +0.0035]` (`0.0776`) | `+0.0698 [+0.0221, +0.1113]` (`0.0062`) | `+0.0310 [+0.0057, +0.0497]` (`0.0205`) |
| 2024 Dec25 (Xmas Eve) vs Feb14 | `all` | `+0.0258 [+0.0002, +0.0490]` (`0.0487`) | `+0.0010 [-0.0286, +0.0460]` (`0.9580`) | `+0.0248 [-0.0319, +0.0634]` (`0.3590`) | `+0.0285 [+0.0105, +0.0482]` (`0.0024`) |
| 2024 Dec25 (Xmas Eve) vs Oct12 | `all` | `+0.0291 [+0.0126, +0.0454]` (`0.0007`) | **`+0.0433 [+0.0226, +0.0769]`** (`0.0001`) | `-0.0141 [-0.0545, +0.0154]` (`0.3770`) | `+0.0354 [+0.0227, +0.0542]` (`0.0001`) |

**Sister-night grids (`b2_sister_summary.csv`, backed by 4,740 pairs in `b2_sister_pairs.parquet`):**
- `2024_aug9_x_feb13` (`117` pairs, `all` and `inland`): `pre1981` and `born_village_le2k` are positive and WCR-significant (`p < 0.05`) in **117/117** pairs (median `all`: `+0.0788` and `+0.0805`; `inland`: `+0.1498` and `+0.1168`, `max p = 0.0002`).
- `2024_aug9_x_oct12` (`108` pairs, `all` and `inland`): `pre1981` and `born_village_le2k` are significant in **108/108** pairs (median `all`: `+0.0834` and `+0.0872`, `max p = 0.0015`; `inland`: `+0.1457` and `+0.1239`, `max p = 0.0001`). By contrast, `recent_2016_2023` is significant in only `71/108` (`all`) and `46/108` (`inland`).
- `2022_aug9_x_feb12` (`108` pairs) and `2023_aug9_x_feb12` (`108` pairs): `pre1981` and `born_village_le2k` are significant in **108/108** pairs in both years (`all` and `inland`), whereas `recent_2016_2023` is significant in only `16/108` (2022 `all`) and `81/108` (2023 `all`).
- **Other holidays (`b2_main_contrasts.csv`):** Good Friday (`Mar29`) reproduces the summer profile at ~70% amplitude (`pre1981` `+0.0566` vs Feb14 / `+0.0599` vs Oct12; `village_le2k` `+0.0671` / `+0.0740`). Christmas Eve (`Dec25`) behaves **qualitatively differently** from August: against Oct12, **all** groups rise — including `recent_2016_2023` (`+0.0433`, `p=0.0001`) and `born_urban_gt10k` (`+0.0745`, `p=0.0001`) — reflecting universal family visits across cities and generations rather than an older-cohort village pull.

*Verdict (B.2):* **Confirmed for `pre1981`, `born_village_le2k`, and `pre1981 − recent` (S6); `recent` alone (S5) is baseline- and year-dependent.** Notice that when pooling all 9 August nights against 12 October nights in the pre-registered stacked estimator (`2024_Aug9_vs_Oct12_stacked`), `recent` reaches `p = 0.029` (`-0.0340 [-0.0580, -0.0034]`), but it stays null in 2022 (`-0.0210`, `p=0.272`) and on single-night `Aug15 vs OctMean12` (`p=0.081`). The paper should anchor on `pre1981`, `born_village_le2k`, and the contrast `pre1981 − recent`, rather than claiming a standalone negative `recent` drop in every year/baseline.

### 1.4 B.3 — Direct `inland × coast` Stacked Interaction (`b3_inland_coast_interaction.csv`)
Instead of comparing separate inland and coastal regressions, we fit one 4-block stacked PPML (`ordinary-inland`, `holiday-inland`, `ordinary-coast`, `holiday-coast`) with block-specific fixed effects and covariates, and directly test the triple-difference parameter $\gamma = \delta_{\text{inland}} - \delta_{\text{coast}}$ via joint 9,999-draw WCR clustered on origin province ($G=44$). We evaluate four supports:
1. **`common_19_provs` (binding strict common support):** the 19 coastal provinces that contain *both* inland and coastal destination zones (same city–province dyads in both strata, $n_{\text{4-block}} \approx 49,290$).
2. **`common_19_provs_exact_cell_intersect`:** further restricts to exact $(c, p)$ pairs active in all 4 blocks.
3. **`common_19_provs_stratum_centroid`:** replaces province-wide centroids with separate inland-zone and coastal-zone population-weighted centroids inside each of the 19 provinces.
4. **`full_47_inland_vs_19_coast`:** compares all 47 inland-zone destination provinces (including the 28 fully landlocked provinces) against the 19 coastal-zone provinces in one joint 4-block model ($n_{\text{4-block}} \approx 86,500$).

| Support | Contrast | Term | $\delta_{\text{inland}}$ | $\delta_{\text{coast}}$ | $\gamma = \delta_{\text{inland}} - \delta_{\text{coast}}$ | 95% WCR CI | $p_{\text{WCR}}$ |
|---|---|---|---:|---:|---:|---|---:|
| **`common_19_provs`** | **Aug15 vs Feb14** | **`pre1981`** | `+0.0971` | `-0.0173` | **`+0.1144`** | `[+0.0631, +0.1630]` | **`0.0002`** |
| **`common_19_provs`** | **Aug15 vs OctMean12** | **`pre1981`** | `+0.1306` | `+0.0030` | **`+0.1276`** | `[+0.0846, +0.1650]` | **`0.0001`** |
| `common_19_provs` | Aug15 vs Feb14 | `pre1981 − recent` | `+0.1739` | `+0.0287` | `+0.1452` | `[+0.0792, +0.2290]` | `0.0001` |
| `common_19_provs` | Aug15 vs OctMean12 | `pre1981 − recent` | `+0.1774` | `+0.0214` | `+0.1560` | `[+0.0948, +0.2067]` | `0.0001` |
| **`common_19_provs`** | **Aug15 vs Feb14** | **`born_village_le2k`** | `+0.0679` | `+0.0183` | **`+0.0496`** | `[-0.0048, +0.1079]` | **`0.0746`** |
| **`common_19_provs`** | **Aug15 vs OctMean12** | **`born_village_le2k`** | `+0.0753` | `+0.0408` | **`+0.0345`** | `[-0.0058, +0.0843]` | **`0.0924`** |
| `common_19_provs` | AugMean9 vs OctMean12 | `born_village_le2k` | `+0.0802` | `+0.0362` | `+0.0440` | `[+0.0062, +0.0929]` | `0.0240` |
| `common_19_provs` | 2022 Aug15 vs Feb16 | `pre1981` | `+0.1018` | `-0.0013` | `+0.1031` | `[+0.0361, +0.1607]` | `0.0048` |
| `common_19_provs` | 2023 Aug15 vs Feb15 | `pre1981` | `+0.0740` | `-0.0024` | `+0.0765` | `[+0.0032, +0.1455]` | `0.0405` |
| `common_19_provs_exact_cell` | Aug15 vs OctMean12 | `pre1981` | `+0.1305` | `+0.0029` | `+0.1275` | `[+0.0847, +0.1648]` | `0.0001` |
| `common_19_provs_centroid` | Aug15 vs OctMean12 | `pre1981` | `+0.1331` | `+0.0068` | `+0.1264` | `[+0.0851, +0.1628]` | `0.0001` |
| **`full_47_inland_vs_19_coast`** | **Aug15 vs Feb14** | **`pre1981`** | `+0.1545` | `-0.0173` | **`+0.1717`** | `[+0.1269, +0.2052]` | **`0.0001`** |
| **`full_47_inland_vs_19_coast`** | **Aug15 vs OctMean12** | **`pre1981`** | `+0.1540` | `+0.0030` | **`+0.1510`** | `[+0.1144, +0.1816]` | **`0.0001`** |
| **`full_47_inland_vs_19_coast`** | **Aug15 vs Feb14** | **`born_village_le2k`** | `+0.1219` | `+0.0183` | **`+0.1036`** | `[+0.0680, +0.1399]` | **`0.0001`** |
| **`full_47_inland_vs_19_coast`** | **Aug15 vs OctMean12** | **`born_village_le2k`** | `+0.1287` | `+0.0408` | **`+0.0879`** | `[+0.0529, +0.1253]` | **`0.0002`** |

**Surprise & interpretation:**
1. **`pre1981` passes the strictest possible `inland × coast` test with flying colours.** Even within the *same* 19 coastal provinces (and invariant to exact cell intersection or stratum-specific centroids), `pre1981` rises `+0.097` to `+0.131` in inland zones while showing zero detected lift on the coast (`-0.017` / `+0.003`), giving $\gamma = +0.1144$ (`p=0.0002`) vs Feb14 and $+0.1276$ (`p=0.0001`) vs Oct12 — replicating in 2022 (`+0.1031`, `p=0.0048`) and 2023 (`+0.0765`, `p=0.0405`).
2. **`born_village_le2k` behaves differently on strict 19-province support.** On the full peninsula (`full_47_inland_vs_19_coast`), its inland premium is huge (`+0.1036`, `p=0.0001` vs Feb14; `+0.0879`, `p=0.0002` vs Oct12) because the 28 fully landlocked provinces carry the strongest village-born summer pull (`+0.1637` / `+0.1693`, `b2_main_contrasts.csv`, `inlandprov`). Inside the 19 coastal provinces alone, however, $\gamma$ drops to `+0.0496` (`p=0.075`) vs Feb14 and `+0.0345` (`p=0.092`) vs Oct12 (passing 5% only on the 9-night average `AugMean9 vs OctMean12`, `+0.0440`, `p=0.024`, and not detected in 2022/2023 on `common_19_provs`), because coastal zones of those 19 provinces also see a positive village-born summer lift (`+0.0408 [+0.0107, +0.0744]`, `p=0.0086` vs Oct12; e.g., migrants born in coastal-province villages returning to coastal provinces).

*Verdict (B.3):* **Strongly confirmed for `pre1981` (P3); confirmed peninsula-wide (`full_47`) but borderline within coastal provinces (`common_19_provs`) for `born_village_le2k`.**

### 1.5 B.4 — Age-Matched Cohort and Birthplace Contrasts (`b4_age_matched.csv`)
In the 2024 census, everyone who arrived before 1981 is aged $\ge 43$ (`59,060` aged 43–64; `184,798` aged 65+). To separate arrival cohort and rural origin from age, we run eight specifications (`A0`–`A7` in `out/b4_age_matched.csv`) with direct within-age-band WCR differences:
- **Audit fix (`A0` vs `A1`, `Aug15 vs Feb14`):** dropping the 3,245 `VARANOP='7777'` records aged 65+ that V1 had folded into `post1981_ge65` shifts `post1981_ge65` from `-0.0058` to `-0.0087` (`all`) and `+0.0118` to `+0.0088` (`inland`), slightly *widening* the within-65+ contrast from `+0.0926` to **`+0.0965 [+0.0490, +0.1460]`** (`p=0.0003`, `all`) and `+0.1170` to **`+0.1210 [+0.0640, +0.1713]`** (`p=0.0003`, `inland`).

| Spec & Within-Band Contrast | Baseline | `all` $\Delta$ (95% WCR CI; $p$) | `inland` $\Delta$ (95% WCR CI; $p$) |
|---|---|---|---|
| **A2 (`43–64`): `pre1981 − post1981`** | Feb14 | `+0.0374 [-0.0002, +0.0786]` (`0.0518`) | `+0.0471 [-0.0002, +0.0975]` (`0.0512`) |
| **A2 (`43–64`): `pre1981 − post1981`** | Oct12 | **`+0.0239 [+0.0010, +0.0477]`** (`0.0406`) | **`+0.0297 [+0.0035, +0.0551]`** (`0.0285`) |
| **A2 (`65+`): `pre1981 − post1981`** | Feb14 | **`+0.0947 [+0.0484, +0.1430]`** (`0.0003`) | **`+0.1187 [+0.0646, +0.1666]`** (`0.0003`) |
| **A2 (`65+`): `pre1981 − post1981`** | Oct12 | **`+0.0565 [+0.0285, +0.0851]`** (`0.0004`) | **`+0.0914 [+0.0558, +0.1212]`** (`0.0002`) |
| **A3 (`43–64`): `pre1981 − recent`** | Feb14 | **`+0.0395 [+0.0119, +0.0710]`** (`0.0086`) | **`+0.0661 [+0.0292, +0.1038]`** (`0.0021`) |
| **A3 (`43–64`): `pre1981 − recent`** | Oct12 | **`+0.0339 [+0.0128, +0.0565]`** (`0.0050`) | **`+0.0557 [+0.0297, +0.0824]`** (`0.0010`) |
| **A3 (`65+`): `pre1981 − recent`** | Feb14 | **`+0.1205 [+0.0596, +0.1775]`** (`0.0001`) | **`+0.1772 [+0.1165, +0.2299]`** (`0.0001`) |
| **A3 (`65+`): `pre1981 − recent`** | Oct12 | **`+0.0724 [+0.0301, +0.1137]`** (`0.0006`) | **`+0.1332 [+0.0898, +0.1705]`** (`0.0001`) |
| **A4 (`43–59` / `60+`):** `pre1981 − post1981` | Feb14 | `43–59`: `+0.0384` (`0.0752`) / `60+`: **`+0.0943`** (`0.0003`) | `43–59`: `+0.0411` (`0.0642`) / `60+`: **`+0.1228`** (`0.0003`) |
| **A5 (`<43`): `village − urban`** | Feb14 / Oct12 | `+0.0266` (`0.355`) / `+0.0091` (`0.621`) | `+0.0442` (`0.085`) / `+0.0174` (`0.402`) |
| **A5 (`43–64`): `village − urban`** | Feb14 / Oct12 | **`+0.0837`** (`0.0001`) / **`+0.0580`** (`0.0002`) | **`+0.0764`** (`0.0008`) / **`+0.0598`** (`0.0016`) |
| **A5 (`65+`): `village − urban`** | Feb14 / Oct12 | **`+0.0748`** (`0.0001`) / **`+0.0618`** (`0.0001`) | **`+0.0742`** (`0.0002`) / **`+0.0637`** (`0.0001`) |
| **A6/A7** (+46 city `share65+ × dest` slopes) | Feb14 / Oct12 | `pre1981`: `+0.0796` / `+0.0824`; `vil`: `+0.0819` / `+0.0914` (`p <= 0.0002`) | `pre1981`: `+0.1561` / `+0.1545`; `vil`: `+0.1233` / `+0.1335` (`p=0.0001`) |

*Verdict (B.4):* **Confirmed.** Age alone does not generate the summer roots lift: within the **same** 65+ band, `pre-1981` arrivals beat `post-1981` arrivals by `+0.057` to `+0.119` (`p <= 0.0004`), and within `43–64`, `pre-1981` beats `recent` by `+0.034` to `+0.066` (`p <= 0.0086`). Notice also that `post1981_ge65` itself shows a small positive inland lift against Oct12 (`+0.0256 [+0.0058, +0.0467]`, `p=0.011`), still less than a quarter of `pre1981_ge65` (`+0.1170`). Meanwhile, `born_village_le2k` dominates `born_urban_gt10k` equally in **both** `43–64` (`+0.058` to `+0.084`) and `65+` (`+0.062` to `+0.075`), while showing no detected effect among `<43` migrants (`village_lt43` `-0.012`, `p >= 0.20`). Controlling for 46 origin-city age-composition × destination slopes (`A6`/`A7`) leaves the headline coefficients virtually unchanged.

### 1.6 B.5 — Zero-Cell Rules on a Common Predicted-Flow Scale (`b5_zero_cells_common_scale.csv`)
Raw PPML coefficients vary by an order of magnitude across zero-cell transformations (`+0.017` under `log(share+1e-5)` to `+0.201` under `asinh(1000*share)`) simply because each transform compresses or stretches small shares differently. To put all seven Table A2 rules on a **common, directly interpretable scale**, for each dyad $d$ we compute the exact nonlinear shift $\Delta g_r(d) = g_r(s_d + \Delta s) - g_r(s_d)$ induced by a **+1 within-FE SD increase in the raw group share** $\Delta s$ (same physical increment in persons per resident across every rule), and report the flow-weighted mean percentage change in the August/ordinary expected-flow ratio, $100(\exp(\hat{\delta}_r \overline{\Delta g}_r^{\text{fw}}) - 1)$, alongside the exact aggregate-flow ratio and `+1 within-FE SD of g_r` (source: `out/b5_zero_cells_common_scale.csv`).

| Rule (`Aug15 vs Feb14`) | `pre1981` raw $\delta$ (`p`) | `pre1981` % flow (`all`, +1 SD raw share) | `pre1981` % flow (`inland`, +1 SD raw share) | `village_le2k` raw $\delta$ (`p`) | `village_le2k` % flow (`all`, +1 SD raw share) | `village_le2k` % flow (`inland`, +1 SD raw share) |
|---|---:|---:|---:|---:|---:|---:|
| **`haldane_0.5_main`** | `+0.0836` (`0.0003`) | **`+7.9% [+4.2, +11.3]`** | **`+13.1% [+9.2, +16.4]`** | `+0.0857` (`0.0001`) | **`+13.6% [+9.5, +17.5]`** | **`+15.5% [+11.6, +19.1]`** |
| `plus_1.0_log1p` | `+0.1125` (`0.0002`) | `+9.9% [+5.9, +13.2]` | `+14.9% [+10.8, +18.2]` | `+0.0994` (`0.0001`) | `+13.9% [+9.9, +17.3]` | `+14.8% [+10.9, +18.8]` |
| `asinh_1000_share` | `+0.1409` (`0.0001`) | `+10.3% [+6.2, +13.5]` | `+14.5% [+10.4, +17.7]` | `+0.2006` (`0.0001`) | `+14.7% [+10.1, +18.7]` | `+14.1% [+10.9, +17.5]` |
| `log_share_plus_1e4` | `+0.0589` (`0.0015`) | `+5.6% [+2.4, +8.5]` | `+10.6% [+6.8, +13.7]` | `+0.0898` (`0.0001`) | `+12.8% [+9.0, +16.5]` | `+14.2% [+10.4, +18.0]` |
| `plus_0.1_pop_plus1` | `+0.0379` (`0.0063`) | `+4.1% [+1.1, +7.0]` | `+8.4% [+4.9, +11.7]` | `+0.0447` (`0.0001`) | `+8.8% [+4.4, +13.7]` | `+11.9% [+8.4, +15.4]` |
| `log_share_plus_1e5` | `+0.0168` (`0.0907`) | `+2.0% [-0.3, +4.1]` | `+5.2% [+2.2, +7.9]` (`p=0.0003`) | `+0.0278` (`0.0025`) | `+5.9% [+2.0, +9.9]` | `+8.1% [+4.9, +11.2]` |
| `R2_112_large_cities` | `+0.0815` (`0.0250`) | `+7.3% [+0.8, +14.3]` | `+12.5% [+5.7, +19.2]` | `+0.0988` (`0.0001`) | `+17.4% [+10.7, +26.4]` | `+18.5% [+10.4, +28.6]` |

On `Aug15 vs OctMean12` (`b5_zero_cells_common_scale.csv`), **all 7 rules** are significant for both `pre1981` and `born_village_le2k` (even `log_share_plus_1e5` gives `+2.1% [+0.1%, +3.8%]`, `p=0.035` `all` and `+4.8% [+2.1%, +7.1%]`, `p=0.0001` `inland` for `pre1981`; `+6.2%` `all` and `+7.8%` `inland`, `p=0.0001` for `village_le2k`). Across the 6 non-extreme rules, the predicted flow increase for `+1 within-FE SD` of raw share clusters tightly around **`+4%` to `+10%`** (`all`) and **`+8%` to `+15%`** (`inland`) for `pre1981`, and **`+9%` to `+17%`** (`all`) and **`+12%` to `+19%`** (`inland`) for `village_le2k` (and on `+1 within-FE SD` of $g_r$, `inland` `pre1981` spans just `+7.1%` to `+10.3%` across **all 7 rules**, including `1e-5`!).

**Complete 3-term partitions (`b5_zero_cells_common_scale.csv`, Haldane 0.5):** adding the omitted middle category confirms a clean monotonic gradient on both baselines:
- `cohort_complete_3term` (`inland`, `Aug15 vs OctMean12`): `pre1981` **`+0.1301`** (`p=0.0001`, `+10.9%`), `mid_1981_2015` **`+0.0508`** (`p=0.0004`, `+5.1%`), `recent_2016_2023` **`-0.0333`** (`p=0.0408`, `-3.9%`). (Vs Feb14: `+0.1381`, `+0.0221` n.s., `-0.0654`.)
- `birthsize_complete_3term` (`inland`, `Aug15 vs OctMean12`): `born_village_le2k` **`+0.1122`** (`p=0.0001`, `+14.0%`), `born_mid_2k_10k` **`+0.0493`** (`p=0.0001`, `+5.6%`), `born_urban_gt10k` **`-0.0071`** (`p=0.712`, `-0.7%`).

*Verdict (B.5):* **Confirmed, with transparent `1e-5` caveat.** Moving from raw elasticities to a common predicted-flow contrast collapses most of the arbitrary scale difference across zero-cell rules; only `log(share+1e-5)` attenuates magnitudes materially (because `1e-5` assigns ~4 log-units of leverage to zero-count cells in a 10% sample), and even there inland `pre1981` and all `village_le2k` contrasts remain `p <= 0.0025`.

### 1.7 B.6 — Province Jackknives and Dropping Madrid & Barcelona (`b6_jackknife_summary.csv`, `b6_dest_jackknife_47provs.csv`, `b6_drop_mad_bcn.csv`)
- **44-origin-province jackknife (`b6_jackknife_summary.csv`):** leaving out each origin province one at a time leaves `pre1981` significant in **44/44** runs (`Aug15 vs Feb14` `all`: range `[+0.0717, +0.0889]`, `max p = 0.0007`; `inland`: `[+0.1356, +0.1591]`, `max p = 0.0001`; `Aug15 vs OctMean12` `all`: `[+0.0786, +0.0931]`, `max p = 0.0002`; `inland`: `[+0.1435, +0.1589]`, `max p = 0.0001`) and `born_village_le2k` significant in **44/44** runs (`max p = 0.0001` everywhere). By contrast, `recent_2016_2023` is significant in `44/44` vs Feb14 (`max p = 0.0265`) but only `3/44` (`all`) / `2/44` (`inland`) vs Oct12.
- **47-destination-province jackknife (`b6_dest_jackknife_47provs.csv`):** leaving out each destination province one at a time preserves **47/47** significant runs for `pre1981` (`max p = 0.0003`) and `born_village_le2k` (`max p = 0.0001`) across both baselines and `all`/`inland`.
- **Dropping Madrid, Barcelona, or both (`b6_drop_mad_bcn.csv`):**

| Drop spec (`G`) | Contrast | `all` `pre1981` (CI; $p$) | `inland` `pre1981` (CI; $p$) | `all` `village_le2k` (CI; $p$) | `common_19` `inland×coast` `pre1981` (CI; $p$) |
|---|---|---|---|---|---|
| `drop_Madrid` (`43`) | Feb14 | `+0.0889 [+0.0487, +0.1211]` (`0.0005`) | `+0.1580 [+0.1101, +0.1948]` (`0.0001`) | `+0.0904 [+0.0639, +0.1121]` (`0.0001`) | `+0.1274 [+0.0818, +0.1667]` (`0.0001`) |
| `drop_Madrid` (`43`) | Oct12 | `+0.0875 [+0.0503, +0.1194]` (`0.0001`) | `+0.1532 [+0.1081, +0.1911]` (`0.0001`) | `+0.0989 [+0.0706, +0.1256]` (`0.0001`) | `+0.1346 [+0.0928, +0.1673]` (`0.0001`) |
| `drop_Barcelona` (`43`) | Feb14 | `+0.0717 [+0.0337, +0.1078]` (`0.0008`) | `+0.1356 [+0.0905, +0.1759]` (`0.0001`) | `+0.0826 [+0.0542, +0.1095]` (`0.0001`) | `+0.0871 [+0.0350, +0.1442]` (`0.0012`) |
| `drop_Barcelona` (`43`) | Oct12 | `+0.0797 [+0.0434, +0.1139]` (`0.0001`) | `+0.1436 [+0.0961, +0.1883]` (`0.0001`) | `+0.0908 [+0.0663, +0.1209]` (`0.0001`) | `+0.1037 [+0.0623, +0.1469]` (`0.0001`) |
| **`drop_Madrid_and_Barcelona` (`42`)** | **Feb14** | **`+0.0745 [+0.0301, +0.1142]`** (`0.0016`) | **`+0.1333 [+0.0799, +0.1802]`** (`0.0001`) | **`+0.0869 [+0.0521, +0.1161]`** (`0.0001`) | **`+0.1010 [+0.0507, +0.1461]`** (`0.0005`) |
| **`drop_Madrid_and_Barcelona` (`42`)** | **Oct12** | **`+0.0769 [+0.0338, +0.1195]`** (`0.0004`) | **`+0.1368 [+0.0837, +0.1918]`** (`0.0001`) | **`+0.0972 [+0.0622, +0.1319]`** (`0.0001`) | **`+0.1137 [+0.0746, +0.1471]`** (`0.0001`) |

Notably, dropping Madrid actually *strengthens* `common_19_provs` `born_village_le2k` to `+0.0678` (`p=0.0069`) vs Feb14 and `+0.0448` (`p=0.0308`) vs Oct12, whereas dropping Barcelona slightly reduces `pre1981` (by ~0.012) while keeping every primary estimate `p <= 0.0016`.

*Verdict (B.6):* **Confirmed.** Neither Madrid, Barcelona, nor any single origin or destination province drives the `pre1981` or `born_village_le2k` results.

---

## 2. Analysis D — The Receiving End: Emptied Villages

### 2.1 Fixed Sample and 306-Night Coverage (`d_night_manifest.csv`, `d_zones_main_nights.csv`, `d_quintiles_all_nights.csv`)
- **Fixed zone sample:** starting from the 2,435 mainland zones (`pipeline/config.py`), excluding the 3 artefact zones (`31905_AM`, `19103_AM`, `40149_AM`) and zones missing historical peak-to-2021 population loss rank or covariates leaves **`n = 2,376` zones** across 47 mainland provinces (`Q1`: 476; `Q2`–`Q5`: 475 each; exact match to coursework R1). Subsets:
  - **`single_mun` (`n = 1,473`):** strictly `n_municipalities == 1` (`Q1`: 432, `Q2`: 403, `Q3`: 320, `Q4`: 235, `Q5`: 84). Note that zone `03015` (*Almoradí*, Alicante) groups 2 municipalities despite lacking an `_AM` suffix and is therefore excluded from `single_mun`.
  - **`inland` (`n = 2,053`):** `coastal_true == False`. Coastal zones (`323` total) are heavily skewed toward low-loss quintiles: `Q1` has 94 (`19.8%`), `Q2` has 102 (`21.5%`), `Q3` has 70 (`14.7%`), `Q4` has 44 (`9.3%`), and `Q5` has only 13 (`2.7%`).
  - **`inland_single_mun` (`n = 1,186`):** `coastal_true == False` and `n_municipalities == 1`.
- **Flow reconstruction & manifest:** per-zone `s_in` (arrivals from outside zone $i$) and `s_out` (residents of zone $i$ sleeping outside $i$) were computed directly from the cached parquets across all **306 cached nights** (65 in 2022, 61 in 2023, 180 in 2024) and verified against `results/p_zone_daily.parquet` (`max |diff| < 1e-6`). Zero cached nights fail the 99% zone-coverage rule or the 10% total-persons deviation rule; the 3 DANA nights (`2024-10-29..31`) are flagged and excluded from `2024-10-mean12`.
- **All 93 cached August nights (2022–2024, `d_quintiles_all_nights.csv`):** Q5 (most emptied fifth) is overwhelmingly net-receiving on **every single cached August night** across all three years:
  - `all` (`n = 475` in Q5): Q5 net-receiver share spans **`96.63%` to `98.95%`** (median **`98.11%`**; minimum `2023-08-29`).
  - `single_mun` (`n = 84` in Q5): Q5 net-receiver share spans **`86.90%` to `97.62%`** (median **`92.86%`**; minimum `2023-08-30`).
- **Baseline surprise across ordinary nights (`d_quintiles_main_boot.csv`):** Q5 net-receiver share on ordinary nights varies substantially by season and year: **`24.00% [20.36%, 27.80%]`** on `2024-02-14` and **`25.47% [21.54%, 29.28%]`** on `2023-02-15`, rising to **`47.37% [43.02%, 51.92%]`** on `2024-10-mean12` and **`55.58% [51.45%, 60.00%]`** on `2022-02-16` (in `single_mun`, Q5 sits near `32%–40%` on all ordinary nights). This reflects lower ordinary departures from rural zones in February 2022 (`0.039/res` vs `0.095/res` in Feb 2024) and slightly higher autumn weekend/weekday arrivals in October 2024 (`0.092/res` vs `0.066/res` in Feb 2024). Even from the highest ordinary baseline, the August shift is massive (`+51.6 pp [+47.1, +56.0]` vs Oct12; `+74.9 pp [+71.0, +78.7]` vs Feb14; `+43.2 pp` in 2022).

### 2.2 D.1 & D.2 — Arrivals vs Departures: Unadjusted Quintiles and Controlled FE-OLS

All uncertainty intervals below use within-province zone resampling (`B = 4,999`, `zone boot`); controlled regressions additionally report province-cluster bootstrap (`B = 4,999`, `prov boot`) and 9,999-draw Webb WCR/WCU clustered by province ($G = 47$).

#### Table 2A — Unadjusted Q5 and Q5 − Q1 Levels & Paired Changes (`d_quintiles_main_boot.csv`, `all` & `single_mun`)

| Sample & Period | Group | Share net receiver (95% zone CI) | Arrivals / res (95% zone CI) | Departures / res (95% zone CI) | Net / res (95% zone CI) |
|---|---|---:|---:|---:|---:|
| **`all` levels** | | | | | |
| 2024-02-14 (Feb ordinary) | Q5 (`n=475`) | `24.0% [20.4, 27.8]` | `0.066 [0.063, 0.068]` | `0.095 [0.092, 0.099]` | `-0.029 [-0.033, -0.025]` |
| 2024-10-mean12 (Oct ordinary) | Q5 (`n=475`) | `47.4% [43.0, 51.9]` | `0.092 [0.089, 0.096]` | `0.090 [0.088, 0.093]` | `+0.002 [-0.001, +0.005]` |
| 2024-08-15 (Summer) | Q5 (`n=475`) | **`98.9% [98.0, 99.8]`** | **`0.710 [0.672, 0.750]`** | **`0.082 [0.080, 0.084]`** | **`+0.628 [+0.590, +0.668]`** |
| 2024-12-25 (Christmas Eve) | Q5 (`n=475`) | **`99.6% [98.9, 100.0]`** | **`0.380 [0.370, 0.390]`** | **`0.143 [0.139, 0.147]`** | **`+0.237 [+0.227, +0.246]`** |
| 2024-03-29 (Good Friday) | Q5 (`n=475`) | `98.5% [97.4, 99.5]` | `0.536 [0.508, 0.567]` | `0.104 [0.101, 0.107]` | `+0.432 [+0.404, +0.462]` |
| 2024-11-01 (All Saints) | Q5 (`n=475`) | `94.1% [92.0, 96.1]` | `0.236 [0.226, 0.248]` | `0.099 [0.094, 0.104]` | `+0.138 [+0.126, +0.150]` |
| **`all` paired changes** | | | | | |
| **Aug15 − Feb14** | Q5 | `+74.9 pp [+71.0, +78.7]` | **`+0.644 [+0.607, +0.683]`** | **`-0.013 [-0.017, -0.009]`** | **`+0.657 [+0.618, +0.699]`** |
| **Aug15 − Feb14** | Q5 − Q1 | `+83.8 pp [+77.1, +90.7]` | **`+0.479 [+0.431, +0.527]`** | **`-0.148 [-0.155, -0.141]`** | **`+0.627 [+0.576, +0.679]`** |
| **Aug15 − OctMean12** | Q5 − Q1 | `+62.5 pp [+55.8, +68.9]` | **`+0.468 [+0.423, +0.514]`** | **`-0.133 [-0.140, -0.127]`** | **`+0.602 [+0.552, +0.650]`** |
| **Dec25 − Feb14** | Q5 | `+75.6 pp [+71.7, +79.2]` | **`+0.314 [+0.305, +0.323]`** | **`+0.048 [+0.045, +0.051]`** | **`+0.266 [+0.257, +0.275]`** |
| **Dec25 − Feb14** | Q5 − Q1 | `+67.4 pp [+60.7, +74.1]` | **`+0.204 [+0.193, +0.216]`** | **`-0.049 [-0.055, -0.044]`** | **`+0.254 [+0.242, +0.265]`** |
| **`single_mun` paired changes** | | | | | |
| **Aug15 − Feb14** | Q5 (`n=84`) | `+53.6 pp [+41.8, +64.5]` | `+0.278 [+0.233, +0.329]` | `+0.020 [+0.014, +0.026]` | `+0.258 [+0.211, +0.311]` |
| **Aug15 − Feb14** | Q5 − Q1 | `+61.9 pp [+48.8, +74.4]` | `+0.113 [+0.057, +0.172]` | `-0.116 [-0.125, -0.107]` | `+0.229 [+0.170, +0.291]` |
| **Dec25 − Feb14** | Q5 (`n=84`) | `+58.3 pp [+47.9, +68.5]` | `+0.227 [+0.213, +0.242]` | `+0.038 [+0.032, +0.045]` | `+0.189 [+0.174, +0.204]` |
| **Dec25 − Feb14** | Q5 − Q1 | `+50.0 pp [+38.0, +61.2]` | `+0.123 [+0.107, +0.138]` | `-0.059 [-0.066, -0.051]` | `+0.181 [+0.165, +0.198]` |

In unadjusted terms, **Q5 itself flips to net-receiving in August overwhelmingly because of visitors arriving**: arrivals per resident surge tenfold (`0.066 -> 0.710`, `+0.644`), accounting for **98%** of Q5's own net gain (`+0.657`), while Q5 residents' own departures barely budge (`0.095 -> 0.082`, `-0.013`). Relative to Q1 (growing/urban zones whose departures jump from `0.054` to `0.189` in August), **76.4%** of Q5's differential net gain (`+0.627`) comes from higher arrivals (`+0.479`) and **23.6%** from smaller departures (`-0.148`; in `single_mun`, the `Q5 − Q1` split is `49.3%` arrivals / `50.7%` departures because multi-municipality `_AM` rural zones receive the largest arrival ratios). On Christmas Eve, Q5 arrivals jump `+0.314` (`0.066 -> 0.380`) while Q5 departures *rise* `+0.048` (`0.095 -> 0.143`, still smaller than Q1's `+0.097` rise), so Q5's `Q5 − Q1` net gain (`+0.254`) splits **80.6%** arrivals / **19.4%** smaller departure rise.

#### Table 2B — Controlled FE-OLS Decomposition (`d_decomp_controlled.csv` and `d_decomp_controlled_inland.csv`)
Regressing each outcome on `loss_rank` (0 to 1) with 47 province FE, $\log(\text{pop})$, $\log(\text{area})$, `median_income`, and `coastal_true` (for `all` and `single_mun`) satisfies the exact linear identity $\hat{\beta}_{\text{net}} \equiv \hat{\beta}_{\text{arr}} - \hat{\beta}_{\text{dep}}$.

- **Exact reproduction of coursework R1 level slopes (`d_decomp_controlled.csv`, `2024-08-15`, `net_importer`):** `+39.2 pp` (`all`, zone boot `[+27.9, +49.8]`, WCU `[+23.0, +55.3]`, `p=0.0001`) and `+43.8 pp` (`single_mun`, zone boot `[+29.2, +56.3]`, WCU `[+21.0, +67.0]`, `p=0.0033`).
- **Major finding — why controlled `all` vs `inland` diverge in August:**

| Sample (`n`) | Contrast | `arrivals_per_res_diff` (zone CI; WCU CI) | `departures_per_res_diff` (zone CI; WCU CI) | `net_per_res_diff` (zone CI; WCU CI) | `log_s_in_ratio` (zone CI) | `net_importer_diff` (zone CI) |
|---|---|---:|---:|---:|---:|---:|
| **`all`** (`2,376`) | **Aug15 − Feb14** | `-0.040 [-0.169, +0.073]` (`[-0.211, +0.132]`) | **`-0.069 [-0.083, -0.055]`** (`[-0.092, -0.044]`) | `+0.029 [-0.105, +0.146]` (`[-0.158, +0.218]`) | `+0.459 [+0.298, +0.613]` | `+51.5 pp [+34.0, +68.5]` |
| **`all`** (`2,376`) | **Aug15 − Oct12** | `-0.047 [-0.176, +0.065]` (`[-0.225, +0.132]`) | **`-0.070 [-0.082, -0.057]`** (`[-0.093, -0.046]`) | `+0.023 [-0.111, +0.138]` (`[-0.174, +0.213]`) | `+0.380 [+0.240, +0.512]` | `+76.8 pp [+60.2, +93.3]` |
| **`all`** (`2,376`) | **Dec25 − Feb14** | **`+0.111 [+0.091, +0.130]`** (`[+0.075, +0.147]`) | **`-0.079 [-0.089, -0.068]`** (`[-0.097, -0.061]`) | **`+0.190 [+0.166, +0.212]`** (`[+0.143, +0.237]`) | `+0.687 [+0.577, +0.796]` | `+34.0 pp [+17.0, +50.4]` |
| `single_mun` (`1,473`) | Aug15 − Feb14 | `-0.176 [-0.337, -0.052]` (`[-0.356, +0.005]`) | `-0.055 [-0.070, -0.036]` (`[-0.085, -0.024]`) | `-0.121 [-0.291, +0.008]` (`[-0.327, +0.087]`) | `+0.141 [-0.061, +0.323]` | `+31.3 pp [+10.0, +50.3]` |
| `single_mun` (`1,473`) | Dec25 − Feb14 | `+0.067 [+0.046, +0.085]` (`[+0.037, +0.097]`) | `-0.079 [-0.090, -0.066]` (`[-0.098, -0.059]`) | `+0.145 [+0.119, +0.167]` (`[+0.104, +0.187]`) | `+0.515 [+0.391, +0.641]` | `+17.2 pp [-4.4, +37.2]` |
| **`inland`** (`2,053`) | **Aug15 − Feb14** | **`+0.144 [+0.058, +0.224]`** (`[+0.005, +0.285]`) | **`-0.076 [-0.090, -0.060]`** (`[-0.103, -0.049]`) | **`+0.221 [+0.130, +0.306]`** (`[+0.073, +0.375]`) | **`+0.713 [+0.553, +0.871]`** | **`+66.9 pp [+48.2, +85.9]`** |
| **`inland`** (`2,053`) | **Aug15 − Oct12** | **`+0.126 [+0.031, +0.211]`** (`[-0.030, +0.275]`) | **`-0.076 [-0.089, -0.063]`** (`[-0.102, -0.050]`) | **`+0.202 [+0.101, +0.292]`** (`[+0.036, +0.365]`) | **`+0.543 [+0.402, +0.683]`** | **`+83.4 pp [+64.9, +101.9]`** |
| **`inland`** (`2,053`) | **Dec25 − Feb14** | **`+0.120 [+0.097, +0.143]`** (`[+0.081, +0.159]`) | **`-0.080 [-0.092, -0.068]`** (`[-0.102, -0.059]`) | **`+0.200 [+0.173, +0.226]`** (`[+0.150, +0.250]`) | **`+0.713 [+0.589, +0.839]`** | **`+36.4 pp [+17.6, +54.1]`** |
| `inland` (`2,053`) | 2022 / 2023 Aug−Feb | `+0.149` / `+0.144` (zone `p <= 0.004`) | `-0.078` / `-0.075` (`p=0.0004`) | `+0.227` / `+0.219` (WCR `p <= 0.012`) | `+0.679` / `+0.772` | `+63.4 pp` / `+66.8 pp` |
| **`inland_single_mun`** (`1,186`) | **Aug15 − Feb14** | **`+0.066 [+0.018, +0.109]`** (`[-0.016, +0.147]`) | **`-0.061 [-0.079, -0.040]`** (`[-0.098, -0.024]`) | **`+0.127 [+0.067, +0.179]`** (`[+0.017, +0.240]`) | `+0.434 [+0.228, +0.624]` | `+50.7 pp [+26.1, +73.8]` |
| `inland_single_mun` (`1,186`) | Dec25 − Feb14 | `+0.070 [+0.045, +0.091]` (`[+0.038, +0.102]`) | `-0.079 [-0.093, -0.064]` (`[-0.103, -0.056]`) | `+0.149 [+0.119, +0.175]` (`[+0.106, +0.192]`) | `+0.527 [+0.383, +0.665]` | `+18.2 pp [-7.6, +40.7]` |

**Why this matters for the paper:** The coursework report noted that in the full-sample controlled regression, August appeared to work through fewer residents leaving (`-0.069`) while linear net gain per resident was not detected (`+0.029`), in contrast to Christmas Eve (`+0.190`). Our decomposition exposes the exact reason: **coastal resort zones** sit overwhelmingly in low-loss `Q1–Q2` (196 of 323 coastal zones, vs 13 in `Q5`) and have huge right-tail August arrivals per resident (`arrivals/res` of 2–6+), which an additive `coastal_true` dummy cannot remove from a linear regression (whereas `log_s_in_ratio` damps the tail and was already `+0.459`!). Once coastal zones are excluded (`inland`, `n = 2,053`), **August and Christmas Eve look remarkably parallel**:
- `Aug15 − Feb14` (`inland`): `arrivals_per_res_diff = +0.144` (**65.5%** of net, zone boot `p=0.0008`, WCR `p=0.049`) + fewer departures `-0.076` (**34.5%**, `p=0.0002`) = **`+0.221`** (`WCR p=0.007`), replicating identically in 2022 (`+0.227`) and 2023 (`+0.219`), and holding on `inland_single_mun` (`+0.127`, `WCR p=0.044`, 51.9% arrivals / 48.1% fewer departures).
- `Dec25 − Feb14` (`inland`): `arrivals_per_res_diff = +0.120` (**59.9%** of net) + fewer departures `-0.080` (**40.1%**) = **`+0.200`** (`WCR p=0.0001`). Even `log_s_in_ratio` is identical to three decimals (`0.713` vs `0.713`)!

*Verdict (D.1–D.2):* **Confirmed, with a major mechanistic refinement.** In the emptied interior, August does **not** work merely by fewer residents leaving: it combines a large surge in relative arrivals (~65% of controlled inland net gain; 98% of raw Q5 net gain) with a smaller departure rate (~35% of controlled inland net gain, because urban residents leave in droves while rural residents stay put).

### 2.3 D.3 — Non-Primary Housing Slope Attenuation: August vs Christmas Eve (`d_housing_attenuation.csv`, `d_housing_aug_vs_xmas_diff.csv`)
Using Census 2021 table 59525 non-primary dwellings (`share_np` = non-primary / total dwellings, plus `log1p_np_per_res` as a scale check), we compare the raw `loss_rank` slope $\hat{\beta}_{\text{raw}}$ against the housing-adjusted slope $\hat{\beta}_{\text{adj}}$ on the **exact same sample** for August (14–15 Aug) and Christmas Eve (24–25 Dec), computing the absolute slope attenuation $A = \hat{\beta}_{\text{raw}} - \hat{\beta}_{\text{adj}}$, the percentage attenuation $A / \hat{\beta}_{\text{raw}}$, and the paired difference $\Delta A = A_{\text{Aug}} - A_{\text{Xmas}}$ with both within-province zone bootstrap and province-cluster bootstrap 95% CIs (`B = 4,999`). Note that baseline (`Feb14` vs `OctMean12`) cancels algebraically from $\Delta A = A_{\text{Aug}} - A_{\text{Xmas}}$.

*(Terminology note: our verification script cross-checked the legacy coursework file `results/v1_c3_mediation.csv` to confirm exact numerical agreement to `< 1e-10`; all new outputs and analysis strictly label this **slope attenuation**, never mediation, because cross-sectional housing stock cannot identify whether visitors sleep in non-primary dwellings or family homes.)*

| Sample (`n`) | Outcome & Proxy | Baseline | $\hat{\beta}_{\text{raw}}$ (Aug / Xmas) | $A_{\text{Aug}}$ (share %) | $A_{\text{Xmas}}$ (share %) | **$\Delta A = A_{\text{Aug}} - A_{\text{Xmas}}$** | Zone boot 95% CI ($p$) | Prov boot 95% CI ($p$) |
|---|---|---|---:|---:|---:|---:|---|---|
| **`all`** (`2,376`) | `log_s_in_ratio` (`share_np`) | `vs_Feb14` | `0.459` / `0.687` | `0.221` (`48.1%`) | `0.012` (`1.8%`) | **`+0.209`** | `[+0.120, +0.292]` (`0.0004`) | `[+0.061, +0.361]` (`0.0052`) |
| `all` (`2,376`) | `log_s_in_ratio` (`share_np`) | `vs_OctMean12` | `0.380` / `0.608` | `0.202` (`53.1%`) | `-0.007` (`-1.1%`) | **`+0.209`** | `[+0.120, +0.292]` (`0.0004`) | `[+0.061, +0.361]` (`0.0052`) |
| `all` (`2,376`) | `log_s_in_ratio` (`log1p_np_per_res`) | `vs_Feb14` | `0.459` / `0.687` | `0.173` (`37.7%`) | `-0.007` (`-1.0%`) | **`+0.180`** | `[+0.095, +0.261]` (`0.0004`) | `[+0.042, +0.306]` (`0.0104`) |
| **`inland`** (`2,053`) | **`log_s_in_ratio` (`share_np`)** | **`vs_Feb14`** | **`0.713` / `0.713`** | **`0.313` (`43.9%`)** | **`0.042` (`6.0%`)** | **`+0.271`** | **`[+0.198, +0.343]`** (`0.0004`) | **`[+0.161, +0.398]`** (`0.0004`) |
| `inland` (`2,053`) | `log_s_in_ratio` (`share_np`) | `vs_OctMean12` | `0.543` / `0.543` | `0.311` (`57.3%`) | `0.040` (`7.4%`) | **`+0.271`** | `[+0.198, +0.343]` (`0.0004`) | `[+0.161, +0.398]` (`0.0004`) |
| `inland` (`2,053`) | `arrivals_per_res_diff` (`share_np`) | `vs_Feb14` | `0.144` / `0.120` | `0.186` (`129%`) | `0.042` (`35.2%`) | **`+0.144`** | `[+0.103, +0.187]` (`0.0004`) | `[+0.083, +0.233]` (`0.0004`) |
| `inland` (`2,053`) | `net_per_res_diff` (`share_np`) | `vs_Feb14` | `0.221` / `0.200` | `0.204` (`92.5%`) | `0.034` (`16.9%`) | **`+0.170`** | `[+0.122, +0.219]` (`0.0004`) | `[+0.100, +0.268]` (`0.0004`) |
| `inland` (`2,053`) | `net_importer_diff` (`share_np`) | `vs_Feb14` | `0.669` / `0.364` | `0.145` (`21.7%`) | `0.018` (`5.0%`) | **`+0.127`** | `[+0.089, +0.169]` (`0.0004`) | `[+0.072, +0.199]` (`0.0004`) |
| **`single_mun`** (`1,473`) | **`log_s_in_ratio` (`share_np`)** | **`vs_Feb14`** | `0.141` / `0.515` | `0.086` (`61.2%`) | `0.004` (`0.7%`) | **`+0.083`** | `[-0.039, +0.184]` (`0.1788`) | `[-0.099, +0.278]` (`0.3212`) |

For completeness across other 2024 holidays (`all`, `vs_Feb14`, `log_s_in_ratio`, `share_np`, `d_housing_attenuation.csv`), attenuation sits cleanly between August and Christmas Eve: Good Friday (`Mar29`) `0.178` (`26.3%` of `0.675`) and All Saints (`Nov01`) `0.069` (`12.8%` of `0.537`).

**Caveat (`single_mun`):** When multi-municipality `_AM` zones (which carry the highest non-primary housing stock: mean `share_np = 0.465` vs `0.312` in `single_mun`, and `0.545` vs `0.454` inside Q5) are dropped, the raw August `log_s_in_ratio` slope shrinks from `0.459` to `0.141` (`vs_Feb14`) and the paired Aug − Xmas attenuation difference (`+0.083`) is **not detected** (`zone p = 0.179`, `prov p = 0.321`). Meanwhile, on `inland` (`n = 2,053`), where raw August and Christmas Eve slopes on `log_s_in_ratio` are identical (`0.713` vs `0.713`), the contrast is sharpest (`+0.271 [+0.198, +0.343]`, `43.9%` vs `6.0%`), showing that non-primary housing stock tracks August's emptied-village inflow far more strongly than Christmas Eve's.

*Verdict (D.3):* **Confirmed on `all` and `inland` across baselines, outcomes, and proxies; not detected in `single_mun` alone (explicitly flagged).**

---

## 3. Deviations and Operational Notes vs `PREREGISTRATION.md`

1. **Discovery-cache scope only (by design of this assignment):** `PREREGISTRATION.md` §3 evaluates P1–P3 on the held-out 2022–2024 confirmation set (`_local/external_journal/mitma_pernoct/`). Per task instructions, this worktree analysed only the 306 cached discovery nights so the confirmation set stays untouched for the dedicated confirmation analyst. All pre-registered stacked multi-night estimators (`b2_stacked_multinight.csv`, `b3_inland_coast_interaction.csv`) are implemented and verified here on the discovery cache.
2. **October clean ordinary window (`2024-10-mean12`):** `PREREGISTRATION.md` §2 drops DANA (`2024-10-29..31`) and regional holidays of the five largest regions (which includes Valencia Day, 9 October). Our main October mean uses the 12 non-DANA Mon/Tue/Wed-evening weekday files (`2024-10-mean12`), and we explicitly verified `2024-10-mean10_noVal` (further excluding files `2024-10-09` and `2024-10-10`, evenings 8–9 Oct): `pre1981` shifts by `+0.0008` (`0.0870 -> 0.0878`), `born_village_le2k` by `+0.0015` (`0.0927 -> 0.0942`), and `recent` by `-0.0009` (`-0.0291 -> -0.0300`, `p=0.084`).
3. **Fast WCR CI inversion in `b_common.py`:** Exact Kline–Santos 9,999-draw Webb WCR inversion was accelerated ~4× using adaptive demeaning tolerance (`1e-6 -> 1e-9`) and bracketed secant root-finding on $|t_0(b)| - c_{0.95}(b)$, benchmarked against exact bisection with `max |diff| <= 2.3e-5`.
4. **Auxiliary audit rows in `b2_main_contrasts.csv`:** The 18 `Aug15_vs_OctMean15` rows (which retain DANA solely to compare directly against coursework V1) omit WCR CI inversion (`ci_wcr_lo/hi = NaN`, `p_wcr` populated); all 300 binding contrast rows have complete 9,999-draw WCR CIs.
5. **Single-municipality definition (`n_municipalities == 1`):** Strict `n_municipalities == 1` (`1,473` zones; `1,186` inland) was used rather than string suffix `~id.endswith("_AM")` (`1,474`), because `03015` (*Almoradí*) groups 2 municipalities without an `_AM` suffix.
6. **Added `inland` controlled decomposition (`d_inland_decomp.py`):** Added beyond the minimal D checklist specifically to diagnose why controlled `all` August `arrivals_per_res_diff` was `-0.040` despite `Q5` raw arrivals jumping `+0.644`.

---

## 4. Output Inventory (Code and Tidy CSVs)

### Analysis B — [`journal/analysis/B_cohorts/`](B_cohorts/)
- [`b_common.py`](B_cohorts/b_common.py) — Data loaders, audit recode (`7777` & artefact zones), `FastStackedMultiBlock` HDFE PPML + WCR bootstrap.
- [`run_b_cohorts.py`](B_cohorts/run_b_cohorts.py) — 14-core staged runner producing all B.1–B.6 outputs and `fig_who_pulls_tidy.csv`.
- [`out/b1_v1_audit_comparison.csv`](B_cohorts/out/b1_v1_audit_comparison.csv) (84 rows) — Exact V1 reproduction and stepwise audit diffs (`7777`, artefacts, DANA).
- [`out/b2_key_single_nights.csv`](B_cohorts/out/b2_key_single_nights.csv) (304 rows) — Single-night PPML slopes across key 2022–2024 dates.
- [`out/b2_main_contrasts.csv`](B_cohorts/out/b2_main_contrasts.csv) (318 rows) — Main holiday vs ordinary contrasts (`Feb14`, `OctMean12`, `FebMean13`, `OctMean10_noVal`, 2022–2023, placebo, other holidays) with WCR CIs.
- [`out/b2_stacked_multinight.csv`](B_cohorts/out/b2_stacked_multinight.csv) (60 rows) — Pre-registered multi-night stacked PPML across 2022–2024 with joint WCR CIs.
- [`out/b2_sister_summary.csv`](B_cohorts/out/b2_sister_summary.csv) (48 rows) — Summaries of 4,740 sister pairs (full pair table in `_local/journal_data/BD/b2_sister_pairs.parquet`).
- [`out/b3_inland_coast_interaction.csv`](B_cohorts/out/b3_inland_coast_interaction.csv) (144 rows) — 4-block direct `roots × holiday × inland` interactions across 4 supports.
- [`out/b4_age_matched.csv`](B_cohorts/out/b4_age_matched.csv) (160 rows) — Age-matched cohort and birthplace models (`A0`–`A7`) with within-band WCR CIs.
- [`out/b5_zero_cells_common_scale.csv`](B_cohorts/out/b5_zero_cells_common_scale.csv) (136 rows) — 7 zero-cell rules on common predicted-flow contrasts + complete 3-term partitions.
- [`out/b6_jackknife_44provs.csv`](B_cohorts/out/b6_jackknife_44provs.csv) (704 rows), [`out/b6_jackknife_summary.csv`](B_cohorts/out/b6_jackknife_summary.csv) (16 rows), [`out/b6_dest_jackknife_47provs.csv`](B_cohorts/out/b6_dest_jackknife_47provs.csv) (752 rows), and [`out/b6_drop_mad_bcn.csv`](B_cohorts/out/b6_drop_mad_bcn.csv) (132 rows).
- [`out/fig_who_pulls_tidy.csv`](B_cohorts/out/fig_who_pulls_tidy.csv) (710 rows) — Plotting-ready tidy table (`panel`, `group`, `model`, `baseline`, `dest_type`, `estimand`, `unit`, `point`, `ci_lo`, `ci_hi`, `p_wcr`, `n`, `G`, `note`) covering Panels A–E.

### Analysis D — [`journal/analysis/D_villages/`](D_villages/)
- [`run_d_villages.py`](D_villages/run_d_villages.py) — Full 306-night zone panel, quintile bootstraps, controlled decomposition (`all` & `single_mun`), housing slope attenuation, and `fig_receiving_end_tidy.csv`.
- [`d_inland_decomp.py`](D_villages/d_inland_decomp.py) — Controlled FE-OLS arrivals/departures decomposition on `inland` (`n=2,053`) and `inland_single_mun` (`n=1,186`) with zone boot, province boot, and Webb WCR/WCU.
- [`out/d_night_manifest.csv`](D_villages/out/d_night_manifest.csv) (306 rows) — Cached night quality flags and calendar classes.
- [`out/d_quintiles_all_nights.csv`](D_villages/out/d_quintiles_all_nights.csv) (3,672 rows) — Nightly quintile series (`share_net_importer`, arrivals/departures/net per resident) across all 306 cached nights (`all` & `single_mun`).
- [`out/d_quintiles_main_boot.csv`](D_villages/out/d_quintiles_main_boot.csv) (1,610 rows) — Main nights and paired changes with within-province zone bootstrap 95% CIs (`B = 4,999`).
- [`out/d_zones_main_nights.csv`](D_villages/out/d_zones_main_nights.csv) (2,376 rows) — Per-zone main-night values for mapping/scatter checks (full 727,056-row daily panel in `_local/journal_data/BD/d_zone_daily_panel.parquet`).
- [`out/d_decomp_controlled.csv`](D_villages/out/d_decomp_controlled.csv) (266 rows) and [`out/d_decomp_controlled_inland.csv`](D_villages/out/d_decomp_controlled_inland.csv) (92 rows) — Controlled FE-OLS decompositions.
- [`out/d_housing_attenuation.csv`](D_villages/out/d_housing_attenuation.csv) (196 rows) and [`out/d_housing_aug_vs_xmas_diff.csv`](D_villages/out/d_housing_aug_vs_xmas_diff.csv) (49 rows) — Common-sample non-primary housing slope attenuation and paired August − Christmas Eve differences.
- [`out/fig_receiving_end_tidy.csv`](D_villages/out/fig_receiving_end_tidy.csv) (832 rows) — Plotting-ready tidy table (`panel`, `sample`, `period`, `group`, `metric`, `stat`, `point`, `ci_lo`, `ci_hi`, `n`).
