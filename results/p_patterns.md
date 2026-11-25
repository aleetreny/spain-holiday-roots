# P: pattern mining for question B ("Roots or beaches?")

*Agent P (wave 2), 2026-11-24. Code: `pipeline/p_patterns.py` (stages `daily zones hubs cities holidays distance panel plots`, run with `.venv/bin/python -I pipeline/p_patterns.py`). Every number below is read from the `results/p_*.csv` file named next to it. Quick plots: `results/p_diag/`.*

**Scope.** These results are descriptive. I used the b00 loaders throughout (N = 2,435, `coastal_true`, w ≥ 10 and equal-edges backbones). I computed one network per night for all 180 nights of 2024 in the data and all 126 nights of 2022–23 (`p_daily_metrics.csv`, 306 rows). R1 owns the node-level k-core test, R2 the roots inference (PPML, QAP) and R3 the bipartite territories, so I did not repeat their tests. Where a pattern touches their ground, I flag it for them.

**Definitions.**
- *Excess* of a night = w(night) − w_base per pair. w_base is the mean over the 13 Feb-2024 Tue/Wed/Thu sister nights, with a missing pair counted as 0 (`b00.window_mean_from_frame`, audit C3).
- *Emptied* = inland zone with `loss_peak_2021` ≥ 0.5 (457 zones).
- *Roots index* = correlation over the 137 mainland INE cities × the other 46 provinces between log(1+flow) and log(1+born-in-p residents), after removing city means and log distance. It is descriptive only.
- **Classes.** ROBUST = same sign and similar size on every sister night and in 2022 and 2023. FRAGILE = depends on a choice (baseline, threshold, extrapolation) or fails one replication. ARTEFACT = produced by the data or the method, not by behaviour.

## Ranked list

| # | Pattern | Class | Bearing on the question | Where |
|---|---|---|---|---|
| 1 | File date d is the night **ending** on d | ROBUST (data convention) | Changes the meaning of every "main night" | Main report (methods), 1 line + footnote |
| 2 | Holidays come in two modes: summer (August, partly Easter) and family/weekend (All Saints, Christmas, puentes, ordinary Sundays) | ROBUST | Changes the story: All Saints is a big weekend, not a small August | Main report (Result 1 or 2 framing) |
| 3 | The emptied-village share of away nights depends on how many people are away, not on the kind of holiday | ROBUST (in range); FRAGILE above 5 M | Strengthens "villages fill up" and threatens "holiday-specific" | Main report (one panel) |
| 4 | Two numbers, volume and mean trip length, predict the network shape of each night (R² 0.95–0.99 for r) | ROBUST | Threat to Result 1: the topology shift is low-dimensional | Main report (one sentence + appendix figure) |
| 5 | In/out decoupling: Spearman(k_in, k_out) 0.69 → 0.09 in August; in-hubs become resorts, out-hubs cities | ROBUST | Supports the bipartite framing (R3); explains disassortativity | Main report (Result 1) |
| 6 | Net balance: 99 % of emptied inland zones are net importers in August (+1.39 M persons), 94 % on All Saints, 22 % on an ordinary night | ROBUST | Strengthens the answer; robust to aggregation in sign | Main report (Result 1, replaces k-core headline) |
| 7 | The roots index peaks at Christmas (0.77) and All Saints (0.69); August (0.58–0.62) is no higher than ordinary weekdays (0.59) | ROBUST (descriptive) | **Threat**: August is the least roots-aligned holiday across provinces | Flag to R2; main report discussion |
| 8 | "Roots" and "beach" are separate axes at city level (Spearman −0.04, n = 137) | ROBUST (beach share); FRAGILE (city roots index across years) | Changes the title framing: the dichotomy is false | Main report (discussion) or journal |
| 9 | During Easter and Christmas the structure switches 2–4 days before the volume peaks | ROBUST (one year) | Nuance: long trips leave first; the main night is not the structural extreme | Appendix |
| 10 | August excess distance is bimodal (25–150 km and 300–500 km); family holidays are local (< 50 km) | ROBUST | Explains #2 and #4 mechanically | Appendix / journal |
| 11 | Weighted PageRank in August moves to the north coast (Gijón, Sanxenxo, A Coruña, Vigo) | ROBUST | Minor; shows sinks, not hubs | Appendix |
| 12 | Cross-year and cross-season level shifts in volume (Feb 2022 weekday 1.17 M vs Feb 2024 2.1 M; Oct–Dec weekdays +27 % vs Feb) | ARTEFACT/FRAGILE | Threat to w ≥ 10 replications and to the Feb baseline | Appendix (robustness) + note to all agents |
| 13 | Zone drop-outs: Beriáin_AM has almost no residents in Feb–Mar 2024; Checa_AM, Olombrada_AM partial | ARTEFACT | Minor; contaminates the node-level gain of 3 zones | Note to R1; appendix |
| 14 | Mid-December 2024 (Dec 12–17) has fewer reported pairs and a higher median flow, a sign of a smaller device panel | ARTEFACT (mild) | None for the design (no main night there) | Appendix footnote |
| 15 | Sister nights Feb 13 (Carnival) and Feb 28 (Andalusia Day eve) | ROBUST, harmless | None: their excess is 0.16–0.21 M | Appendix footnote |

---

## 1. Date convention: file date d = the night from d−1 to d (ROBUST)
**Evidence** (`p_date_convention.csv`):
- In ordinary weeks of Feb, Oct and Nov 2024, the median away volume by file weekday is Mon 2.78 M, Tue 2.62, Wed 2.64, Thu 2.69, Fri 2.76, **Sat 3.73, Sun 4.53**. If the file date were the night *following* d, Sunday night would be the busiest night of the week, which is implausible because people return on Sunday. If the file date is the night *ending* on d, "Sat" is Friday night and "Sun" is Saturday night, which fits.
- 2022: All Saints fell on a Tuesday. File 1 Nov is 4.65 M and file 2 Nov is 2.65 M (ordinary). The puente night (Mon→Tue) is high and the holiday night itself is not, because people work on Wednesday.
- 2024: file 25 Dec (6.56 M) is the December peak, which is Christmas Eve night (Nochebuena). File 10 Dec (2.74 M) is back to normal after the Monday-9 regional holiday.

**Consequences.**
- "15 Aug" is the night of 14→15 Aug and "1 Nov" is the night of 31 Oct → 1 Nov, both holiday eves. The main nights stay valid as holiday nights, but the report must say what they are.
- The All Saints sister set 31 Oct / 1 Nov / 2 Nov covers the nights of 30 Oct (ordinary, 2.80 M), the eve (4.23 M) and the holiday night (5.59 M). The ordinary-looking 31 Oct pulls the sister range down, so the representativeness check for All Saints mixes regimes. Suggestion: replace 31 Oct with 3 Nov (5.70 M).
- Confirm with the MITMA methodology PDF (I could not fetch it: 403). **Main report**, methods.

## 2. Two holiday modes; All Saints is a big weekend, not a small August (ROBUST)
**Evidence.** The test compares the destination excess vectors (2,435 zones) of each night by cosine similarity (`p_holiday_similarity.csv`, heatmap `p_diag/p_holiday_similarity.png`).
- **August** nights agree with each other at 0.96–1.00, and Aug 15 2023 vs 2024 gives 0.99.
- **Nov 1** is close to Sun 17 Nov (0.97), Dec 8 (0.96), Christmas (0.93), a February Sunday (0.92) and Nov 1 2022 (0.98). It is far from Aug 15 (0.36).
- Regressing Nov 1 excess on Aug 15 excess gives R² = −0.02. On a November Sunday it gives R² = 0.93.
- **Easter** lies between the two modes: 0.85 with August and 0.64–0.66 with the family mode.
- Pair-level correlations show the same structure: Aug15–Nov01 0.53, Nov01–Sun17Nov 0.88.

**Mechanism** (`p_holiday_features.csv`). The coastal share of the excess is 0.44–0.49 in August, 0.32–0.36 at Easter and 0.22–0.27 for All Saints, Christmas, puentes and Sundays. The emptied-inland share of the excess is similar in every mode: 0.17–0.22 in August, 0.19–0.23 at Easter, 0.16–0.22 for All Saints, 0.15 at Christmas and 0.14–0.15 on Sundays.

**Bearing.** All Saints works as a "no-beach control" only in the sense of "a long weekend". Its similarity to August's inland excess is not evidence that August is about roots. Frame All Saints, Christmas and weekends as one family/second-home mode, and August as that mode plus a beach layer. **Main report.**

## 3. The emptied-village share depends on volume, not on the holiday type (ROBUST in range; FRAGILE when extrapolated)
**Evidence.**
- Across the 180 nights of 2024, log(away persons) alone explains R² = 0.93 of the share of away nights spent in emptied inland zones (`p_volume_vs_structure.csv`). The scatter in `p_diag/p_metric_vs_volume.png` shows one curve shared by Feb weekdays, weekends, Easter, All Saints and August.
- I fitted that curve on non-holiday nights only (R² 0.78, slope +0.071 per log unit, residual sd 0.009) and measured residuals in sd units (`p_one_curve_residuals.csv`): August +1.7 (2022 +1.3, 2023 +1.6), Easter +1.9, All Saints +2.9, Christmas −1.3.
- The coastal share behaves very differently: August +19.5 sd (2022 and 2023 +20.7), Easter +6.6, All Saints −1.7.
- All August nights lie beyond the largest ordinary volume (5.1 M), so the August residual is an extrapolation.

**Bearing.** Villages fill up whenever Spaniards are away, at a share set by how many are away. August is only slightly above that curve. What makes August special is the beach. This strengthens the claim "emptied villages are a structural holiday destination". It weakens "August is a return to roots" as a summer-specific claim. All Saints sits most above the curve and Christmas below it, because Christmas goes to cities. **Main report** (one panel, possibly inside Result 1).

## 4. Two scalars predict the network shape (ROBUST; a threat to over-reading Result 1)
**Evidence** (`p_volume_vs_structure.csv`, 180 nights of 2024, OLS R²):

| Metric | log(volume) alone | Mean trip length alone | Both |
|---|---:|---:|---:|
| w ≥ 10 assortativity r | 0.77 | 0.84 | 0.99 |
| Equal-edges r | 0.78 | 0.77 | 0.95 |
| Share in own province | | | 0.97 |
| Equal-edges reciprocity | | | 0.76 (0.91 with coastal share) |
| Spearman(k_in, k_out) | | | 0.81 (0.93 with coastal share) |

At matched volume, holidays have longer trips and more negative r than Sundays (`p_weekend_vs_holiday_matched.csv`):
- Good Friday (6.95 M): r = −0.190, trip length 198 km.
- Dec 8 (6.15 M): r = −0.165, 160 km.
- Dec 27 vs Nov 17 (about 4.6 M): r = −0.183 vs −0.144.

**Bearing.** The topology shift (disassortativity, rich-club collapse) follows from "more people, travelling further from a few cities". This agrees with the b02 finding that r is forced by the degree sequence. Report r and the rich club as consequences of volume and trip length, not as separate evidence. **Main report**, one sentence; plot in the appendix.

## 5. In/out decoupling and who becomes a hub (ROBUST)
**Evidence** (`p_degree_summary.csv`, `p_hubs_top20.csv`, `p_hubs_composition.csv`).
- **Decoupling.** Spearman(k_in, k_out) across zones is 0.69 on Feb 14 (sisters 0.69), 0.09 on Aug 15 (sisters 0.10; 2022–23 0.08) and 0.63 on Nov 1. The equal-edges backbone gives 0.69 / 0.14 / 0.63.
- **In-hubs.** The coastal share of the top 50 zones by k_in is 0.36 (Feb), 0.80 (Aug, w ≥ 10), 0.88 (Aug, equal edges) and 0.36 (Nov 1). August in-hubs: Benidorm 676, Torrevieja 524, Oropesa 485, Salou 483, Roquetas 459, Gijón 447, Conil 432. Madrid's k_in drops from 1,324 to 995. Top-50 sets are stable: Jaccard with the main night is 0.89 for sisters and 0.80 for years.
- **Out-hubs.** These are the cities and their suburbs: Madrid 2,203 of 2,434 possible targets, Barcelona 1,893, then l'Hospitalet, Móstoles, Alcalá, Getafe and Alcorcón. Gini(k_out) on equal edges rises from 0.60 to 0.84 while Gini(k_in) stays at 0.63, and the CV of k_in halves (2.85 → 1.40).
- **All Saints in-hubs** stay the inland provincial capitals (Madrid 1,482, Sevilla, Granada, Córdoba, Salamanca). Emptied zones never enter a top 50 (share 0.00 on every night).
- Zone-level gain in log(1+k_in) on Aug 15 (`p_centrality_gain_ols.csv`, HC1): coastal +1.03 (se 0.04), emptied +0.49 (0.03), log(#municipalities) +0.31 (0.02), with log population and the February value as controls (n = 2,379). The emptied coefficient survives the aggregation control and is the same on every sister night (0.44–0.49) and in both years (0.47).

**Bearing.** In August the network splits into senders and receivers, close to a bipartite structure. That supports R3's bipartite design, and it is why disassortativity appears (#4). The "hubs" are resorts, not villages, so village promotion must be measured by gain (#5, #6), not by hub rank. **Main report**, Result 1.

## 6. Net importers: emptied villages flip sign (ROBUST, including single-municipality zones)
**Evidence.** Balance = (s_in − s_out)/(s_in + s_out) over away flows (`p_net_balance_by_class.csv`).

| Night | Emptied inland: net importers | Emptied inland: median balance | Emptied inland: net persons | Coastal: net importers | Coastal: median balance | Coastal: net persons |
|---|---:|---:|---:|---:|---:|---:|
| Feb 14 | 22 % | −0.18 | −70 k | | | |
| Aug 15 | 98.9 % | +0.76 | +1.39 M | 84 % | +0.71 | +1.82 M |
| Nov 1 | 94 % | +0.40 | +309 k | | | −4 k |

- Cities over 100k residents are net exporters of −3.22 M on Aug 15.
- Sisters: August 97–99 %, All Saints 62–97 %, ordinary 19–55 %. 2022: 98.5 %, 2023: 96.9 %.
- Single-municipality zones only (n = 76 emptied): August 94.7 %, median +0.52; Feb 14 38 %.

**Bearing.** This is a simple, ratio-based statement of the answer that the zone-aggregation problem (audit C1) affects much less than coreness does: the sign holds in single-municipality zones. Proposed for R1 as a companion to k_in. **Main report.**

## 7. The roots index is highest at Christmas and All Saints, not in August (ROBUST as description; a THREAT)
**Evidence.**
- Daily index (`p_daily_metrics.csv`, `roots_corr_partial`): ordinary Feb weekdays 0.58–0.60, Feb Sundays 0.63; August 0.56–0.62 (Aug 15 0.61); Easter 0.66–0.68; All Saints 0.68–0.69; Christmas 0.72–0.77 (Dec 25 0.768).
- 2022: Feb weekdays 0.54, Aug 0.59, Nov 1 0.64. 2023: Feb 0.58, Aug 0.58.
- Volume and trip length explain only R² 0.10 of the index. Adding the coastal share gives 0.80: the index falls whenever the beach share rises.
- On the excess, per city (`p_cities_rank_stability.csv`, medians over 137 cities): Aug 15 0.37 (sisters 0.34–0.45), Good Friday 0.56, Nov 1 0.52, Dec 25 0.71.

**Bearing.**
- Across provinces, after distance, August's *total* flow is no more aligned with birthplaces than an ordinary Tuesday. Its *excess* is the least aligned of all holidays.
- This does not contradict R2's PPML (origin-province × destination-province fixed effects, identified within provinces), which answers a different question.
- It does mean the honest headline is: "roots show most clearly at Christmas and All Saints; in summer, roots are diluted by beaches."
- **R2 should test whether the August − February difference survives its own design across sisters. If it does not, the main answer must shift to the family holidays.** Main report discussion.

## 8. Roots and beach are separate axes at city level
**Evidence** (`p_cities_rooted_vs_beach.csv`, `p_cities_correlates.csv`).
- **Beach share of the August excess** is very stable: Spearman vs Aug 15 is 0.95–0.99 for sisters and 0.96–0.98 for 2022–23 (ROBUST).
  - Least beach: Soria 0.13, Zamora 0.15, Cuenca 0.18, Salamanca 0.20, León 0.20, Ávila 0.20. Their emptied share is 0.48–0.71 and their own-province share about 0.5: castellanos go to nearby villages.
  - Most beach: Lorca 0.78, Murcia 0.74, Jerez 0.74, Ferrol 0.74, Utrera 0.73, Dos Hermanas 0.73.
- **City roots index** is highest in Barcelona's working-class suburbs: Sant Boi 0.78, Santa Coloma 0.77, Viladecans 0.76, El Prat 0.76, Badalona 0.75, Mataró 0.75, Hospitalet 0.73. Their top birth provinces are Córdoba, Jaén, Badajoz and Granada, the 1950–75 migration corridors. Their beach share is still about 0.40.
- **Correlations across the 137 cities.**
  - Beach share vs roots index: −0.04.
  - Share of residents born in another province: roots index +0.31, beach share −0.38.
- **Stability of the city roots index.** Across August sisters 0.84–0.94 (ROBUST); across years 0.39 (2022) and 0.63 (2023) (FRAGILE); between August and Nov 1 / Christmas about 0 (FRAGILE).

**Bearing.** "Roots or beaches?" is not a single axis. Migrant-origin cities do both. Castilian capitals go "home" because home is near. Good material for the discussion or the journal paper. The city-level roots ranking should not be a main result.

## 9. The switch: structure leads volume (ROBUST in 2024; a single year)
**Evidence** (`p_switch.csv`; frac = share of the eventual peak shift away from the Feb median).
- **Easter.** On the Holy Week Monday file (night 24→25 Mar), r is at 0.68, the coastal share at 0.77 and the own-province share at 0.60 of their peaks, while volume is only at 0.31. Volume peaks on Mar 30 (Good Friday night) at 7.59 M, when r has already peaked and the coastal share is falling (0.87).
- **Christmas.** On Dec 23–24, r is at 0.77–0.91 while volume is at 0.45–0.53. Nochebuena (file 25 Dec) has peak volume but pulls structure back (r at 0.57) because it adds short family trips: trip length 156 km vs 187 km the day before.
- **August.** The regime is already on from Aug 1 (r at 0.70, coastal at 0.87). Each weekend in August *reduces* the structural shift: r is at 0.58–0.72 on Sunday files because weekend local trips are added. The structural extreme falls on the weekdays of Aug 6–14 (r −0.221 on Aug 14), not on the Aug 15 volume peak.
- **All Saints.** The switch happens only on the long-weekend nights. The day after, everything is back within the ordinary range.

**Bearing.** Main nights chosen by volume are not the structural extremes. That is good for the representativeness rule ("the main night is not the extreme"), but the text should not call Aug 15 "the most switched night". **Appendix.**

## 10. Distance of holiday links (ROBUST)
**Evidence** (`p_distance_bins.csv`, `p_diag/p_distance.png`).
- The August excess has two modes: 25–150 km (44 % of the excess) and 300–500 km (22.6 %, 67 % of it to the coast). Aug 15 2023 has the same shape: 45 % and 21.9 %, 68 % coastal.
- The excess of Christmas and Nov 1 falls < 50 km in 36–37 % of cases (August 22 %) and at 25–150 km in 48–50 %. Its 300–500 km band carries only 13 %, and 32 % of that goes to the coast.
- Ratio night / Feb baseline on Aug 15: 6.4× at 300–400 km, 1.9× at 10–25 km.

**Bearing.** This explains #2 and #4. The long-range beach layer belongs to August and Easter. The village layer sits at 25–150 km in every mode. **Appendix / journal.**

## 11. Centrality shifts (ROBUST, minor)
- **Weighted PageRank** (`p_hubs_top20.csv`). Madrid falls from 0.049 (Feb) to 0.012 (Aug). The August top 10 is mostly coastal: Gijón 0.0082, Málaga, Sanxenxo, Barcelona, A Coruña, Vigo, València, Benidorm. Rank correlation with February is 0.56 (`p_centrality_rank_stability.csv`). The northern coast leads because it is a sink that does not send people back.
- **Betweenness** (w ≥ 10, undirected) keeps the same cities in the same order: Spearman with February 0.81, and Madrid stays first. With L ≈ 2, betweenness measures little.

**Bearing.** Appendix at most.

## 12. Volume level shifts across years and seasons (ARTEFACT/FRAGILE)
**Evidence** (`p_daily_metrics.csv`).
- Median ordinary-weekday away volume: Feb 2022 1.17 M, Feb 2023 1.83 M, Feb 2024 2.11 M. Feb 2022 may reflect the Omicron tail or MITMA panel changes.
- The 2024 weekday level also rises between seasons: Feb–Mar 2.11 M, Oct–Dec 2.65–2.71 M (+27 %). Equal-edges reciprocity moves from 0.52 to 0.58–0.59 and the roots index from 0.59 to 0.63–0.64.
- Total resident counts also jump between March and August (+1.3 %, padrón calibration; audit m3).

**Consequences.**
- w ≥ 10 backbones are not comparable across years. Feb 2022 has E = 22,219 vs 35,812 in 2024 (b01 manifest). Use equal edges and shares for year replication.
- Feb is the *lowest* ordinary baseline. An October baseline would shrink every holiday-vs-ordinary contrast somewhat; All Saints vs October is the fairest test.
- All agents should report one robustness row with an Oct/Nov weekday baseline. Appendix.

## 13. Zone drop-outs (ARTEFACT)
`p_artefact_zone_residents.csv` compares each zone's daily resident count (all destinations, including home) with its own median:
- 31905_AM Beriáin is about 105 residents instead of about 5,700 on 60 days, all of Feb–Mar 2024, including all 13 ordinary sister nights.
- 19103_AM Checa has a median of 440 against a padrón of 2,511, and drops 34 days, 9 of them in analysis nights.
- 40149_AM Olombrada drops 10 days, 4 of them ordinary sisters.
- Five more zones drop on 1–4 days, none in analysis nights. In 2022–23, 16 zones have 349 drop-out zone-days.

Otherwise residents / padrón has median 1.009 (IQR 0.993–1.027).

**Action.** R1 should exclude or flag Beriáin, Checa and Olombrada in any February-to-holiday gain, because their February out-flows are near zero by construction.

## 14. Mid-December sampling dip (ARTEFACT, mild)
`p_artefact_days.csv` and `p_artefact_flow_floor.csv` show the pattern on Dec 12–17:
- The number of reported pairs falls (Sun 15 Dec 279 k vs about 360 k on other Sundays).
- The median reported flow rises from 3.8 to 5.0 and the 1st percentile from 1.15 to 1.45.

This is what a smaller device panel expanded with larger factors looks like. National away totals barely move. No design night is affected; it could matter if the daily panel enters a model. Footnote.

## 15. Doubtful ordinary sisters (ROBUST, harmless)
Feb 13 (Carnival) and Feb 28 (eve of Andalusia Day) have excess of only +0.21 M and +0.16 M over the Feb mean (`p_holiday_features.csv`). Their destination patterns do not resemble any holiday (cosine ≤ 0.27). Keep them.

---

## What I would change in the design
1. Methods: state the date convention (#1). Replace 31 Oct by 3 Nov in the All Saints sister set.
2. Result 1: lead with in/out decoupling (#5) and the net-importer flip of emptied villages (#6). Present r and the rich club as consequences of volume and trip length (#4).
3. Results 1–2 framing: two holiday modes (#2). The emptied share follows volume (#3) and the beach is what is specific to August.
4. R2: test August − February on the roots index across sisters (#7). If it is null, the answer becomes "yes, at Christmas and All Saints; in summer, beaches dilute roots".
5. Robustness: an October baseline (#12) and excluding the drop-out zones (#13).

## Files
- Tables:
  - Daily panel and its analyses: `p_daily_metrics.csv`, `p_date_convention.csv`, `p_volume_vs_structure.csv`, `p_one_curve_residuals.csv`, `p_weekend_vs_holiday_matched.csv`, `p_switch.csv`.
  - Hubs and centrality: `p_degree_summary.csv`, `p_hubs_top20.csv`, `p_hubs_composition.csv`, `p_centrality_rank_stability.csv`, `p_centrality_gain_ols.csv`, `p_net_balance_by_class.csv`.
  - Cities: `p_cities_scores_by_night.csv`, `p_cities_rooted_vs_beach.csv`, `p_cities_rank_stability.csv`, `p_cities_correlates.csv`.
  - Holidays and distance: `p_holiday_features.csv`, `p_holiday_similarity.csv`, `p_holiday_allsaints_vs_aug.csv`, `p_distance_bins.csv`.
  - Artefacts: `p_artefact_days.csv`, `p_artefact_zone_residents.csv`, `p_artefact_flow_floor.csv`.
- Caches: `p_zone_daily.parquet`, `p_node_metrics.parquet`, `p_holiday_dest_excess.parquet`.
- Plots: `p_diag/p_daily_timeseries.png`, `p_metric_vs_volume.png`, `p_degree_ccdf_in_out.png`, `p_distance.png`, `p_holiday_similarity.png`.
