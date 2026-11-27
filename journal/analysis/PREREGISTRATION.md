# Pre-registration of the journal analyses

**Frozen on 27 November 2026, before any of the newly downloaded MITMA nights (2022–2024 files not present in the coursework cache) has been opened or analysed.** The commit that adds this file is the timestamp. Machine-readable twin: [`design.yml`](design.yml). Changes after this commit go to the "Deviations" section at the end, with date and reason; the original text stays.

Author: sole author. Target: EPJ Data Science. Plan: [`../research/epjds_plan.md`](../research/epjds_plan.md).

## 0. What has already been seen

- **MITMA nights already analysed (discovery set):** the 306 dates in the coursework cache: Feb and Aug 2022 and 2023 (full months), a few Oct/Nov/Dec dates of 2022 and 2023, and Feb, Mar, Aug, Oct, Nov (27 days), Dec 2024. The exact list is the set of files in the main checkout's `_local/external/mitma_pernoct/` on 27 Nov 2026 (308 entries including 2 text files).
- **Confirmation set:** every other standard MITMA municipal daily file dated 2022-01-01 to 2024-12-31 that passes the quality rules of §2 (about 780 files listed). None of these has been opened.
- **INE municipal OD (2019-07 to 2026-09) and ECEPOV 2021** have been explored in full (N1 and the design pilot). Analysis C therefore has **no held-out data**; its safeguards are the fixed gate of §6, not a hold-out.
- **Census cohort results** (pre-1981, village-born, recent) were estimated on 2024 nights and on the 2022 replication nights; the confirmation nights test them out of sample.

## 1. Definitions fixed in advance

- **Night.** File dated *d* = night ending on the morning of *d* (evening *d*−1). All calendar rules below refer to the **evening** date *e* = *d* − 1.
- **Nodes and flows.** 2,435 mainland zones (`pipeline/config.py`); artefact zones `31905_AM`, `19103_AM`, `40149_AM` excluded. Flows are full weights (no w ≥ 10 threshold) from residents of zone *i* sleeping in zone *j* ≠ *i*.
- **Roots exposure.** INE Census (Censo anual) table 68545, residents of city *c* by province of birth, **1 January 2024** for all nights (main); year-matched table year as sensitivity. Share = born in *p* / all residents of *c*.
- **City sets.** (i) Pooled model: the 137 mainland census cities, of which 112 in the 22 provinces with ≥ 2 cities identify the slope (as R2). (ii) Cohort models: the 675 single-zone municipalities > 10,000 residents of the 2024 Census 10 % microdata (as V1), +0.5 added to group counts (Table A2 rules as sensitivity).
- **Pooled model (the R2 "all, main" specification).** PPML of persons from city *c* to mainland zones of province *p* ≠ prov(*c*) on log roots share + log km + (log km)² + 46 income × destination slopes (ADRH 2024), with city FE and origin-province × destination-province FE. Distance = city to population-weighted centroid of *p*.
- **Cohort models (the V1 specification).** Same structure on the 675-city set with group-specific shares: (a) pre-1981 arrivals and 2016–2023 arrivals jointly ("cohort_joint"); (b) born in a municipality ≤ 2,000 and born in a municipality > 10,000 jointly. Arrival codes: only explicit codes enter a group; special or missing codes are dropped from both (audit correction).
- **Inland / coast.** Destination zone flag `coastal_true` (never the buffer method).

## 2. Night classes and quality rules

**Quality rules (applied before any estimate is looked at):** a night is excluded if (a) it is on the ministry's incident or outage list, (b) it falls in 29–31 Oct 2024 (DANA), (c) fewer than 99 % of the median number of mainland residence zones report residents, or (d) total mainland persons (home + away) deviate by more than 10 % from the 28-night centred rolling median. Exclusions are listed with reasons in the night manifest.

**Calendar.** Official national public holidays 2022–2024 and the regional holidays of the mainland autonomous communities from the annual BOE resolutions on the labour calendar. A day is a *regional holiday day* if it is a holiday in at least one of Madrid, Cataluña, Andalucía, Comunidad Valenciana or País Vasco (the five largest origin regions).

**Ordinary reference nights (O).** Evenings Monday, Tuesday or Wednesday that are (i) not a national or regional holiday day, nor the day before or after one; (ii) outside 24 June–10 September; (iii) outside 20 December–7 January; (iv) outside Holy Week (Friday before Palm Sunday to Easter Monday); (v) not excluded by the quality rules. Matched to holiday nights **within the same calendar year**.

**Holiday classes.**
- **AUG (primary):** files dated 10–18 August (evenings 9–17 August), 9 nights per year.
- **XMAS:** file 25 December (Christmas Eve).
- **EASTER:** file of Good Friday (evening of Holy Thursday): 2022-04-15, 2023-04-07, 2024-03-29.
- **ALLSAINTS:** file 1 November (evening 31 October), where available (2023 is missing).
- Descriptive only (no test): every other night, including all of July–August, New Year's Eve, long weekends, weekends.

## 3. Primary hypotheses (family P, Holm-corrected at α = 0.05)

- **P1 — August holiday premium (pooled).** The roots elasticity is larger on AUG nights than on ordinary nights of the same year: δ_AUG = β(AUG) − β(O) > 0.
- **P2 — Exodus cohort.** In the cohort_joint model, the pre-1981 term increases on AUG nights relative to O nights: δ_pre1981 > 0.
- **P3 — Inland.** δ_AUG is larger for flows to inland zones than to coastal zones (roots × AUG × inland interaction > 0 on common city–province support).

**Estimation.** Stacked PPML over all AUG and O nights of the confirmation set, with every fixed effect and covariate interacted with the night (night × city FE, night × origin-province × destination-province FE, night-specific distance and income slopes), a common roots slope and a roots × AUG term. Year-specific contrasts are averaged with equal weight across years.

**Inference.** Wild cluster restricted score bootstrap (Kline & Santos 2012), Webb weights, 9,999 draws, clustered by **origin province jointly across all stacked nights** (22 clusters pooled, 44 cohort). Two-sided p-values; Holm across P1–P3.

**Required diagnostic, always reported (not a test).** Placebo windows: the same contrast for every block of 9 consecutive ordinary-eligible nights in the same year, giving the rank of AUG among placebo windows.

**Decision rule.** P1–P3 are judged on the **confirmation set**. Discovery-set and all-nights estimates are reported next to them. A claim enters the abstract only if it holds on the confirmation set after Holm correction.

## 4. Secondary hypotheses (family S, Holm within family)

- **S1–S3.** δ > 0 for XMAS, EASTER, ALLSAINTS (pooled model).
- **S4.** Village-born (≤ 2,000) term increases on AUG nights.
- **S5.** Recent-arrival (2016–2023) term decreases on AUG nights.
- **S6.** δ_pre1981 − δ_recent > 0 on AUG nights.

## 5. Descriptive outputs (no tests)

- Nightly series β_t for every valid night 2022–2024, pooled and cohort terms, with per-night 95 % intervals.
- Village receiving end: in the most emptied fifth of zones (loss rank as in R1), share of zones with inflow > outflow, and arrivals and departures per resident, for every night.
- Network structure (density, in/out degree correlation, reciprocity, assortativity vs degree-preserving rewiring) on the main nights only, for the SI.

## 6. Analysis C (municipal, parents' birthplaces): gate fixed now

Frame: the 137 ECEPOV mainland cities × mainland municipalities with ≤ 2,000 residents (2021) in other provinces. Exposure: survey-weighted share of city-born residents' parent slots born in the village (strict second generation), jointly with first-generation share. Outcome: INE municipal monthly interprovincial trips (published cells ≥ 30).

The parental result enters the main text only if **all** hold:
1. The August − February and August − October changes in publication probability associated with a parental tie are positive with 95 % intervals above zero on the fixed frame.
2. The count-margin change is positive under at least **two** treatments of suppressed cells (e.g. interval/censored model with the documented threshold, and bounds with suppressed = 0 and = 29), not only on always-published pairs.
3. Intervals propagate ECEPOV sampling uncertainty (household resampling) and cluster by origin and destination province.
4. Real municipal geometry (GISCO LAU 2021) is used for distance and coast.

Otherwise it is reported in the SI as a measurement limitation and the paper stays at city–province level.

## 7. Robustness set (reported, not used to choose results)

Year-matched roots and income; road travel time instead of km where available; dropping Madrid and Barcelona origins; leave-one-destination-province-out; inland-province-only destinations; zero-cell rules (Table A2); seasonal reference (June and September evenings only) instead of all ordinary evenings; equal weighting by night vs by flow.

## Deviations

*(none yet)*
