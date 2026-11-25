# Audit of legacy question-B code (`explore/p60`–`p100`, `common.py`) and reported numbers

*Audit W1-audit, 2026-11-24. Read-only on `explore/`. Every number quoted as "re-run" comes from a script in `research_notes/audit_scratch/` (listed at the end), run with `.venv/bin/python -I` on the local data in `data/derived/`.*

Scope of reported numbers checked: `research_notes/Preguntas alternativas desde cero/B_robustness_checks.md` (Checks 1–7), `network_layer_A_B.md` Part 2, `reports/alternativas_desde_cero.md`, and the README/HANDOFF headlines that quote them.

## Executive table

| ID | Severity | Finding | Where | Measured impact | Fix for `pipeline/` |
|---|---|---|---|---|---|
| C1 | **Critical** | k-core "promotion of emptied villages" is confounded by MITMA zone aggregation. Emptied zones are mostly `_AM` aggregates of many municipalities, which mechanically have more partner zones. | `p92:298–360` | Extreme-loss small zones: 84 % `_AM`, median 6 municipalities, vs 34 % and 1 for low-loss. OLS Aug coreness on loss rank (+log pop + Feb coreness): **13.6 (se 1.15) → 1.7 (se 1.2)** after adding log(#municipalities); **−4.8 (se 2.2)** in single-municipality zones. The in-degree version survives (log k_in: 1.13, se 0.07). | Control for log(#municipalities) and log(area), or restrict to single-municipality zones; report k_in/strength alongside coreness; check at municipal level (INE 53001). |
| C2 | **Critical** | "k-core ×4.5" is a ratio of quintile means driven by a low February baseline. Quintiles are not size-matched. | `p92:307–335`, network note §2.1 table, README | Absolute Aug−Feb gain: Q1 **+35.0**, Q5 **+33.6**. Q1 (grew, median pop 9,623) still has a higher August core (49.2) than Q5 (median pop 4,932; 43.1). | Report differences and a regression with size controls; never a ratio of means across groups with different baselines. |
| C3 | **Critical** | Window means are computed as `sum(p)/count(DISTINCT date)` per group. A pair seen on 1 of 12 nights is divided by 1, not 12, so sparse pairs are inflated. | `p92:128`, `p84:181,221`, `p93:102,256,597`, `p62–p69`, `p72`, `p82`, `p86`, `p88` (all multi-night windows) | Pair level, all-Spain Feb Tue–Thu: away total **4.27 M vs 2.38 M** correct; pairs ≥10 **48,303 vs 33,618**. Feb-mean network: density 1.09 % → **0.68 %**, r −0.111 → **−0.140**. Mid-Aug mean: density 3.89 % → **3.01 %**, r −0.213 → **−0.229**. Zone-level roots elasticity (`c+z+ttbin`): Feb **0.269 → 0.401**, mid-Aug 0.505 → **0.568**, so the "1.8–2× holiday jump" becomes **≈1.4×**. BiCM input links (w≥5): Feb 48,237 → **23,470**; Mid-Aug 116,039 → **80,069**. City×province PPML: unaffected (0.721 → 0.722). | Divide by the number of nights in the window, with missing = 0, or build one network per night and summarise (median/range). Add a unit test with a pair present on only one night. |
| C4 | **Critical** | The Result 3 evidence ("12 of 16 territories are diaspora archipelagos", "Q_B z = +1,168") does not beat its own null. The null is Poisson around the fitted gravity mean, which ignores overdispersion. | `p93:663–700`, network note §2.2B | Observed non-contiguous communities **12 vs null 11.8 ± 1.2**. Observed diaspora-destination share **35 % < null 80 %**. Pearson dispersion of flows around the gravity null is **17.8** (Poisson = 1), so the Poisson Q null (0.0089 ± 0.0001) is far too narrow. The null is not refitted per draw, and its optimiser gets 3 restarts against 15 for the observed network. Inputs also carry C3. | Use a null that keeps overdispersion and margins (e.g. refit gravity on each draw from a negative-binomial or a strength- and distance-bin-preserving randomisation). Use the same optimiser effort. Test contiguity against that null. Report partition stability. Only roots enrichment (1.36 vs 1.10 ± 0.04) is currently a defensible signal. |
| C5 | **Critical (claim)** | "Disassortativity is structurally forced by the degree sequence" is presented as a holiday/roots finding. It means the opposite: once degrees are known, r carries no extra information about who links to whom. | Network note §2.1 point 1, HANDOFF result 1 | Re-run: r = −0.131 vs rewired −0.140 (Feb 14); −0.214 vs −0.217 (Aug 15). Hubs reach k = 2,233 of 2,435 nodes, so any simple graph with that degree sequence is disassortative. | Report r only as a descriptive consequence of hub growth. Do not use it as evidence of roots. If mixing is to be tested, use a null that fixes degrees and test r_pop or the JDD. |
| M1 | Major | Rich-club ρ is evaluated at fixed k (k > 100). In August that threshold covers ~300 nodes, against ~80 in February, and the rewired top club is already almost complete. | `p92:90–100, 213–220` | By rank (top-50 nodes): ρ **1.23 (Feb 14) → 1.05 (Aug 15)**. The "collapse" is real but smaller than 1.25 → 0.99. Aug top-25 φ_rew = 0.967, so ρ cannot exceed ≈1.03. "1.25" is the Feb-mean network, which carries C3. | Compare rich clubs at matched rank or degree percentile. Report the ceiling 1/φ_rew and the spread across rewirings. |
| M2 | Major | Buggy coastal flag (`buffer(50).intersects(union.boundary)`) still used. | `p73:18`, `p74:5–6` (fig9), `p87:63`, `p92:57` | 157 mainland false positives, 0 false negatives (41 of them are the phantom zones of M3). Inland refill Aug Q1→Q5 **0.90→1.49 (+0.563)** becomes **0.90→1.55 (+0.584)** with `coastal_true`. k-core table changes little: Q5 ×4.55→×4.62; loss-rank OLS 12.2→13.6. | `coastal_true` from `explore/out_p86_coastal_flags.csv`. Document the 13 estuary/ria zones (Bilbao, Barakaldo, Irún, Betanzos…) where an independent coastline (`p100`) disagrees. |
| M3 | Major | The node set "N = 2,484" includes 49 phantom zones: uninhabited shared territories ("Comunero de Ansó y Hecho", "Comunidad de Campoo…") with 11-digit ids. All get `prov = '34'` (Palencia), population NaN, and no flows. | `common.py:13–19`, `nodes_municipios.parquet`, `p92:36–40` (`fillna(100)`), `p93:44` | True mainland N = **2,435**. Aug 15 density 3.31 % → **3.44 %**, Feb 14 0.86 % → **0.89 %**. ⟨k⟩ and C_avg (mode='zero') are understated by about 2 %. In `p92` the phantoms count as pop-100 inland villages in the "small village k_in" row. | Drop ids that never occur in the MITMA data or have length > 8. Assert N and the isolated-node count. |
| M4 | Major | The refill gradient is also partly an aggregation effect. | `p63`, `p73`, `p86`, `p88` | True-inland `lmult_aug ~ loss_rank + lres | prov`: **0.584 (0.050)**. Adding log(#municipalities): **0.251 (0.049)**. Single-municipality zones: **0.242 (0.033)**. The sign holds; the magnitude roughly halves. | Report both. Use INE 53001 (unaggregated municipalities) as the main refill test or as a co-equal one. |
| M5 | Major | Numbers in the notes that no script in `explore/` produces. | `B_robustness_checks.md` Check 4 table, last row; "Key bonus insight" | "City × 2,148 inland zones, `c+z+pairFE+ttbin+p_incq4`": 0.289/0.365/0.524/0.518/0.484/0.492. The `c+z+ttbin_adj` Saturday 0.365 and Good Friday 0.526. `p84` runs only 4 slices and no income term. `p_incq4` exists only in `p82` at province level. Note links point to `~/Proyectos/Spain-Mobility/...`. | Drop these numbers. Every reported number must have a script + output file (plan rule). |
| M6 | Major | Small-sample p-values reported as decisive. | `p82:223` (OLS HC3, n = 19–23 cities), `p82:312–335` (district OLS clustered on 19/23 cities with the regressor constant within cluster) | "p = 1e-35" (Madrid→Extremadura, n = 23) and "p = 1.3e-23" (n = 23) come from normal-approximation HC3 z-tests. With n = 23 and R² = 0.93 the largest possible OLS t (df 20) is ≈16. District "p = 1.7e-9" uses 19 clusters for a city-level regressor. | Use t with df = n − k, HC3 or wild-cluster bootstrap. Report n. Treat these as illustrations, not tests. |
| M7 | Major | QAP: the holiday-vs-weekday contrast is never tested, and the province-level evidence is weaker than the headline. | `p93:174–350` | Province spatial-null r: **Feb 0.573 > Mid-Aug 0.545**. Province partial QAP, Aug excess vs roots: **r = 0.085, p = 0.02** (z 1.9). City-level neighbouring QAP: Feb 0.330 vs Mid-Aug 0.559, but no test of the difference. ~120 one-sided tests with p-floor 0.0005. | Test Δr (holiday − weekday) with a paired within-province permutation or bootstrap. Fix the family of tests in advance and report it. |
| M8 | Major | Comparisons across networks with different node sets. | Network note §2.1 ("moving toward the Roots network, r = −0.268") | The roots r comes from a 47-node province graph, thresholded at ≥2,000 persons, density 35.8 %. It cannot be compared with r on 2,435 zones. | Compare networks only on the same node set and threshold logic, or say that they are not comparable. |
| m1 | Minor | Rewiring is undirected only. Directed r(k_out,k_in) and reciprocity get no null. L uses 8 of 25 rewirings. The text says "25 degree-preserving swaps (10E steps)"; the code makes 25 rewired graphs with 10E attempted swaps each. | `p92:200–212` | igraph `rewire` keeps degrees exactly (asserted in the re-run). Edge overlap with the original after 10E swaps: 16 % (Feb 14) and 26 % (Aug 15), so mixing is adequate. | Directed degree-preserving rewiring for directed metrics. Fix the wording. Log the acceptance rate. |
| m2 | Minor | Historical-loss construction. | `p73:4–11`, `p86:49–58`, `p92:43–52` | 188 mainland zones have constituent municipalities missing some census years (`sum(min_count=1)` then undercounts the peak). 105 mainland zones have no loss (municipalities created later, codes 9xx). INE 03903 is split across two zones, so its history is counted twice. | Build loss on a constant-territory basis (INE alteration history); merge segregated municipalities back to the parent for 1900–1960. |
| m3 | Minor | Inconsistent date windows and outage handling. | `p63:7`, `p88:79–83` vs `p72`, `p92`; `p100` | Feb Tue–Thu starts on Feb 1 (13 nights) in `p63`/`p88` but Feb 5 (12 nights) elsewhere. Mid-Aug is 10–18 in most scripts and 12–18 in `p100`. "All Saints 2023" is Sun 29 Oct, an ordinary Sunday, not a holiday (outage 30 Oct–5 Nov 2023 locally). 2022 has only 29 Oct–2 Nov. Nov 2024 lacks the 9th, 10th and 26th. Padrón calibration jumps between March and August 2024. | One date table in `pipeline/config.py`. Sister nights as in the design. Drop 2023 All Saints instead of substituting a Sunday. |
| m4 | Minor | City sets differ between results that are compared side by side: 151 (with 14 island/Ceuta/Melilla cities and island destination provinces, `p72`, `p82`, `p86`), 137 mainland (`p84`, `p93`; the `p93` docstring says 142). 25 of 47 mainland provinces have a single city and add nothing to identification under `pairFE`. | `p72:8`, `p82:30`, `p93:241` | Mixes 0.717 (151 cities) with 0.644 (137 cities) in Check 4. | One city set (137 mainland) for the main analysis; report identifying provinces (22 with ≥2 cities, 112 cities). |
| m5 | Minor | `p93` QAP details. "Freedman–Lane" is in fact Dekker's double semi-partialling (permute X residuals). That is valid, but should be named correctly. In the full bipartite QAP, independent row/column permutation moves the masked own-province zeros into the mask. Residual `log((W+1)/(E+1))` depends on scale. The gravity fit has no assertion (it does converge: relative margin error ≤ 5e-7). | `p93:66–96, 193–208, 312–328` | Within-province permutation is valid: the mask is invariant and 112 cities in 22 provinces are exchangeable. | Name the method correctly; use a permutation-invariant mask; assert convergence as `p100` does. |
| m6 | Minor | `p94` tail test: the lognormal tail is fitted with sample moments, not a truncated MLE, and the power law is continuous on discrete degrees, so the Vuong test is not valid. (Mostly question A.) | `p94:118–135` | — | Use `powerlaw` (Clauset et al.) or drop. |
| m7 | Minor | Specification search and wording. The inland/strict/10-km/28-province variants, 6 slices × 8 PPML specs, and "8.5/10" scores read as stronger than they are. "Removing coast strengthens the effect" changes the outcome population, so it does not test the confound. | Notes | — | Pre-specify one main specification and a short robustness list. |

**Checked and found correct:** residence-district → municipal-zone mapping (3,792 districts → 2,618 zones, many-to-one, no district maps to two zones, daily totals preserved to 1e-6); self-loops removed after mapping (`r <> s`); no foreign zones in the overnight data; ADRH income matched for 151/151 cities (no silent sample drop between specs); city×province PPML elasticity reproduces (`c+pairFE` mid-Aug 0.721, se 0.023); IPF gravity null converges in 20–26 iterations; BiCM solver reaches max degree error < 1e-6; `p94` joint-degree-preserving swap logic is correct; within-province permutation QAP is correctly implemented.

## Five fixes that matter most for the clean pipeline
1. **Window means** (C3): divide by nights in the window, with missing = 0, or work per night. Everything pair-level built on Feb Tue–Thu or the mid-Aug mean must be rebuilt: mean-network properties, zone-level PPML (the "1.8×" becomes ≈1.4×), and all bipartite territories.
2. **Zone aggregation** (C1, M4): add log(#municipalities) and log(area) as controls, or analyse single-municipality zones and INE 53001. The k-core promotion result does not survive this control; the k_in and refill results do, at about half the size.
3. **Result 3 null** (C4): replace the Poisson MC with an overdispersion-preserving, refitted null and the same optimiser effort. Until then, "12/16 diaspora archipelagos" is not evidence (null 11.8 ± 1.2).
4. **Node set and coastal flag** (M2, M3): N = 2,435 (drop the 49 phantom zones); `coastal_true` everywhere; assert both in tests.
5. **Claims discipline** (C2, C5, M1, M5–M8): no ratio-of-means headlines; disassortativity is "explained by degrees", not "evidence of roots"; rank-matched rich club with its ceiling; drop untraceable numbers; small-n p-values with t/bootstrap; a direct test of the holiday − weekday difference in QAP.

---

## Detail

### (a) Coastal flag
- **Buggy method** (`buffer(50).intersects(union_all().boundary)`): `p73:18`, `p74:5` (it `exec`s `p73`, so `reports/figs/fig9_exodus_holidays.png` panel a uses it), `p87:63`, `p92:57`. `p86:22` recomputes it only as a comparison column. `p82`, `p82b`, `p86`, `p88` and `p89` use `coastal_true`. `p100` uses its own flag (within 200 m of an external coastline).
- **Counts on the 2,484 "mainland" ids:** `coastal_p73` = 493, `coastal_true` = 336. All 157 differences are p73 false positives; there are no false negatives. 41 of the 157 are the phantom zones of M3 (that is the "42 in Palencia" in Check 6). By province: Palencia 42, Navarra 18, Burgos 12, Ourense 9, Girona 8, Badajoz 7, Cáceres 7.
- **`coastal_true` against `p100`'s coastline:** they agree on 2,420 of 2,435 zones. 13 zones are coastal for `p100` but not `p86`, all estuary or ria zones (Bilbao, Barakaldo, Erandio, Irún, Betanzos, Cambre, Culleredo, Gondomar, Teià…). 2 zones go the other way (Colindres, Bárcena de Cicero). This is a definitional choice to document, not a bug.
- **Numbers that depend on the buggy flag:**
  - `alternativas_desde_cero.md`: "×1,48 agosto / ×1,11 Todos los Santos / ×0,90" and "+0,55 log (ee 0,05)". Corrected (from `out_p86`): 0.90→1.55, 1.00→1.13, +0.584 (0.050).
  - fig9 panel a.
  - Check 7 quintile medians (`p87`).
  - Network note §2.1: k-core quintile table, size-matched block, "small inland village k_in" row.
  - Check 6 row "Inland (p73 flag=0)", which is labelled as such.
- With `coastal_true` the k-core table moves only slightly (Q5 ×4.62, n = 411 per quintile). The table's real problems are C1 and C2.

### (b) Node set
- `nodes()` and `nodes_municipios.parquet` contain 1,528 five-digit municipal zones, 907 `_AM` zones and **49 eleven-digit "comunidad/comunero" zones**. The eleven-digit codes start with `34`, so `prov = '34'` (Palencia) although they lie in Huesca, Cantabria, Burgos, etc. None of them appears in the MITMA overnight data. With them, "N = 2,484" contains 49 isolated nodes. `p92:40` gives them population 100, and the only reason they stay out of the quintile tables is that their loss is NaN.
- Mainland filter: `~foreign & ~island & prov ∉ {51,52}` is used the same way in `p65`, `p86`, `p92`, `p93` and `p94`. `p72` and `p82` use all of `N` for cities and province centroids, which brings in islands and Ceuta/Melilla (14 of 151 cities; island destination provinces). `p100` keeps every zone, islands included.
- `_AM` zones are treated as villages everywhere. They are aggregates: in the extreme-loss tercile of small inland zones the mean is 10 municipalities per zone (C1, M4).

### (c) Residence code → zone and self-loops
`p61` joins the district-level residence code to `municipio_mitma` (3,792 → 2,618, each district maps to one zone). Daily totals match to 1e-6 (Feb 14: 47,433,299; Aug 15: 48,047,683). All `s` codes are municipal zones. After mapping, stays in another district of the same municipality become `r = s` and are dropped, as intended. Inside an `_AM` zone, trips between its own villages are invisible. That makes the refill estimates conservative, but it also means "village" results are really results about zones.

### (d) Dates
- The C3 bug applies to every multi-night window computed at a level where groups are sometimes missing. Distribution of nights present for mid-Aug pairs (9 nights): 396k pairs on 1 night, 157k on all 9. At province or city×province level almost every group is present every night, so the bias is negligible there (PPML 0.721 vs 0.722).
- Visible symptom in the network note: the "Mid-Aug mean" column has more edges and higher density (3.74 %) than the peak night of Aug 15 (3.31 %). The "Feb Wed mean" moves 3.79 M away-nights against 2.11 M on Feb 14.
- 2024 data: Feb 29/29 nights, Mar 31, Aug 31, Oct 31, Nov 27 (9, 10 and 26 Nov missing), Dec 31. 2022–23 data: Feb and Aug complete; 2022 Oct 29–Nov 2; 2023 only Oct 29 and Dec 24.

### (e) Threshold and density
w ≥ 10 expanded persons per directed pair, then the undirected union (`p92:166–176`). With N = 2,435 the densities are: Feb 14 0.89 %, Aug 15 3.44 %, Mid-Aug mean 3.01 % (corrected), Feb mean 0.68 % (corrected). All stay below 5 %. The "backbone captures 64–83 % of flow" figures for the mean networks carry C3. The threshold is applied to directed weights and the graph is then made undirected, so one direction ≥ 10 is enough for an edge. This is fine, but it should be stated.

### (f) Rewiring
`Graph.rewire(n=10E)` in igraph's default simple mode keeps the degree sequence exactly (checked with `assert` in `a07`) and allows no multi-edges or loops. After 10E attempted swaps, 16–26 % of the original edges remain; for graphs with hubs of degree 1,500–2,233 that is close to the expected random overlap. There is no directed null, and the spatial embedding is not preserved.

### (g) Rich club
φ(k) uses `deg > k` and ρ = φ_obs / mean(φ_rew), which is the standard definition. The problems are the fixed-k comparison across networks whose ⟨k⟩ quadruples and the saturation ceiling (M1). Re-run, ρ by top-N:

| Top-N | Feb 14 ρ | Aug 15 ρ | Aug 15 φ_rew |
|---|---:|---:|---:|
| 10 | 1.08 | 1.02 | 0.98 |
| 25 | 1.18 | 1.01 | 0.97 |
| 50 | 1.23 | 1.05 | 0.92 |
| 100 | 1.19 | 1.05 | 0.81 |
| 200 | 1.10 | 1.02 | 0.59 |

### (h) k-core and quintile stratification
- **Loss definition:** 1 − pop2021 / max(pop1900…1960) per zone, after summing constituent municipalities (m2). Quintiles are taken over inland zones with pop ≥ 200 and non-missing loss.
- **Size matching:** quintile medians run from 9,623 (Q1) to 4,932 (Q5), so the quintiles are not matched. The "size-matched" block (pop < 6,000, terciles) does have similar median populations (4,309 / 4,418 / 4,811), but it is matched on population, not on the number of municipalities (1.8 vs 10.2 on average).
- **Label bug:** `p92:352` prints "pop < 2,000" although the code uses < 6,000.

Re-run (`a08`, true coast, N = 2,435), OLS of Aug 15 coreness on loss rank:

| Controls | b(loss_rank) | se |
|---|---:|---:|
| log pop, Feb coreness | 13.64 | 1.15 |
| + log #municipalities | 1.66 | 1.18 |
| + province FE | 2.16 | 1.70 |
| single-municipality zones (n = 1,186), + province FE | −4.78 | 2.17 |
| log(1+k_in) instead of coreness, all controls | 1.127 | 0.069 |

Controlling for the number of municipalities may over-control, because aggregation is itself a consequence of depopulation. The defensible statement is narrower: at MITMA zone resolution, the k-core result cannot be separated from zone aggregation, while the in-degree result can.

### (i) QAP and partial QAP
- **Province QAP** (`p93:174–191`): simultaneous row/column permutation of X with an off-diagonal mask. Correct.
- **"Partial QAP"** (`p93:193–208`): both matrices are residualised on the Feb residual and the X residuals are permuted. That is Dekker's double semi-partialling, not Freedman–Lane. It is valid.
- **Neighbouring-city QAP** (`p93:288–350`):
  - Gravity residuals are demeaned within origin province × destination province.
  - Cities are permuted within their province. Only provinces with ≥ 2 cities enter (22 provinces, 112 cities), and the mask is invariant under these permutations. This is correct.
  - It answers "within a province, do cities with more natives of p send more (gravity-residual) stays to p?"
  - The headline "partial QAP r = 0.49" is the Mid-Aug full flow (not the excess) residualised on Feb; the excess version gives 0.345.
- **Full bipartite QAP:** rows and columns are permuted independently while the mask (own province) moves. This is a minor bias in a secondary number.
- **Missing test:** none of the tests compares holiday r with weekday r (M7).

### (j) Gravity null
Doubly constrained IPF with a distance-bin deterrence function (20 quantile bins for city×province; fixed bins 2–25 km wide for city×zone). Converged in the re-run: relative margin errors are 2e-10 (rows), 2e-16 (columns) and 5e-7 (bins) at city×province, with similar values at city×zone. `p93` checks only the change in f and never asserts the margins; `p100` asserts them. The gravity specification is adequate. What is weak is the stochastic model placed on top of it (C4).

### (k) PPML, fixed effects and income
- **Model:** `w ~ log1p(born) | c + pairFE` with CRV1 by city. pairFE is origin province × destination province, so identification comes from cities sharing a province: 22 mainland provinces plus the islands.
- **Income:** `p[linc]` is a destination-specific slope on city log income, which is correct in pyfixest. The ADRH join is complete for 151/151 cities.
- **Reproduction:** `c+pairFE` re-runs to 0.721/0.722 (C3 negligible at this level) and the Feb value to 0.633/0.644.
- **Zone level** is the exception: correct Feb 0.401, Feb 14 0.386–0.412, mid-Aug 0.555–0.568, Aug 15 0.566–0.574. The holiday excess is about 0.16 (≈1.4×), not 0.23 (1.8–2×).
- **Other problems:** M5 (no code for the zone-level income spec) and M6 (small-n p-values in the metro OLS).

### (l) BiCM, Barber modularity, contiguity
- **BiCM:** fixed-point iteration converges to a max degree error of 9.9e-7 in all three slices. Q_B is normalised by the link count.
- **Optimiser:** the bipartite Louvain is a home-made greedy algorithm. Partition stability across restarts is never reported.
- **Inputs:** w ≥ 5 on window means, so the Feb and Mid-Aug inputs carry C3 (links ×2.06 and ×1.45).
- **BiCM vs weekday:** Mid-Aug BiCM Q_B (0.159) is *lower* than Feb (0.268).
- **Contiguity:** components are computed on exact polygon adjacency, which is fine. The test against the null fails (C4).
- **Roots enrichment:** 1.36× under the spatial null, against 2.8–3.3× under BiCM/Barber. Most raw enrichment is therefore spatial.

### (m) Other claims and leakage
- **No leakage found.** Roots are Census stocks; flows are MITMA nights. The "excess" matrices subtract the Feb mean, which carries C3 at pair level, so Aug/Nov/Xmas excess at city×zone or pair level is biased. At city×province the bias is negligible.
- **Claims stronger than the evidence:**
  - "rich club dissolves completely" (M1)
  - "structurally forced … moving toward the Roots network" (C5, M8)
  - "z = +1,168", "12/16 archipelagos" (C4)
  - "k-core ×4.5" (C1, C2)
  - "p = 1e-35" (M6)
  - "1.8–2× holiday jump" (C3, M5)
  - "Beach tourism diluted the roots signal", based on restricting the outcome to inland destinations (m7)

## Status of headline numbers

| Headline (where quoted) | Status |
|---|---|
| Roots elasticity 0.72 mid-Aug, 0.63 Feb, city×province `c+pairFE` (README, alternativas, Check 1) | **Trustworthy**: reproduced (0.722 / 0.644 with correct means). The weekday baseline is already high; report the difference. |
| Income-robust 0.717 / 0.694 (Check 1A) | Likely trustworthy (same data level, complete ADRH join). Rebuild on one city set. |
| Multi-year 0.733 / 0.733 / 0.717 (Check 2) | Likely trustworthy at city×province. The 2023 "All Saints" column is a Sunday: drop it. |
| Zone-level 0.29 → 0.52, "1.8–2×" (Check 4) | **Wrong**: ≈0.40 → 0.57, ≈1.4×. Income row has no code. |
| Refill ×1.48 / +0.55 (alternativas), +0.584 (Checks 1, 2, 6) | Direction holds. ×1.48 → ×1.55 with the correct coast. Halves (≈0.25) with an aggregation control. |
| Assortativity −0.13 → −0.21 (README) | Single-night values reproduce (−0.131, −0.214). Not evidence of roots (C5). Feb-mean −0.111 is wrong (−0.140). |
| Rich club 1.25 → 0.99 | **Rebuild**: 1.25 is the bugged Feb mean; by rank 1.23 → 1.05 with a ceiling ≈1.03. |
| k-core ×4.5 (README, network note) | **Not supported** as stated (C1, C2). Report k_in instead. |
| Partial QAP r = 0.49, neighbouring cities (README) | Implementation correct. Inputs aggregated to province, so C3 has little effect. Needs a Δ(holiday − weekday) test. |
| 12/16 diaspora archipelagos; Q_B z = 1,168 (README, network note) | **Not supported** (C4). Inputs carry C3. |
| Density 0.86 % / 3.31 % (design) | 0.89 % / 3.44 % with N = 2,435. |

## Scratch scripts (all in `research_notes/audit_scratch/`)
| Script | What it checks |
|---|---|
| `a01_mapping.py` | District → zone mapping, totals, months covered |
| `a02_dates.py` | Missing days, pair-presence distribution, size of the window-mean bug |
| `a03_coastal.py` | p73 vs p86 vs p100 coastal flags |
| `a04_hist.py`, `a05_nodes.py`, `a06_phantom.py` | Historical-loss join, nodes without population/loss, phantom zones |
| `a07_topology.py` | Feb 14 / Aug 15 networks with N = 2,484 vs 2,435; rewiring checks; corrected mean networks; k-core table with both coastal flags |
| `a08_richclub_rank_AM.py` | Rank-based rich club; k-core vs zone aggregation |
| `a09_ppml_qap_means.py` | City×province PPML, buggy vs correct means |
| `a10_zone_level.py` | City×zone PPML, buggy vs correct means and single nights |
| `a11_gravity_conv.py` | Gravity IPF convergence; overdispersion against the Poisson null |
| `a12_misc.py` | ADRH coverage, city sets, 2022–23 date coverage |
| `a13_refill_AM.py` | Refill gradient vs zone aggregation |
| `a14_bip_input.py` | Effect of the window-mean bug on BiCM inputs |
