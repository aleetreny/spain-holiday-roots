# R3: summer territories (mesoscale), with a fair null

*Agent R3, wave 2, 2026-11-24. Code: [pipeline/r3_territories.py](pipeline/r3_territories.py). Every number below comes from a `results/r3_*.csv` file named next to it. Logs: `r3_run.log`, `r3_run_addnights.log`, `r3_summary.log`. Diagnostic plots: `results/r3_diag/`.*

## Verdict: No. The "summer territories / diaspora archipelagos" result does not survive.

1. **Community structure beyond gravity exists on every night, and it is strongest in August.** It is not diaspora-shaped, though. Under a refitted, overdispersion-preserving gravity null, the observed bipartite modularity beats the null on all 27 nights. Aug 15: Q = 0.170 vs null 0.115 ± 0.002. Feb 14: 0.141 vs 0.106 ± 0.002. The August excess is the largest: median over sister nights +0.055 vs +0.030 on ordinary nights, and summer is higher in 117 of 117 summer × ordinary pairs.
2. **Diaspora signatures run the wrong way.** The observed partitions are *more* contiguous than the null and send *less* flow to detached pieces, on every night.
   - Share of territories that are non-contiguous, Aug 15: 0.71 observed vs null 0.99 ± 0.03.
   - Diaspora flow share, Aug 15: 0.48 vs 0.70 ± 0.03.
   - Feb 14: 0.28 vs 0.52 ± 0.03.

   Noise alone scatters a community's destinations, so non-contiguity cannot be evidence of diaspora under any noise null. This confirms audit C4 with a fair null: legacy "12 of 16" vs null 11.8.
3. **Roots enrichment beats the null, but it is not a holiday effect.** Inside gravity-residual communities, city–destination pairs have 1.4–1.8× the birthplace share of pairs across communities (null ≈ 1.07–1.12).
   - Aug 15: 1.41 vs 1.09 ± 0.04.
   - Feb 14: 1.76 vs 1.07 ± 0.04.
   - Sister medians: ordinary 1.48, summer 1.47, All Saints 1.46, Christmas Eve 1.60.

   The summer excess is larger than the ordinary excess in only 47 of 117 pairs. Most of the enrichment is shared by cities of the same province: under the within-province permutation, the August median z is only 2.2 (main night 1.55, p = 0.06).
4. **The BiCM partitions are geography.** They form 5–7 contiguous regional blocks: NMI with province 0.49–0.60, diaspora share 1–8% vs null 93%. Their roots enrichment (2.8–3.4×) disappears under the within-province permutation (z −0.5 to 0.8).

**Answer to R3 as asked.** Holiday communities do not group cities with distant villages beyond what geography and noise produce. They are not more roots-aligned on holiday than on an ordinary night. A mesoscale "summer territories" result should not appear in the report.

**Recommended third result (tested below).** A per-city **two-axis decomposition of the holiday ego-network** against the same refitted NB gravity null.
- On holiday nights, a city's excess flows line up more with its own birthplace profile than on ordinary nights: median roots axis 0.45 (August) vs 0.36 (ordinary).
- This holds in every one of the 117 summer × ordinary pairs, with 59–81% of cities rising.
- It is not a volume artefact: August thinned to February volume gives 0.446.
- The roots axis is independent of the beach axis: Spearman −0.16 to 0.16, never significant at 5%.

## 1. Design

| Item | Choice | Where |
|---|---|---|
| Rows | 137 mainland Census cities (INE 68545 ∩ N = 2,435). Recorded in decisions.md because R2 had not fixed a set yet | `load_inputs` |
| Columns | every mainland zone with inflow > 0 from those cities that night (2,431–2,434); own zone masked | `build_context` |
| Weights | persons on a **single night** (b01 full edges = b00 loaders). No window means | `night_matrix` |
| Nights | 13 ordinary (Feb Tue–Thu), 9 summer (Aug 10–18), All Saints 11-01/02/03 (new convention); appendix: Christmas Eve 12-25, Good Friday 03-29. 10-31 was run under the old convention and dropped | `night_table` |
| Modularity | Barber bipartite Q = (1/W) Σ (O − E) δ, maximised with Leiden. Uses CPM with γ = 0 on a signed 2-layer multiplex, so the quality equals Σ B δ exactly | `leiden_bipartite` |
| Effort | 10 seeded restarts and run to convergence, **identical** for observed and null networks; best Q kept | `N_RESTARTS` |
| Gravity null | 3-way IPF (rows, columns, distance bins), error < 1e-6 asserted (all nights 1e-6) | `fit_gravity` |
| Fair gravity null | NB draws, mean E, **empirical variance function** V(μ) from the binned observed (O − E)², gravity **refitted per draw** | `variance_function`, `draw_nb` |
| BiCM null | binary links w ≥ 10, BiCM with mask and saturated nodes fixed; Bernoulli draws, BiCM **refitted per draw** | `fit_bicm` |
| Draws | 100 per family on the 5 main/appendix nights, 20 per family on each sister night (1,934 clusterings kept; 42 more on the dropped 10-31) | `r3_draws*.csv` |

**Does the null keep the dispersion?** (`r3_null_diagnostics_all.csv`, `r3_diag/variance_function.png`)
- Pearson dispersion, observed vs one draw: Aug 15 23.5 vs 21.8; Feb 14 6.06 vs 5.95; Nov 1 9.10 vs 8.75.
- Binned variance, observed / draw: median 0.99–1.11, IQR within [0.90, 1.83] on all nights.
- Weak point: draws have fewer exact zeros on low-volume nights (Feb 14: 0.84 obs vs 0.77). These extra zeros are tiny cells and carry almost no weight in Q.
- The variance is calibrated on all residual variation, including real structure, so the null is conservative for Q.

**Statistics.**
- A *territory* is a community with ≥ 1 city and ≥ 5 zones.
- A territory is *non-contiguous* if it has ≥ 2 polygon pieces of ≥ 3 zones (gpkg `intersects` adjacency, 7,026 pairs).
- The *diaspora share* is the share of within-territory flow (persons, or links for BiCM) that goes to pieces neither containing nor touching a city of the territory.
- *Roots enrichment* is the mean of r[c, prov(d)] over cross-province pairs in the same community, divided by the same mean over pairs in different communities. r = Census 2024 share of city residents born in province p (other mainland provinces).
- Geography-aware check: 1,000 permutations of cities' roots rows **within province** (22 provinces, 112 cities).

## 2. Main nights, observed vs fair null (mean ± sd)
Source: `r3_null_comparison.csv`. Gravity arm: 100 NB draws. BiCM arm: 100 Bernoulli draws.

| Statistic (gravity null) | Ordinary Feb 14 | Summer Aug 15 | All Saints Nov 1 | Christmas Eve | Good Friday |
|---|---|---|---|---|---|
| Q | 0.141 (0.106 ± 0.002) | **0.170** (0.115 ± 0.002) | 0.141 (0.103 ± 0.001) | 0.149 (0.115 ± 0.002) | 0.158 (0.109 ± 0.002) |
| restart NMI | 0.62 (0.43 ± 0.02) | 0.74 (0.47 ± 0.02) | 0.64 (0.45 ± 0.02) | 0.75 (0.48 ± 0.03) | 0.72 (0.45 ± 0.02) |
| territories | 18 (12.3 ± 0.9) | 17 (10.2 ± 1.1) | 19 (11.1 ± 1.1) | 18 (11.0 ± 1.1) | 14 (10.0 ± 0.8) |
| non-contiguous | 15 (12.2 ± 0.9) | 12 (10.1 ± 1.1) | 12 (11.1 ± 1.0) | 12 (10.8 ± 1.1) | 12 (9.9 ± 0.8) |
| share non-contiguous | 0.83 (0.99 ± 0.02) | 0.71 (0.99 ± 0.03) | 0.63 (0.99 ± 0.03) | 0.67 (0.99 ± 0.03) | 0.86 (0.99 ± 0.03) |
| diaspora flow share | 0.28 (0.52 ± 0.03) | 0.48 (0.70 ± 0.03) | 0.26 (0.60 ± 0.03) | 0.28 (0.58 ± 0.04) | 0.40 (0.72 ± 0.03) |
| roots enrichment | 1.76 (1.07 ± 0.04) | 1.41 (1.09 ± 0.04) | 1.41 (1.09 ± 0.04) | 1.60 (1.12 ± 0.04) | 1.48 (1.09 ± 0.04) |
| within-province perm z of enrichment | 3.7 | 1.6 | 4.5 | 5.6 | 3.0 |
| NMI(destination communities, province) | 0.48 (0.11 ± 0.01) | 0.37 (0.12 ± 0.01) | 0.47 (0.13 ± 0.02) | 0.52 (0.15 ± 0.02) | 0.39 (0.11 ± 0.01) |

| Statistic (BiCM null) | Feb 14 | Aug 15 | Nov 1 |
|---|---|---|---|
| Q | 0.466 (0.194 ± 0.002) | 0.266 (0.085 ± 0.001) | 0.425 (0.131 ± 0.001) |
| territories / non-contiguous | 7 / 4 (15.4 / 12.7) | 5 / 4 (10.6 / 10.2) | 6 / 5 (13.1 / 12.0) |
| diaspora share (links) | 0.035 (0.933 ± 0.011) | 0.078 (0.926 ± 0.010) | 0.028 (0.933 ± 0.011) |
| roots enrichment | 3.21 (1.01 ± 0.01) | 2.75 (1.00 ± 0.01) | 3.37 (1.00 ± 0.01) |
| within-province perm z | 0.8 | 0.7 | −0.5 |

The raw count of non-contiguous territories is slightly above the null (Aug 15: 12 vs 10.1 ± 1.1). That is only because the observed partitions have more territories (17 vs 10.2). The share is far below the null.

## 3. Holiday vs ordinary, each relative to its own null (sister nights)
Source: `r3_sister_summary.csv`. Excess = observed − null mean. Pairs = every (holiday sister, ordinary sister) combination.

| Gravity arm | Ordinary (13) median excess | Summer (9) | All Saints (3) | Summer > ordinary pairs | All Saints > ordinary pairs |
|---|---|---|---|---|---|
| Q | +0.030 | **+0.055** | +0.038 | 117 / 117 | 39 / 39 |
| restart NMI | +0.16 | +0.28 | +0.22 | 117 / 117 | 39 / 39 |
| roots enrichment | +0.40 | +0.39 | +0.35 | **47 / 117** | 11 / 39 |
| within-province perm z (observed) | 2.98 | **2.17** | 4.05 | 25 / 117 | 36 / 39 |
| diaspora share | −0.16 | −0.25 | −0.36 | below null on all nights | below null on all nights |
| share non-contiguous | −0.16 | −0.21 | −0.32 | below null on all nights | below null on all nights |

- The main ordinary night (Feb 14) is the unique maximum of its 13 sisters on roots enrichment (1.76; range 1.34–1.76), so these comparisons use sister medians.
- **Partition stability across nights** (`r3_partition_stability.csv`, gravity arm, NMI of city partitions with the main night):
  - summer sisters 0.81 (destination partitions 0.56);
  - ordinary sisters 0.60 (0.39);
  - summer vs ordinary 0.47 (0.27).
- The BiCM arm is far more stable: 0.94–0.96 for cities.
- The August gravity partitions are stable and specific to August. What makes them stand out is beach demand from specific city groups, not roots.
  - Example: the largest Aug 15 community is 51 cities, mostly Alicante/Valencia, with 930 zones.
  - By contrast, ordinary nights hold the clearest "diaspora" pairings, for example L'Hospitalet, Cornellà and El Prat grouped with Badajoz, Mérida and Seville province on Feb 14.

## 4. Alternatives (step 5)

### 4a. NMI between holiday city communities and roots-network communities (weak)
- Roots communities come from the city × birth-province network under the same gravity null and Leiden routine: 6 city communities, Q = 0.130, restart NMI 0.83 (`r3_alt_roots_partition*.csv`).
- Gravity-arm city partitions vs the roots partition (`r3_alt_roots_partition_nmi.csv`, sister medians):

| Night group | Ordinary | Summer | All Saints | Christmas Eve |
|---|---|---|---|---|
| NMI observed | 0.21 | 0.24 | 0.26 | 0.28 |
| NB null | 0.14 ± 0.02 | 0.12 ± 0.02 | 0.13 ± 0.03 | 0.14 ± 0.02 |
| within-province permutation z | 0.8 | 0.8 | 1.2 | 3.6 |

- The city partitions match province far better than they match roots (NMI 0.61–0.69), and the roots partition itself has NMI 0.40 with province.
- Only Christmas Eve passes the neighbour-city test. Not recommended.

### 4b. Two-axis decomposition of each city's holiday ego-network (recommended)
Each city gets two scores for each night (`_axes`):
- **Roots axis:** Pearson correlation, over the 46 other provinces, between log((O_cp + 1)/(E_cp + 1)) and log(1 + born_cp). E is the zone-level gravity fit summed to province.
- **Beach axis:** log(Σ_coastal O / Σ_coastal E), using `coastal_true`. Gravity keeps destination totals, so this axis measures a city's coastal preference *relative to other cities* at the same distance.

Each score is tested against the same NB gravity null refitted per draw (100 or 20 draws), and the roots axis also against province-mates' roots profiles. Sources: `r3_alt_two_axis_summary.csv`, `_pairs.csv`, `_thinned.csv`, and `_cities.csv` (per city).

| Sister medians (range) | Ordinary | Summer | All Saints | Christmas Eve | Good Friday |
|---|---|---|---|---|---|
| median roots axis | 0.36 (0.31–0.39) | **0.45** (0.43–0.47) | 0.42 (0.42–0.43) | 0.48 | 0.43 |
| null median | 0.067 | 0.044 | 0.047 | 0.049 | 0.044 |
| cities with roots z > 2 | 53% | 73% | 75% | 81% | 73% |
| own − province-mates' roots (median; share > 0) | 0.029; 75% | 0.045; 84% | 0.036; 73% | 0.046; 83% | 0.044; 85% |
| cities with beach z > 2 | 15% | 23% | 18% | 8% | 20% |
| Spearman(roots axis, beach axis) | −0.09 (p ≥ 0.06) | 0.03 (p ≥ 0.37) | 0.14 (p ≥ 0.057) | 0.06 | 0.11 |

- **Paired holiday − ordinary, per city** (`r3_alt_two_axis_pairs.csv`). Summer: median difference +0.083. All 117 pairs are positive (range 0.038–0.128), 59–81% of cities rise, and the largest Wilcoxon p is 1.4e-3. All Saints is positive in 39 of 39 pairs (+0.070), Christmas Eve in 13 of 13 (+0.116), Good Friday in 13 of 13 (+0.076).
- **Volume check** (`r3_alt_two_axis_thinned.csv`). Each holiday main night was binomially thinned to Feb 14 volume (keep probability 0.17–0.40, 20 replicates). The roots axis barely moves: summer 0.446 ± 0.010, All Saints 0.408, Christmas Eve 0.450, Good Friday 0.418, against ordinary 0.363. So noise attenuation does not explain the gap.
- **Reading.** Once each destination's total attraction (beaches included) and distance are absorbed by the gravity null, every holiday shifts each city's excess towards its own birth provinces. This includes August.
  - Christmas Eve shows the largest shift, then August, then Good Friday and All Saints.
  - This matches pattern-mining #7 only for Christmas. Pattern mining found August no higher than ordinary nights because its index does not remove destination attraction.
  - Roots and beach are separate axes across cities, which confirms pattern-mining #8 with a null.
- **Caveats.**
  - The ordinary-night baseline is already high (0.36).
  - The own-vs-neighbour gain is small (+0.03 to +0.05).
  - This is a city-level ego-network statistic, close to R2's city × province alignment. R2 and R3 should agree on which of them owns it.
  - The beach axis is relative by construction. The overall beach pull is in the destination margins (R1 / pattern-mining #5), not in this axis.

## 5. Notes for the orchestrator
- Legacy p93 normalised Q differently (spatial Q 0.0089, on window means). The new Q values are not comparable with it.
- The legacy observed "12 non-contiguous" count is close to ours (12 on Aug 15), but the fair null makes the *share* the meaningful statistic, and the share is below the null.
- The cache `results/r3_diag/cache/` (170 MB, `.npz`) is needed by `--summarise`. I added it to `.gitignore`, because `*.npz` was not ignored.
- I changed no files of other agents. My decisions are appended to `decisions.md`.
