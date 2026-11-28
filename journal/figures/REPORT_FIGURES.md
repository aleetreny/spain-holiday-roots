# EPJ Data Science figure production report

All nine figures are 170 mm wide. Each has a vector PDF, 300-dpi PNG and a plotted-data table in `out/`. The latter identifies panels and retains the input estimates at full precision. File dates are the mornings at the end of the nights.

To rebuild, run this command from this worktree. The two source checkouts are opened read-only by absolute path.

```sh
.venv/bin/python -I -B journal/figures/make_all.py
```

`make_all.py` also accepts a list of figure numbers such as `1 3 S2`. Figure 3 uses the two cached urban-born model fits in `out/fig3_urban_contrast.csv`. Run `urban_contrast.py` with the same interpreter and flags to refresh them. The fit uses the existing stacked PPML solver with the same night sets and model controls as the pre-registered village-born test (S4), with 95% cluster-robust intervals by origin province.

## Main figures

| Figure | What is visible | Trace |
|---|---|---|
| **1 · Design and data** | Net balance on 14–15 August ranges from −24.2% to +71.4% of resident population across the 47 provinces. The 12 largest excess corridors carry approximately 38,500 to 85,700 additional people per night relative to the 92 ordinary nights in 2024. The largest is Madrid → Alicante. | [fig1_data.csv](out/fig1_data.csv) province and corridor rows, computed from A's zone and 137-city panels. Coverage rows come from the A and A2 night manifests and C's monthly inventory. |
| **2 · Calendar** | P1 is +0.100, +0.106 and +0.113 in 2022, 2023 and 2024. The equal-year estimate is +0.106 [0.010, 0.197]. The 2025 hold-out gives +0.134 [0.037, 0.233]. Every yearly point exceeds its highest placebo. The pre-1981 term rises in August and the recent-arrival term rises around Christmas and New Year. | [fig2_data.csv](out/fig2_data.csv) uses A/A2 nightly, pre-registered and placebo tables. |
| **3 · Who pulls** | Pre-1981 and village-born premiums are +0.071 and +0.073 in 2022–2024 and +0.057 and +0.056 in 2025. The early-minus-recent gap is +0.115 and +0.089. The two urban-born contrasts are −0.028 [−0.059, +0.004] and −0.009 [−0.039, +0.021]. In B, the pre-1981 contrast is concentrated inland, with a common-support interaction of +0.114 vs February and +0.128 vs October. The within-65+ early-minus-later gap is +0.095 and +0.057 against the two baselines. | [fig3_data.csv](out/fig3_data.csv) uses A/A2 pre-registered tests, the supplementary urban-born fit, B's tidy table and its age table. |
| **4 · Receiving end** | Net receivers constitute 24.0%, 98.9% and 99.6% of the most emptied fifth on the three mapped nights. On 14–15 August, mean arrivals per resident are 0.710 and mean departures are 0.082. Relative to 13–14 February, this is +0.644 arrivals and −0.013 departures. The adjusted inland August slopes are +0.144 arrivals and −0.076 departures, with corresponding Christmas Eve slopes of +0.120 and −0.080. | [fig4_data.csv](out/fig4_data.csv) uses the D zone, bootstrap and controlled-decomposition outputs and A's daily Q5 series. The single-municipality series is calculated across all valid nights from A's zone panel, using D's 84-zone subset. Both Q5 series agree exactly with D over all 302 overlapping valid nights. |

## Supplementary figures

| Figure | What is visible | Trace |
|---|---|---|
| **S1 · Network** | The in/out-degree Spearman correlation falls from 0.69 to 0.09. Rewiring reproduces 93% of the assortativity shift. August second-order mixing is −0.289 observed and −0.252 under rich-club-preserving rewiring. | [figS1_data.csv](out/figS1_data.csv) uses the coursework node metrics and structural constraints. |
| **S2 · Parents' villages** | August publication premiums recur across years. The outer envelopes for the August count contrast are approximately −65 to +93 trips against February and −65 to +88 against October. | [figS2_data.csv](out/figS2_data.csv) uses C's monthly and August estimates. |
| **S3 · Robustness** | P1 is +0.07 to +0.15 across general specifications, with a coastal-only estimate of +0.01. P2 varies across zero rules (+0.014 to +0.122). The leave-one-destination strips show all 47 fits per hypothesis. | [figS3_data.csv](out/figS3_data.csv) uses A's robustness and destination-omission outputs. |
| **S4 · Bridge** | Mean absolute paired differences are 0.014 for pooled roots and 0.005 for the pre-1981 term. | [figS4_data.csv](out/figS4_data.csv) uses A2's date-level method bridge. |
| **S5 · Zero cells** | With the main +0.5 rule and the February baseline, the predicted-flow contrasts are +7.9% for pre-1981 arrivals and +13.6% for village-born residents. Both contrasts are positive across the displayed zero-cell and sample rules. | [figS5_data.csv](out/figS5_data.csv) uses B's flow-weighted one-raw-SD prediction contrasts. |

## Design choices and limits

- **Fig. 1** arrows show aggregated excess travel from origin cities into whole destination provinces. The endpoints are province representative points. Arrow widths encode magnitude and adjacent corridors can overlap near a shared city or province. The 2025 hatching marks a method era, while the ticks show actual file availability.
- **Fig. 2** annual 2022–2024 rows in panel b use cluster-robust intervals supplied by A. The pooled and 2025 tests use wild-cluster bootstrap intervals. Placebo windows overlap, so their distributions are descriptive. Nightly dots and ribbons have substantial overplotting at a 170 mm width. The 2025 series is sparse after July.
- **Fig. 3** the urban-born row is an additional calculation with 95% CRV1 intervals. The other rows use 95% wild-cluster bootstrap intervals. The two published arrival and birth-size groups use distinct colour roles from the destination panel. The direct B interactions have a different support and estimand from the pooled P3 test. The pooled interaction is inconclusive. The B coastal village-born interval excludes zero against October.
- **Fig. 4** 24–56% describes D's four selected ordinary baselines. A's full ordinary-night class ranges from about 17% to 88%, as the strip shows. The strip shows the August test window and the rest of August separately. Across all August nights the range is about 96–99% in A and about 97–99% in the 93-night D cache. The full-calendar Christmas Eve values cover three years. Inland coefficients describe a controlled gradient across zones. Housing attenuation remains in the source report's supplementary results.
- **Figs. S2 and S3** the count envelope in S2 is based on the suppression bounds and their endpoint uncertainty. The strips in S3 show the spread of leave-one-destination estimates.

Birth shares and overnight flows are measured between places. The figures describe associations at that level.

## Production checks

I inspected all nine PNGs at their intended width and corrected collisions in panel titles, legends and labels. The final rebuild completed in figure order. PyMuPDF confirmed that all nine PDFs are single-page vector files with embedded Arial, no raster layers and all extracted text inside the page. Pillow confirmed 300 dpi and the expected image dimensions.
