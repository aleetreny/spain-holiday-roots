# Figure specification (EPJ Data Science version)

27 November 2026. Numbers and claims: [`../analysis/RESULTS_SUMMARY.md`](../analysis/RESULTS_SUMMARY.md). Journal rules: [`../research/epjds_style.md`](../research/epjds_style.md) §3 and §5. Four main figures, each answering one question; supplementary figures listed at the end.

## Production rules (all figures)
- Output: vector PDF (fonts embedded, `pdf.fonttype = 42`) + PNG preview at 300 dpi, in `journal/figures/out/`. One script per figure, `journal/figures/figN_*.py`, plus a shared `journal/figures/style.py`. Every plotted value is also written to `journal/figures/out/figN_data.csv`.
- Width 170 mm (double column) unless stated; height ≤ 200 mm so figure + legend fit 225 mm.
- Font: Helvetica/Arial-like sans (TeX Gyre Heros or Arial, whichever is installed; same family in all figures). Tick and label text 7–8 pt, panel titles 8.5 pt, panel letters **a b c** bold 10 pt, top-left outside the axes.
- Lines 0.75–1.0 pt for data, 0.5 pt for axes; no line thinner than 0.25 pt. Markers ≥ 3.5 pt. Grid: very light (`#ececec`), y only, behind data. No top/right spines. No chart titles inside the figure beyond short panel titles; the caption carries the explanation.
- Uncertainty: 95 % intervals as thin lines (points) or light ribbons (series), always defined in the caption.
- Text in ink (`#222222`) or secondary (`#666666`), never in a series colour. Direct labels where there are ≤ 4 series; a legend otherwise.
- Never a dual y-axis.

## Colour roles (validated with the dataviz validator, light surface `#ffffff`)
| Role | Hex | Use |
|---|---|---|
| Inland destinations | `#c0623b` (terracotta) | destination type |
| Coastal destinations | `#2b6f9e` (blue) | destination type |
| All destinations / reference | `#444444` | pooled estimate, reference line (neutral, not a series hue) |
| Exodus generation (pre-1981; village-born) | `#7b3fa0` (violet) | cohort |
| Recent arrivals (2016–2023; urban-born) | `#0f8f7a` (teal) | cohort |
| Net senders / net receivers | `#7b3fa0` → `#efefef` → `#1b7837` | diverging maps only |
| Ordinary nights | `#bdbdbd` points / `#e6e6e6` band | background |

Validator results: terracotta/blue adjacent CVD ΔE 15.8, normal 24.2 (pass); violet/teal CVD 15.3, normal 25.9 (pass); violet/green diverging poles CVD 17.7, normal 28.4 (pass). **Never put the destination pair and the cohort pair in the same panel** (all-pairs fails: violet↔blue CVD ΔE 3.9). Period (2022–24 vs 2025 hold-out) is encoded by **marker fill** (filled = 2022–24, open = 2025), never by a new colour.

## Figure 1 — Design and data (170 × ~120 mm)
Question: what is measured and how the two networks are aligned.
- **a** (left, ~60 % width): map of mainland Spain, provinces outlined (thin `#bdbdbd`), zones not drawn. Shading = net overnight balance per resident on 14–15 Aug 2024 by province (diverging map roles) or, if cleaner, the August-minus-ordinary change in overnight visitors per resident by province (sequential, one hue). Overlay the 10–15 largest **excess** city → province corridors (August mean of the 10–18 Aug files minus the mean of same-year ordinary nights), as curved arrows with width ∝ excess persons, neutral dark grey. Label 3–4 big origin cities only. Scale bar, no north arrow needed.
- **b** (right, top): schematic of the dyad: a city node, two destination provinces, the census tie (residents of the city born in province *p*, dashed) and the overnight flow (persons from the city sleeping in *p*, solid), with the model in words: "flow ~ birth share, comparing cities of the same province against the same destination, with distance and income". Monochrome boxes and arrows; no colour except one accent for the tie.
- **c** (right, bottom): coverage strip 2022–2025: one row per year, one tick per available MITMA night; ordinary reference nights grey, August window terracotta-free neutral dark, Christmas Eve / All Saints / Good Friday marked with small symbols; 2025 after 1 July hatched or labelled "new method"; a second thin row per year for INE monthly data (months available). Legend below.
- Data: zone geometry `clean_data/v2/zones/municipios_mitma.gpkg`; night dataset and panel `_local/journal_data/A_calendar/` (+ `A2_confirmation2025/`); manifests `journal/analysis/A_calendar/out/night_manifest.csv`, `journal/analysis/A2_confirmation2025/out/night_manifest.csv`; INE inventory `journal/analysis/C_municipal/out/od_inventory.csv`.

## Figure 2 — The roots calendar (170 × ~150 mm). Central figure.
Question: when do flows follow birthplaces most?
- **a**: nightly pooled roots elasticity β_t for every valid night 2022–2025 (x = date, one long panel or four stacked year rows sharing the y-axis). Ordinary reference nights as small grey points; all other nights as a thin dark line; August window nights as filled dark points; Christmas Eve, All Saints, Good Friday as distinct small symbols with direct labels the first time they appear. Per-night 95 % CI as a very light ribbon. Mark the 1 July 2025 method change with a thin vertical dotted line and label. Horizontal reference = mean of ordinary nights per year (thin dashed).
- **b**: the pre-registered test as a dot-and-interval plot: δ_AUG (P1) for 2022, 2023, 2024 (confirmation against new ordinary nights) and 2025 (hold-out), plus the pooled 2022–24 estimate; filled vs open markers by period. Behind each year, the distribution of placebo-window contrasts as a light strip or small violin, so the reader sees August above all placebos. Y = "Change in roots elasticity, August vs ordinary nights".
- **c**: cohort terms over the calendar: pre-1981 term (violet) and 2016–2023 arrivals term (teal) by night, smoothed only by plotting points + a 7-night rolling median line; August windows and Christmas/New Year/Epiphany shaded very lightly in grey. Shows the exodus generation peaking in August and recent arrivals at Christmas.
- Data: `journal/analysis/A_calendar/out/{nightly_pooled,nightly_cohorts,placebo_windows,preregistered_tests,night_manifest}.csv` and the same files in `journal/analysis/A2_confirmation2025/out/`.

## Figure 3 — Who pulls (170 × ~110 mm)
Question: which residents' birthplaces gain weight in August, and where.
- **a**: forest plot of August contrasts by group: pre-1981 arrivals, arrived 2016–2023, born in a municipality ≤ 2,000, born in a municipality > 10,000, and the pre-1981 minus recent gap. For each group: 2022–24 confirmation (filled) and 2025 hold-out (open), both against ordinary nights, 95 % intervals; a zero line. Cohort colours.
- **b**: the same contrasts for pre-1981 and village-born by destination type on the coursework nights (B): all / inland / coastal, against 13–14 Feb and against October weekdays (two marker shapes), destination colours. Include the direct inland × coast interaction as a separate row.
- **c** (small): age-matched check: pre-1981 minus later arrivals within 43–64 and 65+ (from B), neutral colour.
- Data: `journal/analysis/A_calendar/out/preregistered_tests.csv`, `journal/analysis/A2_confirmation2025/out/preregistered_tests.csv`, `journal/analysis/B_cohorts/out/fig_who_pulls_tidy.csv` and the age table in `journal/analysis/B_cohorts/out/`.

## Figure 4 — The receiving end (170 × ~120 mm)
Question: what happens to the most emptied places.
- **a**: three small maps (ordinary 13–14 Feb 2024, 14–15 Aug 2024, Christmas Eve 2024): zones in the most emptied fifth coloured by net receiver (green) / net sender (violet), other zones very light grey; province outlines.
- **b**: share of net receivers in the most emptied fifth on every valid night 2022–2024 (strip or dot plot by night class: ordinary, weekends, August window, Christmas Eve, other holidays), with n = 475 and the 24–56 % ordinary range visible; second series for single-municipality zones (open markers).
- **c**: arrivals per resident and departures per resident (two bars or dot pairs) for the most emptied fifth: ordinary vs August vs Christmas Eve 2024, with bootstrap intervals; plus the inland controlled decomposition (arrivals vs fewer departures) as text annotation or a small inset.
- Data: `journal/analysis/D_villages/out/` (`fig_receiving_end_tidy.csv`, `d_zones_main_nights.csv`, `d_quintiles_all_nights.csv`, `d_quintiles_main_boot.csv`), `journal/analysis/A_calendar/out/descriptive_q5_villages_daily.csv`, zone geometry gpkg.

## Supplementary figures (85 or 170 mm, same style)
- **S1** Network diagnostic: in/out degree split (ρ 0.69 → 0.09), assortativity observed vs degree-preserving rewiring (null reproduces ~93 % of the shift), second-order mixing. Restyle from coursework `pipeline/f_figures.py` (`fig2`, `figbasic`) and `results/b02_*.csv`.
- **S2** Parents' villages: monthly change in publication probability with a parental tie, 2019–2026, against February and October, 95 % intervals; and the count-margin envelope. `journal/analysis/C_municipal/out/calendar.csv`, `estimates.csv`.
- **S3** Robustness of P1 and P2 (forest): `journal/analysis/A_calendar/out/robustness.csv`, `robustness_lodo_dest.csv`.
- **S4** Method bridge 2025: β old vs new method on the 6 recalculated dates. `journal/analysis/A2_confirmation2025/out/method_bridge_dates.csv`.
- **S5** Zero-cell rules as predicted-flow contrasts. `journal/analysis/B_cohorts/out/`.

## Captions
Draft one caption per figure (≤ 15-word title, ≤ 300-word legend, panel-by-panel, units, sample, interval definition) in `journal/figures/captions.md`.
