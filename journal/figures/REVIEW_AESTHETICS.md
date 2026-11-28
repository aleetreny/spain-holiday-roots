# Figure aesthetics review (EPJ Data Science)

27 November 2026. Scope: the nine figures in `out/` (4 main, 5 SI), all 170 mm wide. No plotted number or claim changed: every `out/fig*_data.csv` column that existed before is identical row for row (checked with pandas). Fig. S1 adds three columns for the new elements (zone type and the structural-constraint range).

Rebuild: `.venv/bin/python -I -B journal/figures/make_all.py`

## One theme for all nine figures (`style.py`)

- **Type scale (pt):** ticks, annotations and legends 7; axis labels 7.5; panel titles 8.5 (regular); panel letters 10 bold. Arial throughout, including maths text (Arial italic replaces the DejaVu fallback in Fig. 1b).
- **Layout in millimetres:** `place(fig, left, top, width, height)` places every axes from the figure's top-left corner, so panels in a row share their top and bottom edges and panels in a column share their left and right edges.
- **Panel letters and titles:** `panel(ax, letter, title, at=mm)` puts the letter at a fixed x position for each column and the title 4.2 mm to its right. Both sit on one baseline 3.2 mm above the axes. Every figure now uses the same letter offset (2 mm from the left edge for the first column).
- **Marks:** estimates use 4 pt markers and 0.85 pt interval lines. Background observations use 3.5 pt markers. Reference lines are 0.55 pt dashed `#444444`. The grid is `#ececec` and sits behind the data.
- **Numbers:** axes use a fixed number of decimals per axis and a true minus sign (−), via `decimals(n)`.
- **Legends:** legends are frameless and use the same handle sizes as the data. A panel's key sits on its title baseline, right-aligned, or directly under its axis label, never far from the data.
- **Colour roles:** unchanged validated hex values. The destination pair (terracotta/blue) and the cohort pair (violet/teal) never share a panel. Text is always ink or grey, never a series colour.

## Per figure: critique, changes, reasons

### Fig. 1 — Design and data
*Before:* the "Madrid" label covered Ávila, the dark-green neighbouring province, so it read as if Ávila were Madrid. Arrowheads were tiny and the arrows lost contrast where they crossed the green provinces. The width of an arrow had no key. Panel b's two-line key crowded its title, and the right column started lower than the map, leaving a block of empty space. In panel c the source labels ("Daily ticks MITMA / Monthly bars INE") floated under the rows and did not line up with the month ticks. The legend sat far below, and "New method" was loose text.
*After:*
- **a:** Madrid's label now sits to the right of the city dot, with a white halo, away from Ávila. The arrows are ink with a white casing, and their heads scale with their width. A key for arrow width (40,000 and 80,000 excess people per night) sits beside the colour bar, both on one row under the map. The colour bar has five ticks (−80 % to 80 %). Provinces have white borders.
- **b:** the panel's top now lines up with panel a. The census and MITMA key moved under the diagram, and *p* and *q* are italic.
- **c:** each month has a one-letter label centred on its span, with tick marks at month boundaries. INE months appear as bars under each year. The legend sits directly under the axis, in two columns of four, and covers every mark: ordinary, August window, other MITMA nights, the three holidays, INE month and the hatched new MITMA method.
*Considered and rejected:* swapping the shading to the sequential August-minus-ordinary change. Rejected because the diverging net balance is validated and already in the caption. Once the arrows were cased and darker, the shading and the arrows read as two separate layers.

### Fig. 2 — The roots calendar
*Before:* in panel a the year labels were on the left and the y-axis was drawn by hand on the right, with only 0.4 and 0.7 as ticks. No key explained the grey dots, the dark points, the line or the ribbon. Panel b had a 7.2 pt two-line y-label and a "Grey dots: …" note floating inside the data area. In panel c the key crowded the title, the shaded windows were almost invisible (alpha on `#e6e6e6`), and the method-change line had no label.
*After:*
- **a:** each row has a normal left y-axis (0.4/0.6/0.8) and one shared label, "Roots elasticity". Year labels are bold, in the right margin, so the left edge is shared by all three panels. A key on the title line covers ordinary night, August window, other nights and the 95 % CI. Holiday labels are centred above their first marker, with a halo. The 2025 row reads "1 July 2025: new MITMA method".
- **b:** a single "Change in roots elasticity" label at the same x position as panel c's label. Two-decimal ticks with true minus signs. A faint divider separates the yearly tests from the pooled and hold-out tests. The key, including the placebo dots, sits on the title line.
- **c:** the panel is taller (24 mm). The cohort lines are 1.0 pt. The windows are solid `#e6e6e6` and named in the key. The method-change line is labelled, and the axis gains a "Cohort term" label.

### Fig. 3 — Who pulls
*Before:* the letter "b" sat higher than "a". The x labels were partial phrases, and a bold "Change in roots elasticity" floated at the bottom. The legends sat at the very bottom, far from panels a and c.
*After:* the a/b tops and the c/b bottoms now line up on the mm grid. Each axis has a full two-line label. The confirmation/hold-out key sits under panel a, and the February/October key under the figure, shared by b and c. The group headers in b are placed relative to the axes, so they stay aligned with their rows.

### Fig. 4 — The receiving end
*Before:* panel b's seven vertical category labels collided ("Christmas Eve", "Other holidays"). The 24–56 % bracket was an unlabelled vertical line at the left edge, with a separate "Bracket: …" footnote. Panel c mixed a dot plot with a three-line prose block, and its legend was a text glyph inside the axes. The map titles crowded the shares.
*After:*
- **a:** each map has its night as the title and "xx.x % net receivers" in grey on a second line. The legend is centred under the three maps.
- **b:** the strip is now horizontal, with categories as single-line row labels, so nothing collides. All zones are filled dots and the 84 single-municipality zones are open dots, offset within each row. The bracket is drawn horizontally above the Ordinary row and labelled "Ordinary baselines 24–56%".
- **c:** "Arrivals" and "Departures" are labelled directly beside the August points, with no legend.
- **c (table):** the inland adjusted changes now form a small aligned table under the panel (rows 14–15 Aug and Christmas Eve, columns Arrivals and Departures), with the same numbers. It is clearer than prose and adds no new claim; I chose it over plotting the coefficients with intervals for that reason.

### Fig. S1 — Network structure
*Before:* the a/b scatters were monochrome, so the reader could not see *which* zones break the in/out symmetry. Panels c and d joined three categorical nights with lines, which implies a trend between nights. The 93 % sentence sat as a loose footer.
*After:*
- **a/b:** coloured with the project's zone-type roles: terracotta inland, blue coastal, ink cities over 100,000. This is taken from the coursework Fig. 2 (see below). The legend sits inside a, and ρ is placed bottom-right in both panels.
- **c/d:** now dot-and-tick rows, one per night. Dot = observed; tick = mean of 50 rewirings. In c, grey bars show the range reachable by degree-preserving rewiring (Zhou & Mondragón structural constraints; values from the same coursework CSV, `r_min_maxdisassortative`/`r_max_maxassortative`).
- **Shared:** one legend for c and d, with the 93 % note under it. The caption is updated.

### Fig. S2 — Parents' villages
*Before:* panel b's title nearly touched a's axis. Panel c's "Outer 95% envelope" label sat inside the data area and the estimated bounds had no key. The y ticks used hyphens.
*After:* equal 11 mm gaps between panels. A key on the title line of a (monthly change, August, 95 % CI) and of c (estimated bounds, outer envelope). The y ticks run −4/0/4/8 with a true minus. Row labels read "vs February / vs October".

### Fig. S3 — Robustness
*Before:* the row spacing differed between a (13 rows) and b (15 rows). The "b" letter sat mid-figure, away from its labels. A grey bold footer stood in for the axis labels.
*After:* both panels use an identical 4.6 mm row pitch and are top-aligned. The leave-one-out strips line up on a shared baseline, each carrying a full two-line x label. Main-specification diamonds and all estimates are in ink.

### Fig. S4 — Method bridge
*Before:* the letters sat at different offsets and the decimals were uneven (0.095 … 0.115 spaced oddly). The legend was text glyphs.
*After:* both panels have the same width and the same letter offset. Each panel has five ticks (2 or 3 decimals, as its scale needs). The legend uses real markers, and the night labels repeat in panel b so each panel reads on its own.

### Fig. S5 — Zero-cell rules
*Before:* the letters for b and d sat close to their axes while a and c sat far left, and a long single-line x label spanned both columns.
*After:* the four panels sit on one mm grid with identical letter offsets. Each bottom panel carries the full two-line x label. The legend sits just under the figure.

## Coursework figures: reuse decision
- **Reused (an element):** the zone-type colouring of the degree scatter and the structural-constraint range bar from `reports/figures/fig2_senders_receivers.png`, now in Fig. S1. They make the SI figure answer "which zones change role" and "how far could rewiring push r" at a glance. They use only the validated destination roles plus ink, and they add no claim beyond the coursework results already in the paper (ρ 0.69 → 0.09; range [−0.22, −0.09] → [−0.42, −0.11]).
- **Not reused:**
  - The **two-Spains flow map** (`figures_drafts/F1_A_two_spains.png`) is striking. But it shows total ≥100 km flows by destination type on two nights, a descriptive claim ("4.5 million people…") that the journal paper does not make. Its serif styling and thousands of semi-transparent arcs would also break the uniform theme and cannot be read at 85 mm. Fig. 1a's excess-corridor arrows tell the paper's own story (excess August flows) with far less ink.
  - The **inland/coast calendar** (`fig1_option2_calendar.png`) counts people sleeping ≥100 km away. That is a different quantity from the nightly roots elasticity in Fig. 2a, and putting both in the main text would duplicate the calendar idea with a second, unvalidated measure. It would fit an SI "descriptive context" figure if the authors want one later.
  - `fig3_roots_generation.png` and `fig4_villages.png` are superseded by the journal's Figs 3 and 4, which use the newer pre-registered estimates.

## Production checks (final build)
- All nine PDFs: single page, vector only (0 raster images), fonts embedded (Arial, Arial Bold, Arial Italic), no text outside the page.
- All nine PNGs: 300 dpi, 2007 px wide (170 mm).
- Heights (mm): Fig. 1 106, Fig. 2 182, Fig. 3 119, Fig. 4 134, S1 124, S2 113, S3 117, S4 60, S5 116. All are ≤ 200 mm.

## Remaining issues
- **Fig. 2a:** ordinary-night dots still overlap into grey bands, because consecutive weekdays sit 1 day apart at about 0.37 mm per day. Smaller markers would break the 3.5 pt minimum. The caption carries the meaning.
- **Fig. 4b:** the dense dark clusters near 100 % (August, other August) overplot. The two series stay separated by their vertical offset.
- **Fig. 1a:** three Madrid arrows overlap near the origin. With 12 corridors this is inherent to a city-centred star.
- **Fig. S3:** panel a has 2 fewer rows than b, leaving white space above its strip. I kept it on purpose, so both strips share a baseline and the row pitch is equal.
- `captions.md` was edited only for Fig. 4c ("table below the panel") and Fig. S1 (zone colours, ticks and range bars replace the connecting curves). `REPORT_FIGURES.md` still describes the earlier layout under "Design choices" for Fig. 1 and Fig. 4. The numbers in it are unchanged.
