# Plan for the EPJ Data Science paper

27 November 2026. Synthesis of four research reports in this folder:
[`literature_epjds.md`](literature_epjds.md) (novelty map, 81 must-cite works),
[`epjds_style.md`](epjds_style.md) (house style, 10 papers read in full, submission rules),
[`audit_figures_pipeline.md`](audit_figures_pipeline.md) (inventory, referee report, headline checks),
[`design_first_principles.md`](design_first_principles.md) (designs, data, pilots).
Status: **proposal, awaiting author approval.** Nothing below is a final result.

## 1. The verdict in one paragraph

The paper is not a complex-networks paper. It is a data-science paper about **whether the geography of past internal migration shows up in where people sleep on holiday**, measured by linking mobile-phone overnight flows to census and survey birthplaces. EPJDS states that "purely descriptive statistics" or standard methods on mainstream data are unlikely to be considered ([scope](https://link.springer.com/journal/13688/aims-and-scope)), so the paper must sell a *new measurement* (birthplace ties as a layer of holiday mobility) and a *test that survives the obvious objections* (distance, beach, baseline choice, ecological inference). The coursework already has most of the ingredients. What it lacks, and what a referee will hit first, is (a) a holiday effect that does not depend on picking one February or one October baseline, and (b) evidence below the city-to-province level. The plan fixes (a) with a night-by-night series on ~1,000 nights that include never-analysed nights, and attempts (b) with parents' birth municipalities, under an explicit go/no-go gate.

## 2. Where the novelty sits (corrected statement)

**Closest works** ([`literature_epjds.md`](literature_epjds.md) §1): Li et al. 2016 (Spring Festival vs census, 31 provinces, raw QAP, no distance control, no ordinary baseline); Wen et al. 2026 *J. Transp. Geogr.* (holiday vs regular intercity networks in China, cultural proximity and housing, no birthplace matrix, no degree null); García-González et al. 2024 and Garzón-Paredes & Royo-Vela 2026 (Spanish phone data, return tourism and summer-receiving areas, node-level, no birthplaces). The style report adds He et al. 2026 in EPJDS itself (Nanjing Spring Festival), which shows that "holiday + mobility network" alone is not novel.

**Defensible novelty.** No published study links national mobile-phone overnight flows, on holiday and ordinary nights, to census birthplace stocks by arrival cohort and to survey-reported parents' birth municipalities, with distance, destination and corridor fixed effects. The contribution is the *dyadic* question (which cities send visitors to which places) and its *calendar* (when the alignment strengthens), not the observation that rural Spain fills up in August.

**Claims we must not make** (consolidated from all four reports):
- First to compare holiday travel with census migration (Li et al. 2016); first to show holiday network reconfiguration (Wesolowski 2017, Li 2020, Wen 2026); first to show rural Spain refills in summer (several Spanish phone studies); first to put migration or social ties in gravity models (Bailey et al., Santana-Gallego & Paniagua 2022, Prieto-Curiel et al. 2024).
- That individuals return to *their* family village; that second homes *mediate* anything (it is coefficient attenuation); that the coast effect is zero (it is "not detected"); that roots grow in August without naming the baseline.
- Note: the novelty paragraph in `literature_epjds.md` §1.2 over-reaches. "~93 %" applies to the assortativity shift only, not to rich club, k-core and communities, and "via inherited non-primary homes" is a causal reading the data do not support. Use the statement above instead.

## 3. Research design

Four analyses, in order of priority. Each says what it answers, what already exists, what is new, and how it can fail.

### A. The roots calendar (core, new)
- **Question.** On every available night, how strongly does a city's overnight flow to each other province follow the share of its residents born there? Is the elasticity higher on holiday nights than on the whole distribution of ordinary nights, and on which holidays?
- **Why it is the centre.** It turns the biggest weakness (the August increase is +0.096, p = 0.18 against 13–14 Feb and +0.137, p = 0.006 against October) into a result: instead of one baseline, every ordinary night is a baseline. It also answers "is 15 August special or is it all of August".
- **Data.** MITMA lists ~1,000 standard daily files for 2022–2024; the project holds 306 of them. The ~780 missing nights of 2022–2024 cost about 1.9 GB at ~2.4 MiB per night ([`design_first_principles.md`](design_first_principles.md), "MITMA nights available"). 2025 is excluded because of the 1 July 2025 method break.
- **Model.** The existing R2 PPML (city FE + origin × destination-province FE, distance, income × destination), estimated night by night on a fixed city set; holiday contrasts against calendar-matched ordinary nights; uncertainty by calendar-block resampling plus the wild cluster bootstrap. Same series for the cohort and village-born terms (675-city design).
- **Hold-out.** The ~780 nights never analysed are a genuine confirmation sample. Freeze the primary hypotheses and contrasts in `journal/analysis/design.yml` and commit it **before** running them. This is the answer to the multiple-testing and specification-search objection (audit risk 4).
- **Can fail if** the holiday elasticity sits inside the range of ordinary nights. Then that is the finding (roots shape travel every night; holidays add volume, not alignment), and the cohort split carries the holiday story.

### B. Who pulls: arrival cohort and birthplace size (core, exists, needs hardening)
- Existing: pre-1981 arrivals +0.084, village-born +0.086, recent −0.072; inland +0.154 / +0.123 (`results/v1_verify_n1.md`, re-fitted independently by the audit).
- To do: the age recode found by the audit (3,245 older people with a special arrival code; older post-1981 term moves from −0.0055 to −0.0085); direct inland-minus-coast interaction on common support; Table A2 zero-cell rules shown as predicted-flow contrasts; the nightly series from A. Recent-arrival inversion stays secondary.

### C. Parents' villages (novelty upgrade, gated)
- **Question.** Within one destination province, do a city's flows go more to the specific small municipalities where its city-born residents' parents were born (ECEPOV 2021), and does this strengthen in August?
- **Pilot** (`design_first_principles.md`, Pilot 1, 137 cities × 5,849 villages, INE municipal monthly trips 2022–2026): a parental tie raises the August-minus-February probability that INE publishes the flow by **+3.90 pp [3.04, 4.75]** (vs October +3.15 pp); all five yearly estimates lie between 3.44 and 4.26 pp. Among pairs published every month, the count growth is **+0.056 [−0.027, 0.140]** vs February and **+0.050 [0.015, 0.084]** vs October.
- **Main risk.** INE suppresses cells under 30 trips, and suppression is strongly linked to having a parental tie (24.6 % vs 1.8 % of dyads visible in both months). The earlier N1 results filled suppressed cells with zeros.
- **Gate** (audit §2): proceed to the main text only if the result is positive on a fixed eligible frame, under at least two treatments of suppression, with ECEPOV sampling uncertainty propagated and a held-out year. Otherwise it goes to the SI as a measurement limitation and the paper stays at city-to-province level. That is still publishable, but with less novelty.

### D. The receiving end: emptied villages (keep, decompose)
- Existing: in the most emptied fifth of zones, net receivers 114/475 (24 %) on the ordinary night, 470/475 (99 %) on 15 August, 473/475 on Christmas Eve (re-computed by the audit).
- To do: split into arrivals and departures per resident (August works mainly because fewer residents leave); housing as attenuation only, in the SI.

### What leaves the main text
- **To the SI, as one compact diagnostic:** the network-structure part (density, in/out split, assortativity vs degree-preserving rewiring, second-order mixing). It earns one sentence in the results and one SI figure, because Wen et al. 2026 and the Chinese rich-club papers read holiday structure as reorganisation, and our null shows the assortativity shift is mostly forced by the degrees.
- **Cut:** BOSAM, power laws, small world, k-core story, community "territories" (one line in the SI as a failed test), PageRank and betweenness.
- **SI only:** out-of-sample prediction (roots cut publication log loss by 2.1 % over geography alone but only 0.4–0.6 % once lagged flows are included); SCI comparison; ECEPOV second homes in the family's birth region (77.4 % of inland cross-region second homes).
- **Dropped for now:** local festivals (fiestas). No national machine-readable calendar for 2022–2025, and the Castilla y León pilot found one usable event and none in villages of up to 2,000 people.

## 4. Figures (5 main + SI)

EPJDS papers in the sample average 5–6 main figures with 35–75-word captions and a Figure 1 that explains the design ([`epjds_style.md`](epjds_style.md) §3).

| # | Figure | Panels | Status |
|---|---|---|---|
| 1 | **Design and data** | a, signed excess holiday overnight flows on a map of mainland Spain; b, schematic of how census birthplace ties and phone flows are aligned on city → province (and city → village) dyads; c, calendar strip of nights and months per source | New; reuses `f01_hero_drafts` geometry |
| 2 | **The roots calendar** | a, nightly roots elasticity 2022–2024 with holidays marked; b, holiday nights against the distribution of matched ordinary nights; c, cohort and village-born terms on the same calendar | New, central |
| 3 | **Who pulls** | pre-1981 / village-born / recent, inland vs coast, against several baselines, common scale with 95 % intervals | Restyle of Fig. 3 |
| 4 | **Parents' villages** | publication probability and count margins by month and year; robustness to thresholds | New; only if gate C passes, else SI |
| 5 | **The receiving end** | a, net-balance maps; b, receiver share by loss quintile; c, arrivals and departures per resident | Restyle of Fig. 4 |
| SI | Network diagnostic, Table 1, zero-cell rules, sister nights, province jackknife, SCI, prediction, housing, data coverage | — | Mostly existing |

## 5. Data to add

| Data | Size | Use | Where |
|---|---|---|---|
| ~780 MITMA nights 2022–2024 not yet local | ~1.9 GB | Analysis A (and hold-out) | `_local/external_journal/mitma_pernoct/` in the journal worktree; never in git |
| GISCO LAU 2021 municipal geometry | ~123 MB | Real municipal distances and coast for C | same |
| Everything else (INE OD 87 months, ECEPOV, census microdata, housing, SCI) | already local | B, C, D, SI | main checkout `_local/external/` (read-only) |

## 6. Work plan and effort

Effort from the audit and design estimates (hours of analysis; writing excluded).

| Step | Content | Effort |
|---|---|---|
| 0 | `design.yml`: freeze primary hypotheses, contrasts and hold-out nights; commit before touching new nights | 0.5–1 day |
| 1 | Download and validate MITMA nights; date and coverage manifest (DANA window, missing days) | 1 day |
| 2 | Analysis A: nightly roots and cohort series, matched baselines, block bootstrap | 2–3 days |
| 3 | Analysis B hardening (age recode, inland × coast interaction, zero-cell contrasts) | 1–2 days |
| 4 | Analysis C: municipal panel with suppression-aware models, ECEPOV uncertainty, held-out year; apply the gate | 3–5 days |
| 5 | Analysis D decomposition; SI diagnostics (network, SCI, prediction) | 1–2 days |
| 6 | Figures in house style; tidy result tables feeding text and figures | 2–3 days |
| 7 | Manuscript in the Springer Nature template (`sn-jnl`, pdfLaTeX, numbered citations), SI, declarations, Zenodo release | 4–6 days |

Roughly **3–4 weeks** of focused work. New code lives in `journal/analysis/` and new files only; every change against the coursework is logged in `journal/CHANGES_vs_coursework.md`.

## 7. Format decisions already fixed by the journal

Regular Article; abstract 150–250 words (coursework: 264), 3–10 keywords; numbered square-bracket citations; pdfLaTeX with the Springer Nature template, double spacing, line numbers; figures 85 or 170 mm wide, titles ≤ 15 words, legends ≤ 300; declarations: availability of data and materials, competing interests, funding, authors' contributions, acknowledgements, LLM-use statement in Methods; single-anonymous review; suggested reviewers optional; median 10 days to first decision ([`epjds_style.md`](epjds_style.md) §5). Target length ~6,000–7,000 words, 50–80 references (coursework: 20).

## 8. Decisions for the author

1. Approve the framing (birthplace ties as a layer of holiday mobility; network part to the SI) and the design A–D with gate C.
2. Approve downloading ~780 MITMA nights (~1.9 GB) and the municipal geometry (~123 MB).
3. Authorship: sole author ("I") or with co-authors ("we"), e.g. a supervisor.
4. Title. Options: "Historical internal migration and holiday mobility in Spain" (style report's first choice); "Roots or beaches? Historical internal migration and holiday mobility in Spain" (keeps the coursework hook).
5. Write to the UCL Library open access team before submission: EPJDS is in the Springer Nature fully OA agreement and UCL is listed as "Selective Centralised", but UCL says SpringerOpen titles follow its normal eligibility rules and its public text for unfunded research only covers August 2025–July 2026.
