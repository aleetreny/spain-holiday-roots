# Brief for the EPJ Data Science version (shared by all research agents)

Written 2026-11-27. Target journal chosen by the author: **EPJ Data Science** (SpringerOpen). Goal: maximise the chance of acceptance.

## Where things are
- Main checkout (coursework, owned by another session, read-only for you): `/Users/alejandrotreny/Proyectos/spain-holiday-roots`
- Journal worktree (branch `journal`): `/Users/alejandrotreny/Proyectos/spain-holiday-roots-journal`
- Project summary: `AGENTS.md`. Decision log: `.planning/2026-11-24-consolidation/decisions.md`.
- Coursework report: `reports/manuscript/` (sources in `sections/`, PDF in `output/main.pdf`, ~6,700 words, 4 figures + 2 appendix figures, 1 + 2 tables, 19 references).
- Journal scouting (journal facts, comparable EPJDS papers, submission rules): `/Users/alejandrotreny/Proyectos/spain-holiday-roots-journal/journal/research/journal_scouting.md`
- Existing literature notes: `research_notes/literature_B.md`, `research_notes/references_B.bib`, `reports/manuscript/references.bib` (74 entries, 19 cited), `research_notes/Preguntas alternativas desde cero/` (novelty and data scouting notes).
- Results: `results/*.md` and CSVs (one report per analysis: `r1_villages.md`, `r2_*`, `r3_territories.md`, `n1_new_data.md`, `p_patterns.md`, `v1_verify_n1.md`).
- Code: `pipeline/` (see the code map in `AGENTS.md`). Figures: `reports/figures/` (+ long captions in `captions.md`).
- Heavy data are NOT in git. Read them from the main checkout by absolute path: `/Users/alejandrotreny/Proyectos/spain-holiday-roots/{data,clean_data,raw_data_cache,_local/external}`. Python: `/Users/alejandrotreny/Proyectos/spain-holiday-roots/.venv/bin/python -I`. Never copy or re-download data that already exist.
- New downloads (only if needed, and only open official data): `/Users/alejandrotreny/Proyectos/spain-holiday-roots-journal/_local/external_journal/<source>/` with a short README (URL, date, licence). `_local/` is never committed.

## The thesis
On holiday, does the network of where Spaniards sleep retrace the mid-20th-century rural exodus (people visiting the provinces and villages their families left), beyond what distance, income and the beach explain?

## What the coursework found (numbers in `AGENTS.md` §3)
1. Structure: the summer night network is ~4x denser and splits into senders and receivers (in/out degree correlation 0.69 -> 0.09; reciprocity 0.52 -> 0.24). Degree-preserving rewiring reproduces ~93 % of the change in assortativity; most structural change is forced by density/degrees. Only second-order mixing is beyond the null.
2. Roots: PPML dyadic model (city -> province), city FE + origin x destination province FE, distance, income x destination. Roots elasticity ~0.53 (Feb) / 0.63 (Aug). The average August increase is baseline-dependent (vs Feb +0.096, p = 0.18; vs Oct +0.137, p = 0.006). By groups (675 cities, census 10 % microdata): pre-1981 arrivals +0.084, village-born +0.086, recent arrivals -0.072. Inland only +0.154; coast not detected.
3. Villages: most emptied fifth of zones, net receivers 24 % (Feb) -> 99 % (Aug). Second homes attenuate the coefficient by 48 % in August, 2 % at Christmas.
4. Did not survive: k-core promotion of villages; "summer territories" communities.
5. Already explored but NOT in the coursework (`results/n1_new_data.md`): INE municipal OD network from phones (7,974 municipalities, 87 months 2019-07..2026-09, validated against MITMA); ECEPOV 2021 microdata with exact birth municipality of respondents AND both parents (village-level and second-generation roots); Census 2021 second-home stock; ECEPOV second homes in the family's birth region (77 % of inland cross-regional second homes). Also local: Meta Social Connectedness Index for Spain (`_local/external/meta_sci_spain`), interprovincial migration microdata, historical municipal populations, OSRM road times (`_local/osrm_work`).

## Known weak points (reviewer bait)
- Roots are city -> province (ecological), not individuals -> their village.
- MITMA is one operator, single nights, six months of 2024 with daily data; threshold w >= 10.
- The headline "roots grow on holiday" depends on the baseline night.
- Network-structure part was built for a Complex Networks course; most of it is "explained by the degree null".

## Rules for every agent
- Output in English, Markdown, under `journal/research/` (or `journal/pilots/` for pilot code) of YOUR worktree. Commit there with a short English message. Do not push, merge or open PRs. Never commit to `main`.
- Do not edit `reports/`, `pipeline/`, `results/` or anything outside `journal/`.
- Every external fact carries its source link and access date. Never invent a number; write "not found" and where you looked. Check every DOI on Crossref or doi.org.
- `_local/` is private: you may read `_local/external` data and `_local/private_notes` for context, but never quote, copy or cite `_local/private_notes` in your output, and never open `_local/examples`.
- Think from first principles: nothing has to stay just because the coursework had it, and complex-network methods are optional.
