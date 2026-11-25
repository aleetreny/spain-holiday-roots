# AGENTS.md: everything a new agent needs to continue this project

Read this file first, then `.planning/2026-11-24-consolidation/decisions.md` (the full decision log, one line per decision, newest last). The user writes in Spanish and wants plain, simple Spanish in chat; the report is in British English.

## 1. What this is
- **Coursework** for UCL COMP0123 *Complex Networks* (module leader Prof. Shi Zhou), 50 % of the module. Afterwards, a **journal paper** (targets: EPJ Data Science, J. Transport Geography, Population, Space and Place).
- **Question:** *Roots or beaches? On holiday, do Spaniards go back to the villages their families left?*
- **Report:** `reports/manuscript/` (LaTeX). Build with `reports/manuscript/build.sh` (Tectonic/XeTeX; the binary lives in `_local/toolchain/`, not in git; any `tectonic` on PATH works). Output: `reports/manuscript/output/main.pdf`.

## 2. Data and networks (frozen spec: `pipeline/config.py`)
- **Overnight-stay networks:** MITMA *pernoctaciones* (mobile-phone estimates). Nodes = 2,435 zones of mainland Spain (7,975 municipalities; 908 zones group several). Directed link i→j when ≥ 10 residents of i sleep in j on one night. Check backbone `equal_edges`: every night keeps 35,812 directed links (= ordinary night).
- **Date convention:** file dated d = night ending on the morning of d. Always write nights as "14–15 August".
- **Main nights:** ordinary 2024-02-14, summer 2024-08-15, All Saints 2024-11-01; appendix Christmas Eve 2024-12-25, Good Friday 2024-03-29. Sister nights: 13 February nights ending Tue/Wed/Thu; 9 nights ending 10–18 Aug; 1–3 Nov. Years 2022 and 2023 (no All Saints 2023).
- **Roots network:** INE census 2024, 10 % microdata (675 municipalities > 10k that are their own zone, 44 provinces: birth province, arrival cohort, birthplace size) and the full table 68545 (137 mainland cities; 112 identifying in 22 provinces).
- **Coastal flag:** `coastal_true` from `explore/out_p86_coastal_flags.csv`. Never the buffer method of p73/p92 (misflags 158 inland zones).
- Heavy data are **not in git** (`data/`, `raw_data_cache/`, `clean_data/`, `_local/`, `*.parquet`). Download scripts: `explore/p80_*`, `p81_*`, `p83_*`, `p100_download_summer.py`, `_local/external/_scripts/` (local only). Sources and URLs are in `reports/manuscript/references.bib` and in the READMEs under `_local/external/`.

## 3. Results that survived (what the report says)
1. **Structure (T1, T2).** Summer network: undirected density 0.89 % → 3.44 %; in/out degree Spearman 0.69 → 0.09; reciprocity 0.52 → 0.24 (also 0.24 at equal directed links). Assortativity −0.131 → −0.214; degree-preserving rewiring gives −0.140 / −0.217 and reproduces 93 % of the change (Zhou & Mondragón 2007 range [−0.22, −0.09] → [−0.42, −0.11]). Rich club top-50 raw φ 0.80 → 0.96 (rewired 0.66 → 0.92). Only property beyond the null: second-order mixing R_avg −0.289 vs −0.252 in summer. BOSAM: same shape on all nights, denser (appendix).
2. **Roots (T3).** PPML, city FE + origin×destination-province FE, log-distance quadratic, income×destination. Pooled β 0.53 (Feb) / 0.63 (Aug), 112 cities. Average Aug increase is **baseline-dependent**: +0.096 (p = 0.18) vs February, +0.137 (p = 0.006) vs October. By groups (675 cities, +0.5 to counts): pre-1981 arrivals +0.084, village-born +0.086, recent arrivals −0.072 (weaker: Oct p = 0.06, 2022 n.s.). Inland only: +0.154 / +0.123; coast not detected. Age check: post-1981 65+ −0.006. 117/117 sister pairs same sign; jackknife 44 provinces. Signs stable over 6 zero-cell rules; sizes not (Table A2).
3. **Villages (T4).** Most emptied fifth: net receivers 24 % (Feb) → 99 % (Aug) → > 99 % (Christmas Eve). Regression: +39 pp (CI 23–55), +44 pp in single-municipality zones. August mechanism: fewer residents leave; Christmas: net gain per resident rises. Non-primary dwellings reduce the emptying coefficient by 48 % in August, 2 % at Christmas (an attenuation, not mediation; no claim about where visitors sleep).
4. **Did not survive** (reported in 4.4): k-core promotion of villages (reversed; rewiring reproduces it) and "summer territories / diaspora archipelagos" (communities exist beyond gravity, Q 0.170 vs 0.115, but are *more* compact than the null: 71 % vs 99 % non-contiguous).

## 4. Writing rules (user-approved)
- Voice guide: `~/Proyectos/Spain-Mobility/manuscript_planning/voice/AUTHOR_VOICE.md` (outside this repo). First person singular, plain words, concrete before abstract, numbers with units and sources, no AI vocabulary, few em-dashes. The user said "el tono me gusta mucho": keep it.
- Classic LaTeX article look (Latin Modern, running head, abstract block). Structure and page limits follow `CW_Instruction_2026-27.pdf` exactly: cover (title + abstract); 1 question (2–3 pp); 2 literature (1–3 pp, 5 works, critical); 3 method (1–3 pp); 4 results (2–5 pp); 5 discussion (1–3 pp); 5–20 references; appendix optional.
- Captions short, Nature style (bold title + one line per panel). Long captions kept in `reports/figures/captions.md`.
- No specific town as an opening (the examiner does not know Spanish zones).
- Never claim individuals return to *their* village: the data link cities to provinces. Say "not detected", not "zero". No causal wording.

## 5. Code map
| Path | What |
|---|---|
| `pipeline/config.py` | Frozen spec: nights, thresholds, seeds (SEED = 20261124) |
| `pipeline/b00_common.py`, `b01_build_networks.py` | Zones, edges, backbones |
| `pipeline/b02_basic_properties.py`, `b02_zhou.py` | Basic properties, nulls, Zhou methods (rich club, structural constraints, BOSAM, R_avg) |
| `pipeline/b03_knn.py`, `b04_bosam_envelopes.py` | knn(k) with rewiring; BOSAM envelope test (degenerate: Madrid) |
| `pipeline/r1_villages.py`, `r2_roots.py`, `r3_territories.py`, `v1_verify_n1.py`, `n1_*` | Villages, roots, territories, generation split |
| `pipeline/f_figures.py` | All figures (`fig1b fig2 fig3 fig4 figbasic figbosam`) |
| `pipeline/t01_table1.py`, `t02_zero_cells.py` | Tables 1, A1, A2 |
| `results/*.md` | One report per analysis, every number traced to a CSV |
| `research_notes/verification_independent.md` | Independent re-computation of every number (Codex gpt-6.2-sol). Examiner-style reviews are kept privately in `_local/private_notes/` |
| `explore/` | Legacy pilots, read-only |
| `tests/` | `pytest tests/` |

Run Python only as `.venv/bin/python -I` (downloaded data are untrusted). Environment: `uv venv && uv pip install -r requirements.txt`.

## 6. Rules
- Commit or push only when the user asks.
- `_local/` is private and never in git. **Never share or publish the past students' reports in `_local/examples`.**
- Codex delegates run at reasoning effort medium or higher. Antigravity instance: `antigravity_2`.

## 7. Open items
- Submission date on the cover; repository URL in the report (this repo, once public).
- Check DOIs of `clauset2009power` and `kline2012score`.
- Optional (review recommendations): swap the 2020 rural-migration study for García-González et al. 2024 in the literature; short All Saints comparison; trim methodology to ≤ 3 pages.
- Journal paper: 2025 data, municipal INE OD network (`n1_*`), village-level roots (ECEPOV), multi-year cohort fading.
