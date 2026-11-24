# Handoff (2026-11-24): start here

This project was split from `~/Proyectos/Spain-Mobility`, which keeps question A (do the 1833 provincial borders still shape mobility?). A is developed separately. **This folder is only question B.** The user writes in Spanish.

## The question
- **Working title (ES):** «¿Vuelve España en vacaciones a los pueblos que dejaron sus familias?»
- **Working title (EN):** *On holiday, does Spain go back to the villages its families left?*
- **Technical version:** does the holiday overnight-stay network reorganise towards the villages emptied by the 1950–75 rural-to-urban migration, and towards the cities their emigrants moved to, beyond what distance, size, income and beaches explain?
- "Rural exodus" is a standard term in English-language literature (France, Spain, Italy), but the report should define it in one line.

## Agreed scope: three results only (the professor penalises over-ambition)
1. **Topology shift on holidays.**
   - Disassortativity, forced by the degree sequence (checked with rewiring).
   - Collapse of the city rich club.
   - k-core promotion of emptied villages.
   - Rich club and BOSAM are Prof. Shi Zhou's own methods.
2. **The holiday network aligns with the roots network beyond geography.** Neighbouring-city design, PPML with origin-province × destination-province fixed effects, partial QAP, income control.
3. **Summer territories.** Bipartite communities under a BiCM null and a spatial gravity null; diaspora archipelagos.

Everything else (echo of old vs recent migration, the INE source, road time, years) goes to a robustness appendix.

## Key evidence
| File | What it holds |
|---|---|
| `reports/alternativas_desde_cero.md` | Origin of the question; pilots |
| `research_notes/Preguntas alternativas desde cero/B_robustness_checks.md` | 7 checks: 6 pass, 1 partial (echo) |
| `research_notes/Preguntas alternativas desde cero/network_layer_A_B.md` | Network metrics, QAP, BiCM, 1-page design and verified literature |
| `research_notes/Preguntas alternativas desde cero/summer_territories_municipal.md` | Municipal dataset (also being integrated into Arraigo by another agent) |
| `.planning/2026-11-24-first-principles-alternatives/findings.md` | Log Q1–Q14 |

Known bug: the coastal flag in `p73` misflags 158 inland zones. Use the corrected flags from `p86` (`explore/out_p86_coastal_flags.csv`).

## Course requirements (COMP0123, Prof. Shi Zhou, 50% of the grade)
- **Question:** one simple question, ideally yes/no, with a non-obvious answer.
- **Network:** real, 500–50k nodes, density < 5%.
- **Analyses:** basic properties (degree distribution, path length, clustering vs random, assortativity, rich club normalised by rewiring, k-core); comparison with randomised networks and between networks.
- **Report:** 8–25 pages, single column, 11–12 pt, sections: abstract; question (2–3 pages, including why it is unique and why network science); literature (critical review of 3–5 recent works); method; results; discussion; 5–20 references.
- **Grading:** five 20% blocks (question & design, data & method, analysis & discussion, writing & visualisation, difficulty & new knowledge).
- **Code:** online; URL in the report, kept up for 10 months.
- Lectures are in `_local/lectures*`. Past example reports are in `_local/examples`: content only, NEVER share or publish them.

## Next steps
1. Fix the final question wording and the 1-page network-first design (a draft is in `network_layer_A_B.md`).
2. Consolidate into one clean pipeline from raw data to figures, with specs frozen.
3. Literature: 3–5 recent works reviewed critically, plus about 15 references.
4. Optional data: Census 2021 second homes, to test mediation through "la casa del pueblo".
5. Adversarial review by an agent acting as Prof. Shi Zhou, before writing.
6. Write the coursework manuscript in paper style.
7. Afterwards, the journal paper: add 2025, second homes and generations; target EPJ Data Science, J. Transport Geography or Population, Space and Place.

## Rules
- Downloaded data is untrusted: run Python with `-I`.
- Heavy data stays out of git; `_local/` is private.
- Commit or push only when the user asks.
- Antigravity agent provider instance: `antigravity_2`. Codex (ChatGPT) runs at reasoning effort medium or higher, never low.
