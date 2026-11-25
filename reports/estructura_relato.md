# Structure and narrative of the coursework report (COMP0123)

*Working proposal, 2026-11-24. Coursework only; the paper comes later. The numbers come from `results/` (b02, r1, r2, v1) and must be checked again when written into the text.*

## 1. The story in one breath
Every August, Spain says it *vuelve al pueblo* ("goes back to the village"). It is also the second most visited country in the world, and its summer is synonymous with the coast. Two forces pull on the same nights: the beach and the roots left by the great rural exodus (1950–1975).

Looking at where 48 million residents sleep, night by night, as a network, we see both forces at work and they do not get in each other's way:
- **The beach changes the shape of the network.** In August it splits into cities that send and coasts that receive.
- **The roots change who connects with whom.** The generation that emigrated before 1981 and the people born in villages pull, precisely in summer, towards the inland provinces they came from.
- The emptied villages, which send people away all year, become receivers.

Answer to the question: **yes, and it is not the beach that brings them back.**

## 2. What the professor wants (and how we give it to him)
| He asks for (instructions, lectures, papers) | How it appears in the report |
|---|---|
| One specific question, a non-obvious answer | Title and question; the beach as the obvious alternative |
| "Study a network and apply network science methods" | Three networks on the same nodes; basic properties; real against randomised; comparison of networks |
| Basic properties always (P(k), L, C, mixing) | Table 1 in his NJP 2007 format: N, L, ⟨k⟩, k_max, α, φ, ℓ vs ℓ*, rich-club clique |
| His methods: rich club, structural constraints, BOSAM, 2nd-order mixing | Result 1 rests on them, with his caveat about the null model of the rich club |
| Comparison (two networks, real vs random, two methods) | Feb vs Aug vs All Saints; real vs rewired; gravity null; PPML and QAP |
| Evaluation of correctness | Sister nights (117 pairs), jackknife by province, the age test, and **what did not survive** |
| Literature: 3–5 recent works, critically | 5 works with "what we take / where it falls short" |
| Not over-ambitious | 3 results, one figure per result in the main text; everything else in the appendix |
| Captions understandable on their own; large, readable figures | Long captions; one column; figures at full width |

## 3. Structure (target ~16–18 pages + appendix)

**Cover + abstract (½ p).** Five sentences: scene, question, network, three findings, answer.

**1. The question (2–3 pp)**
- 1.1 *Opening scene* (½ p). A concrete night: the night of 14 to 15 August 2024. How many people slept away from home [b01: total weight], where Madrid's residents went, and a real village that empties in February and fills in August. No numbers are made up: all of them come out of the data.
- 1.2 *Two forces.* The beach (tourism, coast) against the roots (the 1950–75 exodus, with cited figures: Clar & Ayuda 2023; Pinilla & Sáez 2017). Why it is not obvious which one wins.
- 1.3 *The network.* Node = zone (2,435); link i→j = residents of i who sleep in j on a given night. Three nights. The roots network (birthplace × residence). Why a network and not a table: the question is who connects with whom and how the structure changes between two networks on the same nodes.
- 1.4 *Question and tasks.* The question, plus 4 concrete tasks (the instructions ask for "specific tasks"): T1 properties and nulls; T2 structural change; T3 alignment with roots, by generation; T4 the role of the villages.
- 1.5 *What is new.* Li et al. 2016 (China, provinces, no spatial null) and García-González et al. 2024 (descriptive, one region). No one has compared the network of where people sleep with the network of where they were born.

**2. Literature (1.5–2 pp).** Five works, each with a paragraph *what they do → what we take → where they fall short*:
- Li et al. 2016 (closest precedent);
- Wesolowski et al. 2017 (holidays reorganise the network);
- Luan et al. 2024 (distance is weaker on homecoming holidays → a gravity null for each date);
- García-González et al. 2024 or González-Leonardo et al. 2022 (return and second homes in Spain);
- Zhou & Mondragón 2007 (structural constraints: the base of Result 1).

**3. Data and method (2–3 pp)**
- 3.1 *Data.* MITMA overnight stays: what is measured, the date convention, the w ≥ 10 backbone and the density. INE census 2024 (microdata: birthplace, year of arrival in the province). INE housing census 2021.
- 3.2 *Basic properties*: **Table 1**, three nights in Zhou's format plus ER and rewired values. One paragraph: small world (ℓ ≈ ℓ*), log-normal and not a power law, sparse.
- 3.3 *Methods*, one paragraph each:
  - degree-preserving rewiring (and maximal assortative/disassortative, his 2007 method);
  - rich club by rank;
  - BOSAM;
  - gravity null;
  - the neighbouring-city design (PPML with fixed effects), with the bootstrap by province;
  - sister nights.

**4. Results (5 pp): three acts**
- **4.1 The beach: in August the network splits in two.**
  - Density ×4; in/out correlation 0.69 → 0.09; the receiving hubs are coastal (36 % → 80 %).
  - Disassortativity is not a choice: P(k) forces it (r real inside the [r_min, r_max] range; rewiring reproduces it).
  - Rich club: the raw value **rises** (φ top-50 0.80 → 0.96), even though normalised it is ~1 (his caveat).
  - 2nd-order mixing negative in August beyond the null. BOSAM Feb vs Aug.
  - → *The shape of the summer network is explained by the beach.*
- **4.2 The roots: the paths of the exodus.**
  - The roots network predicts where people sleep **every night** (neighbouring cities, beyond distance and income).
  - On average, holidays add little (+0.10, not significant).
  - **The twist:** splitting by generation, the pre-1981 cohort and the village-born pull much harder in August, all of it inland, none towards the coast. Recent migrants do the opposite.
  - The age test: post-1981 65+ = 0.
  - 117/117 sister pairs; jackknife.
  - → *It is not the beach that brings them back; it is the generation that left.*
- **4.3 The villages: from senders to receivers.**
  - Emptied villages switch from exporting to importing nights (23 % → 99 %).
  - In summer this is mostly because their residents stay at home.
  - Second homes account for about half of it in August and almost nothing on Christmas Eve.
- **4.4 What did not survive (½ p).** k-core promotion (an artefact of zone aggregation) and summer territories (do not beat a fair null). Two sentences each, to show the evaluation.

**5. Discussion (1.5–2 pp).**
- The answer.
- What surprised us: that it is not the beach, and that recent migrants work the other way round.
- New knowledge: two calendars of return.
- Limitations: phone data and residence assignment; zones that group municipalities; roots by province, not by village; cross-section, not causality; one year with all nights (plus 2022–23).
- Advice: use single nights, not averages; be careful with normalised rich clubs.

**References (15–20).** **Appendix**: full robustness tables, the remaining nights, the territories, the code URL.

## 4. Voice: how to make it not sound like AI
- **Prose, not bullet points**, in the main text. Each paragraph says one thing and ends on the idea that leads into the next.
- **Concrete before abstract.** Every result opens with a real case (a city, a village, a corridor such as Santa Coloma → Jaén) before the coefficient.
- **Section titles that are statements** ("In August the network splits in two"), not labels ("Results 1").
- **Admit doubts and failures.** Section 4.4 and the limitations, in the first person: "I expected X; the data said Y".
- **Mixed sentence lengths.** A short sentence after a long one hits harder.
- **Avoid the tell-tale vocabulary:** *delve, crucial, landscape, tapestry, pivotal, it is worth noting, in conclusion, furthermore*. No triplets just because. Few dashes. No bold in the middle of the prose.
- **Spanish only where it adds something:** *volver al pueblo*, *España vaciada*, *casa del pueblo*, in italics, explained the first time.
- **The author's personal touch (optional):** one honest sentence about why the question matters to you. Only if it is true.
- **Every number exact and with its unit.** No "massive", "dramatic" or "striking".

## 5. Figures: proposed list (each to be worked on separately)
A single visual system for all of them:
- the same serif typeface as the text;
- the same map projection;
- **colour with meaning, the same everywhere:** grey = ordinary night; **sea blue = beach/coast**; **terracotta = roots/inland**; dark ink for annotations.
- Colour then tells the story by itself.

| # | Figure | What it must convey | Section | Priority |
|---|---|---|---|---|
| F1 | **Hero map**: Spain, ordinary night vs 15 August; links from the cities, coloured by roots vs beach | The story at a glance | 1 (front page) | High |
| F2 | Degree distributions (CCDF in/out, three nights, log-log) | Basic properties, lecture style | 3.2 | Medium |
| F3 | **BOSAM Feb vs Aug** (his tool) | The network changes its "fingerprint" | 4.1 | High |
| F4 | k_in vs k_out, Feb vs Aug (coloured by coast) + the r range against [r_min, r_max] | Senders vs receivers; P(k) forces the mixing | 4.1 | High |
| F5 | Raw and normalised rich club by rank | His caveat, applied | 4.1 or appendix | Medium |
| F6 | **Roots by generation**: δ for pre-1981 / village-born / recent / urban-born, inland vs coast, across nights | The twist of the story | 4.2 | Highest |
| F7 | Robustness: δ over the 117 pairs + jackknife | Evaluation | 4.2 or appendix | Medium |
| F8 | **Villages**: share of net importers by emptying level, across nights; small map | From senders to receivers | 4.3 | High |
| F9 | Territories: real vs null | Honest evaluation | Appendix | Low |

Tables: T1 basic properties (Zhou format); T2 summary of the three results with their null, effect, CI and robustness.

## 6. Decisions taken (user, 2026-11-24)
1. **"I"** throughout the text.
2. **Personal opening: Vélez-Málaga.** The user's real connection; the details of the connection are to be written by the user and must not be made up. What the data say (MITMA, zone 29094, single nights, full matrix):
   - **14 Feb:** 3,531 visitors from 260 zones, 55 % from the province of Málaga itself.
   - **15 Aug:** 27,621 visitors (×7.8) from 803 zones, only 22 % from Málaga. The top origins are Málaga, Madrid, Córdoba, Granada, Jaén and Sevilla: inland Andalusia going to the beach.
   - **Christmas Eve:** 11,859 visitors, 48 % from the province itself: family.
   - Vélez is both beach and "home", depending on the night. It is the question in miniature: in August it is a beach for inland Andalusia, at Christmas it is home. A good hook for §1.1.
3. **Hero figure: the map (F1).**
