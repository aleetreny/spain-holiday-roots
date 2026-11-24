# Review: Santander Cycles Network vs. Spain Mobility Candidates (A & B) for UCL COMP0123

*Independent review for UCL COMP0123 Complex Networks & Online Social Networks (Prof. Shi Zhou, 50% of module grade). Date: 2026-11-24.*

---

## Executive Verdict

**No — pivoting to `santander-cycles-network` is NOT a safer path to a 10 than Candidate A (1833 borders), and it carries a significantly lower ceiling on novelty and publishability than both A and B.**

- **What the Santander repo gives you:** A genuinely impressive, already-completed 14-chapter interactive portfolio essay (*"London, as the bikes see it"*), built on clean, fast, laptop-reproducible code (6 min runtime), covering $N = 796$ stations and 857,151 May 2026 journeys. Turning one focused slice of it — specifically **S1: *"Do e-bikes rewire London's bike-share network?"*** — into a single-question COMP0123 report would be **purer network science** than A's econometric core and very feasible (low engineering risk).
- **Why it is NOT a safer 10:**
  1. **Crowded UCL backyard:** London Santander Cycles / Barclays Cycle Hire is one of the most heavily studied urban mobility datasets in the world — including by UCL's own CASA faculty (O'Brien, Cheshire & Batty 2014) and Netzschleuder's benchmark `london_bikes` dataset (Munoz-Mendez et al. 2018). Passing the *"Difficulty & Achievement: new knowledge"* 20% rubric bar is structurally harder when the examiner has almost certainly seen past London bike-share reports.
  2. **The density trap:** The raw monthly bike graph has **33.6% density** ($\langle k \rangle = 267$), violating Prof. Zhou's explicit `<5%` guideline by $6.7\times$ (requiring backbone filtering or sparse sub-layers to justify). Spain-Mobility naturally sits at **$N = 2{,}567$, density 3.8%**.
  3. **Unverified core risks in the bike candidates:** The two strongest bike findings currently in the repo — e-bike backbone Jaccard = 8% (**S1**) and rider-vs-van net balance $\rho = -0.97$ (**S2**) — have **major unverified confounds**: S1 has a $4.1\times$ sample-size disparity (167k e-bike vs 691k classic trips) that mechanically depresses reciprocity and disparity-filter overlap; S2's $\rho = -0.97$ is largely a bike-conservation accounting identity and mixes van rebalancing with 50 m staff shuffles and workshop repairs.
  4. **Existing public portfolio optics:** The repo is already live on GitHub Pages (Nov 2026) as a 14-chapter general showcase. Reusing it requires throwing away the 14-chapter tour (which would trigger Prof. Zhou's *"over-ambitious"* penalty) and building fresh null models (volume downsampling, spatial nulls, rewiring constraints, 2nd-order mixing, BOSAM) around a single question.

**Bottom line:** Keep **Candidate A (1833 borders)** — ruthlessly trimmed to one clean COMP0123 network narrative — as the primary recommendation (already computed, highest new-knowledge ceiling, publishable), keep **Candidate B (Rural exodus → holiday map)** as the highest-clarity challenger pending its income-control test, and hold **S1 (E-bikes rewiring London)** strictly as a **low-stress, network-pure Plan C** if you want to avoid econometric/historical framing altogether (contingent on a 1-hour volume-downsampling check).

---

## 1. Audit of `aleetreny/santander-cycles-network`

Cloned read-only (`HEAD` `e501160`, committed `2026-11-04`). No code was executed; all metrics below come from inspecting the source code (`cycles/`), documentation (`README.md`, `PRODUCT.md`, `docs/STUDY.md`, `docs/DATA_DICTIONARY.md`), and published output JSON (`site/data/study.json`).

### 1.1 Data & Network Profile

| Attribute | Value / Detail | Fit with COMP0123 Guidelines |
|---|---|---|
| **Source** | TfL Cycle Hire extracts 443–444 (1–31 May 2026: 19 working days, 12 weekend/bank-holiday days); BikePoint snapshot (20 Oct 2026); OSM Thames/parks; Sep 2026 BikePoint occupancy (`feregrino/london-cycles`, CC BY-SA 3.0) | Clean, open, fully reproducible; **caveat:** only 1 month of journeys (May 2026); occupancy validation is cross-month (Sep vs May). |
| **Scale** | 863,014 valid hires → **857,151 matched** (99.32%); 11,531 unique bikes; **166,628 e-bike hires (19.4%)** | Manageable on a laptop (~6 min full pipeline). |
| **Nodes ($N$)** | 801 BikePoint docks (**796 active**; 6 retired/renamed docks dropped) | **Ideal** (inside Zhou's 500–50,000 range). |
| **Edges & Density** | **212,532 directed** (142,446 undirected); **density = 33.59%** ($\langle k \rangle = 267$, diameter = 3, $L = 1.55$, 1 giant SCC) | **FAILS raw `<5%` guideline** ($6.7\times$ too dense). Binary metrics saturate (`structure.py` explicitly acknowledges this). |
| **Sparse Backbones / Sub-layers** | Disparity filter $\alpha=0.05$: 16,228 directed edges (**2.56% density**, 796 nodes); $\alpha=0.01$: 7,011 edges (**1.11%**, 796 nodes); $\alpha=0.001$: 2,578 edges (**0.41%**, 767 nodes). Inferred relocation graph $R$: 37,139 moves across 796 docks (naturally sparse). | **Passes `<5%`** once properly filtered or focused on sparse layers, but requires explicit justification. |

### 1.2 State of Completion & Methods Already Implemented

The repository is **100% complete as a portfolio interactive essay** (with unit tests in `tests/` and a live D3/Canvas web app), structured as 14 short chapters rather than a single research question:

1. **Graph basics & weighted reciprocity:** Binary reciprocity 0.660, weighted 0.620 vs Poisson strength null 0.570 ($\rho = 0.117$). Self-loops 3.75% (excluded from topology).
2. **Heterogeneity:** CSN MLE + Vuong tests show station strength is **not** scale-free ($p = 0.055$ vs log-normal; dock slot count caps hubs; Gini = 0.293), whereas link weights are log-normal (Gini = 0.540; 13.4% of links carry 50% of trips; 39.6% are single-trip links). Superlinear $s \propto k^{1.256}$ ($R^2 = 0.795$).
3. **Disparity backbone (Serrano et al. 2009):** At $\alpha = 0.001$, keeps 767 stations vs 660 for an equal-edge weight threshold. At $\alpha = 0.05$, small-world ($C = 0.248$ vs $C_{\text{ER}} = 0.037$, $L = 2.56$ vs $L_{\text{ER}} = 1.97$).
4. **Spatial models:** Doubly constrained gravity via IPF + Poisson MLE: exponential $\beta = 0.569\text{/km}$ gives CPC = 54.6% (power law 49.6%; radiation fails at 14.0%). **River Thames term:** $\theta = 0.265$ ($\approx 0.479\text{ km}$ penalty); below 1 km, cross-river pairs get 70.6% of gravity prediction vs 110.7% same-bank. Standardised residuals flag corridors (e.g., Great Dover St halls → LSE/Kingsway, Hyde Park loops).
5. **Centrality:** PageRank, HITS, eigenvector, and effective-distance ($1 - \ln P_{ij}$) betweenness/closeness; strength $\approx$ PageRank ($\tau = 0.929$), betweenness diverges ($\tau = 0.537$, picking bridges like St James's Square, Bank, Waterloo).
6. **Communities & roles:** Leiden sweep ($\gamma \in [0.25, 4]$, 24 seeds) gives **4 regions** at $\gamma=1$ ($Q = 0.374$, $z = 1565$ vs Chung–Lu null $0.0168$), corroborated by Infomap (4 modules, NMI = 0.812). At $\gamma=0.5$ and by normalised-Laplacian eigengap, London splits **East vs West** (not North vs South across the Thames). Guimerà–Amaral roles: 20 connector hubs, 3 provincial hubs.
7. **Higher-order structure:** Opsahl weighted rich club reaches **$\rho = 7.26$ for top 10 stations** (weight reshuffling on fixed topology); neutral mixing ($r_{\text{deg}} = -0.013$, $r_{\text{str}} = +0.046$); directed triad census on $\alpha=0.01$ backbone vs 40 rewirings finds fully reciprocal cliques (`300`) massively over-represented ($z = 1072.8$).
8. **Temporal layers & archetypes:** 48 hourly layers (work/free), 5 regimes, commute reversal (station balance $r = -0.933$, reversed-evening overlap 36.2% vs same-direction 16.5%); 6 k-means rhythm archetypes (weak silhouette 0.065; rail termini act like residential sources).
9. **Markov dynamics & out-of-sample validation:** $|\lambda_2| = 0.790$ (mixing time 4.25 rides). May weekday hourly net flows predict September weekday BikePoint dock changes ($r = 0.665$ pooled across 18,285 station-hours; median station $r = 0.804$; slope 0.476; stress vs empty/full $\rho = 0.662$).
10. **Hidden rebalancing network:** From 845,620 consecutive bike-ID handoffs across 11,531 bikes, **37,139 spatial discontinuities** (4.39%) are inferred as operator relocations (median 1.92 km, median gap 14.8 h). Rider vs operator station net flows correlate at $\rho = -0.968$, cancelling 87.4% of rider imbalance.
11. **Robustness & weak ties:** Targeted strength removal shatters the $\alpha = 0.001$ backbone at ~25% closed; Onnela topological overlap rises with weight ($0.066 \to 0.182$) and weak-first percolation breaks the backbone sooner than strong-first.
12. **Multiplex layers:** Classic (690,523 trips) vs e-bike (166,628 trips): e-bikes travel farther (median 2.93 vs 1.94 km), decay half as fast ($\beta = 0.337$ vs $0.645\text{/km}$), show collapsed reciprocity (0.369 vs 0.609) and lower modularity ($Q = 0.252$, 3 communities, NMI to global = 0.434 vs 0.931 for classic), share only **8.0%** of $\alpha=0.05$ backbone links (Jaccard = 0.0795) and 20.3% weighted edge overlap, despite high station-strength correlation ($r = 0.842$).
13. **Embeddings:** t-SNE flow space, force-directed, spectral (Procrustes-aligned to geography; flow-vs-geo Spearman 0.939).

### 1.3 What Is MISSING for Prof. Shi Zhou's COMP0123

Despite the impressive breadth, the repo currently **lacks** several hallmarks of a top COMP0123 report:
- **No single focused research question:** A 14-chapter tour is the exact "broad exploratory sweep" Prof. Zhou cautions against (*"single specific question"*, *"not over-ambitious"*). (Compare past example `2016_TRIPON_S.pdf` on London Underground, which covered basic metrics + day/time comparisons at MSc level, vs focused hypothesis tests like `2014_GRIGAT_D.pdf`.)
- **Missing Prof. Zhou's signature methods:** No **degree-preserving rewiring structural-constraint bounds** (max-assortative / max-disassortative / $P(k)+\phi(k)$ rewiring; Zhou & Mondragón 2007), no **second-order assortative mixing** ($R_{\text{avg}}, R_{\text{max}}$; Zhou et al.), and no **BOSAM** visual adjacency fingerprint.
- **Methodological gaps in key findings:**
  - *Rich-club:* `structure.py` uses Opsahl weight reshuffling on fixed topology, **ignoring geography** — because top-10 docks cluster in Central London, $\rho = 7.26$ is likely heavily inflated by spatial proximity.
  - *E-bike vs classic:* No **volume-matched downsampling** (classic has $4.1\times$ more trips!), no elevation/cycling-route distance control, and no pricing/availability confound analysis.
  - *Rebalancing:* No topological comparison of the operator graph $R$ against the rider graph $A$ (only node-level net-balance correlation), and no filtering of workshop transfers / 50 m adjacent-dock shuffles / retired-dock gaps.
  - *Single month:* May 2026 only (easy to expand to 3–12 months via existing `cycles/acquire.py`).

---

## 2. The 3 Best Simple Yes/No Research Questions on Santander Cycles & Novelty Check

### Background: How Crowded Is London Bike-Share Network Science?

Very crowded. Classic network-science studies on London's cycle hire scheme (formerly Barclays Cycle Hire, now Santander Cycles) already cover spatial communities, centrality, temporal volatility, transit-strike disruptions, gender/group cycling, and the River Thames barrier:

- **Austwick, O'Brien, Strano & Viana (2013)**, *PLOS ONE*, [10.1371/journal.pone.0074685](https://doi.org/10.1371/journal.pone.0074685) — spatial networks, scaling, and community detection in London + 4 global cities.
- **O'Brien, Cheshire & Batty (2014)**, *J. Transport Geography*, [10.1016/j.jtrangeo.2013.06.007](https://doi.org/10.1016/j.jtrangeo.2013.06.007) — station temporal archetypes & global classification including London (**UCL CASA**).
- **Beecham & Wood (2014)**, *Transp. Plan. Technol.* [10.1080/03081060.2013.844903](https://doi.org/10.1080/03081060.2013.844903) & *Transp. Res. Part C* [10.1016/j.trc.2014.03.007](https://doi.org/10.1016/j.trc.2014.03.007) — gendered and group cycling flows in London.
- **Munoz-Mendez, Klemmer, Han & Jarvis (2018)**, *UbiComp '18* (arXiv:1804.05584), [10.1145/3267305.3274156](https://doi.org/10.1145/3267305.3274156) — Louvain & Infomap community structures and hourly temporal dynamics on London Santander Cycles (archived as Netzschleuder's `london_bikes`).
- **Saberi et al. (2018)**, *J. Transport Geography*, [10.1016/j.jtrangeo.2017.11.018](https://doi.org/10.1016/j.jtrangeo.2017.11.018) & **Yang, Beecham & Heppenstall (2021)**, *J. Transport Geography*, [10.1016/j.jtrangeo.2021.103255](https://doi.org/10.1016/j.jtrangeo.2021.103255) — graph-based network changes during London Tube strikes.
- **Morton, Kelley, Monsuur & Hui (2021)**, *J. Transport Geography*, [10.1016/j.jtrangeo.2021.103125](https://doi.org/10.1016/j.jtrangeo.2021.103125) — spatial demand in London Santander Cycles, explicitly highlighting the **River Thames** south-bank deficit.
- **Beecham, Yang & Tait (2023)**, *Environment and Planning B*, [10.1177/23998083231165122](https://doi.org/10.1177/23998083231165122) — connected bikeability across London localities.
- **Corneck, Cohen & Martin (2025)**, *Statistics and Computing*, [10.1007/s11222-025-10606-w](https://doi.org/10.1007/s11222-025-10606-w) — network Poisson processes with community structure on Santander Cycles.
- Broader bike-share topology classics: **Borgnat et al. (2011)**, *Adv. Complex Syst.*, [10.1142/S0219525911002950](https://doi.org/10.1142/S0219525911002950); **Kou & Cai (2019)**, *PLOS ONE*, [10.1371/journal.pone.0213106](https://doi.org/10.1371/journal.pone.0213106) (local/non-local dynamics); **Xiao & Sha (2022)**, *Design Science*, [10.1017/dsj.2021.27](https://doi.org/10.1017/dsj.2021.27) (directed network motifs in bike-share).

Because generic "communities / centrality / temporal regimes of London bike-share" is thoroughly saturated, **only questions that exploit post-2022 fleet changes (e-bikes), hidden inferred layers (rebalancing), or sharp null-model tests (Zhou's structural constraints)** remain viable.

---

### Candidate S1 (Strongest Bike Question): *"Do e-bikes rewire London's bike-share network?"*

- **Technical RQ:** On the exact same 796 Santander Cycles docking stations, do electric bikes (19.4% of May 2026 hires; introduced late 2022 and scaled through 2024–2026) generate a structurally different mobility network — backbone links, communities, reciprocity, Thames permeability, rich-club, and assortativity — from classic bikes, **beyond sample-size sparsity and longer trip distances**?
- **Preliminary signal already in repo:** Classic (690.5k trips) vs e-bike (166.6k): median distance $+51\%$ (1.94 → 2.93 km), exponential decay halves ($\beta = 0.645 \to 0.337\text{/km}$), weighted reciprocity collapses ($0.609 \to 0.369$), Leiden modularity drops ($0.409 \to 0.252$, 4 → 3 regions, NMI to global = 0.434 vs 0.931), $\alpha=0.05$ backbone Jaccard = **8.0%** and edge overlap = 20.3%, yet station strength correlates at $r = 0.842$ ("same docks, different destinations").
- **Why the answer is NOT obvious (the real test):** E-bikes have $4.1\times$ fewer trips, which mechanically lowers both observed reciprocity and disparity-filter overlap ($\alpha = (1-p_{ij})^{k-1}$). Comparing the e-bike network against **100 volume-matched downsamples** of classic trips (166,628 trips each) AND a **distance-matched / spatial-gravity null** isolates genuine topological rewiring from sampling and range stretching. Adding Prof. Zhou's toolkit (degree-preserving rewiring, rich-club $\rho(k)$, 2nd-order mixing, k-core, and side-by-side BOSAM fingerprints of the two backbones) gives a textbook COMP0123 comparative study on two sparse networks (<2.6% density on $\alpha=0.05$ backbones).
- **Closest literature & novelty verdict (GAP in network topology; MODERATE overall novelty):**
  - **Schnieder (2023)**, *Sustainability*, [10.3390/su151813898](https://doi.org/10.3390/su151813898) — compares e-bike vs classic bike demand in London Santander Cycles, Philadelphia, and LA using Random Forests & DNNs (demand forecasting, **not** complex-network topology).
  - **Xie, Smart & Noland (2024)**, *Transp. Plan. Technol.*, [10.1080/03081060.2024.2341301](https://doi.org/10.1080/03081060.2024.2341301) — NYC Citi Bike e-bikes vs regular bikes via spatial CAR regression (station-level determinants, **not** multiplex/community/backbone structure).
  - **Grau, Bassolas & Vicens (2024)**, *EPJ Data Science*, [10.1140/epjds/s13688-024-00486-x](https://doi.org/10.1140/epjds/s13688-024-00486-x) — Barcelona Bicing e-bike vs m-bike mobility patterns for predictive maintenance.
  - **Li, Zhang, Luca & Fuerst (2024)**, *Transp. Res. Part D*, [10.1016/j.trd.2024.104179](https://doi.org/10.1016/j.trd.2024.104179) — dockless e-bike vs regular bike spatial clustering (finds e-bikes are more spatially dispersed).
  - **Yang et al. (2024)**, *J. Cycling Micromob. Res.*, [10.1016/j.jcmr.2024.100023](https://doi.org/10.1016/j.jcmr.2024.100023) & **Bernieri, Schweizer & Rupi (2025)**, *Sustainability*, [10.3390/su172411153](https://doi.org/10.3390/su172411153) — trip-level speed/range comparisons.
  - **Zhang & Cui (2024)**, *Physica A*, [10.1016/j.physa.2024.129962](https://doi.org/10.1016/j.physa.2024.129962) — multilayer dockless bike + metro (not classic vs e-bike).
  - *Verdict:* **No verified paper** performs a volume-controlled multiplex/backbone/community/rich-club/BOSAM network comparison of classic vs e-bikes on the same docked system. However, the prior intuition ("e-bikes go farther and smooth barriers") lowers the "wow", and e-bike pricing surcharges (£1/min-style tariffs vs classic) + fleet availability confound rider self-selection (must be discussed openly).

---

### Candidate S2 (Second Bike Question): *"Is London's hidden van-rebalancing network a structural mirror of its rider network?"*

- **Technical RQ:** Reconstructing the unobserved operator relocation network $R$ from 11,531 individual bike trajectories (37,139 spatial discontinuities across 845,620 handoffs), does the operator's rebalancing network merely reverse the rider flow network $A$, or does logistics produce a fundamentally different network architecture (hub-and-spoke disassortativity, rich-club, directed motifs, k-core, BOSAM)?
- **Preliminary signal already in repo:** Operator net deliveries correlate with rider net arrivals at Spearman $\rho = -0.968$, cancelling 87.4% of rider imbalance; median jump 1.92 km, median gap 14.8 h, peaking at morning rush hour; top moves are ultra-short dock-to-dock shuttles (Queen Street / Waterloo). Out-of-sample September occupancy validates rider drift ($r = 0.665$ pooled, median station $r = 0.804$, slope 0.476 showing ~52% damping).
- **Why it fits COMP0123:** Comparing two directed weighted networks on identical nodes where $R$ is **naturally sparse** (37k relocations across 796 docks $\implies$ naturally satisfies `<5%` density!). Uses Prof. Zhou's full suite (assortativity, rich-club vs rewiring, 2nd-order mixing, k-core, directed motifs, BOSAM) to test whether a human-engineered logistical mirror has the same topology as self-organised commuter flows.
- **Closest literature & novelty verdict (PARTIAL novelty; high measurement-noise risk):**
  - **Zhu, Zhang, Kondor, Santi & Ratti (2020)**, *Computers, Environment and Urban Systems*, [10.1016/j.compenvurbsys.2020.101483](https://doi.org/10.1016/j.compenvurbsys.2020.101483) — explicitly groups trips by `bike_id` and reconstructs hidden operator repositioning from spatial discontinuities between consecutive trips (dockless bikes & scooters).
  - Rebalancing optimization & inferred demand is a massive OR literature: **Regué & Recker (2014)**, *Transp. Res. Part E*, [10.1016/j.tre.2014.10.005](https://doi.org/10.1016/j.tre.2014.10.005); **Shui & Szeto (2020)**, *Transp. Res. Part C*, [10.1016/j.trc.2020.102648](https://doi.org/10.1016/j.trc.2020.102648).
  - *Verdict:* Comparing the *complex-network topology* of $R$ vs $A$ on Santander Cycles is fairly fresh, **but two major vulnerabilities weaken it:** (1) station net-balance $\rho = -0.97$ is nearly a **conservation-of-mass tautology** (a 30-dock station cannot accumulate +500 bikes over a month without vans removing ~500), so a sharp grader will dismiss it unless the paper proves non-trivial *edge-level* divergence; (2) discontinuities conflate van rebalancing, workshop maintenance (`Mechanical Workshop Penton`), 50 m staff wheeling, and the 6 unmatched/retired stations.

---

### Candidate S3 (Third Bike Question): *"Is London's bike-share rich club real, or an illusion of Central London geography?"*

- **Technical RQ:** The Santander Cycles network exhibits a striking $7.3\times$ weighted rich club ($\rho = 7.26$ for top 10 hubs) under standard non-spatial weight reshuffling alongside neutral degree assortativity ($r = -0.013$) and extreme reciprocal-triangle excess ($z = 1073$). Applying Prof. Zhou's **degree-preserving rewiring structural constraints** (Zhou & Mondragón 2007) alongside a **spatial configuration null** (doubly constrained gravity + Thames barrier), does hub-to-hub concentration survive once Central London dock proximity and dock-capacity caps are controlled for?
- **Why it fits COMP0123:** **Direct bullseye on Prof. Shi Zhou's own research** — Zhou & Mondragón (2004, *IEEE Comm. Lett.*), Zhou & Mondragón (2007, *NJP* structural constraints), Colizza et al. (2006, *Nat. Phys.*), Opsahl et al. (2008, *PRL*), 2nd-order mixing, and BOSAM. Also exposes a genuine methodological error in the existing repo (and many spatial network papers): non-spatial rich-club nulls generate false positives whenever hubs cluster in a city centre.
- **Closest literature & novelty verdict (LOW–MODERATE novelty; lower "wow"):**
  - **Zhang & Ng (2021)**, *Physica A*, [10.1016/j.physa.2021.126377](https://doi.org/10.1016/j.physa.2021.126377) — rich-club phenomenon in urban mobility passenger-flow networks.
  - **Ramasco, Colizza & Panzarasa (2009)** & **Colizza et al. (2006)** — normalised rich-club in transportation/mobility networks; **Cazabet, Borgnat & Jensen (2017)** (*CompleNet*) — degree-constrained spatial null models.
  - *Verdict:* Methodologically pristine for COMP0123, **but** conceptually narrow ("null-model critique of one metric on an already-studied London dataset") with weak general-audience wow and low publishability.
- *(Note on runner-up S4: "Does the River Thames split London's bike-share network?" — interesting tension between the 0.48 km link-level penalty and the East–West community split, but strictly dominated by Candidate A, which tests 86 historical borders with OSM road travel times and a 1789/1833 natural experiment, whereas the Thames split is already noted in Morton et al. 2021 and heavily confounded by bridge detours without OSM cycling routes.)*

---

## 3. Head-to-Head Scoring Matrix: S1–S3 vs. Candidate A (1833) & Candidate B (Rural Exodus)

Scores on a 1–10 scale calibrated strictly to **UCL COMP0123 grading criteria** ($5 \times 20\%$: Question & design; Data & method; Analysis & discussion; Writing & visualisation; Difficulty & new knowledge) plus publishability.

| Criterion | **A · 1833 Borders** (`propuesta_final.md`) | **B · Rural Exodus → Holidays** (`alternativas_desde_cero.md`) | **S1 · E-bikes Rewire London** (Bike #1) | **S2 · Hidden Rebalancing Mirror** (Bike #2) | **S3 · Spatial Rich-Club Illusion** (Bike #3) |
|---|---|---|---|---|---|
| **Clarity** (1-sentence yes/no, zero jargon) | 8.5 *(clear Q; proof needs 1789/1833 history)* | **9.5** *(instantly clear & human)* | **9.0** *(very clear, intuitive to UCL grader)* | 8.5 *(clear, needs 1 sentence on bike-ID inference)* | 7.5 *(requires explaining "rich club" & null models)* |
| **Uniqueness / Novelty** (lit gap) | **9.0** *(no 1789 vs 1833 new/old/erased mobility test)* | 8.5 *(China Chunyun is contemporary migrants; here 60-yr echo)* | 6.5 *(crowded London + e-bike lit, though topological gap exists)* | 6.0 *(Zhu et al. 2020 did bike-ID inference; conservation tautology risk)* | 5.5 *(rich-club + spatial nulls well trodden)* |
| **Wow / Memorability** | **9.0** *(new 1833 lines = −41%, erased 1789 lines vanished)* | **9.0** *(ghost towns refill ×1.48; r=0.85 in BCN suburbs)* | 7.0 *(if rewiring survives downsampling; drops to 5 if artifact)* | 7.5 *("invisible van network" is cool, if edges diverge)* | 6.0 *(methodological "gotcha", low general wow)* |
| **Network-Nativeness** (fit with Shi Zhou's COMP0123) | 7.5 *(core is PPML gravity; network methods support it; density 3.8% ideal)* | 8.5 *(3 coupled networks, bipartite modularity, spatial null; density ideal)* | **9.0** *(pure 2-layer network comparison; Zhou's rewiring/BOSAM/2nd-order fit naturally; must use backbone for <5%)* | **9.0** *(pure 2-network mirror test; naturally <5% sparse; ideal for BOSAM/motifs/rich-club)* | **9.5** *(tailor-made to Zhou & Mondragón 2004/2007, 2nd-order mixing, BOSAM)* |
| **Robustness (Today)** | **9.0** *(OSM travel time, 1789/1822, placebos, 12m, 2025 method, LORO)* | 7.0 *(pending income control — confounded in Madrid! — & multi-year check)* | 6.0 *(raw only; $4.1\times$ volume downsampling & pricing/hill confounds untested)* | 6.0 *(workshop/shuffle noise & edge-level topology untested)* | 6.5 *(needs spatial-null rich-club + rewiring bounds implemented)* |
| **Feasibility** (time to polished COMP0123 report) | **9.0** *(95% computed; main work is ruthless trimming & writing)* | 7.5 *(4–5 days of decisive robustness checks still needed)* | **8.5** *(pipeline done; 2–3 days for downsampling, Zhou nulls, extra months)* | **8.5** *(pipeline done; 2–3 days for noise filtering & edge-topology analysis)* | **9.0** *(fastest to finish: 1–2 days of null-model scripts on existing data)* |
| **Publishability** (*EPJ DS* / *JTG* / *RSOS* / *PNAS Nexus*) | **8.5** *(strong journal paper; top-tier if 2023+GAL/CAT added)* | 8.0 *(strong journal paper if income control holds)* | 6.5 *(solid *Applied Network Science* / *JTG* note if rewiring survives)* | 6.0 *(workshop / short note; needs ground-truth van validation)* | 5.0 *(coursework / methods note only)* |
| **COMP0123 Weighted Fit (5×20% proxy)** | **8.7 / 10** *(IF strictly trimmed to 1 thread)* | **8.3 / 10** *(9.0 if income check passes)* | **7.7 / 10** *(up to 8.2 if downsampling holds)* | **7.4 / 10** | **7.3 / 10** |

---

## 4. Honest Recommendation: Is a Smaller, Purer Bike Project a Safer 10?

### The Tempting Trap — and Why It's Partially True
It is completely understandable why `santander-cycles-network` looks tempting right now:
- Candidate A has grown sprawling (1833, 1789, 1822, CCAA half-hour rule, language laws, partidos judiciales, dioceses, purpose graduation bias, OSM travel times), and Prof. Zhou **explicitly penalises over-ambition**. Its core estimator is PPML gravity regression, not pure graph theory.
- Candidate B is emotionally sharper and cleaner, but still has an open **kill-switch** (income confounding in Barcelona/Madrid).
- Meanwhile, `santander-cycles-network` sits finished on GitHub, runs in 6 minutes, and speaks pure complex-networks language (backbones, rich-club, motifs, Leiden/Infomap, multiplex layers).

**However, a smaller, purer network project is ONLY a "safer 10" if it doesn't sacrifice the 20% "Difficulty & Achievement (new knowledge)" band — and Santander Cycles DOES sacrifice it.**

Look closely at the risk asymmetry:
- **Risk with Candidate A (1833 borders):** *Presentation / scope risk* (trying to tell five stories instead of one, or sounding like an econometrics paper rather than a Complex Networks paper). **This is a writing/framing risk — 100% controllable by editing!** The hard empirical work is already done and battle-tested.
- **Risk with Candidate B (Rural exodus):** *Empirical kill-switch risk* (if income controls collapse the neighbour-city elasticity). Checkable in 1–2 days.
- **Risk with Santander Cycles (S1/S2/S3):** *Structural ceiling + unverified empirical risk.*
  1. You cannot edit your way out of **"London bike-share is a classic, heavily-mined UCL dataset"** (especially with UCL CASA's O'Brien/Cheshire/Batty and Netzschleuder's `london_bikes` right there).
  2. You have to defend the **33.6% raw density** (via disparity backbone or sparse layer) right on page 1.
  3. You haven't yet tested whether **S1** survives $4.1\times$ volume downsampling — if the 8% backbone overlap and reciprocity collapse mostly evaporate once classic trips are downsampled to 166k, the headline collapses.
  4. **Self-recycling optics:** Because `aleetreny/santander-cycles-network` is already published as a 14-chapter portfolio site, reusing it for COMP0123 adds awkward boundary management (even though it's your own non-assessed personal work, an examiner visiting the required GitHub URL will see a pre-existing 14-chapter showcase unless you carve out a clean, distinct single-question study with genuinely new methodology).

### Final Actionable Ranking

1. **First Choice — Candidate A (1833 Borders), "Disciplined Network Edition" (Recommended for both a 10 and a paper):**
   - **Why:** Highest combined score on **new knowledge + robustness + natural `<5%` density (3.8%, $N=2{,}567$)**, and already computed.
   - **How to make it a safe 10 for Prof. Zhou:** Cut the sprawl ruthlessly. Make the report **15–18 pages, ONE single thread**: *"Do the 1833 provincial borders still partition Spain's daily mobility network?"* Lead with **COMP0123 core network properties** (degree/strength distributions, small world, neutral assortativity + strong **2nd-order mixing** $R_{\text{avg}}=+0.60$, rich-club vs degree-preserving rewiring, k-core, BOSAM), then frame the doubly-constrained travel-time model explicitly as a **spatial configuration null model** (network science, not econometrics) and show the **modularity partition race** (1833 vs 1789 vs 1822 vs CCAA vs Voronoi vs rotated placebos) + the killer 1789→1833 new/old/erased border test. Relegate language zones and purpose/routine to a crisp 1-page robustness/mechanism subsection (or appendix).
2. **Co-First Choice (Conditional) — Candidate B (Rural Exodus → Holiday Map):**
   - Run the **income-control test** immediately (1 day). If the neighbour-city elasticity and village refill gradient survive cleanly with income fixed effects, **B actually edges out A for COMP0123** because it is **more network-native** (3 coupled networks on the same nodes + bipartite modularity with BiCM & spatial nulls), **needs zero Spanish history lecture**, and has equal wow + publishability.
3. **Fallback Only (Plan C) — Santander Cycles S1 (*"Do e-bikes rewire London's bike-share network?"*):**
   - Switch to S1 **only if** (a) you decide you strictly want a 100% pure-graph, zero-econometrics project on a small UK dataset, **and** (b) a 1-hour pre-test (downsampling classic trips to 166,628 × 50 draws, testing reciprocity, Leiden NMI, backbone Jaccard, and Thames penalty at matched volume + distance) confirms genuine topological rewiring beyond range extension.
   - If you do pick S1, expand from 1 month (May 2026) to 3–6 months (seasonal stability), add OSM cycling elevation/bridges, and add Prof. Zhou's **rewiring structural constraints, 2nd-order mixing, and BOSAM** to lift it clearly above the existing portfolio site and past literature.
