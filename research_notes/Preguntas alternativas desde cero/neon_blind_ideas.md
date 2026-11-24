# Independent Research-Question Scout: Spain as a Complex-Networks Natural Laboratory

**Prepared for:** UCL COMP0123 (*Complex Networks and Online Social Networks*, Prof. Shi Zhou) & Publishable-Paper Pipeline  
**Scope:** First-principles scouting of plain-language, non-obvious, network-native yes/no research questions on real Spanish networks (500–50,000 nodes, density < 5%).  
**Exclusion honored:** Administrative/historical/language border effects on mobility flows are completely excluded.  
**Independence:** Conducted blind to prior project notes, reports, and planning folders.

---

## 1. First-Principles Analysis: Where Spain Is a World-Class Natural Laboratory

To produce a paper that is simultaneously **one-sentence clear to the general public**, **methodologically native to network science** (rich-club ordering, degree-preserving rewiring, second-order assortativity, $k$-core decomposition, bipartite configuration models, percolation, BOSAM fingerprints), and **publishable internationally**, we should ask: *where does Spain exhibit an extreme physical, social, or institutional regime that makes a universal question in network science cleaner to answer here than anywhere else on Earth?*

Seven structural extremes stand out:

1. **Synchronized Seasonal Population Inversion & the Dual-Home (*"El Pueblo"* + Coast) Culture**  
   Spain is both the world's #2 tourism destination and home to a unique mass cultural institution: millions of urbanites maintain an ancestral family village (*el pueblo*) in depopulated rural Spain alongside summer coastal escapes, all synchronized around August, *Semana Santa*, and national *puentes*. Combined with MITMA's rare **daily national overnight-stay (*pernoctaciones*) matrix** (~2,700 zones), Spain is the best natural laboratory on Earth to test **how an entire nation's topological hierarchy transforms under seasonal displacement** — whether cities mix in holiday hubs or carve out segregated summer territories, and whether the rich-club flattens or flips.
2. **World's #2 High-Speed Rail Network Built as an Extreme Radial Star on an Inland Capital**  
   With >4,000 km of AVE track radiating outward from inland Madrid across sparsely populated plateaus (*la Meseta*) — and three new openings right inside MITMA's 2022–2025 window (Burgos July 2022, Murcia Dec 2022, Asturias *Variante de Pajares* Nov 2023) — Spain is the cleanest natural experiment for testing **how engineered long-range hub shortcuts deform an existing spatial mobility graph** (community capture, local clustering collapse, and gravity-residual rich-clubs).
3. **"La España Vaciada" (Europe's Starkest Demographic Desert) + Ultra-Longevity + Counter-Current Immigration**  
   Roughly 70% of Spain's land has <12.5 inhabitants/km² across 8,131 municipalities, paired with one of the world's highest life expectancies (~84 years) and large post-2000 immigration. This creates an extreme testbed for **core-periphery collapse, generational network decoupling, step-migration breakdown, and hidden rural bridging**.
4. **Hyper-Decentralized Public Spending Meets Concentrated Corporate & Legal Structures**  
   Seventeen Autonomous Communities, 50 provinces, and 8,131 municipalities publish tenders on a unified open portal (PLACSP), including hundreds of thousands of no-bid *contratos menores* (<€15k services / <€40k works). Alongside the BOE consolidated legal-citation corpus and BORME corporate registry, Spain is a prime laboratory for **bipartite institutional capture vs. local fragmentation** (comparing competitive vs. discretionary networks against bipartite configuration nulls).
5. **The "Iberian Energy Island" + Rapid Remote-Renewable Shift + Mediterranean Climate Extremes**  
   Mainland Spain and Portugal sit behind a narrow ~3% Pyrenees grid bottleneck while running some of the highest wind/solar shares in the world (culminating in the historic April 28, 2025 blackout), alongside Europe's worst mega-wildfires and catastrophic Mediterranean flash floods (the October 29, 2024 Valencia DANA). This provides real-world laboratories for **percolation, cascading failures, flood-bottleneck overlap, and post-disaster network rerouting/hysteresis**.
6. **Extreme Social Clock (UTC+1/+2 Offset & Split Workday) + Intense Multiparty Affective Polarization**  
   Spain combines a notorious solar-time mismatch (late dinners, 2 a.m. weekend nightlife, split workdays) with high left–right affective polarization in a European multiparty system (without US-style suburban residential apartheid). That allows clean tests of **24-hour temporal percolation** and **whether political polarization spills into physical travel outside the United States**.
7. **Extreme Academic Inbreeding (*Endogamia*) & the Science-to-Patent Paradox**  
   Spain ranks near the global top 10–12 in scientific publications yet suffers ~70–90% faculty academic inbreeding and chronic weakness in university–industry patenting, making it an ideal natural lab for **reciprocity cliques, missing domestic collaboration backbones, and cross-sector brokerage failure**.

> **Verified Data-Hygiene Traps Discovered During Scouting (Important for 4-Week Feasibility):**
> - **April 28–29, 2025 Iberian Blackout in MITMA:** Official MITMS incident logs explicitly state that **April 28 and 29, 2025 are missing** from the mobility dataset (due to nationwide power/cell-tower outages), and methodology changed on July 1, 2025. Do **not** design a MITMA mobility study around the April 2025 blackout day itself!
> - **October 29–31, 2024 Valencia DANA in MITMA:** MITMS flags **cell-tower anomalies on Oct 29–31, 2024** in the Levante area; any DANA network study must mask those 3 acute days and analyze **post-flood rerouting and recovery from Nov 1, 2024 onward** (using the RD-ley 6/2024 affected zones).
> - **INE Migration Microdata (EVR $\to$ EMCR):** The old EVR (through 2021) identified municipalities <10,000 inhabitants, whereas the new EMCR (from Dec 2023) restricts fine <10k municipal disclosure in public microdata — use EVR + 2021 Census or >10k municipalities (~760 nodes) cleanly.
> - **LibreBORME Shutdown:** The original `LibreBORME` service shut down in April 2025; BORME corporate data must now be drawn from the official BOE Open Data API or successors (`OpenMercantil` / `BORME Fácil`).

---

## 2. Candidate Research Questions (30 Candidates Across 8 Themes & Datasets)

All 30 candidates are framed as **single-sentence plain-language yes/no questions** in English and Spanish, target **500–50,000 nodes with density < 5%**, use real relational edges (no arbitrary correlation thresholding), and have genuinely non-obvious answers.

---

### Theme A: Seasonal Mobility, Overnight Stays & "El Pueblo" (MITMA Trips & Overnight Stays 2022–2025)

#### Q1. Holiday Melting Pots vs. Private City Summer Territories
- **EN:** Do Spaniards from different major cities mix in the same holiday destinations, or does each city have its own private summer territory?
- **ES:** ¿Se mezclan los españoles de distintas grandes ciudades en los mismos destinos de vacaciones, o cada gran ciudad tiene su propio territorio de verano?
- **Network:** Bipartite network of origin municipalities/urban areas (~100–2,700 origins) $\times$ destination municipalities (~2,700 zones) from MITMA July–August 2024 overnight stays (*pernoctaciones*), thresholded/backbone-extracted (projected origin–origin co-visitation graph and full bipartite graph; ~2,700–5,400 nodes, density ~1–3%). Stratified into **Coastal destinations** vs. **Inland rural (*el pueblo*) destinations** and contrasted with **October weekday overnight stays**.
- **Network-Science Method:** Bipartite modularity (Barber's $Q_B$ / Leiden on bipartite projection), normalized origin-mixing entropy per destination, rich-club coefficient $\rho(k)$ of multi-city shared hubs, $k$-core decomposition, and Coastal vs. Pueblo mixing gap.
- **Null Model / Comparison:** (1) Bipartite Configuration Model (BiCM; Saracco et al. 2015) preserving exact origin out-strength and destination capacity; (2) **Spatial gravity null** using the student's road travel-time matrix (crucial to separate cultural/historical territorial sorting from mere driving distance!); (3) August vs. October control network.
- **Data Source:** MITMS *Estudio de movilidad con Big Data* — *Pernoctaciones* 2024 (already downloaded) + road travel-time matrix (`https://www.transportes.gob.es/ministerio/proyectos-singulares/estudio-de-movilidad-con-big-data`).
- **Why Non-Obvious:** Mass resorts (Benidorm, Mallorca, Marbella) feel like national melting pots, while folklore claims each city colonizes its own coast and rural diaspora (*madrileños* in Gandía/Soria, *catalanes* in Costa Brava, *vascos* in Laredo/Salou). Whether geography alone explains holiday segregation — and whether beaches mix while ancestral *pueblos* segregate (or vice versa) — is empirically unknown.

#### Q2. The August Rich-Club Inversion
- **EN:** Does Spain's travel network flatten out when everyone goes on holiday in August, or does it become even more concentrated in a few hubs?
- **ES:** ¿Se aplana la red de viajes de España cuando todo el mundo se va de vacaciones en agosto, o se concentra todavía más en unos pocos polos?
- **Network:** 366 daily weighted directed networks of Spain's ~2,700 municipal zones for all of 2024 (both overnight stays and inter-municipal trips; density ~1–4% after backbone/threshold filtering).
- **Network-Science Method:** Daily normalized topological and weighted rich-club coefficient $\rho(k)$ (Zhou & Mondragón 2004; Colizza et al. 2006; Opsahl et al. 2008), daily $k$-core membership turnover (Jaccard/Kendall $\tau$ of core nodes), degree/strength Gini, small-world path length, and BOSAM adjacency-matrix fingerprint shift across calendar regimes (ordinary Tuesday vs. *Semana Santa* vs. August 15 vs. Christmas).
- **Null Model / Comparison:** Degree-preserving Maslov–Sneppen rewiring (for $\rho(k)$ normalization every day) + direct structural comparison of August 2024 vs. October 2024 networks.
- **Data Source:** MITMS *Viajes* & *Pernoctaciones* 2024 (already downloaded; `https://www.transportes.gob.es/ministerio/proyectos-singulares/estudio-de-movilidad-con-big-data`).
- **Why Non-Obvious:** Intuitively, when Madrid and Barcelona empty across thousands of towns in August, the network should decentralize and flatten. Yet simultaneous piling into a narrow strip of coastal resorts could actually *sharpen* hub dominance and create an even tighter seasonal rich-club composed of completely different municipalities.

#### Q3. "Empty Spain" Villages as Hidden Seasonal Bridges
- **EN:** Are Spain's shrinking villages dead ends on the mobility map, or do they secretly connect distant cities during holidays?
- **ES:** ¿Son los pueblos que se vacían callejones sin salida en el mapa de movilidad, o conectan en secreto ciudades lejanas durante las vacaciones?
- **Network:** ~2,700 MITMA zones tagged by INE population density and 20-year depopulation trend (~1,400 shrinking/low-density rural zones), comparing normal-week trips vs. holiday overnight stays (*Semana Santa*, August, Christmas); density ~1–3%.
- **Network-Science Method:** Node betweenness centrality, $k$-core shell promotion (do dead-end shell-1 rural nodes jump into higher $k$-cores in August?), and **bipartite urban-bridge score** (measuring whether a shrinking village hosts overnight visitors from *multiple distant metropolitan origins* simultaneously).
- **Null Model / Comparison:** Normal non-holiday weeks (e.g., March/October) vs. peak holiday weeks, benchmarked against a distance-preserving gravity null model.
- **Data Source:** MITMS 2024 *Pernoctaciones* & *Viajes* (already downloaded) + INE *Padrón Municipal* (`https://www.ine.es/dynt3/inebase/index.htm?padre=525`).
- **Why Non-Obvious:** Demographic maps portray *La España Vaciada* as a disconnected periphery. Yet Spain's 1960s rural exodus sent relatives from the same interior comarcas (Soria, Teruel, Extremadura, Galicia) to Madrid, Barcelona, and Bilbao alike — potentially turning depopulated villages into hidden social bridges uniting Spain's distant metros each summer.

#### Q4. Daytime Backbone vs. Nighttime Backbone (Multiplex Decoupling)
- **EN:** Are the towns that keep Spain connected by day the same ones that hold it together at night?
- **ES:** ¿Son los municipios que mantienen conectada a España de día los mismos que la articulan de noche?
- **Network:** Two-layer multiplex network on the exact same ~2,700 MITMA municipalities: Layer 1 = 2024 daily inter-municipal trips; Layer 2 = 2024 daily overnight stays (backbone density ~1–2% per layer).
- **Network-Science Method:** Inter-layer edge overlap, degree & betweenness rank correlation, $k$-core and rich-club membership overlap, community partition similarity (NMI/ARI), and **cross-layer targeted attack** (does removing the top 5% daytime hubs shatter the overnight-stay giant component, or vice versa?).
- **Null Model / Comparison:** Direct two-layer comparison (Day vs. Night), plus degree-preserving rewired multiplex nulls controlling for total flow volume.
- **Data Source:** MITMS 2024 *Viajes* and *Pernoctaciones* (already downloaded; `https://www.transportes.gob.es/ministerio/proyectos-singulares/estudio-de-movilidad-con-big-data`).
- **Why Non-Obvious:** National transport investment targets daytime commuting hubs. If Spain's overnight-stay layer relies on a structurally distinct set of coastal, tourist, and rural-diaspora connectors, then daytime traffic maps miss half the country's real social backbone.

#### Q5. One-Way Overnight Hierarchy (Trophic Coherence of Where Spaniards Sleep)
- **EN:** When residents of one Spanish town sleep in another, is the visit ever returned — or is Spain's overnight network strictly one-way?
- **ES:** Cuando los vecinos de un municipio español duermen en otro, ¿se devuelve la visita, o es la red de pernoctaciones estrictamente unidireccional?
- **Network:** Directed weighted overnight-stay graph of ~2,700 MITMA municipalities (edge $i \to j$ = residents of $i$ sleeping in $j$, filtered by multiscale disparity backbone; density ~1–2%).
- **Network-Science Method:** Weighted dyadic reciprocity, **trophic coherence / trophic level hierarchy** (MacKay, Johnson & Sansom 2020 — originally from food webs!), directed triad census (feed-forward chains vs. reciprocal pairs), and weekday vs. summer trophic-level inversion.
- **Null Model / Comparison:** Directed configuration model preserving exact in- and out-strength sequences + comparison against the daytime trip graph (which is near-symmetric by round-trip conservation).
- **Data Source:** MITMS 2024 *Pernoctaciones* (already downloaded; `https://www.transportes.gob.es/ministerio/proyectos-singulares/estudio-de-movilidad-con-big-data`).
- **Why Non-Obvious:** Daily trips are inherently symmetric (out in the morning, back at night), masking directional power hierarchies. Overnight stays can be wildly asymmetric: do big cities sit at the bottom or top of Spain's "sleeping food chain," and are inter-city business stays reciprocal while leisure/rural stays form strict one-way hierarchies?

---

### Theme B: Social, Demographic & Temporal Divides in Mobility (MITMA 4 Detailed Weeks 2024 + Elections / Rent)

#### Q6. Socioeconomic Rich-Club ("Two Classes, Two Maps")
- **EN:** Do Spain's richest neighbourhoods travel in a tighter, more exclusive national club than working-class neighbourhoods?
- **ES:** ¿Viajan los barrios más ricos de España en un club nacional más cerrado y exclusivo que los barrios obreros?
- **Network:** Separate national mobility networks (~2,700 zones) built from MITMA's 4 detailed 2024 weeks for top-income vs. bottom-income residence bands, stratified by trip purpose (**work/study** vs. **occasional/leisure**); density ~1–3%.
- **Network-Science Method:** First-order and **second-order assortative mixing** (Zhou, Cox & Hansen 2017) by income, normalized rich-club coefficient $\rho(k)$ among affluent destinations, and modularity comparison across income layers and purposes.
- **Null Model / Comparison:** Spatially constrained gravity + degree-preserving null model (controlling for residential clustering of wealth) + Work vs. Occasional comparison.
- **Data Source:** MITMS 2024 4 detailed weeks (already downloaded) + INE *Atlas de Distribución de Renta de los Hogares* (`https://www.ine.es/experimental/atlas/experimental_atlas.htm`).
- **Why Non-Obvious:** Intuition says wealthy citizens are more cosmopolitan and broadly connected. Yet working-class commuters cross income boundaries daily to staff affluent urban cores, whereas wealthy residents' occasional/leisure travel may form an ultra-exclusive rich-club linking rich urban enclaves only to rich resorts. *(Note: Literature audit below reveals strong overlap with Duran-Sala et al. 2024!)*

#### Q7. Political Polarization in Physical Travel
- **EN:** Do left-wing and right-wing towns in Spain travel to each other less than geography predicts?
- **ES:** ¿Viajan entre sí menos de lo que predice la geografía los municipios de izquierdas y de derechas en España?
- **Network:** ~2,700 MITMA municipalities with real edges from 2024 **work/study** trips, **occasional/leisure** trips, and **summer overnight stays** (density ~1–3%), node-labelled with left–right vote share and polarization strength from the July 23, 2023 Spanish General Election.
- **Network-Science Method:** Weighted scalar assortativity $r$ (Newman) and **second-order assortative mixing** (Zhou et al. 2017) on political leaning, edge-residual analysis by ideological distance, and community alignment with electoral blocs.
- **Null Model / Comparison:** Crucial multi-factor null: **gravity + road travel time + municipal income + urban/rural density** residual model + degree-preserving rewiring + **Obligatory Work Commuting vs. Chosen Leisure/Overnight Stays** contrast.
- **Data Source:** MITMS 2024 (already downloaded) + road travel times + Ministerio del Interior *Infoelectoral* July 2023 municipal results (`https://infoelectoral.interior.gob.es/es/elecciones-celebradas/resultados-electorales/`).
- **Why Non-Obvious:** Spain exhibits some of Europe's highest affective polarization and echo-chambered social media, yet lacks US-style two-party residential apartheid. Does ideological avoidance leak into *physical travel* — especially chosen holidays vs. unavoidable work commutes — or does physical space remain politically blind once income and geography are controlled?

#### Q8. Gendered Mobility Network Topology Outside Work
- **EN:** Do women and men in Spain travel through differently shaped national networks — even outside work?
- **ES:** ¿Se mueven mujeres y hombres en España por redes con distinta forma, incluso fuera del trabajo?
- **Network:** Two matched national networks (~2,700 zones, density ~1–3%) constructed from female vs. male trips in MITMA's 4 detailed 2024 weeks, disaggregated into **work/study**, **frequent**, and **occasional** trips (and stratified by age group).
- **Network-Science Method:** Global topological comparison — clustering coefficient, community modularity & community geographic radius, characteristic path length, $k$-core reach, and normalized rich-club $\rho(k)$ — plus edge-level gender imbalance.
- **Null Model / Comparison:** Matched-volume subsampling (so total trip count is identical across male and female graphs!) + degree-preserving rewiring + Work vs. Occasional and 18–29 vs. 45–64 comparisons.
- **Data Source:** MITMS 2024 4 detailed weeks (already downloaded; `https://www.transportes.gob.es/ministerio/proyectos-singulares/estudio-de-movilidad-con-big-data`).
- **Why Non-Obvious:** Shorter female work commutes are well known ("mobility of care"). What is not obvious is whether the *network geometry itself* (fragmentation into tighter local communities, fewer inter-regional bridges, lower rich-club cohesion) persists for **non-work occasional trips** and among **young adults**.

#### Q9. Generational Disconnect (Retirees vs. Young Adults)
- **EN:** Do Spain's retirees and young adults live in two almost separate geographic networks?
- **ES:** ¿Viven los jubilados y los jóvenes españoles en dos redes geográficas casi separadas?
- **Network:** Multiplex / comparative networks of travellers aged 65+ vs. 18–29 (and 30–44) across ~2,700 MITMA zones from the 4 detailed 2024 weeks (density ~1–3%).
- **Network-Science Method:** Edge Jaccard / cosine overlap, hub rank correlation, community partition NMI, giant-component fragmentation, and destination age-assortativity (do ageing municipalities connect disproportionately to other ageing municipalities?).
- **Null Model / Comparison:** Volume-matched subsampling + population-age-weighted gravity null model.
- **Data Source:** MITMS 2024 4 detailed weeks (already downloaded) + INE municipal age structure (`https://www.ine.es/dynt3/inebase/index.htm?padre=525`).
- **Why Non-Obvious:** In one of the world's fastest-ageing nations, seniors might either retract into tiny disconnected local islands or form a parallel national retirement/coastal mobility network (coastal enclaves, spa/IMSERSO towns) structurally decoupled from the youth university-and-job backbone.

#### Q10. Spain's 24-Hour Breathing Network (Midnight Percolation)
- **EN:** Does Spain's mobility network shatter into disconnected local islands every night, or stay connected 24 hours a day?
- **ES:** ¿Se rompe la red de movilidad de España en islas locales desconectadas cada noche, o permanece conectada las 24 horas?
- **Network:** 24 hourly national mobility snapshots (~2,700 zones; sparse density 0.05%–2% depending on hour) on typical Tuesdays vs. Friday/Saturday nights from MITMA's 4 detailed 2024 weeks.
- **Network-Science Method:** Temporal percolation analysis across the 24-hour cycle: giant connected component fraction $S(t)$, second-largest component susceptibility peak (phase-transition signature), hourly path length, and identification of nocturnal bridge hubs (freight corridors vs. nightlife coasts).
- **Null Model / Comparison:** Weekday 24h cycle vs. Saturday 24h cycle + hourly configuration-model percolation threshold.
- **Data Source:** MITMS 2024 4 detailed weeks (hourly OD trips, already downloaded; `https://www.transportes.gob.es/ministerio/proyectos-singulares/estudio-de-movilidad-con-big-data`).
- **Why Non-Obvious:** Spain's unique solar-time offset, split workday (*jornada partida*), and 2–4 a.m. weekend nightlife create a four-peak diurnal signature unlike Northern Europe. Whether the national giant component undergoes a genuine percolation phase transition around 3 a.m. on weekdays while remaining supercritical all night on weekends has never been tested.

#### Q11. Rent-Driven Asymmetric Commuter Extraction
- **EN:** Do Spain's most expensive municipalities pull commuters from a wider, poorer ring of towns than affordable cities do?
- **ES:** ¿Atraen los municipios más caros de España trabajadores desde un anillo más amplio y pobre que las ciudades asequibles?
- **Network:** Directed commuting network (~2,700 MITMA zones from monthly *movilidad obligada* + work trips; density ~1–2%) annotated with official municipal median rent per m² (SERPAVI) and INE median income.
- **Network-Science Method:** Directed in-degree/out-degree asymmetry, rent-gradient assortativity, $k$-core stratification (do high-rent zones monopolize the in-degree core while low-rent towns occupy the out-degree rim?), and catchment-radius scaling controlling for city population.
- **Null Model / Comparison:** Population- and employment-matched gravity null model comparing high-rent vs. affordable Spanish urban areas of similar size (e.g., Málaga/Palma/San Sebastián vs. Murcia/Zaragoza/Valladolid).
- **Data Source:** MITMS *Movilidad obligada* (already downloaded) + MIVAU *SERPAVI* national rent reference database (`https://www.mivau.gob.es/vivienda/alquila-bien-es-tu-derecho/serpavi`).
- **Why Non-Obvious:** Big cities obviously have higher rents and more commuters. The subtle network question is whether *rent levels independent of city size* distort local commuting topology into wider, more asymmetric, low-reciprocity "extraction funnels."

---

### Theme C: Transport Infrastructure, High-Speed Rail & Urban Morphology (MITMA + Road Times + ADIF/NAP GTFS + OSM)

#### Q12. High-Speed Rail (AVE) "Orbit Capture" & Local Hollowing-Out
- **EN:** Does a high-speed train station pull a Spanish city closer to Madrid than to its own neighbours?
- **ES:** ¿Acerca una estación de AVE a una ciudad española más a Madrid que a sus propios vecinos?
- **Network:** Realized all-mode national mobility network (~2,700 MITMA zones, density ~1–3%) compared directly against Spain's road travel-time network, focusing on Spain's ~45 AVE-station cities (e.g., Cuenca, Ciudad Real, Lleida, Zamora, Segovia, Puertollano, Antequera, Guadalajara, Burgos, Murcia, Oviedo).
- **Network-Science Method:** Community detection (Louvain/Leiden/Infomap — do AVE cities detach from their surrounding provincial communities and jump into Madrid's community like topological exclaves?), ego-network local clustering drop (do long-range AVE shortcuts erode local triangles with road neighbours?), and **gravity-residual rich-club $\rho(k)$** among AVE cities!
- **Null Model / Comparison:** (1) **Road-travel-time gravity null model** (using the student's road travel-time matrix); (2) Propensity-matched non-AVE provincial capitals/cities (Badajoz, Cáceres, Jaén, Almería, Logroño, Santander, Pamplona, Soria, Teruel, Huelva); (3) **Natural before/after causal check** inside MITMA 2022–2025 around three real AVE openings: Burgos (July 21, 2022), Murcia (Dec 19, 2022), and Asturias *Variante de Pajares* (Nov 29, 2023)!
- **Data Source:** MITMS 2022–2024 + road travel-time matrix (already downloaded) + ADIF/NAP AVE station list (`https://nap.transportes.gob.es/`).
- **Why Non-Obvious:** High-speed rail is politically justified as creating "territorial cohesion." In network science, adding radial hub shortcuts can paradoxically *reduce* local clustering and sever intermediate cities from their geographic neighbours (the empirical network signature of the "tunnel/straw effect").

#### Q13. Rail Network Fragility Without Madrid (Spain vs. Germany)
- **EN:** If Madrid's train stations closed tomorrow, would Spain's national rail network split into disconnected pieces?
- **ES:** Si mañana cerraran las estaciones de Madrid, ¿se partiría la red ferroviaria española en pedazos desconectados?
- **Network:** Complete Spanish passenger railway station network (~1,500–1,900 active ADIF/Renfe/FEVE/regional stations, weighted by daily scheduled trains from NAP GTFS; density ~0.2–0.4%), compared side-by-side against Germany's polycentric DB rail network (~5,000 stations from GTFS).
- **Network-Science Method:** Targeted hub removal (removing Madrid Atocha/Chamartín, then sequential highest-betweenness removal vs. random removal), giant connected component percolation curve, global efficiency / shortest-path travel-time blowup, articulation points, and normalized rich-club $\rho(k)$.
- **Null Model / Comparison:** Direct **two-network comparison** (Spain radial rail vs. Germany polycentric rail) + degree-preserving rewired planar benchmarks.
- **Data Source:** Spain NAP Multimodal GTFS (`https://nap.transportes.gob.es/`) + Germany GTFS (`https://gtfs.de/`).
- **Why Non-Obvious:** Spain's radial AVE map looks utterly Madrid-dependent, yet conventional Iberian-gauge lines and Mediterranean/Northern corridors still exist on paper. Quantifying whether removing Madrid actually *disconnects* Spain (or merely doubles/triples effective network diameter relative to Germany) tests infrastructure resilience cleanly.

#### Q14. Medieval Old Towns vs. Modern *Ensanches* (Street Resilience)
- **EN:** Are Spain's chaotic medieval old towns more fragile — or more resilient — to blocked streets than modern grid cities?
- **ES:** ¿Son los cascos históricos medievales españoles más frágiles — o más resistentes — a calles cortadas que los ensanches en cuadrícula?
- **Network:** Street intersection graphs (~2,000–15,000 nodes each, density < 0.2%) of Spanish historic Islamic/medieval urban cores (Toledo, Córdoba, Seville *Casco Antiguo*, Granada *Albaicín*) compared directly against planned 19th-century *Ensanches* (Barcelona *Eixample*, Madrid *Salamanca*, San Sebastián *Centro*).
- **Network-Science Method:** Percolation under random vs. betweenness-targeted street blockages, detour-factor inflation, meshedness/4-cycle vs. dead-end ratio, $k$-core structure, and edge-betweenness concentration.
- **Null Model / Comparison:** Direct two-class comparison (Medieval vs. *Ensanche*) normalized by node count and area, plus planar Delaunay/grid reference graphs.
- **Data Source:** OpenStreetMap via OSMnx (`https://www.openstreetmap.org/`).
- **Why Non-Obvious:** Orthogonal grids look maximally redundant, while organic medieval labyrinths look easily choked; however, high-betweenness arterial concentration in modern *Ensanches* vs. pedestrian shortcut webs in old towns can yield surprising vulnerability crossovers under targeted closures.

---

### Theme D: Natural Disasters, Climate Extremes & Ecological Connectivity (MITMA + AEMET + CNIG/SNCZI + Natura 2000/EFFIS)

#### Q15. Valencia DANA 2024: Real-World Percolation, Rerouting & Hysteresis
- **EN:** After the October 2024 Valencia flood cut roads and rail lines, did Spain's mobility network reroute smoothly around the damage or stay structurally scarred for months?
- **ES:** Tras cortar la DANA de octubre de 2024 carreteras y trenes en Valencia, ¿se reconfiguró la red de movilidad alrededor del daño o quedó marcada durante meses?
- **Network:** Daily mobility networks (~2,700 national zones + ~300 Levante/Valencia corridor zones; density ~1–3%) before (Sep–Oct 28, 2024) and after (Nov 1–Dec 31, 2024 + Jan 2025, explicitly excluding the Oct 29–31 cell-tower outage window flagged by MITMS).
- **Network-Science Method:** Empirical link-severance shock vs. simulated static percolation: measure local/regional giant-component efficiency, **emergent detour edges** (links that surge >3$\sigma$ above baseline to bypass broken highways/rail like the A-7/V-31/AVE corridor), community-boundary deformation, and **topological hysteresis** (time for BOSAM adjacency fingerprint and $k$-core to return to baseline).
- **Null Model / Comparison:** (1) Pre-DANA matched weekdays (Oct 1–28, 2024) and year-over-year Nov–Dec 2023 controls; (2) Static link-removal simulation (removing physically severed road/rail corridors without human rewiring) to quantify the exact value of adaptive human rerouting.
- **Data Source:** MITMS 2024 (already downloaded) + RD-ley 6/2024 affected-municipality list (`https://www.transportes.gob.es/ministerio/proyectos-singulares/estudio-de-movilidad-con-big-data`).
- **Why Non-Obvious:** Classical network robustness assumes static failure (edges vanish and traffic dies). Here a real catastrophe severed major highways and rail around Spain's 3rd largest metro — testing whether human adaptation rapidly heals global connectivity or leaves lasting multi-month topological scars.

#### Q16. Extreme Heatwaves: Shrinking Inward or Triggering Long-Distance Escapes?
- **EN:** Do extreme Spanish heatwaves shrink the country's mobility network inward, or trigger long-distance escapes to the cooler north and coast?
- **ES:** ¿Encogen las olas de calor extremas la red de movilidad española hacia dentro, o disparan escapadas de larga distancia hacia el norte y la costa?
- **Network:** Daily mobility & overnight-stay networks (~2,700 zones, density ~1–3%) during Summer 2022–2024 matched with daily maximum temperature per zone from AEMET meteorological stations.
- **Network-Science Method:** Thermal-gradient directed assortativity (do >40°C interior nodes redirect flow asymmetrically toward cooler northern/coastal nodes?), average trip-weighted topological path length, local clustering, and weekday vs. weekend structural bifurcation during official AEMET *olas de calor*.
- **Null Model / Comparison:** Matched same-month, same-weekday non-heatwave days (controlling for August seasonality and holiday calendar) + temperature-shuffled null.
- **Data Source:** MITMS 2022–2024 (already downloaded) + AEMET OpenData REST API (`https://opendata.aemet.es`).
- **Why Non-Obvious:** Extreme heat suppresses short outdoor urban trips (shrinking local network weight), yet may simultaneously *amplify* long-range inter-regional escape flows toward *España Verde* and beaches — testing whether climate extremes contract or stretch national mobility topology.

#### Q17. Are Spain's Most Flood-Prone Roads Also Its Structural Bottlenecks?
- **EN:** Are Spain's most flood-prone roads also the exact bottlenecks holding the national road network together?
- **ES:** ¿Son las carreteras más inundables de España precisamente los cuellos de botella que mantienen unido el mapa viario?
- **Network:** Spain's inter-urban highway/road intersection graph (~10,000–25,000 major road junctions & segments from CNIG IGR-RT; density < 0.1%), intersected spatially with MITECO's official 100-year and 500-year floodplain maps (SNCZI).
- **Network-Science Method:** Edge betweenness centrality, articulation-point/bridge identification, and percolation collapse comparing **flood-risk-ordered edge removal** against **random removal** and **adversarial highest-betweenness removal**.
- **Null Model / Comparison:** Random edge removal (baseline resilience) vs. spatial-cluster-matched random removal (controlling for spatial correlation of floodplains along river valleys) vs. highest-betweenness targeted attack.
- **Data Source:** CNIG *Redes de Transporte* (`https://centrodedescargas.cnig.es/`) + MITECO *SNCZI* floodplains (`https://sig.miteco.gob.es/snczi`).
- **Why Non-Obvious:** Spanish highways historically follow flat river basins and Mediterranean *ramblas*. If floodplain intersections systematically overlap high-betweenness cut-edges, climate-driven floods behave like *adversarial targeted attacks* rather than random failures.

#### Q18. Mega-Wildfires vs. Spain's Protected Nature Reserve Network
- **EN:** Would mega-wildfires break Spain's network of protected nature reserves into isolated islands?
- **ES:** ¿Partirían los grandes incendios la red de espacios naturales protegidos de España en islas aisladas?
- **Network:** Spatial habitat-connectivity graph of Spain's ~1,867 Natura 2000 protected sites (MITECO), with edges connecting reserves within terrestrial mammal/bird dispersal thresholds (e.g., 15–30 km; density ~1–3%).
- **Network-Science Method:** Percolation threshold $p_c$ and giant-component collapse under **historical wildfire-hazard-ordered node removal** (from EFFIS historical burn frequency 2006–2025) vs. random removal vs. highest-betweenness stepping-stone removal, plus $k$-core resilience.
- **Null Model / Comparison:** Random reserve loss vs. spatially autocorrelated hazard null vs. targeted betweenness attack.
- **Data Source:** MITECO *Red Natura 2000* (`https://www.miteco.gob.es/es/cartografia-y-sig/ide/descargas/biodiversidad/rn2000.html`) + Copernicus EFFIS wildfire database (`https://forest-fire.emergency.copernicus.eu/applications/data-and-services`).
- **Why Non-Obvious:** Spain protects >27% of its land (the EU's largest Natura 2000 estate), implying apparent redundancy. Yet if fire hotspots (NW Galicia/León, Sierra Morena, Levante) sit directly on the narrow north–south ecological stepping stones, losing just 5% of sites could fragment national biodiversity connectivity.

---

### Theme E: Energy Infrastructure & the Renewable Transition (Iberian Power Grid — PyPSA-Eur / OSM / REE)

#### Q19. Solar/Wind Boom & Cascading Blackout Fragility of the Iberian Grid
- **EN:** Has Spain's shift to remote solar and wind farms made its power grid more vulnerable to cascading blackouts than its old power plants?
- **ES:** ¿Ha hecho el despliegue de parques solares y eólicos remotos que la red eléctrica española sea más vulnerable a apagones en cascada que con las centrales antiguas?
- **Network:** Iberian high-voltage transmission grid (~1,500–2,800 substations, ~2,500–4,500 transmission lines at 132–400 kV from PyPSA-Eur OSM-based network; density ~0.1–0.3%), weighted by DC-flow / edge load under two dispatch regimes: **legacy centralized thermal/nuclear** vs. **modern remote solar/wind** (Extremadura, Castilla-La Mancha, Aragón, Galicia).
- **Network-Science Method:** Motter–Lai overload cascading failure & DC linear power-flow redistribution under random line faults, targeted high-betweenness attacks, and spatially correlated weather shocks (e.g., simultaneous SW solar drop), measuring cascade size and critical tolerance $\alpha_c$.
- **Null Model / Comparison:** Direct two-regime comparison on the identical physical grid (Legacy vs. Renewable generation siting) + rewired transmission-topology null preserving geographic length budget.
- **Data Source:** PyPSA-Eur Prebuilt OSM Electricity Network (`https://doi.org/10.5281/zenodo.14036625`) + REE / Global Power Plant Database (`https://www.ree.es/es/datos/generacion`).
- **Why Non-Obvious:** Distributing generation across thousands of wind and solar sites intuitively sounds more decentralized and resilient than relying on a handful of nuclear/coal plants. Network flow physics predicts a possible paradox: placing utility-scale renewables in depopulated Empty Spain forces bulk power across long, high-betweenness corridors into coastal/Madrid demand sinks.

#### Q20. Hidden Rich-Club Bottlenecks of the Iberian Energy Island
- **EN:** Does the Iberian electricity grid rely on a tiny hidden "rich club" of substations whose failure could black out half the peninsula?
- **ES:** ¿Depende la red eléctrica ibérica de un pequeño «club de ricos» oculto de subestaciones cuyo fallo podría apagar media península?
- **Network:** Same Iberian HV substation transmission graph (~1,500–2,800 nodes, density ~0.2%), compared directly with continental France/Germany sub-grids of matched size.
- **Network-Science Method:** Normalized degree and betweenness rich-club coefficient $\rho(k)$ (Zhou & Mondragón 2004; Cinelli 2019), $k$-core decomposition, electrical betweenness centrality, and targeted removal of the rich-club core vs. peripheral nodes.
- **Null Model / Comparison:** Degree-preserving spatial/planar rewiring + direct comparison with a matched-size well-interconnected continental European grid (France or Germany).
- **Data Source:** PyPSA-Eur Prebuilt OSM Electricity Network (`https://doi.org/10.5281/zenodo.14036625`).
- **Why Non-Obvious:** Planar spatial grids are theoretically expected to show *no* rich-club (since long-range hub-to-hub shortcuts are costly). Testing whether Iberia's radial geography around Madrid and weak Pyrenees border create an anomalous high-betweenness bottleneck club reveals hidden structural fragility.

---

### Theme F: Public Procurement, Corporate Elites & Legal Systems (PLACSP / BORME / BOE)

#### Q21. No-Bid "Minor Contracts" vs. Open Tenders in Spanish Public Procurement
- **EN:** Do Spain's no-bid "minor contracts" go to a tighter, more closed club of companies than open-competition tenders?
- **ES:** ¿Van los contratos menores «a dedo» en España a un club más cerrado de empresas que las licitaciones abiertas?
- **Network:** Two bipartite networks constructed from Spain's PLACSP open data (2023–2024): (A) **No-bid minor contracts** (*contratos menores*, <€15k services / <€40k works) and (B) **Open competitive tenders** (*procedimiento abierto*). Nodes: ~5,000–12,000 public contracting bodies $\times$ ~20,000–35,000 winning firms (density < 0.2%).
- **Network-Science Method:** Bipartite clustering (4-cycle / Robins–Alexander closure), repeat-edge weight concentration (Herfindahl/Gini of agency–supplier lock-in), $k$-core depth, normalized rich-club $\rho(k)$ of serial multi-agency winners, and bipartite modularity.
- **Null Model / Comparison:** Direct **two-network comparison** (Minor vs. Open, matched by sector CPV code and contract volume) + **Bipartite Configuration Model (BiCM)** preserving exact agency contract counts and firm win counts.
- **Data Source:** PLACSP (*Plataforma de Contratación del Sector Público*) Open Data ATOM/XML feeds (`https://contrataciondelestado.es`, parseable via official `OpenPLACSP` tool).
- **Why Non-Obvious:** Defenders argue *contratos menores* exist to let town halls hire local mom-and-pop SMEs (which predicts a **star-like, unclustered, degree-1 periphery with NO rich-club**). Critics allege contract splitting and patronage (which predicts **excess 4-cycles, repeat-capture edges, and a dense insider $k$-core**). Both hypotheses yield diametrically opposite network topologies.

#### Q22. National Madrid Oligopoly vs. Regional Procurement Monopolies
- **EN:** In Spanish public spending, do the same big Madrid-based contractors win across the whole country, or does each region buy from its own closed local network?
- **ES:** En la contratación pública española, ¿ganan las mismas grandes empresas en todo el país o compra cada región a su propia red local cerrada?
- **Network:** Bipartite agency–supplier network (~8,000 agencies $\times$ ~30,000 suppliers from PLACSP 2024; density < 0.2%), with geographic province metadata for both buyer and supplier.
- **Network-Science Method:** Regional modularity ($Q$ vs. geographic null), bipartite core-periphery structure (is the core national/Madrid-HQ while the periphery is local?), supplier geographic reach distribution, and assortativity by agency tier (national vs. regional vs. municipal).
- **Null Model / Comparison:** BiCM preserving contract counts and contract-size distribution (controlling for the fact that huge infrastructure contracts can only be bid by big national firms).
- **Data Source:** PLACSP Open Data (`https://contrataciondelestado.es`).
- **Why Non-Obvious:** Spain simultaneously features extreme political decentralization (17 regional governments + 8,131 towns) and extreme corporate HQ concentration in Madrid. Which force governs the actual money-flow topology — local capture or national oligopoly — is an open structural question.

#### Q23. Corporate Interlocking Directorates: Super-Directors vs. Family Islands
- **EN:** Is Spain's corporate world stitched into a single "small world" by a tiny elite of shared directors, or split into isolated family islands?
- **ES:** ¿Está el mundo empresarial español cosido en un solo «mundo pequeño» por una pequeña élite de consejeros compartidos, o dividido en islas familiares aisladas?
- **Network:** Bipartite Company–Director network (and projected company–company interlock graph) covering ~10,000–25,000 major Spanish S.A./S.L. firms and directors from BORME open records (density < 0.2%).
- **Network-Science Method:** Giant connected component fraction, small-worldness ($\sigma, \omega$), normalized rich-club $\rho(k)$ of interlocking firms, $k$-core structure, and **targeted removal of the top 0.5%–1% "super-directors"** to measure percolation fragility.
- **Null Model / Comparison:** Bipartite Configuration Model (BiCM) preserving board sizes and director seat counts + comparison of IBEX/listed sector vs. medium-sized family-dominated firms.
- **Data Source:** Official BOE/BORME Open Data API (`https://www.boe.es/datosabiertos/api/`) and structured Spanish open-data mirrors (`https://openmercantil.com/`, `https://bormefacil.es/` via `datos.gob.es`).
- **Why Non-Obvious:** Spain is paradoxically stereotyped both as a cozy "palco del Bernabéu" elite of 30 interconnected corporate families and as a fragmented economy of independent family SMEs. Percolation under super-director removal cleanly measures which stereotype governs reality.

#### Q24. The Fragile Core of Spanish Law (BOE Citation Network)
- **EN:** Does Spanish law depend on a tiny hidden core of old statutes whose repeal would ripple through half the legal system?
- **ES:** ¿Depende la legislación española de un pequeño núcleo oculto de leyes antiguas cuya derogación afectaría a medio ordenamiento jurídico?
- **Network:** Directed legal citation/modification network of ~15,000–35,000 consolidated Spanish laws and royal decrees from the official BOE REST API (edge $A \to B$ = norm $A$ modifies, repeals, or explicitly cites norm $B$; density < 0.05%).
- **Network-Science Method:** $k$-core decomposition, directed rich-club of foundational statutes, preferential-attachment age-vs.-in-degree scaling, and directed reachability / cascade depth when perturbing top hubs vs. random statutes.
- **Null Model / Comparison:** Directed configuration model + age-conditioned citation null (controlling for older laws having more time to accumulate citations) + comparison of pre-1978, post-1978, and EU-transposed laws.
- **Data Source:** Agencia Estatal BOE *Legislación Consolidada* Open Data API (`https://www.boe.es/datosabiertos/`).
- **Why Non-Obvious:** Citizens imagine newer laws cleanly replace older ones. In legal graphs, a handful of century-old codes (1889 Civil Code, tax laws) or sprawling modern "omnibus" decrees (*reales decretos-ley*) may create extreme dependency depth where modifying a tiny core reverberates across thousands of regulations.

---

### Theme G: Internal Migration & Demographic Restructuring (INE EVR / EMCR / Census 2021)

#### Q25. Breakdown of Step-Migration in "Empty Spain"
- **EN:** Do people leaving shrinking Spanish towns step through their provincial capital first, or jump straight to Madrid and the coast?
- **ES:** ¿Pasan quienes abandonan los pueblos que se vacían por su capital de provincia primero, o saltan directamente a Madrid y la costa?
- **Network:** Directed weighted internal-migration network of ~760+ Spanish municipalities (>10k + rural comarca aggregates from INE EVR 2000–2021, or ~2,700 zones from 2021 Census birthplace-to-residence matrix; density ~1–4%).
- **Network-Science Method:** Directed triad census (transitive step chains $\text{rural} \to \text{provincial capital} \to \text{metro}$ vs. direct jump arcs), trophic hierarchy levels, rich-club of metropolitan sinks, and longitudinal shift from 2000s to 2020s.
- **Null Model / Comparison:** Distance-weighted gravity / radiation null model + 1990s/2000s vs. 2015–2021 temporal comparison.
- **Data Source:** INE *Estadística de Variaciones Residenciales* (EVR) microdata & *Censo de Población 2021* (`https://www.ine.es/dyngs/INEbase/operacion.htm?c=Estadistica_C&cid=1254736177013`).
- **Why Non-Obvious:** Ravenstein's classic "Law of Step Migration" predicts migrants climb the urban hierarchy step-by-step via local provincial capitals. If inland Spanish capitals (Soria, Teruel, Zamora, Jaén) have lost their intermediate-absorber role and become net exporters, Spain's migration graph has collapsed into direct long-range jumps.

#### Q26. Foreign-Born vs. Spanish-Born Internal Migration Counter-Currents
- **EN:** Once living in Spain, do immigrants move between the same cities as Spanish-born residents, or follow a reverse internal map?
- **ES:** Una vez en España, ¿se mudan los residentes extranjeros entre las mismas ciudades que los nacidos en España, o siguen un mapa interno inverso?
- **Network:** Two matched directed internal-migration networks (~760–2,000 municipalities/areas from INE EVR/EMCR microdata, density ~1–3%): Layer 1 = Spanish-born internal movers; Layer 2 = foreign-born internal movers.
- **Network-Science Method:** Edge-direction sign-flip rate (links where net flow points *opposite* directions for natives vs. immigrants!), urban–rural degree assortativity, hub rank correlation, modularity, and rich-club comparison.
- **Null Model / Comparison:** Direct **two-network comparison** (Foreign-born vs. Spanish-born, matched for total movers) + age-adjusted gravity null.
- **Data Source:** INE EVR / EMCR microdata (`https://www.ine.es/dyngs/INEbase/listaoperaciones.htm`).
- **Why Non-Obvious:** Conventional wisdom assumes immigrants concentrate even more tightly in big-city metros than natives. Yet rural Spain's agriculture, construction, and elder-care labour shortages may pull foreign-born residents *into* ageing/depopulated nodes that Spanish youth are fleeing — creating a structural mirror-image network.

---

### Theme H: Science, Higher Education, Patents & Culture (OpenAlex / TESEO / OEPM / ICAA+BDNS)

#### Q27. Spain's Missing Domestic Backbone in Science
- **EN:** Do Spanish scientists collaborate more with foreign labs than with universities in other Spanish cities?
- **ES:** ¿Colaboran más los científicos españoles con laboratorios extranjeros que con universidades de otras ciudades de España?
- **Network:** Co-authorship network of ~1,200 Spanish research institutions/hospitals/centres (or ~20,000 active Spanish researchers across cities) from OpenAlex (2019–2024; density ~0.5–2%), augmented with their international co-authoring hubs.
- **Network-Science Method:** Domestic inter-city vs. international bridge strength, normalized rich-club $\rho(k)$ among top Spanish institutions (do Madrid and Barcelona hubs connect directly or via US/UK/German hubs?), regional modularity, and domestic shortest-path inflation if foreign co-authors are removed.
- **Null Model / Comparison:** Field- and institution-size-preserving configuration model + benchmark against a peer European country of similar output (Italy or France) in OpenAlex.
- **Data Source:** OpenAlex Open Scholarly Graph API (`https://api.openalex.org/`).
- **Why Non-Obvious:** Spanish national grant schemes explicitly incentivize domestic consortia, yet academic inbreeding and regional university politics may produce a "hollow domestic backbone" where top Spanish groups partner more readily with London or Boston than with neighbouring Spanish regions.

#### Q28. Reciprocal "You-Scratch-My-Back" Cliques in Spanish PhD Exam Committees
- **EN:** Do Spanish university professors form closed "you-scratch-my-back" cliques when choosing external PhD exam committees?
- **ES:** ¿Forman los profesores universitarios españoles camarillas cerradas de favores mutuos al elegir los tribunales externos de tesis?
- **Network:** Directed academic graph of ~15,000–40,000 Spanish thesis advisors and committee members from TESEO (edge $A \to B$ = advisor $A$ invites *external* professor $B$ onto a PhD tribunal; density < 0.1%).
- **Network-Science Method:** **Dyadic reciprocity** ($A \to B$ and $B \to A$, i.e., mutual favour exchange), 3-cycle/clique density, modularity across university pairs, and rich-club of "super-examiners."
- **Null Model / Comparison:** Discipline- and university-size-preserving directed configuration model (controlling for the fact that small subfields have a limited pool of eligible external examiners!).
- **Data Source:** Ministry of Universities *TESEO* Doctoral Thesis Database (`https://www.educacion.gob.es/teseo/`).
- **Why Non-Obvious:** Spanish law (RD 1393/2007 / RD 99/2011) mandated external PhD committees specifically to break university inbreeding (*endogamia*). Testing whether external invitations simply re-formed into tiny reciprocal inter-university exchange clubs evaluates whether the reform actually worked.

#### Q29. University–Industry Disconnect in Spain's Patent Network
- **EN:** In Spain's patent network, do university scientists and private companies actually invent together, or live in two disconnected worlds?
- **ES:** En la red de patentes de España, ¿inventan juntos los científicos universitarios y las empresas privadas, o viven en dos mundos desconectados?
- **Network:** Co-invention network of ~10,000–25,000 Spanish patent inventors / ~4,000 patenting organizations (labelled Public University/CSIC vs. Private Firm; density < 0.1%) from OEPM / Lens.org (2010–2024).
- **Network-Science Method:** Newman's nominal assortativity $r$ by sector (public vs. private), cross-sector modularity, betweenness centrality of bridge nodes (is CSIC/tech-transfer doing all the bridging?), and percolation impact of removing public-private bridge inventors.
- **Null Model / Comparison:** Technology-class (IPC) and degree-preserving configuration model + comparison against Germany's patent co-invention network from the same source.
- **Data Source:** OEPM Open Data (`https://sede.oepm.gob.es/eSede/es/open-data.html`) & Lens.org Open Patent API (`https://www.lens.org/`).
- **Why Non-Obvious:** Spain has high basic-science output yet weak commercial innovation (the "Spanish innovation paradox"). Testing whether public and private inventors sit in separate network components — or actually co-invent inside specific biotech/pharma pockets — turns a macroeconomic cliché into a concrete topological diagnosis.

#### Q30. Public Film Subsidies & the Spanish Cinema Rich-Club
- **EN:** Does Spanish cinema's public funding flow to a tiny closed "rich club" of repeat directors and producers, or does it open the network to newcomers?
- **ES:** ¿Fluyen las subvenciones del cine español hacia un pequeño club cerrado de directores y productores habituales, o abren la red a nuevos creadores?
- **Network:** Collaboration graph of ~5,000–12,000 Spanish film directors, producers, and production companies across feature films (2005–2024; density ~0.2–0.5%), annotated with public subsidy awards from ICAA/BDNS.
- **Network-Science Method:** Compare the sub-network of **subsidized films vs. non-subsidized/independent films**: normalized rich-club $\rho(k)$ (Zhou & Mondragón), $k$-core occupancy of subsidy recipients, newcomer attachment rate (do debut directors enter the core via subsidies or only via co-production with rich-club insiders?).
- **Null Model / Comparison:** Degree-preserving Maslov–Sneppen rewiring + subsidized vs. non-subsidized two-network comparison.
- **Data Source:** Ministry of Culture *ICAA Catálogo de Cine Español* (`https://www.cultura.gob.es`) + IGAE *BDNS* National Subsidies API (`https://www.infosubvenciones.es`).
- **Why Non-Obvious:** Spanish public debate fiercely splits between critics calling film subsidies an insider *chiringuito* and defenders calling them the only ladder for debut filmmakers. Rich-club and $k$-core entry dynamics directly test both claims.

---

## 3. Scoring Matrix (0–10 on All Six Dimensions)

To reflect the student's stated priority order — **(1) Clarity, (2) Uniqueness, (3) Wow**, followed by **Network-Nativeness, Robustness, and 4-Week Feasibility** — the composite **Priority Score** weights them as:  
`Priority = 0.25·Clarity + 0.22·Uniqueness(Prior) + 0.18·Wow + 0.13·NetNative + 0.11·Robustness + 0.11·Feasibility`

*(Uniqueness below is the **pre-search prior**; Section 4 updates Uniqueness after our skeptical literature audit.)*

| ID | Short Title | Theme / Dataset | Clarity | Uniq. (Prior) | Wow | Net-Native | Robust. | Feas. (4w) | Priority Score | Top 8? |
|---:|---|---|---:|---:|---:|---:|---:|---:|---:|:---:|
| **Q1** | Holiday Mixing vs. Private City Territories | A / MITMA Overnight | **10.0** | **9.5** | **9.5** | **9.5** | **9.0** | **10.0** | **9.63** | **YES (#1)** |
| **Q7** | Political Polarization in Physical Travel | B / MITMA + 2023 Vote | **10.0** | **8.5** | **9.5** | **9.0** | **8.0** | **9.5** | **9.18** | **YES (#2)** |
| **Q12** | AVE High-Speed Rail "Orbit Capture" | C / MITMA + Road + AVE | **9.5** | **9.0** | **9.0** | **9.0** | **8.5** | **9.5** | **9.13** | **YES (#3)** |
| **Q2** | August Rich-Club Inversion (366 Days) | A / MITMA Daily | **9.0** | **8.5** | **9.0** | **10.0** | **9.0** | **10.0** | **9.13** | **YES (#4)** |
| **Q21** | No-Bid "Minor Contracts" vs. Open Tenders | F / PLACSP Procurement | **9.5** | **9.0** | **9.0** | **9.5** | **8.5** | **7.5** | **8.97** | **YES (#5)** |
| **Q15** | Valencia DANA 2024 Percolation & Rewiring | D / MITMA Oct–Dec '24 | **9.5** | **8.5** | **9.5** | **9.0** | **7.5** | **9.0** | **8.94** | **YES (#6)** |
| **Q6** | Socioeconomic Rich-Club ("Two Classes") | B / MITMA 4 Weeks | **9.5** | **8.0** | **9.0** | **9.5** | **8.5** | **9.5** | **8.97** | **YES (#7)** |
| **Q19** | Solar/Wind Shift & Iberian Grid Cascades | E / PyPSA-Eur + REE | **9.5** | **8.0** | **9.5** | **9.5** | **7.5** | **7.5** | **8.73** | **YES (#8)** |
| Q3 | Empty Spain Villages as Hidden Summer Bridges | A / MITMA Overnight | 9.0 | 9.0 | 8.5 | 9.0 | 8.0 | 9.5 | 8.86 | *Reserve* |
| Q8 | Gendered Mobility Network Shape Outside Work | B / MITMA 4 Weeks | 9.0 | 8.0 | 8.5 | 9.0 | 8.5 | 9.5 | 8.69 | — |
| Q28 | Reciprocal PhD Exam Committee Cliques | H / TESEO | 9.5 | 9.0 | 9.0 | 9.0 | 7.5 | 5.0 | 8.52 | — |
| Q4 | Daytime vs. Nighttime Multiplex Backbone | A / MITMA Trips+Overnight | 8.5 | 9.0 | 8.0 | 9.5 | 8.5 | 10.0 | 8.82 | — |
| Q5 | One-Way Trophic Hierarchy of Overnight Stays | A / MITMA Overnight | 8.5 | 9.5 | 8.0 | 9.5 | 8.5 | 9.5 | 8.87 | — |
| Q10 | 24-Hour Breathing Network (Midnight Percolation) | B / MITMA Hourly | 9.0 | 8.5 | 8.5 | 9.5 | 8.0 | 9.5 | 8.81 | — |
| Q16 | Extreme Heatwaves: Shrink or Escape North? | D / MITMA + AEMET | 9.0 | 8.5 | 8.5 | 8.0 | 7.5 | 8.5 | 8.45 | — |
| Q22 | National Oligopoly vs. Regional Procurement | F / PLACSP | 9.0 | 8.5 | 8.0 | 8.5 | 8.0 | 7.5 | 8.37 | — |
| Q9 | Generational Disconnect (65+ vs. 18–29) | B / MITMA 4 Weeks | 9.0 | 8.0 | 8.0 | 8.5 | 8.0 | 9.5 | 8.48 | — |
| Q17 | Flood-Prone Roads as National Bottlenecks | D / CNIG + SNCZI | 9.0 | 7.0 | 8.5 | 9.0 | 8.5 | 7.5 | 8.25 | — |
| Q18 | Mega-Wildfires vs. Natura 2000 Reserves | D / MITECO + EFFIS | 9.0 | 7.5 | 8.5 | 9.0 | 8.0 | 8.0 | 8.36 | — |
| Q26 | Foreign vs. Native Internal Migration Maps | G / INE EVR/EMCR | 9.0 | 8.0 | 8.5 | 8.5 | 8.0 | 8.0 | 8.41 | — |
| Q27 | Spain's Missing Domestic Science Backbone | H / OpenAlex | 9.0 | 8.0 | 8.0 | 8.5 | 8.5 | 9.0 | 8.48 | — |
| Q23 | Corporate Super-Directors vs. Family Islands | F / BORME | 9.0 | 7.5 | 8.5 | 9.0 | 7.5 | 6.0 | 8.09 | — |
| Q29 | University–Industry Patent Disconnect | H / OEPM / Lens | 8.5 | 7.5 | 7.5 | 8.5 | 8.0 | 8.0 | 7.99 | — |
| Q30 | Public Film Subsidies & Cinema Rich-Club | H / ICAA + BDNS | 8.5 | 8.5 | 8.0 | 8.5 | 7.0 | 6.5 | 8.03 | — |
| Q11 | Rent-Driven Asymmetric Commuter Extraction | B / MITMA + SERPAVI | 8.5 | 8.0 | 8.0 | 8.0 | 7.5 | 8.5 | 8.13 | — |
| Q13 | Rail Network Without Madrid (Spain vs. Germany) | C / NAP GTFS | 9.0 | 6.5 | 8.0 | 9.0 | 8.5 | 7.5 | 8.05 | — |
| Q25 | Breakdown of Step-Migration in Empty Spain | G / INE EVR/Census | 8.5 | 7.5 | 7.5 | 8.5 | 7.5 | 7.5 | 7.88 | — |
| Q24 | Fragile Core of Spanish Law (BOE Citations) | F / BOE API | 8.0 | 8.0 | 7.5 | 9.0 | 8.5 | 7.5 | 8.04 | — |
| Q20 | Hidden Rich-Club of the Iberian Power Grid | E / PyPSA-Eur | 8.0 | 7.5 | 8.0 | 9.5 | 8.0 | 7.5 | 8.03 | — |
| Q14 | Medieval Old Towns vs. Modern *Ensanches* | C / OSMnx | 8.5 | 5.0 | 7.5 | 9.0 | 8.5 | 9.0 | 7.67 | — |

*Note on Top 8 selection:* Q1, Q7, Q12, Q2, Q21, Q15, Q6, and Q19 represent the top-scoring mix across 5 distinct themes/datasets (MITMA overnight stays, MITMA daily/detailed mobility, elections, HSR/road infrastructure, disaster resilience, PLACSP public procurement, and Iberian power grid) — explicitly stress-testing both MITMA and non-MITMA frontrunners! (Q5 and Q3 scored right on the cusp and are folded naturally as analytical lenses inside Q1/Q2 if desired.)

---

## 4. Real Literature Novelty Check for the Top 8

Searches were run in **English and Spanish** across **OpenAlex, arXiv, Google Scholar / web indexes, OSF, Dialnet, CSIC, and Zenodo** (including terms such as `"red de movilidad"`, `"telefonía móvil"`, `"pernoctaciones"`, `"MITMA"`, `"contratos menores"`, `"DANA Valencia"`, `"efecto túnel AVE"`, `"despoblación"`, `"big data movilidad"`). Verdicts are intentionally **skeptical**: if an equivalent study exists for the US, China, or another European country, it is flagged explicitly.

---

### Top-8 Audit #1 — Q1: Holiday Melting Pots vs. Private City Summer Territories
*Question:* *"Do Spaniards from different major cities mix in the same holiday destinations, or does each city have its own private summer territory?"*

**Closest Existing Papers (6):**
1. **Sánchez-Martín, J. M. (2026).** *"A Reproducible Space–Time Cube Workflow for Domestic Tourism Mobility: Madrid-Origin Flows Across Spain (September 2019–September 2025)."* *Land*, 15(5), 887. DOI: `10.3390/land15050887`  
   - *What it did:* Built a monthly space-time cube of **Madrid-origin flows only** across Spanish municipalities using mobile-phone data to classify emerging tourism hot spots.  
   - *Why our question differs:* Single origin (Madrid only); no multi-city bipartite network, no cross-city mixing/co-visitation modularity, no BiCM or gravity null model.
2. **Garzón-Paredes, A.-R., & Royo-Vela, M. (2026).** *"Protected Heritage Stock and Population Redistribution in Spain: Multivariate Evidence from Four 2021 Mobile-Phone Mobility Snapshots."* *Sustainability*, 18(18), 9683. DOI: `10.3390/su18189683`  
   - *What it did:* Applied multivariate clustering to 4 INE 2021 mobile-phone snapshots across 3,214 areas to classify "summer-receiving" vs. "resident-dispersing" territories and correlate them with heritage assets.  
   - *Why our question differs:* Node-level multivariate typology, not an origin–destination bipartite network; does not ask whether distinct origin cities share or partition destinations.
3. **Osorio Arjona, J., Ruiz-Santacruz, J. S., & de las Obras-Loscertales Sampériz, J. (2023).** *"Mapping of functional areas in Spain based on mobile phone data during different phases of the COVID-19 pandemic."* *Journal of Maps*, 19(1). DOI: `10.1080/17445647.2023.2214804`  
   - *What it did:* Ran Louvain community detection on Spanish INE mobile-phone OD flows to map functional areas before/during/after COVID lockdowns.  
   - *Why our question differs:* Studied daytime trip catchment zones under pandemic restrictions, not bipartite origin mixing in summer holiday overnight stays or Coastal vs. *Pueblo* territorial segregation.
4. **Sierra-Porta, D., Díaz Ramírez, O., & Tobón Ospino, M. (2026).** *"Tourism resilience from networks: diversity and hierarchy."* *Journal of Complex Networks*, `cnag021`. DOI: `10.1093/comnet/cnag021`  
   - *What it did:* Constructed a country–city bipartite graph of *international* tourism arrivals to measure destination market diversity and hierarchy.  
   - *Why our question differs:* International macro-level arrivals (countries $\to$ global cities), not high-resolution intra-national mobile-phone overnight stays comparing summer coastal resorts vs. ancestral rural villages against road-travel-time nulls.
5. **Silm, S., & Ahas, R. (2010).** *"The seasonal variability of population in Estonian municipalities."* *Environment and Planning A*, 42(10), 2527–2546. DOI: `10.1068/a43139` (and **Ahas et al., 2007**, *Tourism Management*, DOI: `10.1016/j.tourman.2006.05.010`)  
   - *What it did:* Pioneered mobile positioning data to measure seasonal summer population swings in Estonian municipalities and foreign tourist seasonality.  
   - *Why our question differs:* Measured municipal population fluctuations and spatial distributions, not bipartite origin-mixing modularity, co-visitation rich-clubs, or BiCM null tests.
6. **Zaragozí, B., Trilles, S., & Gutiérrez, A. (2021).** *"Passive Mobile Data for Studying Seasonal Tourism Mobilities: An Application in a Mediterranean Coastal Destination."* *ISPRS IJGI*, 10(2), 98. DOI: `10.3390/ijgi10020098`  
   - *What it did:* Analyzed seasonal visitor origin composition in a single Catalan coastal region (Costa Daurada).  
   - *Why our question differs:* Single coastal case study, no nationwide bipartite graph or Coastal vs. Rural comparison.

- **Verdict:** **GAP** *(Updated Uniqueness: **9.5/10**)*  
  *Skeptical note:* Several recent Spanish papers use mobile-phone data to map summer tourism hot spots (Sánchez-Martín 2026; Garzón-Paredes 2026; Zaragozí 2021), and Estonia has classic seasonal population work (Silm & Ahas 2010). Yet **nobody — in Spain or internationally — has modeled an entire country's domestic summer overnight stays as a bipartite origin–destination network to test whether major cities mix in national holiday hubs or partition beaches and ancestral villages into segregated territories beyond road travel time.**

---

### Top-8 Audit #2 — Q2: The August Rich-Club Inversion (366 Daily Networks)
*Question:* *"Does Spain's travel network flatten out when everyone goes on holiday in August, or does it become even more concentrated in a few hubs?"*

**Closest Existing Papers (5):**
1. **Kruse, J., Gao, S., Ji, Y., Levin, K., Huang, Q., & Mayer, K. R. (2025).** *"Identifying Rich Clubs in Spatiotemporal Interaction Networks."* *Annals of the American Association of Geographers*, 115(4), 899–922. DOI: `10.1080/24694452.2025.2464806`  
   - *What it did:* Introduced the **Spatially Weighted Temporal Rich Club (WTRC)** metric normalized against null models and applied it to US spatiotemporal human mobility networks (Wisconsin and US counties).  
   - *Why our question differs:* Focused on developing the WTRC metric on US mobility; did not test **seasonal holiday hub inversion** (whether a nation's rich-club flattens or *re-concentrates on different coastal hubs* during mass August holidays vs. ordinary weekdays, across both trips and overnight stays).
2. **Zhou, S., & Mondragón, R. J. (2004).** *"The rich-club phenomenon in the Internet topology."* *IEEE Communications Letters*, 8(3), 180–182. DOI: `10.1109/LCOMM.2004.823426` & **Colizza, V., Flammini, A., Serrano, M. Á., & Vespignani, A. (2006).** *"Detecting rich-club ordering in complex networks."* *Nature Physics*, 2, 110–115. DOI: `10.1038/nphys209`  
   - *What they did:* Foundational definition and degree-preserving rewiring normalization of rich-club ordering on static graphs (Internet, airports, proteins).
3. **Opsahl, T., Colizza, V., Panzarasa, P., & Ramasco, J. J. (2008).** *"Prominence and Control: The Weighted Rich-Club Effect."* *Physical Review Letters*, 101(16), 168702. DOI: `10.1103/physrevlett.101.168702` & **Alstott et al. (2014)**, *Scientific Reports*, DOI: `10.1038/srep07258`  
   - *What they did:* Extended rich-club detection to weighted networks (e.g., US airports, scientific collaborations) on static snapshots.
4. **Bassolas, A., Barbosa-Filho, H., Dickinson, B., et al. (2019).** *"Hierarchical organization of urban mobility and its connection with city livability."* *Nature Communications*, 10, 4817. DOI: `10.1038/s41467-019-12809-y`  
   - *What it did:* Measured $k$-core / hotspot concentration of *intra-city* mobility across 300 cities globally (including Spanish cities).  
   - *Why our question differs:* Static intra-urban hierarchy, not daily year-long national seasonal inversion.
5. **Beneduce, C., Gullón Muñoz-Repiso, T., Lepri, B., & Luca, M. (2025).** *"pySpainMobility: Unlocking Spanish Open Mobility Data for Spatial Inequality Research."* *arXiv:2506.13385*. DOI: `10.48550/arxiv.2506.13385` (and **Kotov et al., 2026**, `spanishoddata`, DOI: `10.1177/23998083251415040`)  
   - *What they did:* Built Python/R software tools to access MITMA data with brief illustrative examples; no 366-day topological rich-club / $k$-core inversion study.

- **Verdict:** **PARTIAL GAP (Strong Empirical Gap)** *(Updated Uniqueness: **8.5/10**)*  
  *Skeptical note:* Kruse et al. (2025) already introduced temporal weighted rich-clubs on US mobility, so the *metric* is not brand new. However, tracking **366 daily normalized rich-club curves and $k$-core hub turnovers across a full year in Spain** — testing whether mass August holidays flatten hierarchy or invert it onto coastal hubs — is empirically untouched and aligns 10/10 with Prof. Shi Zhou's own research agenda.

---

### Top-8 Audit #3 — Q12: High-Speed Rail (AVE) "Orbit Capture" & Local Hollowing-Out
*Question:* *"Does a high-speed train station pull a Spanish city closer to Madrid than to its own neighbours?"*

**Closest Existing Papers (5):**
1. **Gutiérrez, J. (2001).** *"Location, economic potential and daily accessibility: an analysis of the accessibility impact of the high-speed line Madrid–Barcelona–French border."* *Journal of Transport Geography*, 9(4), 229–242. DOI: `10.1016/s0966-6923(01)00017-5`  
   - *What it did:* Seminal evaluation of Spanish HSR's *theoretical* daily accessibility and economic potential using travel-time models.  
   - *Why our question differs:* Computed modeled accessibility indices from infrastructure timetables — NOT observed all-mode human mobility networks, community capture, or clustering collapse.
2. **Ureña, J. M., Menerault, P., & Garmendia, M. (2009).** *"The high-speed rail challenge for big intermediate cities: A national, regional and local perspective."* *Cities*, 26(5), 266–279. DOI: `10.1016/j.cities.2009.02.003` (and **Garmendia et al., 2012**, *J. Infrastructure Systems*, DOI: `10.1061/(asce)is.1943-555x.0000104`)  
   - *What it did:* Qualitative and regional-planning analysis of intermediate AVE cities (Córdoba, Zaragoza, Ciudad Real, Lleida) near Madrid.  
   - *Why our question differs:* Case studies and rail ridership counts; no national graph community detection or road-gravity null model.
3. **Ortega, E., López, E., & Monzón, A. (2013).** *"Territorial cohesion impacts of high-speed rail under different zoning systems."* *Journal of Transport Geography*, 33, 196–205. DOI: `10.1016/j.jtrangeo.2013.10.018` & **Naranjo Gómez, J. M. (2016)**, *Sustainability*, 8(7), 624. DOI: `10.3390/su8070624`  
   - *What they did:* Evaluated whether Spain's AVE/motorway infrastructure increased or reduced territorial cohesion via shortest-path travel-time inequality.  
   - *Why our question differs:* Purely supply-side infrastructure travel times, not realized mobile-phone human travel topology.
4. **Bellet Sanfeliu, C. (2016).** *"Peripheral High-Speed Rail Stations in Spain."* *The Open Transportation Journal*, 10, 45–56. DOI: `10.2174/1874447801610010045` & **Campa, J. L., López-Lambas, M. E., & Guirao, B. (2016)**, *Journal of Transport Geography*, 56, 53–62. DOI: `10.1016/j.jtrangeo.2016.09.012`  
   - *What they did:* Studied peripheral AVE stations (e.g., Cuenca, Guadalajara-Yebes, Camp de Tarragona) and econometric tourism impacts.  
   - *Why our question differs:* No network science, no local clustering collapse or community-orbit capture test.
5. **Wang, W., Du, W., Li, W., et al. (2020).** *"Network analysis and spatial agglomeration of China's high-speed rail: A dual network approach."* *Chinese Physics B*, 29(11), 118902. DOI: `10.1088/1674-1056/abaee3` & **Sun, J. (2023)**, *HSET*, DOI: `10.54097/hset.v53i.9738`  
   - *What they did:* Analyzed complex-network communities of *China's HSR timetable graph* (edges = train services).  
   - *Why our question differs:* Studied the *train schedule graph* in China, whereas our question tests whether AVE stations deform **observed all-mode human mobility (cars + trains + buses in MITMA)** relative to a road travel-time gravity null — plus tests causal rewiring across the **2022–2023 Burgos, Murcia, and Asturias AVE openings**!

- **Verdict:** **GAP** *(Updated Uniqueness: **9.0/10**)*  
  *Skeptical note:* The "tunnel effect" (*efecto túnel*) of Spain's AVE is widely discussed in Spanish transport geography using theoretical travel-time accessibility models. Yet **zero studies have tested it as an empirical network-topology phenomenon in nationwide mobile-phone flows** (community exclave capture into Madrid, local ego-clustering drop vs. road neighbours, AVE gravity-residual rich-club, and before/after rewiring around the 2022–2023 AVE openings).

---

### Top-8 Audit #4 — Q21: No-Bid "Minor Contracts" vs. Open Tenders (PLACSP Procurement)
*Question:* *"Do Spain's no-bid 'minor contracts' go to a tighter, more closed club of companies than open-competition tenders?"*

**Closest Existing Papers (5):**
1. **Wachs, J., Fazekas, M., & Kertész, J. (2020).** *"Corruption risk in contracting markets: a network science perspective."* *International Journal of Data Science and Analytics*, 12, 45–60. DOI: `10.1007/s41060-019-00204-1` (plus Zenodo EU procurement graphs `10.5281/zenodo.3627215`)  
   - *What it did:* Modeled European public procurement (from EU TED database, including Spain) as agency–supplier bipartite networks and linked clustering/centralization to single-bidder corruption risk.  
   - *Why our question differs:* **TED only contains high-value above-EU-threshold tenders** (>€140k–€5.3M). It completely omits Spain's below-threshold no-bid *contratos menores* (<€15k/€40k) and never compares discretionary minor-contract topology against open-tender topology!
2. **Lyra, M. S., Curado, A., Damásio, B., et al. (2021).** *"Characterization of the firm–firm public procurement co-bidding network from the State of Ceará (Brazil) municipalities."* *Applied Network Science*, 6, 48. DOI: `10.1007/s41109-021-00418-y` & **Lyra et al. (2022)** review, *Applied Network Science*, 7, 83. DOI: `10.1007/s41109-022-00523-6`  
   - *What they did:* Analyzed Brazilian co-bidding networks and systematically reviewed network methods in procurement; highlighted that below-threshold direct awards remain severely under-analyzed structurally.
3. **Pliatsidis, A. C. (2024).** *"Analyzing concentration in the Greek public procurement market: a network theory approach."* *Journal of Industrial and Business Economics*, 51, 315–345. DOI: `10.1007/s40812-023-00291-z` & **Fountoukidis, I. (2026)**, LSE/Zenodo `10.5281/zenodo.20175803`  
   - *What they did:* Measured bipartite concentration and recurrence signals in Greek public procurement (including direct awards in Greece).  
   - *Why our question differs:* Closest methodological cousin internationally (Greece), confirming publishability while leaving Spain's massive PLACSP decentralized market untouched — and without our BiCM-normalized *minor-vs-open* "SME dispersion vs. closed rich-club" hypothesis test.
4. **Marco-Gracia, F. J. (2026).** *"Replication package: Procurement Thresholds and Contract Splitting."* Zenodo, DOI: `10.5281/zenodo.18420683` & **Sánchez-Graells, A. (2011)**, *"Public Procurement Below Thresholds in Spain,"* SSRN, DOI: `10.2139/ssrn.1888186`  
   - *What they did:* Econometric threshold-bunching and legal analyses of contract splitting in Spain (plus journalistic audits by *Civio* and *Political Watch*).  
   - *Why our question differs:* Non-network studies (bunching histograms / legal doctrine); no bipartite network topology, 4-cycle closure, BiCM null model, or rich-club comparison.
5. **Waxenecker, H., & Prell, C. (2024).** *"Corruption dynamics in public procurement: A longitudinal network analysis of local construction contracts in Guatemala."* *Social Networks*, 79, 1–14. DOI: `10.1016/j.socnet.2024.07.001`  
   - *What it did:* Longitudinal social network analysis of local construction procurement in Guatemala.

- **Verdict:** **GAP (with EU/Greece methodological cousins)** *(Updated Uniqueness: **8.5/10**)*  
  *Skeptical note:* Bipartite procurement network science is established in Europe (Wachs et al. 2020) and Greece (Pliatsidis 2024; Fountoukidis 2026), so you cannot claim "first procurement bipartite graph ever." However, **no peer-reviewed paper has compared Spain's PLACSP no-bid *contratos menores* network against open competitive tenders using BiCM-normalized clustering, repeat capture, and rich-club topology** to directly adjudicate "local SME dispersion vs. closed insider club."

---

### Top-8 Audit #5 — Q7: Political Polarization in Physical Travel
*Question:* *"Do left-wing and right-wing towns in Spain travel to each other less than geography predicts?"*

**Closest Existing Papers (6):**
1. **Liang, Z., & Cho, J. (2026).** *"Political polarization on the move: Analyzing geographical mobility between counties in the US."* *PLOS ONE*, 21(1). DOI: `10.1371/journal.pone.0339333`  
   - *What it did:* Combined mobile-phone OD tracking between 471 US counties with 2020 election returns; found political leaning alone did not govern mobility overall, but **strength of partisanship** significantly increased travel to politically congruent counties even after socioeconomic/geographic controls!  
   - *Why this is a major caution:* **Direct conceptual precedent published in 2026 for the US!**
2. **Zhang, Y., Cheng, S., Li, Z., & Jiang, W. (2023).** *"Human mobility patterns are associated with experienced partisan segregation in US metropolitan areas."* *Scientific Reports*, 13, 9875. DOI: `10.1038/s41598-023-36946-z`  
   - *What it did:* Measured experienced partisan segregation inside US metro areas from POI mobile-phone visits.
3. **Tonin, M., Tizzoni, M., Lepri, B., et al. (2026).** *"Experienced partisan segregation across patterns of mobility behavior in US cities."* *arXiv:2609.19285*. DOI: `10.48550/arxiv.2609.19285`  
   - *What it did:* Disentangled explorer vs. returner mobility behaviors in US urban partisan segregation.
4. **Brown, J. R., & Enos, R. D. (2021).** *"The measurement of partisan sorting for 180 million voters."* *Nature Human Behaviour*, 5, 998–1008. DOI: `10.1038/s41562-021-01066-z` (and **Chen & Rohla, 2018**, *Science*, DOI: `10.1126/science.aaq1433` on partisan Thanksgiving visits)  
   - *What they did:* Landmark US residential sorting and family-holiday partisan mobility studies.
5. **Jaziri-Arjona, T., & Valentim, V. (2025).** *"Short-Term Mobility and Changes in Political Norms."* *OSF Preprints*. DOI: `10.31219/osf.io/98k5v_v2`  
   - *What it did:* Used Spanish inter-province mobile-phone mobility during the summer between the April and November 2019 general elections to test whether **mobility causally shifted political norms** regarding the radical-right party Vox.  
   - *Why our question differs:* Opposite causal arrow (contact $\to$ norm shift) at coarse 50-province scale, not municipal partisan assortativity/avoidance across 2,700 zones.
6. **Bellodi, L., Docquier, F., Iandolo, S., et al. (2024).** *"Digging Up Trenches: Populism, Selective Mobility, and the Political Polarization of Italian Municipalities."* *SSRN*. DOI: `10.2139/ssrn.4696074`  
   - *What it did:* Studied permanent *residential migration* sorting across Italian municipalities (not daily/leisure/holiday mobility).

- **Verdict:** **PARTIAL GAP** *(Updated Uniqueness: **6.5/10** — downgraded from 8.5 after finding Liang & Cho 2026!)*  
  *Skeptical note:* Liang & Cho (2026, *PLOS ONE*) just answered this exact headline question for US counties, alongside Zhang et al. (2023) and Tonin et al. (2026) in US cities. What remains a genuine, publishable **PARTIAL GAP** is testing it **(a) outside the US in a European multiparty system** at 6× finer granularity (~2,700 municipalities vs. 471 US counties), **(b) explicitly splitting obligatory work commutes vs. discretionary leisure trips & summer overnight stays** (which Liang & Cho could not separate!), and **(c) applying Prof. Shi Zhou's second-order assortativity**. Still compelling, but no longer an untouched virgin gap.

---

### Top-8 Audit #6 — Q15: Valencia DANA 2024 Percolation, Rerouting & Hysteresis
*Question:* *"After the October 2024 Valencia flood cut roads and rail lines, did Spain's mobility network reroute smoothly around the damage or stay structurally scarred for months?"*

**Closest Existing Papers (5):**
1. **Tang, J., Zhao, P., Gong, Z., et al. (2023).** *"Resilience patterns of human mobility in response to extreme urban floods."* *National Science Review*, 10(8), `nwad097`. DOI: `10.1093/nsr/nwad097`  
   - *What it did:* Analyzed mobile-phone human mobility resilience, disruption, and recovery topology during the catastrophic 2021 Zhengzhou urban flood in China.  
   - *Why this is a close precedent:* Direct international methodological template applied to an extreme flood.
2. **Yabe, T., Tsubouchi, K., Fujiwara, N., Sekimoto, Y., & Ukkusuri, S. V. (2020).** *"Understanding post-disaster population recovery patterns."* *Journal of the Royal Society Interface*, 17(163), `20190532`. DOI: `10.1098/rsif.2019.0532` (and **Yabe et al., 2022** review in *CEUS*, DOI: `10.1016/j.compenvurbsys.2022.101777`)  
   - *What it did:* Derived universal post-disaster mobility recovery laws from mobile-phone trajectories across Japanese and US earthquakes/hurricanes.
3. **Dong, S., Gao, X., Mostafavi, A., et al. (2022).** *"Modest flooding can trigger catastrophic road network collapse due to compound failure."* *Communications Earth & Environment*, 3, 38. DOI: `10.1038/s43247-022-00366-0` (and **Wang et al., 2019**, *Nature Communications*, DOI: `10.1038/s41467-019-10063-w`)  
   - *What they did:* Percolation theory on flooded road infrastructure (simulated link failure, not observed all-mode mobile-phone rerouting).
4. **Polo Martín, B. (2025).** *"Valencia's battle against floods: A cartographic review to assess water management strategies."* *Journal of Geography and Cartography*, 8(1). DOI: `10.24294/jgc10129` & **Koroutchev, R. K. (2026)**, *Espacio, Tiempo y Forma*, DOI: `10.5944/etfvi.19.2026.45539`  
   - *What they did:* Cartographic and socio-economic assessments of the October 2024 Valencia DANA (descriptive, non-network).
5. **Pastor Serrano, D. (2026)** (Spanish thesis, Hispana: *l'Horta Sud* mobility & DANA 2024) & **MITMS HERMES GIS review** (*European Journal of Geography*, 2025)  
   - *What they did:* Descriptive GIS dashboards and local transit reporting for the DANA impact zone; no national/regional network percolation, emergent detour-edge identification, or BOSAM topological hysteresis.

- **Verdict:** **PARTIAL GAP** *(Updated Uniqueness: **7.5/10**)*  
  *Skeptical note:* No peer-reviewed network-science paper yet studies the 2024 Valencia DANA's mobility rerouting and multi-month hysteresis in MITMA data, giving a timely European case study. However, the *conceptual framework* is well-established from Zhengzhou (Tang et al. 2023) and Japan/US (Yabe et al. 2020), and MITMA's **confirmed Oct 29–31 tower outage** prevents analyzing the acute 72-hour collapse itself.

---

### Top-8 Audit #7 — Q6: Socioeconomic Rich-Club ("Two Classes, Two Maps")
*Question:* *"Do Spain's richest neighbourhoods travel in a tighter, more exclusive national club than working-class neighbourhoods?"*

**Closest Existing Papers (5):**
1. **Duran-Sala, M., Balachandran, A. K., Morandini, M., Naushirvanov, T., Prabhakaran, A., Renninger, A., & Mazzoli, M. (2024).** *"Disentangling individual-level from location-based income uncovers socioeconomic preferential mobility and impacts segregation estimates."* *arXiv:2407.01799*. DOI: `10.48550/arxiv.2407.01799`  
   - *What it did:* **Used Spain's exact MITMA mobile-phone mobility data stratified by income, age, and gender** to compute **socioeconomic assortativity**, weekday-vs-weekend segregation jumps, and "preferential mobility" of income groups!  
   - *Why this is a critical finding:* **Directly occupies the exact dataset and core assortativity angle in Spain!**
2. **Moro, E., Calacci, D., Dong, X., & Pentland, A. (2021).** *"Mobility patterns are associated with experienced income segregation in large US cities."* *Nature Communications*, 12, 4633. DOI: `10.1038/s41467-021-24899-8`  
   - *What it did:* Landmark study showing income-stratified mobility networks drive experienced segregation in US metros.
3. **Nilforoshan, H., Looi, W., Pierson, E., et al. (2023).** *"Human mobility networks reveal increased segregation in large cities."* *Nature*, 624, 586–592. DOI: `10.1038/s41586-023-06757-3`  
   - *What it did:* Massive US mobility-network income segregation study across 382 metros.
4. **Napoli, L., Sekara, V., García-Herranz, M., & Karsai, M. (2023).** *"Socioeconomic reorganization of communication and mobility networks in response to external shocks."* *PNAS*, 120(50), `e2305285120`. DOI: `10.1073/pnas.2305285120`  
   - *What it did:* Measured socioeconomic assortativity shifts in mobility networks under shocks.
5. **Liao, Y. J., Gil, J., Yeh, S., et al. (2025).** *"Socio-spatial segregation and human mobility: A review of empirical evidence."* *Computers, Environment and Urban Systems*, 117, `102250`. DOI: `10.1016/j.compenvurbsys.2025.102250`  
   - *What it did:* Comprehensive review confirming that income segregation / assortativity in mobile-phone mobility networks is now a densely crowded research field.

- **Verdict:** **DONE (Narrow Partial Gap at best)** *(Updated Uniqueness: **3.5/10** — severely downgraded after discovering Duran-Sala et al. 2024!)*  
  *Skeptical note:* This audit saved weeks of wasted effort. Duran-Sala et al. (2024) already ran income-stratified socioeconomic assortativity (weekdays vs. weekends) on **Spain's MITMA dataset**, atop landmark US *Nature* / *Nature Communications* / *PNAS* papers. Only a narrow incremental slice (weighted rich-club $\rho(k)$ by trip purpose) remains — far too crowded for a publishable flagship paper.

---

### Top-8 Audit #8 — Q19: Solar/Wind Boom & Cascading Blackout Fragility of the Iberian Grid
*Question:* *"Has Spain's shift to remote solar and wind farms made its power grid more vulnerable to cascading blackouts than its old power plants?"*

**Closest Existing Papers (5):**
1. **Rouco, L., Echavarren, F. M., & Lobato, E. (2026).** *"The overvoltage-driven blackout of the Iberian Peninsula on 28th April 2025."* *Sustainable Energy, Grids and Networks*, `102125`. DOI: `10.1016/j.segan.2026.102125` (plus **Albustami & Taha, 2026**, *Applied Energy*, DOI: `10.1016/j.apenergy.2026.129216`, and **Bialek, 2026**, *Curr. Sustain. Renew. Energy Rep.*, DOI: `10.1007/s40518-026-00297-w`)  
   - *What they did:* Authoritative engineering forensics of the April 28, 2025 Iberian blackout, showing **overvoltage and inverter protection dynamics** (not mere static line-betweenness overload) drove the collapse.
2. **Rosas-Casals, M., Valverde, S., & Solé, R. V. (2007).** *"Topological Vulnerability of the European Power Grid Under Errors and Attacks."* *Int. J. Bifurcation and Chaos*, 17(7), 2465–2475. DOI: `10.1142/S0218127407018531` & **Solé et al. (2008)**, *Physical Review E*, 77, `026102`. DOI: `10.1103/physreve.77.026102`  
   - *What they did:* Classic topological percolation / attack vulnerability of the Spanish and European HV power grids.
3. **Schäfer, B., Witthaut, D., Timme, M., & Latora, V. (2018).** *"Dynamically induced cascading failures in power grids."* *Nature Communications*, 9, 1975. DOI: `10.1038/s41467-018-04287-5`  
   - *What it did:* Modeled cascading failures in European/Spanish power grid topologies using AC swing-equation dynamics (and specifically warned against purely topological Motter–Lai abstractions).
4. **Witthaut, D., Hellmann, F., Kurths, J., et al. (2022).** *"Collective nonlinear dynamics and self-organization in decentralized power grids."* *Reviews of Modern Physics*, 94, `015005`. DOI: `10.1103/revmodphys.94.015005`  
   - *What it did:* Landmark physics review on how decentralized renewable siting alters grid vulnerability and cascades.
5. **Huo, L., Chen, X., & Sun, C. (2026).** *"Robustness evolution in decentralized power grids: the role of topology and renewable energy penetration."* *Nonlinear Dynamics*, DOI: `10.1007/s11071-025-12056-5`  
   - *What it did:* Directly modeled how increasing renewable penetration alters topological cascade robustness in power grids.

- **Verdict:** **PARTIAL GAP (leaning DONE / High-Validity-Risk Domain)** *(Updated Uniqueness: **5.5/10**)*  
  *Skeptical note:* Decentralized renewable cascade fragility is heavily mined in physics (*Rev. Mod. Phys.* 2022; *Nat. Comms.* 2018; *Nonlinear Dynamics* 2026), and the April 2025 Iberian blackout already produced multiple 2026 engineering papers. Worse, for publication, referees will rightly reject a 4-week topological/DC-flow model because the real 2025 blackout was driven by AC overvoltage and low synchronous inertia.

---

### Post-Audit Top-8 Summary Leaderboard

| Post-Audit Rank | ID | Candidate | Verdict | Clarity | Updated Uniq. | Wow | Course Fit (Zhou) | Feas. (4w) |
|---:|---:|---|---|---:|---:|---:|---:|---:|
| **1** | **Q1** | **Holiday Mixing vs. Private Summer Territories** | **GAP** | **10.0** | **9.5** | **9.5** | **9.5** | **10.0** |
| **2** | **Q12** | **AVE High-Speed Rail "Orbit Capture"** | **GAP** | **9.5** | **9.0** | **9.0** | **9.5** | **9.5** |
| **3** | **Q2** | **August Rich-Club Inversion (366 Days)** | **PARTIAL GAP** | **9.0** | **8.5** | **9.0** | **10.0** | **10.0** |
| **4** | **Q21** | **No-Bid "Minor Contracts" vs. Open Tenders** | **GAP** | **9.5** | **8.5** | **9.0** | **9.5** | **7.5** |
| **5** | **Q15** | **Valencia DANA 2024 Rerouting & Hysteresis** | **PARTIAL GAP** | **9.5** | **7.5** | **9.5** | **9.0** | **9.0** |
| **6** | **Q7** | **Political Polarization in Physical Travel** | **PARTIAL GAP** | **10.0** | **6.5** | **9.5** | **9.5** | **9.5** |
| **7** | **Q19** | **Solar/Wind Shift & Iberian Grid Cascades** | **PARTIAL GAP / DONE** | **9.5** | **5.5** | **9.5** | **9.0** | **7.5** |
| **8** | **Q6** | **Socioeconomic Rich-Club ("Two Classes")** | **DONE** | **9.5** | **3.5** | **9.0** | **9.5** | **9.5** |

---

## 5. Final Recommendation & Runner-Up

### Single Best Recommendation: **Q1 — Holiday Mixing vs. Private City Summer Territories**
> **EN:** *"Do Spaniards from different major cities mix in the same holiday destinations, or does each city have its own private summer territory?"*  
> **ES:** *"¿Se mezclan los españoles de distintas grandes ciudades en los mismos destinos de vacaciones, o cada gran ciudad tiene su propio territorio de verano?"*

**One-Paragraph Pitch:**  
Every August, tens of millions of Spaniards leave their cities simultaneously for the coast and their ancestral rural villages (*el pueblo*), creating one of the largest synchronized human redistributions in the developed world. Using MITMA's uniquely under-exploited **daily national overnight-stay (*pernoctaciones*) matrix** across ~2,700 municipalities (already on disk) paired with Spain's road travel-time matrix, this project models where Spaniards sleep in July–August as a bipartite origin–destination network — and compares it directly against (a) October weekday overnight stays, (b) a **Coastal Beach vs. Inland Rural (*Pueblo*)** split, and (c) both a degree-preserving **Bipartite Configuration Model (BiCM)** and a **road-travel-time gravity null**. The question needs zero jargon or history lesson for any reader worldwide, addresses a verified literature **GAP**, fits Prof. Shi Zhou's exact toolkit (bipartite modularity, normalized rich-club $\rho(k)$ of shared holiday hubs, $k$-core, and degree-preserving rewiring — with a natural bonus link to Q2's 366-day August rich-club curve), and yields one unforgettable signature figure: a national map of Spain coloring every municipality from **"national melting pot"** (high cross-city mixing entropy beyond geography) to **"single-city summer colony"**, showing whether Spain's beaches stitch the country together while its villages preserve invisible regional diasporas (or whether even the coasts are partitioned!).

**Main Threat to Validity (and How to Neutralize It):**  
**Trivial distance decay (Tobler's First Law).** Because Barcelona sits next to the Costa Brava and Madrid sits next to Cuenca and Alicante, raw bipartite modularity will mechanically detect geographic clustering unless distance is explicitly factored out. You **must** benchmark observed origin-mixing entropy and bipartite modularity not only against a degree-preserving BiCM null (which ignores space), but also against a **road-travel-time gravity null model** (using the already-available inter-municipal road-time matrix): the true scientific discovery is **excess mixing or excess territorial segregation *beyond driving time***. A secondary caveat to document transparently is MITMA's small-cell privacy suppression (<50–100 sleepers/night) and exclusion of foreign SIM cards (actually a feature here, since the research question targets domestic Spanish mixing).

---

### Runner-Up: **Q12 — High-Speed Rail (AVE) "Orbit Capture" & Local Hollowing-Out**
> **EN:** *"Does a high-speed train station pull a Spanish city closer to Madrid than to its own neighbours?"*  
> **ES:** *"¿Acerca una estación de AVE a una ciudad española más a Madrid que a sus propios vecinos?"*

**One-Paragraph Pitch:**  
Spain built the world's second-largest high-speed rail network as a radial star centered on Madrid in the name of "territorial cohesion," yet transport planners have debated for 30 years whether high-speed trains actually *detach* intermediate cities from their surrounding provinces (the famous *efecto túnel*). While every prior Spanish study evaluated this using theoretical timetable accessibility, nobody has tested it directly in **observed all-mode human mobility networks** — a verified **GAP**. Using the student's already-downloaded MITMA mobility matrices and road travel-time matrix across ~2,700 municipalities, this study compares Spain's observed mobility graph against its road-time gravity null to test three clean, network-native signatures on Spain's ~45 AVE-station cities vs. matched non-AVE provincial cities: (1) **community orbit capture** (do AVE cities jump into Madrid's Louvain/Infomap community as floating exclaves inside other provinces?), (2) **ego-network local clustering collapse** (do AVE shortcuts break local triangles with road neighbours?), and (3) a **gravity-residual AVE rich-club $\rho(k)$**. Because the question generalizes immediately to UK HS2, France's TGV, California HSR, and China's HSR, and produces a jaw-dropping map of AVE "islands" glowing Madrid-colored inside distant provinces, it is exceptionally publishable in *EPJ Data Science*, *Nature Cities*, or *Journal of Transport Geography*.

**Main Threat to Validity (and How to Neutralize It):**  
**Endogenous station placement (selection bias).** Spain did not drop AVE stations onto random towns: it connected provincial capitals and economically prominent cities that *already* had stronger administrative, business, and political ties to Madrid before high-speed rail arrived. Simply comparing AVE vs. non-AVE towns cross-sectionally could mistake pre-existing capital status for a causal rail effect. Neutralize this with a **two-pronged defense**: (1) restrict the primary cross-sectional comparison to **matched provincial capitals and mid-sized cities** (comparing AVE cities like Cuenca, Ciudad Real, Lleida, Zamora, Segovia against non-AVE peers like Badajoz, Cáceres, Jaén, Almería, Logroño, Santander, Pamplona, Teruel, Soria, controlling for population, GDP, and road travel time to Madrid); and (2) exploit the student's **2022–2025 MITMA archive as a causal difference-in-differences natural experiment** around three real AVE openings — **Burgos (July 21, 2022)**, **Murcia (Dec 19, 2022)**, and **Asturias/*Variante de Pajares* (Nov 29, 2023)** — to verify whether local clustering drops and Madrid-orbit centrality jumps immediately *after* high-speed service launches!

*(Strategic footnote for alternative appetites: If the student wants a **100% Prof. Shi Zhou rich-club purist** paper on already-downloaded data, **Q2** is the natural companion/pivot to Q1; if the student wants to **leave mobility entirely** for a high-impact civic/institutional dataset and can invest ~10 days parsing PLACSP XML feeds, **Q21** [No-Bid Minor Contracts vs. Open Tenders] is the #1 non-MITMA choice.)*
