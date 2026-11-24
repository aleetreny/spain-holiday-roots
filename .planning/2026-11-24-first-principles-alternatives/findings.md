# Findings — first-principles alternatives (2026-11-24)

## Diagnosis of candidate A (why the user feels it is "rebuscada")
- The QUESTION is clear; the EVIDENCE chain is not (1789 vs 1833 decomposition, spatial configuration null, travel-time bins).
- Its core estimator is a border-effect gravity regression; the network part (modularity race, nulls) is supporting.
  For COMP0123 a question whose answer IS a network property (hubs, clubs, communities, bridges, robustness, mixing) scores higher on "why network science".
- It needs a history lesson (what was 1833? 1789?) before the reader cares.

## Spain as a natural laboratory (where Spain is extreme)
1. Rural exodus 1950-75 + "España vaciada" (interior densities below Lapland) and the cultural institution of "el pueblo" (summer/holiday return).
2. Highest share of second homes in the EU (many inherited village houses).
3. Very fragmented municipalities (8,131; thousands < 100 inhabitants).
4. Radial infrastructure (Km 0; radial motorways; 2nd-largest high-speed rail network in the world, low use).
5. Strong regional identities / languages (tested in A: weak in mobility).
6. Mass tourism, coast-interior split.
7. Extreme summer heat; DANA floods; 2025 blackout (missing in MITMA).

## Candidate pool v0 (before scouting)
MITMA-based:
- M1 Homecoming: does holiday travel (August, All Saints 1-Nov, Christmas, Easter) follow the routes of the 1950-75 rural exodus? Two networks: exodus/birthplace network vs holiday-excess network; non-contiguous "diaspora communities".
- M2 Part-time villages: temporal network where depopulated villages switch on as nodes on weekends/holidays.
- M3 Depopulation & connectivity: are the dying villages the disconnected ones?
- M4 Assortative mixing by demographic trajectory: are dying villages connected to dying villages?
- M5 High-speed rail club / tunnel effect: do AVE cities form a rich club that size-matched non-AVE cities cannot join?
- M6 Pajares tunnel (Nov 2023) natural experiment: did Asturias move into Madrid's orbit?
- M7 Achilles' heel of the road network weighted by real demand (DANA validation).
- M8 Fault lines: hierarchical percolation of flows — where does Spain crack first?
- M9 Weak ties hold Spain together (Onnela-style edge removal).
- M10 Worst place to seed an epidemic; M11 effective distance.
- M12 Weekday vs weekend: two network architectures (assortative commuting r=+0.29 vs national hub network).
- M13 Second-order mixing (Zhou) as a signature of spatial embedding.
- M14 Is Catalonia an island? (P42: its border is among the thinnest).
- M15 Rich and poor municipalities mix? M16 Towns that vote differently mix?
- M17 Camino de Santiago fingerprint. M18 Heatwaves push Spain north/up?
- M19 "Disconnected Spain votes for Empty-Spain parties" (Teruel Existe, Soria Ya) — network peripherality vs vote.
Other datasets:
- O1 Internal migration network (INE EVR 1988-2024): is Madrid becoming a black hole?
- O2 Multiplex: do people move house along their commuting links?
- O3 Meta Social Connectedness Index (provinces): friendships follow the exodus?
- O4 Corporate elite rich club (BORME / IBEX boards; Zhou's "oligarchies as rich clubs").
- O5 Public procurement bipartite network (PLACSP).
- O6 Rail infrastructure (GTFS) star vs demand mesh.
- O7 Power grid & 2025 blackout.
- O8 University student flows. O9 Flights. O10 Surnames (isonymy) vs mobility.

## Pilots (MITMA pernoctaciones = residents of district X sleeping in municipal zone Y, daily; 180 days of 2024 downloaded
## to _local/external/mitma_pernoct; parquet data/derived/pernoct_mun_2024_6m.parquet; scripts explore/p60-p64)
Q1. Residence is coded at MITMA DISTRICT level (3,792 codes) -> mapped to municipal zones via relacion_ine_zonificacionMitma (0 unmapped).
Q2. Share of residents sleeping outside their home municipal zone: Feb Tue-Thu 4.7% (2.4% in another province);
    median Sat 8.8%/4.9%, Sun 10.8%/6.0%; mid-August peak 21.8% (10.5M people) with 14.8% (~7.1M) in ANOTHER PROVINCE (Sun 18-Aug);
    Easter Sat 17.3%/12.0%.
Q3. Province-level holiday excess (vs Feb Tue-Thu) of big-city residents: Madrid Aug -> Alicante 11%, Toledo 6, Ávila 6, Cádiz 5...;
    Nov 1 (All Saints) Madrid -> Toledo 13, Ávila 8, Segovia 5, C.Real 5, Cáceres 5 (interior/roots); Bizkaia -> Cantabria, Burgos, La Rioja, Palencia, Zamora, León;
    Barcelona -> mostly Catalan provinces; Andalusian provinces only ~2% each in August (distance dominates raw shares).
Q4. Village refill (present/residents): median zone x1.11 mid-Aug (p90 x1.8, p99 x3.0); Nov 1 median x1.04 (p99 x1.5).
    Interior aggregated small-municipality zones (_AM in Burgos, Segovia, Zamora, Huesca, La Rioja, Lleida) x3.5-5 in August AND x1.6-2.1 on All Saints,
    while beach towns (Oropesa, Peñíscola, Alcalà de Xivert) x4-5 in August but only x1.0-1.2 on All Saints
    -> All Saints' Day isolates the "roots" signal from tourism. Share of population sleeping in zones <5k residents: 6.9% Feb -> 9.6% Aug -> 7.5% Nov 1.
Q5. Holiday-excess networks (vs ordinary Fridays): Aug 15 6.6M excess sleepers, 74% other province, mean 245 km; All Saints 2.16M, 60% other prov, 182 km;
    size correlation residence-destination slightly negative on All Saints (-0.08). Leiden communities of excess networks: Q 0.60 (All Saints) / 0.47 (Aug),
    NMI(prov) 0.69/0.63, many communities spatially non-contiguous (pieces 21,14,10... — needs proper contiguity test).
Lit (2026-11-24): Exa 8-candidate check all PARTIAL (research_notes/.../exa_novelty_8_candidates.md). Homecoming analogues: China chunyun (many),
    Wesolowski 2017 Nat Commun (holiday seasonal asymmetry Kenya/Namibia/Pakistan), INE EM-1 2019 seasonal overnight population,
    Osorio Arjona 2024 Sci Rep (Andalusia weekends). OpenAlex ES+EN query phone x (second homes|rural exodus|seasonal population|despoblación): 6 hits, none relevant.
Q6. ROOTS TEST (explore/p65; INE Padrón 2021 province-of-birth x province-of-residence, 0ccaa001; mainland 47 provinces, PPML rp+sp FE + 20 distance bins):
    elasticity of overnight flows to "born in destination, living in residence" stock: ordinary weekday 0.20 (se .04), ordinary Sat 0.12, Good Fri 0.13,
    mid-Aug 0.15, All Saints 0.17, Christmas Eve 0.28 (excess-over-ordinary: 0.07 Sat, 0.09 GoodFri, 0.14 Aug, 0.17 AllSaints, 0.30 Xmas).
    Reverse stock (born in residence prov, living in destination) 0.22-0.50. Roots explain 33-45% of residual deviance for ALL slices, incl. ordinary weekdays
    -> roots tie is real but NOT specific to holidays (except Christmas). Risk: "people visit family on holidays" = obvious answer. Village-level exodus test pending.
Q7. ACHILLES' HEEL pilot (explore/p66; persistent weekday net, 2,423 mainland nodes, 47,929 edges): removing up to 400 nodes by degree/betweenness/strength
    or at random never fragments the network (giant = all remaining nodes). Spatial lattice-like robustness -> no Achilles' heel; drop as headline.
Q8. NEON BLIND IDEAS (research_notes/.../neon_blind_ideas.md; 30 ideas, 8 themes): #1 Q1 summer territories (bipartite city x destination, BiCM + gravity null) GAP;
    #2 Q12 AVE orbit capture GAP; #3 Q2 August rich-club inversion PARTIAL (Kruse et al. 2025 AAAG WTRC); #4 Q21 PLACSP minor contracts GAP;
    downgraded: Q7 political mixing (Liang & Cho 2026 PLOS ONE), Q6 income rich-club DONE (Duran-Sala 2024 arXiv 2407.01799), Q19 grid.
Q9. SUMMER TERRITORIES pilot (p67/p68; top-100 cities x 2,400 destinations, IPF gravity null with 2-25 km bins): Aug visitor-weighted HHI 0.234 vs null 0.178;
    28% vs 24% destinations dominated by one city; colonies beyond distance: Vitoria -> Cáceres villages x23 (500 km), Zaragoza -> Vila-seca/Salou x7, Logroño -> Salou x7,
    Pamplona -> Tarragona coast x8.6, Córdoba -> Torrox x7.8, Valencia -> Cuenca villages x8.5. Spatial-null bipartite communities noisy (Q_B 0.10-0.16) but reveal
    roots groupings: L'Hospitalet+Cornellà+Granada -> Almuñécar/Motril/Salobreña (Aug); Hospitalet+Almería+El Ejido+Sant Boi -> Almería coast (All Saints).
    City-level roots elasticity (prov-level birth matrix): Feb 0.28 -> Aug 0.42 -> All Saints 0.46.
Q10. RICH-CLUB by day (p70): rho(k) 1.0-1.5 on all days; hubs unchanged (Madrid, Barcelona, Valencia, Zaragoza...; Benidorm enters top-12 in Aug) -> NO August inversion. Drop D.
Q11. MUNICIPAL ROOTS (p72; Census 2024 birthplace for 151 big cities; other-province nights; PPML): elasticity to residents born in destination province
    with c + p + dist FE: 0.46 Feb, 0.57 Sat, 0.61 GoodFri/Aug/AllSaints/Xmas; WITHIN origin-prov x dest-prov (neighbouring cities): 0.63 Feb, 0.69 Sat, 0.72 Aug, 0.71 AllSaints, 0.70 Xmas.
    Barcelona province (19 cities): % born in Andalusia vs % Aug out-of-province nights in Andalusia r = 0.85 (Santa Coloma 14.8 -> 19.6; Sant Cugat 3.0 -> 5.7; Barcelona 5.2 -> 6.8).
    Madrid suburbs: NO relation (Andalusia = beach; rich suburbs Pozuelo/Boadilla 25% vs born 4%) -> income confound must be controlled.
    District level (p69): Barcelona district 05 (Sarrià) 3.7% of Aug away-nights in Andalusia vs Santa Coloma/Sant Boi/Badalona/Mataró districts 13-16%.
Q12. EXODUS -> REFILL (p73; INE municipal population 1900-2021 via scout, MITMA zones; coastal flag = touches outer land boundary):
    inland zones by change since 1900-1960 peak, median present/residents: Q1 (+406%) Aug 0.90 / AllSaints 1.00 / Feb 1.00 ... Q5 (-65%) Aug 1.48 / AllSaints 1.11 / Feb 0.98.
    Holds within size terciles (small 0.98->1.45, mid 0.97->1.52). OLS inland, prov FE + log residents: log mult Aug on loss rank +0.55 (se .05); All Saints +0.11 (se .04).
    Figure: reports/figs/fig9_exodus_holidays.png.
Q13. AVE natural experiments (p71; Oct Tue-Thu 2022/2023/2024, holiday eves dropped): Madrid–Asturias +24.4% 2024 vs 2023 (Pajares tunnel Nov-2023) vs median Madrid pair +10.3%
    -> 3rd of 51 Madrid pairs (behind Ceuta/Melilla), placebo percentile 0.94; León–Asturias flat (+0.6%); Murcia AVE (Dec-2022): +2.1% vs median +2.1% (no effect).
    Thin, one-corridor evidence; network-level 'orbit capture' implausible. Downloaded data/od_matrices_daily/2022 (10 days) and 2023 (9 days).
Q14. DATA SCOUT done (research_notes/.../neon_data_scout.md): birth x residence provinces 1998-2025 (12 snapshots) + 151 cities 2021-2025; municipal population 1842/1900-2025 (8,132 munis);
    interprovincial migration 1988-2024 (37 yrs); Meta SCI 2026 NUTS3 (ES, ES-PT/FR); 59 HSR stations with opening dates. FAILED: 2023/2019 municipal election results
    (infoelectoral behind Cloudflare challenge; manual browser download needed). Idea for B: do OLD migration flows (1988-95) predict holiday nights better than recent ones ("echo")?
