# Literature for question B: "Roots or beaches?"

Research note, 2026-11-24. BibTeX keys in brackets refer to `research_notes/references_B.bib`.

**How references were checked.** Every DOI below was resolved against the Crossref API (`api.crossref.org/works/<doi>`) on 2026-11-24, and authors, year, title, venue, volume and pages were compared with the Crossref record. arXiv items were checked on their arXiv abstract pages. Institutional documents without a DOI (INE, MITMS, CEDDAR, Banco de España working papers) were opened and read; the page numbers quoted are the ones we read. Anything we could not open or confirm is marked **UNVERIFIED**.

**Search limits.** OpenAlex refused requests (the shared daily quota for this IP was used up), and Dialnet timed out from this machine. Instead we used Crossref bibliographic search, general web search (English and Spanish), arXiv, SSRN pages found through web search, and the full texts of open PDFs. Section C should therefore be read as "nothing found in these sources", not as proof that nothing exists. Before submission, re-run two queries on Dialnet by hand (see C.4).

**DOI errors in the earlier notes** (`network_layer_A_B.md`, Part 4):
- Wesolowski et al. 2017: the correct DOI is `10.1038/s41467-017-02064-4`. The notes give `...-02070-z`, which does not resolve.
- Ratti et al. 2010: the correct DOI is `10.1371/journal.pone.0014248`. The notes drop a zero.
- Sarzynska et al.: Crossref dates the online version 2015. The journal issue is 2016, vol. 4(3).

---

## A. Critical review set (5 works)

These five cover the three claims of the coursework: holiday networks look different, holiday flows follow migration ties, and Spain's emptied interior fills up in summer. Each is reviewed for what it did, why it matters to us, what we borrow, and where it falls short.

### A1. Li, Ye, Deng, Liu & Liu (2016). Spring Festival travel rush in China [`li2016spring`]
Li, J., Ye, Q., Deng, X., Liu, Y., Liu, Y. (2016). Spatial-temporal analysis on Spring Festival travel rush in China based on multisource big data. *Sustainability*, 8(11), 1184. https://doi.org/10.3390/su8111184

- **What they did.** They built 31×31 province-level travel networks for the 2015 Spring Festival from three platforms (Baidu and Tencent location data, Qihoo train-ticket searches) and compared them with the 2010 census migration network (registered *hukou* province vs. current residence). They described the networks and correlated them with QAP. Table 2 (p. 7) reports census–Baidu r = 0.88 linear / 0.81 log, and census–Qihoo r = 0.82 log.
- **Why it matters to us.** This is the closest design precedent we found. It compares a holiday travel network with a network of where people come from, on the same nodes, using QAP.
- **What we borrow.** QAP as the test that two weighted networks on the same nodes are aligned, and the idea of a migration-stock network as the reference layer.
- **Limitations.**
  1. **No spatial null.** Both networks decay with distance, so a raw QAP correlation is inflated by geography alone. Our gravity-null residual QAP and within-province permutations address this.
  2. **Too few nodes.** With 31 provinces the matrix is nearly complete, so topology (degree, rich club, cores) cannot be studied. We have 2,484 zones and densities of 0.9–3.3 %.
  3. **The reference is recent migration.** *Hukou*-vs-residence captures current migrant workers, not a 50-year-old exodus.
  4. **Truncated data.** Baidu and Tencent publish only the top-10 inflows and outflows per city per collection (p. 3), which biases toward big links.
  5. **The aim is validation.** The census is used to check data quality, not as the hypothesis being tested.
  6. **Mismatched years.** Census 2010 against travel 2015.

### A2. Luan et al. (2024). Distance decay on homecoming trips [`luan2024different`]
Luan, X., Paryzat, H., Chu, J., Shu, X., Gu, H., Tong, D., Li, B. (2024). Different roads take me home: the nonlinear relationship between distance and flows during China's Spring Festival. *Humanities and Social Sciences Communications*, 11, 1356. https://doi.org/10.1057/s41599-024-03779-8

- **What they did.** They modelled intercity Tencent flows during the Spring Festival with gradient-boosting trees. They found three shapes for the flow–distance relationship: a plateau with almost no distance decay, a drop, and a rebound. The shapes differ by province, depending on relative location and the urban system within each province (Crossref abstract).
- **Why it matters to us.** Homecoming holidays weaken or reshape distance decay. A gravity model fitted to ordinary days is therefore the wrong null for holiday nights.
- **What we borrow.** Fit the doubly-constrained gravity null separately for each date, and test for non-linear distance terms before reading residuals as "roots".
- **Limitations.**
  1. Migrant origin is never measured. "Home" is inferred from the holiday itself, so the roots mechanism is assumed rather than tested.
  2. Gradient-boosting partial dependence describes the data but does not identify causes.
  3. The flows come from a platform with the same truncation problems as A1.
  4. There is no network-topology analysis.

### A3. Garzón-Paredes & Royo-Vela (2026). Seasonal population redistribution in Spain [`garzon2026heritage`]
Garzón-Paredes, A., Royo-Vela, M. (2026). Protected heritage stock and population redistribution in Spain: multivariate evidence from four 2021 mobile-phone mobility snapshots. *Sustainability*, 18(18), 9683. https://doi.org/10.3390/su18189683 (preprint: 10.20944/preprints202608.1368.v1)

- **What they did.** They used INE's aggregated phone indicators for the 3,214 Spanish mobility areas on four dates in 2021 and computed five indicators per area: resident retention, outbound mobility, inbound pressure, and diversity of connections. PCA and NMDS reduce these to three profiles: resident-dispersing (873 areas), broadly connected (992) and summer-receiving (1,344). Protected-heritage stock adds an incremental R² of 3.6–5.0 % (Crossref abstract).
- **Why it matters to us.** It is the most recent work using Spanish phone data on seasonal redistribution, and its "summer-receiving" interior areas are our emptied villages.
- **What we borrow.**
  1. The INE 3,214-area source as an independent replication. Our README already reports an INE three-operator check.
  2. Their caution that travel purpose is not observed, so dates are "snapshots", which matches our robustness design with sibling nights.
  3. Their honest handling of a weak covariate as a model of reporting.
- **Limitations.**
  1. Everything is at the node level. The origin–destination structure is reduced to a diversity score, and nothing is said about where visitors come from.
  2. There is no roots or birthplace variable, which is exactly our gap.
  3. 2021 was a pandemic year with only four dates.
  4. The heritage crosswalk covers 63.8 % of areas, skewed toward small ones, and an earlier "expert/AI-assisted" heritage score was withdrawn after failing validation.
  5. The authors say spatial dependence is unresolved.

### A4. González-Leonardo, Rowe & Fresolone-Caparrós (2022). Rural revival during COVID-19 [`gonzalezleonardo2022rural`]
González-Leonardo, M., Rowe, F., Fresolone-Caparrós, A. (2022). Rural revival? The rise in internal migration to rural areas during the COVID-19 pandemic. Who moved and where? *Journal of Rural Studies*, 96, 332–342. https://doi.org/10.1016/j.jrurstud.2022.11.006

- **What they did.** They used Spanish population-register (*padrón*) flows for 2020 to classify which rural areas gained internal migrants. In-migration rose most in rural areas close to cities and in areas with many second homes. Inflows covered all ages and included high-income foreign-born people (abstract, OSF preprint 10.31219/osf.io/g4wvd).
- **Why it matters to us.** It is direct evidence that second homes in rural Spain channel urban residents back to the countryside, which is the "casa del pueblo" mechanism we want to separate from beach tourism.
- **What we borrow.** Second-home prevalence as a mediator. The handoff lists the 2021 Census second-home data as the optional test. They also show that 2020 is an anomalous year, which supports our use of 2022–2024.
- **Limitations.**
  1. It measures registered changes of residence, not temporary overnight stays.
  2. Birthplace is coded only as Spain vs. abroad, so the municipality of family origin is unknown.
  3. It uses an area typology rather than a network.
  4. It covers one exceptional year.

### A5. Wesolowski et al. (2017). Seasonal asymmetry in phone-measured travel [`wesolowski2017seasonal`]
Wesolowski, A., zu Erbach-Schoenberg, E., Tatem, A. J., Lourenço, C., Viboud, C., Charu, V., Eagle, N., Engø-Monsen, K., Qureshi, T., Buckee, C. O., Metcalf, C. J. E. (2017). Multinational patterns of seasonal asymmetry in human movement influence infectious disease dynamics. *Nature Communications*, 8, 2069. https://doi.org/10.1038/s41467-017-02064-4

- **What they did.** They used phone records in Kenya, Namibia and Pakistan to measure seasonal travel and its direction. National holidays shift the scope of travel, routes gain and lose importance over the year, and rural-to-urban travel rises after holidays. They then show how these swings change the speed of pathogen spread.
- **Why it matters to us.** It is the standard network-science reference for holidays reorganising a national mobility network, and it frames direction (who goes where) as the interesting quantity.
- **What we borrow.** Treat holiday nights as a different network on the same nodes. Report in/out asymmetry: city residents sleeping in villages in August, then the reverse.
- **Limitations.**
  1. It does not link travel to migration history or origin.
  2. The data come from one operator per country, with unknown representativeness.
  3. It uses region-level aggregation with no topology metrics.
  4. The setting is very different (nomadic and low-income), so transfer to Spain is about method only.

**Reserve picks** if a reviewer wants different weights:
- Osorio Arjona (2024) [`osorio2024andalusia`]: Spanish OD phone data; rural population grows at weekends and population moves to the coast in summer, in Andalusia.
- Alonso-Pardo, Oso & Santaballa (2023) [`alonsopardo2023roots`]: 109 interviews on "roots migrants" in shrinking Spain.
- Sánchez-Oro Sánchez (2020) [`sanchezoro2020backhome`]: survey of back-home tourism in Extremadura.
- Kruse et al. (2025, arXiv) [`kruse2025richclub`]: a rich-club coefficient for spatial, weighted, temporal mobility networks, close to our rich-club result.

---

## B. Wider bibliography (verified)

### B.1 Network methods, including Prof. Zhou's
| Key | Reference | Use in the report |
|---|---|---|
| `zhou2004richclub` | Zhou, S., Mondragón, R. J. (2004). The rich-club phenomenon in the Internet topology. *IEEE Commun. Lett.* 8(3), 180–182. doi:10.1109/LCOMM.2004.823426 | Rich-club definition |
| `zhou2007structural` | Zhou, S., Mondragón, R. J. (2007). Structural constraints in complex networks. *New J. Phys.* 9, 173. doi:10.1088/1367-2630/9/6/173 | Assortativity forced by the degree sequence; rewiring logic |
| `guo2009fingerprint` | Guo, Y., Chen, C., Zhou, S. (2009). Fingerprint for network topologies. In *Complex Sciences* (Complex 2009), LNICST (vol. number not confirmed), 1666–1677. doi:10.1007/978-3-642-02469-6_45 | BOSAM plots |
| `zhou2017second` | Zhou, S., Cox, I. J., Hansen, L. K. (2017). Second-order assortative mixing in social networks. *Complex Networks VIII*, Springer Proc. Complexity, 3–15. doi:10.1007/978-3-319-54241-6_1 | Optional second-order mixing |
| `colizza2006richclub` | Colizza, V., Flammini, A., Serrano, M. A., Vespignani, A. (2006). Detecting rich-club ordering in complex networks. *Nat. Phys.* 2, 110–115. doi:10.1038/nphys209 | Normalised ρ(k) |
| `opsahl2008weighted` | Opsahl, T., Colizza, V., Panzarasa, P., Ramasco, J. J. (2008). Prominence and control: the weighted rich-club effect. *PRL* 101, 168702. doi:10.1103/PhysRevLett.101.168702 | Weighted rich club (appendix) |
| `maslov2002specificity` | Maslov, S., Sneppen, K. (2002). Specificity and stability in topology of protein networks. *Science* 296, 910–913. doi:10.1126/science.1065103 | Degree-preserving rewiring |
| `newman2002assortative` | Newman, M. E. J. (2002). Assortative mixing in networks. *PRL* 89, 208701. doi:10.1103/PhysRevLett.89.208701 | Degree assortativity r |
| `newman2003mixing` | Newman, M. E. J. (2003). Mixing patterns in networks. *PRE* 67, 026126. doi:10.1103/PhysRevE.67.026126 | Directed/weighted mixing |
| `batagelj2003cores` | Batagelj, V., Zaveršnik, M. (2003). An O(m) algorithm for cores decomposition of networks. arXiv:cs/0310049 (verified on arXiv) | k-core algorithm |
| `batagelj2011generalized` | Batagelj, V., Zaveršnik, M. (2011). Fast algorithms for determining (generalized) core groups in social networks. *Adv. Data Anal. Classif.* 5(2), 129–145. doi:10.1007/s11634-010-0079-y | Peer-reviewed k-core reference |
| `barber2007bipartite` | Barber, M. J. (2007). Modularity and community detection in bipartite networks. *PRE* 76, 066102. doi:10.1103/PhysRevE.76.066102 | Bipartite modularity Q_B |
| `saracco2015bicm` | Saracco, F., Di Clemente, R., Gabrielli, A., Squartini, T. (2015). Randomizing bipartite networks: the case of the World Trade Web. *Sci. Rep.* 5, 10595. doi:10.1038/srep10595 | BiCM null |
| `squartini2011analytical` | Squartini, T., Garlaschelli, D. (2011). Analytical maximum-likelihood method to detect patterns in real networks. *New J. Phys.* 13, 083001. doi:10.1088/1367-2630/13/8/083001 | Max-entropy nulls (optional) |
| `vallarano2021fast` | Vallarano, N. et al. (2021). Fast and scalable likelihood maximization for exponential random graph models with local constraints. *Sci. Rep.* 11, 15227. doi:10.1038/s41598-021-93830-4 | BiCM solver (cite only if NEMtropy is used) |
| `krackhardt1988predicting` | Krackhardt, D. (1988). Predicting with networks: nonparametric multiple regression analysis of dyadic data. *Social Networks* 10(4), 359–381. doi:10.1016/0378-8733(88)90004-4 | QAP / MRQAP |
| `dekker2007sensitivity` | Dekker, D., Krackhardt, D., Snijders, T. A. B. (2007). Sensitivity of MRQAP tests to collinearity and autocorrelation conditions. *Psychometrika* 72(4), 563–581. doi:10.1007/s11336-007-9016-1 | Double semi-partialling QAP |
| `barthelemy2011spatial` | Barthélemy, M. (2011). Spatial networks. *Phys. Rep.* 499(1–3), 1–101. doi:10.1016/j.physrep.2010.11.002 | Spatial embedding, gravity |
| `expert2011uncovering` | Expert, P., Evans, T. S., Blondel, V. D., Lambiotte, R. (2011). Uncovering space-independent communities in spatial networks. *PNAS* 108(19), 7663–7668. doi:10.1073/pnas.1018962108 | Spatial (gravity) modularity null |
| `sarzynska2016null` | Sarzynska, M., Leicht, E. A., Chowell, G., Porter, M. A. (2016). Null models for community detection in spatially embedded, temporal networks. *J. Complex Netw.* 4(3), 363–406. doi:10.1093/comnet/cnv027 | Gravity vs. radiation nulls |
| `simini2012universal` | Simini, F., González, M. C., Maritan, A., Barabási, A.-L. (2012). A universal model for mobility and migration patterns. *Nature* 484, 96–100. doi:10.1038/nature10856 | Radiation alternative |
| `santossilva2006log` | Santos Silva, J. M. C., Tenreyro, S. (2006). The log of gravity. *Rev. Econ. Stat.* 88(4), 641–658. doi:10.1162/rest.88.4.641 | PPML |
| `barbosa2018human` | Barbosa, H. et al. (2018). Human mobility: models and applications. *Phys. Rep.* 734, 1–74. doi:10.1016/j.physrep.2018.01.001 | Mobility review |

### B.2 Phone-based mobility, holidays and tourism
| Key | Reference | Relevance |
|---|---|---|
| `louail2015uncovering` | Louail, T. et al. (2015). Uncovering the spatial structure of mobility networks. *Nat. Commun.* 6, 6007. doi:10.1038/ncomms7007 | Spanish phone OD networks |
| `deville2014dynamic` | Deville, P. et al. (2014). Dynamic population mapping using mobile phone data. *PNAS* 111(45), 15888–15893. doi:10.1073/pnas.1408439111 | Seasonal population (France, Portugal) |
| `silm2010seasonal` | Silm, S., Ahas, R. (2010). The seasonal variability of population in Estonian municipalities. *Environ. Plan. A* 42(10), 2527–2546. doi:10.1068/a43139 | About 5 % of Estonians change residence seasonally, mostly from cities to coast and dacha areas |
| `ahas2008evaluating` | Ahas, R. et al. (2008). Evaluating passive mobile positioning data for tourism surveys. *Tourism Manage.* 29(3), 469–486. doi:10.1016/j.tourman.2007.05.014 | Phone data for tourism |
| `saluveer2020methodological` | Saluveer, E. et al. (2020). Methodological framework for producing national tourism statistics from mobile positioning data. *Ann. Tourism Res.* 81, 102895. doi:10.1016/j.annals.2020.102895 | Official-statistics framework |
| `lai2022global` | Lai, S. et al. (2022). Global holiday datasets for understanding seasonal human mobility and population dynamics. *Sci. Data* 9, 17. doi:10.1038/s41597-022-01120-z | Holiday calendars |
| `wang2014tracing` | Wang, X., Liu, C., Mao, W., Hu, Z., Gu, L. (2014). Tracing the largest seasonal migration on Earth. arXiv:1411.0983 (verified on arXiv) | Chunyun with smartphone data |
| `zaragozi2021passive` | Zaragozí, B., Trilles, S., Gutiérrez, A. (2021). Passive mobile data for studying seasonal tourism mobilities. *ISPRS IJGI* 10(2), 98. doi:10.3390/ijgi10020098 | Spanish coastal seasonality |
| `osorio2024andalusia` | Osorio Arjona, J. (2024). Analyzing post-COVID-19 demographic and mobility changes in Andalusia using mobile phone data. *Sci. Rep.* 14, 14828. doi:10.1038/s41598-024-65843-2 | Rural weekend gains |
| `osorio2026merindades` | Osorio Arjona, J., Pardo Abad, C. (2026). Los datos de telefonía móvil como instrumento de análisis del turismo de interior en el entorno rural de Las Merindades. *BAGE* 107. doi:10.21138/bage.3750 | Inland rural tourism, phone data |
| `kruse2025richclub` | Kruse, J., Gao, S., Ji, Y., Levin, K., Huang, Q., Mayer, K. R. (2025). Identifying rich clubs in spatiotemporal interaction networks. arXiv:2501.05636 (verified on arXiv) | Rich club in mobility networks |

### B.3 Data documentation
| Key | Reference | Status |
|---|---|---|
| `mitms2024methodology` | Ministerio de Transportes y Movilidad Sostenible (2024). *Estudio de movilidad con Big Data. Informe metodológico* (v8). https://www.transportes.gob.es/recursos_mfom/paginabasica/recursos/a3_informe_metodologico_estudio_movilidad_mitms_v8.pdf | URL listed by the ministry and cited as v8 by `kotov2026spanishoddata`. **Content not opened**: the server blocked scripted downloads. Check the version number by hand. A newer summary v11 is linked on the ministry's methodology page. |
| `kotov2026spanishoddata` | Kotov, E., Vidal-Tortosa, E., Cantú-Ros, O. G., Burrieza-Galán, J., Herranz, R., Gullón Muñoz-Repiso, T., Lovelace, R. (2026). spanishoddata: a package for accessing and working with Spanish open mobility big data. *EPB: Urban Analytics and City Science* 53(6), 1398–1410. doi:10.1177/23998083251415040 | Verified. Peer-reviewed description of the MITMA data (Nommon processing of Orange raw data). |
| `ine2020em1` | INE (2020). *Estudio de movilidad de la población a partir de la telefonía móvil (EM-1). Año 2019*. Nota de prensa, 20 July 2020. https://www.ine.es/prensa/experimental_em1.pdf | Opened. 3,214 areas, four reference days in 2019. Also lists the most frequent origin areas of summer populations. |
| `gonzalezolmos2024measuring` | González-Olmos, B., Velasco-Gimeno, M. (2024). Measuring tourism using mobile network data. *Spanish J. Stat.* 6, 7–19. doi:10.37830/sjs.2024.1.02 | Verified. INE's three-operator tourism statistic; each operator has about 25 % market share (p. 10). |
| `ine2023ecepov` | INE (2023). *ECEPOV 2021. Hogares y viviendas*. Nota de prensa, 22 Feb 2023. https://www.ine.es/prensa/ecepov_2021_feb.pdf | Opened (see D.2). |
| `ine2024etr3t` | INE (2024). *Encuesta de Turismo de Residentes (ETR/FAMILITUR). Tercer trimestre 2024*. Nota de prensa, 20 Dec 2024. https://www.ine.es/dyngs/Prensa/es/ETR3T24.htm | Opened (see D.2). |

### B.4 Rural exodus, depopulation, second homes and VFR (context)
| Key | Reference |
|---|---|
| `collantes2011peaceful` | Collantes, F., Pinilla, V. (2011). *Peaceful Surrender: The Depopulation of Rural Spain in the Twentieth Century*. Newcastle upon Tyne: Cambridge Scholars. ISBN 978-1-4438-2838-3. The book exists (checked on the publisher's page and Dialnet). **Contents not read.** |
| `pinilla2017rural` | Pinilla, V., Sáez, L. A. (2017). *Rural depopulation in Spain: genesis of a problem and innovative policies*. Informes CEDDAR 2017-2, Zaragoza. Opened. No DOI. |
| `recano2017demographic` | Recaño, J. (2017). The demographic sustainability of empty Spain. *Perspectives Demogràfiques* 7, CED, Barcelona. Opened. No DOI. |
| `clar2023rural` | Clar, E., Ayuda, M.-I. (2023). Rural migration and agricultural modernization: provincial Spain during its great rural exodus, 1960–1981. *Historia Agraria* 90, 223–255. doi:10.26882/histagrar.090e07c |
| `collantes2014reducing` | Collantes, F., Pinilla, V., Sáez, L. A., Silvestre, J. (2014). Reducing depopulation in rural Spain: the impact of immigration. *Popul. Space Place* 20(7), 606–621. doi:10.1002/psp.1797 |
| `ayuda2010locational` | Ayuda, M.-I., Collantes, F., Pinilla, V. (2010). From locational fundamentals to increasing returns: the spatial concentration of population in Spain, 1787–2000. *J. Geogr. Syst.* 12(1), 25–50. doi:10.1007/s10109-009-0092-x |
| `bover1999migrations` | Bover, O., Velilla, P. (1999). *Migrations in Spain: historical background and current trends*. Banco de España, Documento de Trabajo 9909. Opened. |
| `bentolila1991mismatch` | Bentolila, S., Dolado, J. J. (1991). Mismatch and internal migration in Spain, 1962–86. In F. Padoa Schioppa (ed.), *Mismatch and Labour Mobility*, CUP, 182–234. doi:10.1017/CBO9780511599316.010. Working-paper version (Banco de España DT 9006, 1990) opened. |
| `garciabarrero2023circular` | García-Barrero, J. (2023/2024). From circular to permanent: the economic assimilation of migrants during Spain's rural exodus, 1955–73. *Econ. Hist. Rev.* 77(3), 765–795. doi:10.1111/ehr.13294 |
| `alonsopardo2023roots` | Alonso-Pardo, P., Oso, L., Santaballa, L. (2023). Newcomers and 'roots migrants': chain migrations and the revitalization of rural shrinking areas in Spain. *J. Int. Migr. Integr.* 24(S5), 979–999. doi:10.1007/s12134-023-01036-6 |
| `gracia2022noemigraron` | Gracia Bernal, A., Gil-Albarova, A., Gracia Gil, A. (2022). Los que no emigraron del todo. Segunda residencia rural y población vinculada. *RIEM* 11(2), 36–55. doi:10.25115/riem.v11i2.4754 |
| `gracia2013retorno` | Gracia Bernal, A. (2013). Retorno rural turístico en la provincia de Zaragoza. *ROTUR* 6(1), 34–50. doi:10.17979/rotur.2013.6.1.1271 |
| `sanchezoro2020backhome` | Sánchez-Oro Sánchez, M. (2020). Back-home tourism in Extremadura, Spain: spending patterns and motivation. *Cuadernos de Turismo* 46, 71–91. doi:10.6018/turismo.451501 |
| `pino2014segunda` | Pino Artacho, J. (2014). Segunda residencia y residencia móvil en España. *PASOS* 12(1), 199–207. doi:10.25145/j.pasos.2014.12.015 |
| `modenes2007second` | Módenes, J. A., López-Colás, J. (2007). Second homes and compact cities in Spain: two elements of the same system? *TESG* 98(3), 325–335. doi:10.1111/j.1467-9663.2007.00400.x |
| `janta2015rethinking` | Janta, H., Cohen, S. A., Williams, A. M. (2015). Rethinking visiting friends and relatives mobilities. *Popul. Space Place* 21(7), 585–598. doi:10.1002/psp.1914. Lists "affirmations of identities and roots" among VFR practices. |
| `backer2012vfr` | Backer, E. (2012). VFR travel: it is underestimated. *Tourism Manage.* 33(1), 74–79. doi:10.1016/j.tourman.2011.01.027 |
| `iorio2013diaspora` | Iorio, M., Corsale, A. (2013). Diaspora and tourism: Transylvanian Saxons visiting the homeland. *Tourism Geogr.* 15(2), 198–232. doi:10.1080/14616688.2012.647327 |
| `bailey2018social` | Bailey, M., Cao, R., Kuchler, T., Stroebel, J., Wong, A. (2018). Social connectedness: measurement, determinants, and effects. *J. Econ. Perspect.* 32(3), 259–280. doi:10.1257/jep.32.3.259 |
| `camarero2019thinking` | Camarero, L., Oliva, J. (2019). Thinking in rural gap: mobility and social inequalities. *Palgrave Commun.* 5, 95. doi:10.1057/s41599-019-0306-x |

Suggested shortlist of about 15 for the coursework (cap of 20):
- Methods: Zhou & Mondragón 2004 and 2007, Guo–Chen–Zhou 2009, Colizza 2006, Maslov–Sneppen 2002, Newman 2002, Batagelj–Zaveršnik 2003, Barber 2007, Saracco 2015, Dekker 2007, Expert 2011 or Barthélemy 2011.
- Context and data: Kotov 2026, Clar & Ayuda 2023, Pinilla & Sáez 2017.
- The five works in section A.

---

## C. Novelty check

### C.1 Closest hits
| # | Work | What it links | Gap relative to us |
|---|---|---|---|
| 1 | Li et al. 2016 [`li2016spring`] | Spring Festival travel network vs. census migration network, QAP r ≈ 0.8 | China, 31 provinces, recent migrants, no spatial null, no topology |
| 2 | INE EM-1 (2020), EM-3 and EM-4 (2021) [`ine2020em1`] | Summer seasonal population by 3,214 areas; lists the main origin areas of summer visitors | Descriptive; no birthplace link, no network analysis |
| 3 | Garzón-Paredes & Royo-Vela 2026 [`garzon2026heritage`] | Seasonal redistribution across 3,214 Spanish areas, phone data | Node indicators, heritage covariate; no roots |
| 4 | *El turismo de retorno como factor de desarrollo en los espacios rurales de Castilla-La Mancha: una oportunidad para la España vaciada* (Dialnet record 9884317, conference paper, c. 2024) [`dialnet9884317`] | Quantifies "return tourism" with INE experimental phone statistics against municipal population and regulated tourism, in Castilla-La Mancha | **UNVERIFIED** beyond the Dialnet title and abstract snippet (authors, venue and pages not retrieved because Dialnet timed out). Appears to be regional and descriptive, with no birthplace network |
| 5 | Sánchez-Oro Sánchez 2020; Gracia Bernal 2013 and 2022 | Back-home / return tourism and "población vinculada" in Extremadura and Aragón | Surveys and secondary statistics; no national network |
| 6 | Alonso-Pardo et al. 2023 | "Roots migrants" returning to shrinking Spain | Qualitative, permanent moves |
| 7 | González-Leonardo et al. 2022 | Rural in-migration pulled by second homes, 2020 | Residential moves; Spain vs. foreign birthplace only |
| 8 | Bailey et al. 2018 | Facebook friendship networks echo historical US migration | Friendships, not travel; US |
| 9 | Luan et al. 2024; Wang et al. 2014 | Homecoming holiday flows change distance decay | No origin data |
| 10 | Silm & Ahas 2010; Wesolowski et al. 2017; Deville et al. 2014 | Seasonal phone-based population and travel | No link to migration origin |

### C.2 What is not new
- Holiday and seasonal mobility measured with phones, including in Spain (INE since 2020; Osorio Arjona; Garzón-Paredes & Royo-Vela).
- Correlating a holiday travel network with a migration-origin network by QAP (Li et al. 2016, China).
- The idea that Spanish rural second homes and summer "return tourism" come from the 1950–75 exodus. This is established in Spanish geography: García González 2009 (UNVERIFIED, book chapter, see C.4), Gracia Bernal 2013, Sánchez-Oro 2020, and the Pinilla & Sáez and Recaño reports.
- Rich club, k-core and assortativity in mobility networks in general (e.g. Kruse et al. 2025 on rich clubs).

### C.3 What appears to be new (no precedent found)
1. **National, fine-grained comparison in Spain.** Phone overnight-stay origin–destination networks on 2,484 mainland zones, compared directly with an INE birthplace-by-residence roots network.
2. **Alignment beyond geography.** It is tested against a doubly-constrained gravity null fitted for each date, with partial QAP controlling for the ordinary-night network, within-province permutations and a neighbouring-city design. Li et al. 2016 had none of these controls.
3. **Topology shift on holidays.** Assortativity against degree-preserving rewiring, normalised rich club, k-core promotion of emptied villages and BOSAM, applied to holiday vs. ordinary nights. We found no paper applying Zhou–Mondragón structural-constraint analysis or BOSAM to holiday mobility.
4. **Bipartite "summer territories".** Communities under both BiCM and a spatial gravity null, showing non-contiguous "diaspora archipelagos". Not found anywhere.
5. **An explicit link to the 1950–75 rural exodus** as the generating process, using a network-science design rather than surveys.

**Verdict: partially new.** The ingredients exist separately, and one Chinese study combines holiday travel with migration by QAP. The combination of a national phone overnight network, a birthplace roots network, spatial nulls and Zhou-style structural analysis, applied to the legacy of a historical rural exodus, has no precedent in what we searched. The report should state this carefully: "to our knowledge, the first … in Spain" and "extends Li et al. (2016) with spatial nulls and topology". It should not claim to be the first to compare holiday travel with migration.

### C.4 Checks to do by hand before submission
- Dialnet: search `"turismo de retorno" telefonía` and `"población vinculada" movilidad`. Retrieve the authors and venue of record 9884317.
- Google Scholar: `"MITMA" "lugar de nacimiento"` and `"mobile phone" "place of birth" holiday Spain`.
- García González, J. A. (2009). El turismo de retorno: modalidad oculta del turismo residencial. In T. Mazón, R. Huete, A. Mantecón (eds.), *Turismo, urbanización y estilos de vida*, Icaria, 351–365 [`garciagonzalez2009retorno`]. **UNVERIFIED** (seen only in web-search snippets; no DOI). Check before citing.

---

## D. Facts with sources

### D.1 Size of the 1950–75 rural exodus
No single primary figure for "people who left villages, 1950–75" could be opened and checked. These are the figures we did verify, with page references:

1. **Small municipalities lost about 4 million people.** In the 1950s just over **11 million** people lived in municipalities under 2,000 inhabitants, against about **7 million** in 2017. Their share of Spain's population fell from **39 % to 18 %**. In 1961–65 these municipalities lost about **100,000 people per year**. The rural population fell about **40 %** in the second half of the 20th century, by more than 50 % in the interior. Source: Pinilla & Sáez 2017, CEDDAR Informe 2017-2, p. 5 (English version) [`pinilla2017rural`].
2. **About two million people in the 1960s alone.** The migratory flow away from the countryside was "some two million people" in the decade before 1970, and 1.7 million agricultural jobs were lost in the 1960s. Source: Clar & Ayuda 2023, *Historia Agraria* 90, p. 225 (flow) and p. 229 (jobs, citing García Barbancho 1975: 84) [`clar2023rural`].
3. **Villages with 1,000 or fewer inhabitants lost almost 2.5 million.** Municipalities with ≤ 1,000 inhabitants in 2016 lost almost 2.5 million inhabitants between 1950 and 1991, and were the main contributors to the exodus toward Madrid, Barcelona and Valencia. Source: Recaño 2017, *Perspectives Demogràfiques* 7, p. 2 [`recano2017demographic`].
4. **Interprovincial migration rates.** These were 1.22 % of the population per year in 1962–69 and 1.09 % in 1970–75. Counting only interregional moves (excluding moves within a region), the rates were 0.65 % and 0.50 %. Source: Bentolila & Dolado, Banco de España DT 9006 (1990), printed p. 23; published as [`bentolila1991mismatch`].
   - *Our own arithmetic, label it as such:* applied to an average population of about 32 million (1962–69) and about 34.6 million (1970–75), these rates imply about **5.4 million registered interprovincial moves in 1962–75** (about 2.7 million interregional).
   - These are gross moves, so repeat and return moves count, and the series only starts in 1962.
5. **Origins and destinations.**
   - Net emigration regions in 1960–73: Andalusia, Castilla-La Mancha, Castilla y León and Extremadura. Net immigration regions: Catalonia, Madrid and the Basque Country. Source: Bover & Velilla 1999, BdE DT 9909, printed p. 15 [`bover1999migrations`].
   - Clar & Ayuda 2023, p. 232, find the highest provincial rural out-migration in the dry inland provinces of Castile, Andalusia and Extremadura. The gainers were Barcelona, Madrid, the Basque provinces and the Mediterranean and island tourism provinces.
   - Galicia's emigration went mainly abroad, at roughly double the per-capita rate of Andalusia or Extremadura (Bover & Velilla, p. 14).
6. **"5 million in the sixties".** A 2025 book review in *Historia y Política* (CEPC), by A. Fernández Asperilla, of M. Díaz Sánchez, *Fronteras de papel* (2024) [`fernandezasperilla2025review`], says about five million people changed residence during the 1960s and about two million emigrated abroad. This is a secondary source with no DOI found: https://www.cepc.gob.es/sites/default/files/2025-12/40977hyp-5420ana-fernandez-asperilla.html
7. **UNVERIFIED.** The round figures "5–6 million people, 1950–75" and "5.23 million interprovincial migrants, 1951–75" appear only on teaching sites and press pages. Do not cite them. The likely primary sources are García Barbancho, A. (1967), *Las migraciones interiores españolas* (Instituto de Estudios Económicos), and Nicolau, R. (2005), "Población, salud y actividad", in Carreras & Tafunell (eds.), *Estadísticas históricas de España*, Fundación BBVA. Neither was accessible here, so their page and table numbers are not given.

**Suggested sentence for the report:**
> Between 1950 and 1975 Spain's rural exodus moved millions of people from the interior, mainly Andalusia, Extremadura and the two Castiles, to Madrid, Barcelona and the Basque Country. About two million left the countryside in the 1960s alone (Clar & Ayuda 2023, p. 225), and municipalities under 2,000 inhabitants went from over 11 million residents in the 1950s to about 7 million today (Pinilla & Sáez 2017, p. 5).

### D.2 Second homes and "casa del pueblo"
- **15.5 % of households in Spain had a second home in 2021** (INE ECEPOV 2021, press release of 22 Feb 2023, p. 4) [`ine2023ecepov`]. The share is highest among households resident in the Comunidad de Madrid (22.1 %), Aragón (21.8 %), País Vasco (19.3 %), La Rioja (18.5 %) and Castilla y León (17.6 %); it is 12.6 % in Andalucía and 12.8 % in Extremadura (annex table, p. 10 of the same release). The ranking fits the exodus geography, with destination regions owning most second homes, but this is our reading, not INE's.
- **No reliable national figure was found for the share of second homes in the family's village of origin.** ECEPOV does not publish where the second home is relative to birthplace. Aragón has 139,109 households with a second residence (IAEST census figure, quoted by Gracia Bernal et al. 2022, abstract); that is regional and does not record origin. Press claims about "tener pueblo" are not survey-based. Treat any national "% with a casa del pueblo" as **unknown**.
- **Visits to family in summer.** In Q3 2024, 28.3 % of trips by Spanish residents (17.08 million of 60.4 million) were "visits to relatives or friends", averaging 6.7 nights (INE ETR/FAMILITUR Q3 2024, p. 2) [`ine2024etr3t`]. In Q4 2024 the share was 39.6 %, which includes All Saints (ETR Q4 2024, p. 2) [`ine2025etr4t`]. This is a useful framing for the All Saints control night.

---
## Orchestrator verification (2026-11-24)
- Crossref checked by orchestrator: Garzón-Paredes & Royo-Vela 2026 (Sustainability, 10.3390/su18189683) exists with that title/authors/year.
- The "turismo de retorno" hit is now VERIFIED: García-González, J. A., Escudero-Gómez, L. A. & Martínez-Navarro, J. M. (2024). "El turismo de retorno como factor de desarrollo en los espacios rurales de Castilla-La Mancha: una oportunidad para la España vaciada". In A. Ferrandis Martínez, C. Zornoza Gallego & J. V. Sánchez Cabrera (eds.), *Repensando los destinos turísticos en tiempos de cambio global*, pp. 304–320. Valencia: FGUV. ISBN 978-84-9133-737-9. Source: UIB research portal / Dialnet record.
  - What it does: quantifies "return/nostalgic tourism" (people spending leisure time in their municipality of origin) with INE experimental phone statistics, relates it to municipal population and regulated tourism, maps it; one region (Castilla-La Mancha); descriptive.
  - What it does not do: no network, no birthplace/roots matrix, no null model, no national comparison of holiday vs ordinary nights. Novelty of our combination holds; cite it as the closest Spanish precedent and as support for "return tourism is hidden in official statistics".
- Data lead for R1 / journal paper: INE now publishes annual "Número de turistas por municipio de origen y destino" (2019–2026) in the same experimental series (INEbase capsel=8581), i.e. a true municipality-level OD table without MITMA '_AM' aggregation.
- Crossref (orchestrator, 2026-11-24): all 5 review works resolve with matching authors/year/venue: Li 2016 su8111184; Luan 2024 s41599-024-03779-8 ('Different roads take me home'); Garzón-Paredes 2026 su18189683 (title: 'Protected Heritage Stock and Population Redistribution in Spain…' — check its phone-data framing when writing the review); González-Leonardo 2022 jrurstud.2022.11.006; Wesolowski 2017 s41467-017-02064-4.
