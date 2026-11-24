# Open Network-Science Datasets for Spain — Data Scout Report

**Date:** 2026-11-24  
**Workspace:** `/Users/alejandrotreny/Proyectos/Spain-Mobility`  
**External data root:** `/Users/alejandrotreny/Proyectos/Spain-Mobility/_local/external/`  
**Processing scripts:** `/Users/alejandrotreny/Proyectos/Spain-Mobility/_local/external/_scripts/`

---

## Summary Table — Priority Datasets

| # | Dataset | Status | Clean CSV Path(s) | Rows × Cols | Coverage |
|---|---|---|---|---|---|
| 1a | INE Province of Birth × Province of Residence | **Downloaded & built** | [`prov_birth_x_prov_residence.csv`](file:///Users/alejandrotreny/Proyectos/Spain-Mobility/_local/external/ine_birth_residence/prov_birth_x_prov_residence.csv) | 33,696 × 5 | 54 birth origins × 52 provinces × 12 snapshots (1998–2025) |
| 1b | INE Municipality of Residence × Province of Birth | **Downloaded & built** | [`mun_residence_x_prov_birth_censo_2021_2025.csv`](file:///Users/alejandrotreny/Proyectos/Spain-Mobility/_local/external/ine_birth_residence/mun_residence_x_prov_birth_censo_2021_2025.csv) | 40,770 × 5 | 151 major cities/capitals × 54 birth origins × 5 years (2021–2025) |
| 2a | INE Historical Municipal Population (1900–2025) | **Downloaded & built** | [`municipal_population_1900_2025.csv`](file:///Users/alejandrotreny/Proyectos/Spain-Mobility/_local/external/ine_historical_pop/municipal_population_1900_2025.csv) | 8,132 × 57 | All 8,132 current INE municipalities: Censos 1900–2021 + Padrón 1996–2025 |
| 2b | INE Intercensal Full Archive (1842–2021, incl. vanished) | **Downloaded & built** | [`municipal_censos_1842_2021_all_entities.csv`](file:///Users/alejandrotreny/Proyectos/Spain-Mobility/_local/external/ine_historical_pop/municipal_censos_1842_2021_all_entities.csv) | 11,927 × 40 | All current (8,131) + historical vanished/merged (3,796) entities, 1842–2021 |
| 3 | INE Interprovincial Migration Flows (1988–2024) | **Downloaded & built** | [`interprovincial_migration_1988_2024.csv`](file:///Users/alejandrotreny/Proyectos/Spain-Mobility/_local/external/ine_interprovincial_migration/interprovincial_migration_1988_2024.csv) | 100,776 × 5 | 52 × 51 directed pairs × 37 consecutive years (38 year-source slices) |
| 4a | Meta/Facebook SCI — Spain (NUTS3) | **Downloaded & built** | [`sci_spain_nuts3.csv`](file:///Users/alejandrotreny/Proyectos/Spain-Mobility/_local/external/meta_sci_spain/sci_spain_nuts3.csv) | 3,481 × 7 | 59 × 59 Spanish NUTS3 pairs (2026 release, NUTS 2024) |
| 4b | Meta/Facebook SCI — Spain ↔ PT/FR (NUTS3) | **Downloaded & built** | [`sci_spain_pt_fr_nuts3.csv`](file:///Users/alejandrotreny/Proyectos/Spain-Mobility/_local/external/meta_sci_spain/sci_spain_pt_fr_nuts3.csv) | 14,986 × 7 | Spain ↔ Portugal (26) & France (101) NUTS3 pairs (both directions) |
| 5 | Spanish General Elections (July 2023 & Nov 2019, municipal) | **Blocked (Cloudflare)** | *Not downloaded* | — | `infoelectoral.interior.gob.es` returns HTTP 403 (`cf-mitigated: challenge`) |
| 6 | Spanish High-Speed Rail (AVE / AV) Stations | **Compiled & verified** | [`spain_hsr_stations.csv`](file:///Users/alejandrotreny/Proyectos/Spain-Mobility/_local/external/spain_hsr_stations/spain_hsr_stations.csv) | 59 × 13 | 59 HSR/high-perf passenger stations across 32 provinces (1992–2023) |

---

## Priority Datasets — Detailed Notes, Previews & Sanity Checks

### 1. INE Population by Province of Birth × Province of Residence (1998–2025)

- **Directory:** [`_local/external/ine_birth_residence/`](file:///Users/alejandrotreny/Proyectos/Spain-Mobility/_local/external/ine_birth_residence/README.md)
- **Script:** [`build_birth_residence.py`](file:///Users/alejandrotreny/Proyectos/Spain-Mobility/_local/external/_scripts/build_birth_residence.py)
- **Clean CSVs:**
  - `prov_birth_x_prov_residence.csv`: **33,696 rows × 5 columns** (`year, source, prov_birth, prov_residence, persons`), 1.14 MB.
  - `mun_residence_x_prov_birth_censo_2021_2025.csv`: **40,770 rows × 5 columns** (`year, ine_code, municipality_name, prov_birth, persons`), 1.29 MB.
- **Sources:**
  - *Estadística del Padrón Continuo* (`0ccaa001.px`, 1998, 2000, 2004, 2009, 2014, 2019, and 2022 provisional).
  - *Censo Anual de Población* (2021–2025): INEbase table `68523` (provincial) and table `68545` (151 provincial capitals & major cities).

**5-row preview (`prov_birth_x_prov_residence.csv`):**
```csv
year,source,prov_birth,prov_residence,persons
1998,padron_continuo,01,01,160752
1998,padron_continuo,02,01,284
1998,padron_continuo,03,01,147
1998,padron_continuo,04,01,108
1998,padron_continuo,05,01,1380
```

**5-row preview (`mun_residence_x_prov_birth_censo_2021_2025.csv`):**
```csv
year,ine_code,municipality_name,prov_birth,persons
2021,01059,Vitoria-Gasteiz,01,140764
2021,01059,Vitoria-Gasteiz,02,162
2021,01059,Vitoria-Gasteiz,03,280
2021,01059,Vitoria-Gasteiz,04,125
2021,01059,Vitoria-Gasteiz,05,719
```

**Sanity checks:**
- Every single one of the 12 provincial snapshots has exactly **54 birth categories** (52 INE provinces + `ABROAD` + `FORMER_ES`) × **52 residence provinces** = **2,808 rows**, summing to the **exact** official INE total population of Spain down to the single person:
  - `1998` (`padron_continuo`): **39,852,601**
  - `2000` (`padron_continuo`): **40,499,751**
  - `2004` (`padron_continuo`): **43,197,684**
  - `2009` (`padron_continuo`): **46,745,807**
  - `2014` (`padron_continuo`): **46,771,341**
  - `2019` (`padron_continuo`): **47,026,208**
  - `2021` (`censo_anual_68523`): **47,400,798**
  - `2022` (`padron_continuo_provi`): **47,435,597** vs `2022` (`censo_anual_68523`): **47,486,727** (+0.11% methodological change)
  - `2023` (`censo_anual_68523`): **48,085,361**
  - `2024` (`censo_anual_68523`): **48,619,695**
  - `2025` (`censo_anual_68523`): **49,128,297**
- Municipal table (`68545`) covers 151 major cities (2025 sum = 26,117,396 = 53.2% of Spain; Madrid `28079` = 3,477,497).

**Problems / caveats:**
- INE retired *Estadística del Padrón Continuo* after 1 Jan 2022 (provisional) in favour of the annual register-based *Censo de Población* (2021–2025). Filter on `source` when taking a single 2022 cross-section.
- INE does not publish all 8,131 municipalities × 52 provinces of birth due to disclosure control; table `68545` (151 cities) is the finest official municipal cross-tabulation.

---

### 2. INE Historical Municipal Population (1900–2021 Censos + 1996–2025 Padrón)

- **Directory:** [`_local/external/ine_historical_pop/`](file:///Users/alejandrotreny/Proyectos/Spain-Mobility/_local/external/ine_historical_pop/README.md)
- **Scripts:** [`download_intercensal_all.py`](file:///Users/alejandrotreny/Proyectos/Spain-Mobility/_local/external/_scripts/download_intercensal_all.py) and [`build_historical_pop.py`](file:///Users/alejandrotreny/Proyectos/Spain-Mobility/_local/external/_scripts/build_historical_pop.py)
- **Clean CSVs:**
  - `municipal_population_1900_2025.csv`: **8,132 rows × 57 columns** (`ine_code, province_code, name, pop_1900..pop_2021 [Derecho/Residente], hecho_1900..hecho_1991, padron_1996..padron_2025, alterations_note`), 2.60 MB.
  - `municipal_censos_1842_2021_all_entities.csv`: **11,927 rows × 40 columns** (`ine_code, province_code, name, is_historical_6digit, in_ine_2026, derecho_1842..derecho_2021, hecho_1857..hecho_1991, alterations_note`), 2.96 MB.
- **Sources:**
  - INE *Alteraciones de los municipios en los Censos de Población desde 1842* (`https://www.ine.es/intercensal/`), downloaded in full across all 52 provinces (11,927 municipal entities total: 9,839 with 5-digit codes and 2,088 historical vanished entities with 6-digit codes).
  - INE *Cifras Oficiales de Población de los Municipios Españoles: Revisión del Padrón Municipal*, table `29005` (`1996`–`2025`, 8,138 municipal rows covering all 8,132 current municipalities plus 6 post-1996 merged/reassigned ones).

**5-row preview (`municipal_population_1900_2025.csv`, selected columns):**
```csv
ine_code,province_code,name,pop_1900,pop_1950,pop_1981,pop_2001,pop_2021,padron_2024,padron_2025,alterations_note
01001,01,Alegría-Dulantzi,650,692,1037,1533,2924,2971,2955,[1] En estos Censos se denominaba Alegría
01002,01,Amurrio,1174,2152,9032,9460,10295,10330,10364,< Entre el Censo de 1981 y el anterior, crece el término del municipio porque incorpora a 01007 (Arrastaria) y a 01035 (Lezama)
01003,01,Aramaio,2109,1791,1416,1441,1434,1381,1354,[1] En estos Censos se denominaba Aramayona
01004,01,Artziniega,1118,1255,1174,1339,1810,1856,1868,[1] En este Censo se denominaba Arceniega | [2] En estos Censos se denominaba Arciniega
01006,01,Armiñón,379,344,129,167,225,247,236,
```

**Sanity checks:**
- **Zero bad 5-digit `ine_code`s** across all 8,132 rows; **0** current INE 2026 municipalities missing from Padrón `29005`, and only **1** missing from 1842–2021 Intercensal (`48916 Usansolo`, created by secession from Galdakao in 2023 after the 2021 Census).
- **National Census totals (*Población de Derecho / Residente*), all 11,927 Intercensal entities vs current 8,132 municipalities only:**
  - `1900`: **18,830,649** (all entities — exact official INE 1900 total) vs `17,967,633` (current 2026 codes only; 863,016 = **4.58%** lived in municipalities later merged/extinguished).
  - `1950`: **28,174,268** (all) vs `27,502,070` (current; 2.39% in later-vanished municipalities).
  - `1981`: **37,682,355** (all) vs `37,658,912` (current; 0.06%).
  - `2001`: **40,847,371** (all — exact official INE Census 2001) vs `40,834,682` (current; 0.03%).
  - `2011`: **46,816,010** (all) vs `46,804,187` (current; 0.03%).
  - `2021`: **47,400,798** (all = current — exact to the single person matching Censo 2021 in Dataset 1!).
- **Official Padrón totals (`29005`, summing all 8,132 current municipalities):** `1996` = 39,655,367; `2001` = 41,103,825; `2011` = 47,178,310; `2021` = 47,385,107; `2024` = 48,587,813; `2025` = **49,114,494** (exact sum of `29005.csv`; differs by 0.028% from Censo Anual 2025 `49,128,297` due to official Padrón vs register-Census definitions).

**How merged/split municipalities are handled & caveats:**
- Following official INE methodology, `municipal_population_1900_2025.csv` preserves each municipality's **contemporaneous boundaries** (no synthetic back-casting), leaving pre-secession cells empty and recording every merger/secession/rename in `alterations_note` (`@` = created by secession/merger, `<` = absorbed another municipality, `#` = disappeared into successor).
- For **constant-boundary** 1900–2025 panels (recovering the 4.58% of 1900 population in vanished municipalities), use `municipal_censos_1842_2021_all_entities.csv`, which preserves all 3,796 vanished historical entities and their exact `#` successor INE codes.
- Note that `padron_1997` is empty across all rows because INE did not publish an official 1997 Padrón revision, and `hecho_*` ends at 1991 because INE retired *Población de Hecho* from 2001 onward.

---

### 3. INE Interprovincial Migration Flows (1988–2024)

- **Directory:** [`_local/external/ine_interprovincial_migration/`](file:///Users/alejandrotreny/Proyectos/Spain-Mobility/_local/external/ine_interprovincial_migration/README.md)
- **Script:** [`build_interprovincial_migration.py`](file:///Users/alejandrotreny/Proyectos/Spain-Mobility/_local/external/_scripts/build_interprovincial_migration.py)
- **Clean CSV:** `interprovincial_migration_1988_2024.csv`: **100,776 rows × 5 columns** (`year, source, prov_origin, prov_dest, migrations`), 2.69 MB.
- **Sources:**
  - `EVR_microdata` (1988–2005): Official INE EVR individual microdata (`datos_1988.zip`–`datos_2005.zip`, fixed-width 37-byte layout from `disreg_vr.xlsx`), aggregated to `PROVBAJA` → `PROVALTA` for domestic interprovincial moves (`01`–`52`, `origin != dest`).
  - `EVR_px` (2006–2007): Official INE EVR PC-Axis tables `0101311.px` (Spanish nationals) + `0101312.px` (foreign nationals).
  - `EM_24379` (2008–2021): *Estadística de Migraciones* table `24379`.
  - `EMCR_69742` (2021–2024): *Estadística de Migraciones y Cambios de Residencia* table `69742`.

**5-row preview (`interprovincial_migration_1988_2024.csv`):**
```csv
year,source,prov_origin,prov_dest,migrations
1988,EVR_microdata,01,02,12
1988,EVR_microdata,01,03,93
1988,EVR_microdata,01,04,16
1988,EVR_microdata,01,05,2
1988,EVR_microdata,01,06,20
```

**Sanity checks:**
- **Exact microdata-vs-published verification (2005):** Aggregating `datos_2005.zip` yielded **626,324** interprovincial moves, matching INE's published 2005 PC-Axis tables (`0101311.px` + `0101312.px`) with **0 differing pairs out of 2,652** (exact to the single move).
- **Selected annual totals:** 1988 = 280,378; 1995 = 318,581; 2000 = 408,578; 2005 = 626,324; 2007 (EVR peak) = 744,716; 2008 (EM start) = 524,445; 2019 (EM) = 475,859; 2020 (COVID dip, EM) = 416,902; 2021 (EM) = 482,346 vs 2021 (EMCR) = 653,128; 2024 (EMCR) = 700,640.
- **Origin/destination directionality verified in 2021 overlap:** Madrid (`28`) net balance is −16,054 in `EM_24379` (in 62,546 / out 78,600) and −15,965 in `EMCR_69742` (in 89,383 / out 105,348); suburbanization `28 → 45` (Madrid→Toledo) exceeds `45 → 28` in both (`13,056` vs `8,690` in EM; `17,659` vs `12,503` in EMCR).

**Problems / caveats:**
- **Two major methodological breaks:** (1) `EVR` (1988–2007) → `EM` (2008–2021) drops total volume by ~30% because EM filters multiple/administrative moves; (2) `EM` (2008–2021) → `EMCR` (2021–2024) raises total volume by +35.4% in 2021 due to the new register-based census methodology. Both 2021 series are retained so users can splice or normalize across the break.

---

### 4. Meta / Facebook Social Connectedness Index (SCI) — Spain (NUTS3, 2026 Release)

- **Directory:** [`_local/external/meta_sci_spain/`](file:///Users/alejandrotreny/Proyectos/Spain-Mobility/_local/external/meta_sci_spain/README.md)
- **Script:** [`build_meta_sci.py`](file:///Users/alejandrotreny/Proyectos/Spain-Mobility/_local/external/_scripts/build_meta_sci.py)
- **Clean CSVs:**
  - `sci_spain_nuts3.csv`: **3,481 rows × 7 columns** (`user_nuts3, friend_nuts3, user_prov_ine, friend_prov_ine, user_nuts3_name, friend_nuts3_name, scaled_sci`), 0.15 MB.
  - `sci_spain_pt_fr_nuts3.csv`: **14,986 rows × 7 columns** (`user_country, friend_country, user_nuts3, friend_nuts3, user_prov_ine, friend_prov_ine, scaled_sci`), 0.41 MB.
- **Source & Release:** Meta Data for Good via HDX (`https://data.humdata.org/dataset/social-connectedness-index`, resource `nuts_2024.zip`), **2026 release** (snapshot window `2025-12-26` to `2026-01-25`, HDX modified `2026-11-22`, NUTS 2024 boundaries, CC0 Public Domain).

**5-row preview (`sci_spain_nuts3.csv`):**
```csv
user_nuts3,friend_nuts3,user_prov_ine,friend_prov_ine,user_nuts3_name,friend_nuts3_name,scaled_sci
ES111,ES111,15,15,A Coruña,A Coruña,5530009
ES111,ES112,15,27,A Coruña,Lugo,1046845
ES111,ES113,15,32,A Coruña,Ourense,563347
ES111,ES114,15,36,A Coruña,Pontevedra,751420
ES111,ES120,15,33,A Coruña,Asturias,109355
```

**Sanity checks:**
- Complete 59 × 59 = **3,481** Spanish NUTS3 ordered pairs (47 peninsular provinces + Ceuta + Melilla + 3 Balearic islands + 7 Canary islands).
- **0 asymmetric pairs** (`scaled_sci(i, j) == scaled_sci(j, i)` for all pairs).
- Diagonal range: 356,389 to 543,435,109; off-diagonal range: 13,235 to 5,839,419 (highest off-diagonal ties are between neighbouring Canary Islands: El Hierro–La Gomera 5,839,419, La Gomera–Tenerife 3,253,319).
- Cross-border file contains all 59 × 26 × 2 = 3,068 ES↔PT pairs and 59 × 101 × 2 = 11,918 ES↔FR pairs (14,986 total).

**Problems / caveats:**
- NUTS3 is finer than GADM2 for Spain (59 units vs 52 provinces), because Balearic (`07`) and Canary (`35`, `38`) provinces are split by island. Joining to 52-province INE matrices requires population-weighted averaging across island NUTS3 units.

---

### 5. Spanish General Elections (July 2023 & Nov 2019) by Municipality — **BLOCKED**

- **Target official source:** Ministerio del Interior — *Infoelectoral* (`https://infoelectoral.interior.gob.es/estaticos/docxl/apliextr/02202307_MUNI.zip` for July 2023 Congress by municipality, `02201911_MUNI.zip` for Nov 2019 Congress by municipality).
- **Status:** **Blocked by anti-bot challenge.** Both URLs return `HTTP/2 403` with header `cf-mitigated: challenge` and CSP pointing to `https://challenges.cloudflare.com` (Cloudflare Turnstile / JS anti-bot gate).
- **Action taken:** Stopped immediately without attempting to bypass the Cloudflare challenge or fetch cached/mirrored copies to circumvent it.
- **How the user can retrieve it directly:** Open the URLs above in a normal browser, complete the Cloudflare check, and save `02202307_MUNI.zip` and `02201911_MUNI.zip` into `_local/external/elections_spain_municipal/raw/`. Inside each ZIP are fixed-width Ministry of Interior files (`03` = candidacies/parties, `05` = municipal census/turnout/blank/null votes, `06` = votes per party per municipality) directly parseable via Python or the `rOpenSpain/infoelectoral` R package. Key "España Vaciada" candidacies to flag in July 2023: *Teruel Existe – Coalición Existe* (`TERUEL EXISTE` / `EXISTE`), *Soria ¡YA!* (`SORIA ¡YA!`), *¡Vamos Palencia!* (`¡VP!`), *Aragón Existe*, *Cuenca Ahora*, *España Vaciada* (plus regionalist provincial lists such as *UPL*, *Por Ávila*, *PRC*, *CC*, *BNG*).

---

### 6. Spanish High-Speed Rail (AVE / Alta Velocidad) Stations

- **Directory:** [`_local/external/spain_hsr_stations/`](file:///Users/alejandrotreny/Proyectos/Spain-Mobility/_local/external/spain_hsr_stations/README.md)
- **Script:** [`build_hsr_stations.py`](file:///Users/alejandrotreny/Proyectos/Spain-Mobility/_local/external/_scripts/build_hsr_stations.py)
- **Clean CSV:** `spain_hsr_stations.csv`: **59 rows × 13 columns** (`station_name, adif_code, ine_code, municipality_name, province_code, province_name, line, hsr_opening_date, gauge_type, lat, lon, uncertain, notes`), 0.01 MB.
- **Sources:** Wikipedia (es) *Alta velocidad ferroviaria en España* (Adif *Declaración sobre la Red* tables & opening timeline) + 59 individual station infoboxes (ADIF station codes, WGS84 coordinates), verified against the official INE 2026 municipality register (`_local/external/historical/ine_municipios_2026.csv`).

**5-row preview (`spain_hsr_stations.csv`):**
```csv
station_name,adif_code,ine_code,municipality_name,province_code,province_name,line,hsr_opening_date,gauge_type,lat,lon,uncertain,notes
Ciudad Real,37200,13034,Ciudad Real,13,Ciudad Real,LAV Madrid-Sevilla,1992-04-21,standard_1435mm,,,false,Opened with LAV Madrid-Sevilla.
Córdoba-Julio Anguita,50500,14021,Córdoba,14,Córdoba,LAV Madrid-Sevilla; LAV Córdoba-Málaga,1992-04-21,standard_1435mm,37.888482,-4.789578,false,AVE service started 1992-04-21 at provisional/old site; new underground station inaugurated 1994-09-09.
Madrid-Puerta de Atocha-Almudena Grandes,60000,28079,Madrid,28,Madrid,LAV Madrid-Sevilla; LAV Madrid-Barcelona; LAV Madrid-Levante,1992-04-21,standard_1435mm,,,false,First HSR line in Spain; commercial AVE service began 1992-04-21 (inaugurated 1992-04-14).
Puertollano,37700,13071,Puertollano,13,Ciudad Real,LAV Madrid-Sevilla,1992-04-21,standard_1435mm,,,false,Opened with LAV Madrid-Sevilla.
Sevilla-Santa Justa,51003,41091,Sevilla,41,Sevilla,LAV Madrid-Sevilla,1992-04-21,standard_1435mm,37.392470,-5.974809,false,Station building opened 1991-05-02; AVE service began 1992-04-21.
```

**Sanity checks:**
- All **59** `ine_code` values pass exact assertion against the official INE 2026 5-digit municipality table (covering **32 provinces**; opening dates span `1992-04-21` to `2023-11-30`).
- **45** stations are on standard-gauge (`1,435 mm`) or mixed-gauge LAVs; **14** stations are on Iberian-gauge HSR / high-performance lines (Galicia Eje Atlántico, Asturias via Variante de Pajares, Extremadura Plasencia–Badajoz).
- **19** rows are explicitly flagged `uncertain=true` — either for stations straddling municipal borders (`Guadalajara-Yebes` in Yebes `19326`, `Camp de Tarragona` in La Secuita `43144`, `Figueres-Vilafant` in Vilafant `17221`, `Callosa de Segura-Cox` in Cox `03058`) or for Iberian/mixed-gauge corridors where "AVE opening" depends on whether one counts initial Avant/Alvia HSR or variable-gauge AVE S-106 (May 2024).

**Problems / caveats:**
- Wikimedia's `coordinates` API property is populated for 27 of the 59 station articles (32 have empty `lat`/`lon`, though exact coordinates can be joined from Renfe/ADIF GTFS `stops.txt` via `adif_code`).

---

## Catalogue of Other Spanish Open Datasets with Network Structure

| Dataset | Official URL(s) | Network / Granularity | Years | Licence | Approx. Size | Difficulty / Notes |
|---|---|---|---|---|---|---|
| **Renfe / ADIF GTFS (AVE, Larga Distancia, Media Distancia, Cercanías, Rodalies)** | https://data.renfe.com/ (CKAN portal; static GTFS ZIPs + GTFS-RT) | Station-to-station scheduled rail graph (`stops.txt` with ADIF codes, `stop_times.txt`, `trips.txt`); national + all Cercanías nuclei | Current schedule (rolling daily updates; archive snapshots via Mobility Database / Transitland) | CC BY 4.0 (Renfe Open Data) | ~20–60 MB (all feeds zipped) | **Easy.** Standard GTFS; note separate ZIPs for Cercanías vs Larga/Media Distancia and separate regional operators (FGC, Euskotren, FGV, SFM) on NAP (`nap.transportes.gob.es`). |
| **AENA Airport-to-Airport Passengers & Routes** | https://www.aena.es/es/estadisticas/inicio.html (and `datos.gob.es` / Dataestur) | Airport-to-airport OD passenger, operation, and cargo counts (46 AENA airports × domestic & international destinations) | 2004–present (monthly & annual) | Reuse allowed with attribution (AENA / Spanish Law 37/2007) | ~5–50 MB/year (XLSX / PowerBI export) | **Medium.** Monthly summaries are clean XLSX; full airport-to-airport OD matrices require querying AENA's Estadísticas interactive portal or Eurostat `avia_par_es` (easy TSV via Eurostat API for airport-pair flows!). |
| **INE Census 2021 Commuting OD (Movilidad Cotidiana)** | https://www.ine.es/censos2021/ (INEbase *Censo 2021* + *ECEPOV* microdata) | Municipality-of-residence × municipality-of-work/study commuting flows; transport mode & commute duration (full OD for municipalities ≥50k & capitals; microdata ~10% sample / register linkage) | 2021 (also 2011 and 2001 Censuses for longitudinal comparison) | INE Open Data (free reuse with attribution) | ~100 MB–2 GB (microdata) | **Medium–Hard.** Aggregate tables cover capitals & >50k cities; small-municipality OD pairs are censored/grouped in microdata for disclosure control. |
| **BORME Company Officers & Interlocking Directorates (and LibreBORME)** | Official daily PDFs/XML: https://www.boe.es/diario_borme/ ; Open parser: https://github.com/PabloCastellano/bormeparser & https://libreborme.net/ | Bipartite network of individuals (officers/directors) ↔ Spanish companies (*sociedades mercantiles*), with province of mercantile registry (`prov`) and appointment/resignation dates | 2009–present (daily PDFs since 2009) | Spanish public-domain / BOE Open Reuse; LibreBORME parser AGPL/GPL | ~10–30 GB raw PDFs; ~1–3 GB parsed DB | **Hard (from raw PDFs) / Medium (via `bormeparser`).** Official BORME publishes legal acts as unstructured Spanish text inside daily PDFs; entity resolution of person names (no DNI in BORME) is the main research caveat. |
| **PLACSP Public Procurement (Contracting Authorities ↔ Suppliers)** | https://contrataciondelestado.es/wps/portal/datosAbiertos (bulk ZIPs of ATOM/CODICE XML feeds + OpenPLACSP tool) | Directed weighted network of Spanish public contracting bodies (municipalities, provinces, ministries, hospitals) → winning firms (with NIF/CIF, CPV sector code, amount, procedure type) | 2012–present (main tenders), 2016–present (aggregated regional platforms), 2018–present (minor contracts) | Spanish Public Sector Open Data (Law 37/2007) | ~1–5 GB/year (zipped XML) | **Medium–Hard.** Deeply nested CODICE/UBL XML inside paginated ATOM feeds; official `OpenPLACSP` or community parsers (`contratacionabierta`) help flatten to Parquet/CSV. |
| **REE / ENTSO-E Transmission Grid Topology (+ PyPSA-Eur / OSM)** | REE REData API: https://www.ree.es/es/apidatos ; ENTSO-E Grid Map: https://www.entsoe.eu/data/map/ ; PyPSA-Eur: https://github.com/PyPSA/pypsa-eur | High-voltage substations (220 kV / 400 kV) and AC/HVDC transmission lines (incl. Portugal, France, Morocco, Balearic interconnectors); REE API gives regional generation/demand/exchange time series | Static/versioned topology + hourly/daily flows 2014–present | PyPSA-Eur MIT / CC BY 4.0; OSM ODbL; REE Open Data | ~10–100 MB (grid topology CSV/GeoJSON) | **Easy (PyPSA-Eur / OSM Overpass) to Medium (REE API).** `PyPSA-Eur` provides clean, pre-validated substation & line CSVs (`entsoegridkit` / OSM-based) ready for network analysis. |
| **University Student Mobility Flows (SIIU / UNIVbase)** | https://estadisticas.universidades.gob.es/ and https://www.ciencia.gob.es/Ministerio/Estadisticas/SIIU.html | Province of family residence × University/Province of study (domestic student migration matrix) + Erasmus/international incoming/outgoing flows | ~2015–present (annual academic years) | Spanish Government Open Data | ~5–50 MB (PC-Axis / CSV / XLSX) | **Easy–Medium.** UNIVbase serves PC-Axis/CSV tables; province-of-origin × university-of-destination matriculation matrices are directly usable for youth talent-drain networks. |
| **INE Tourism by Mobile Phone (Estadística Experimental)** | https://www.ine.es/experimental/turismo_moviles/experimental_turismo_moviles.htm | Province/municipality-of-residence × province/municipality-of-destination domestic tourists (and country-of-origin × destination for inbound tourism), monthly | July 2019–present (monthly) | INE Open Data | ~50–300 MB (CSV/XLSX) | **Easy.** Clean tabular CSV/XLSX downloads directly from INE; cells below anonymity threshold (<15 or <30 tourists) are suppressed. |
| **MITMA / Transportes Big Data Mobility ("Movilidad obligada" & Daily OD Matrices)** | https://www.transportes.gob.es/ministerio/proyectos-singulares/estudios-de-movilidad-con-big-data/opendata-movilidad (also S3 bucket & `spanishoddata` R package) | ~3,200 mobility areas (municipalities / districts / GAU) × daily/hourly OD trips, trip distance, residence-work-study ("movilidad obligada" / recurrent trips), overnight stays (*pernoctaciones*) | Feb 2020–May 2021 (v1) & Jan 2022–present (v2, daily) | CC BY 4.0 (MITMS) | Tens to hundreds of GB (multi-year daily CSV.gz) | **Easy–Medium.** Note the website may return 403 to plain `curl` (WAF/UA check), while the underlying S3 open-data files and the `spanishoddata` R package provide clean bulk access (note: project already has a sample in `_local/external/mitma_pernoct/`). |
| **INE Experimental Mobility Studies (2019–2021 Mobile Phone)** | https://www.ine.es/experimental/movilidad/experimental_em.htm | ~3,214 *áreas de movilidad* (municipalities/aggregates) OD commuter & seasonal flows (Wednesdays & Sundays, plus Nov 2019 pre-COVID baseline & summer 2020/2021) | Nov 2019 & Mar 2020–Dec 2021 | INE Open Data | ~1–3 GB total (CSV/XLSX) | **Easy.** Static completed archive; useful as an independent mobile-operator benchmark against MITMA v1/v2. |
| **INE Surname & First-Name Frequency by Province** | https://www.ine.es/widgets/nombApell/index.shtml (and INE *Apellidos y nombres más frecuentes* bulk XLSX tables) | Province-of-residence × Province-of-birth surname frequency vectors (top surnames & all surnames ≥5 occurrences per province), enabling isonymy / genetic-cultural distance matrices (Lasker/Nei distance across all 52 provinces!) | Annual Padrón / Census snapshots (latest 2023–2025) | INE Open Data | ~10–50 MB | **Easy–Medium.** Bulk province-level top-surname Excel files are downloadable from INE; full surname distribution for all rare surnames requires querying the widget endpoint or INE custom extract. |
| **Spanish Wikipedia Municipality Hyperlink Network** | https://dumps.wikimedia.org/eswiki/latest/ (`eswiki-latest-pagelinks.sql.gz`, `page.sql.gz`, `linktarget.sql.gz`) + Wikidata SPARQL (`P772` INE code) | Directed hyperlink graph among all 8,131 Spanish municipality articles on `es.wikipedia.org` (mapped 1-to-1 via Wikidata property `P772`) | Monthly snapshots (plus full revision history) | CC BY-SA 4.0 / CC0 (Wikidata) | ~500 MB–1.5 GB (dump) | **Medium.** Mapping `page_id` ↔ `ine_code` via Wikidata `P772` takes one SPARQL query; filtering `pagelinks` gives a clean cultural/geographic salience graph across all 8,131 municipalities. |
| **IGN / CNIG Road & Rail Network (RTN — Red de Transporte)** | https://centrodedescargas.cnig.es/CentroDescargas/ (BTN / IGR-RT *Redes de Transporte*) | Full vector topological graph of Spanish roads (motorways, national, regional), railways, ports, and airports with speed/class attributes | Continuously updated (annual snapshots) | CC BY 4.0 (IGN Spain) | ~1–3 GB (GeoPackage / Shapefile) | **Easy–Medium.** Official national GIS graph; ideal for computing road travel-time matrices between all 8,131 municipal centroids (e.g., via OSRM/Valhalla/NetworkX). |
| **Eurostat / INE Interregional Freight & Road Transport (`road_go_ta_rl`, EPTMC)** | https://www.transportes.gob.es/estadisticas (EPTMC — *Encuesta Permanente del Transporte de Mercancías por Carretera*) | Province-of-origin × province-of-destination road freight tonnage (52 × 52 matrix by commodity group) | 1996–present (quarterly/annual) | Spanish Government Open Data | ~10–50 MB (XLSX/CSV) | **Easy–Medium.** Gives physical economic/supply-chain linkages between Spanish provinces to complement human mobility and migration networks. |
| **C-INTEREG Interregional Trade in Spain** | https://www.c-intereg.es/ | CCAA-to-CCAA (and provincial estimates) monetary & physical trade flows in goods across Spain | 1995–present (annual) | Academic open data (attribution required) | ~5–20 MB (XLSX) | **Easy.** Rare subnational input-output / interregional trade matrix for Spain built from transport & customs records. |
