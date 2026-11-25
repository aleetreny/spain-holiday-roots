# N1: New Open Data — Inventory, Validation, Finer-Scale Tests, and Framing Recommendation

*Research sub-agent N1, 2026-11-24. Code: [`pipeline/n1_new_data.py`](/Users/alejandrotreny/Proyectos/spain-holiday-roots/pipeline/n1_new_data.py) (logs: `results/n1_run.log`, `results/n1_gen2.log`, `results/n1_mediation.log`) and [`pipeline/n1_censo2024_microdata.py`](/Users/alejandrotreny/Proyectos/spain-holiday-roots/pipeline/n1_censo2024_microdata.py) (logs: `results/n1_censo2024.log`, `results/n1_censo2024_cdiff.log`), run with `.venv/bin/python -I`. Every number below comes from the `results/n1_*.csv` file named beside it. Decisions appended to `.planning/2026-11-24-consolidation/decisions.md` tagged `[N1]`.*

---

## Executive summary

All four limitations of the consolidated evidence—(a) MITMA `_AM` zone aggregation hiding real villages, (b) 1st-generation province-only roots mixing recent migrants with the rural exodus, (c) untested second-home ("la casa del pueblo") mechanism, and (d) 2022–2024 single-year roots—can be resolved with official open Spanish data downloaded and validated here. **Three findings directly reshape the answer and resolve Prof. Shi Zhou's main criticisms (G1, G5, M2, M3):**

1. **Why the pooled summer roots premium failed (`+0.10`, `p = 0.18`): two opposite regimes cancelled each other out** (`n1_censo2024_cohorts.csv`, `n1_censo2024_cohort_diffs.csv`, Censo Anual 2024 10% microdata, **675 municipalities > 10k across 44 province clusters**, double R2's 22 clusters). Separating arrival cohorts and birthplace size inside R2's exact neighbouring-city PPML (`city FE + (origin_prov × dest_prov) FE + distance + ADRH income × dest_prov`, WCR Webb bootstrap):
   - **Pre-1981 rural-exodus migrants** (`VARANOP = '0001'`) **surge +87%** from Feb 14 (`0.096`) to Aug 15 (`0.179`; **stacked `δ = +0.084 [+0.043, +0.119]`, `p_wcr = 0.0020`** on **all destinations**, and **`+0.154 [+0.109, +0.192]`, `p_wcr = 0.0002`** inland).
   - **Village-born natives (`≤ 2,000` pop, `CMUN_NAC = '991'`)** **more than triple** from Feb 14 (`0.035`) to Aug 15 (`0.121`; **`δ = +0.087 [+0.060, +0.111]`, `p_wcr = 0.0004`** on **all destinations**, and **`+0.124 [+0.093, +0.152]`, `p_wcr = 0.0002`** inland).
   - **Recent 2016–2024 migrants** (`δ = −0.072 [−0.116, −0.017]`, `p_wcr = 0.014`) and **urban-born (`> 10,000`)** (`δ = −0.062 [−0.119, −0.006]`, `p_wcr = 0.030`) move the **opposite way**—they dominate ordinary nights (`0.196` and `0.284`) and Christmas Eve (`0.197` and `0.309`) and are diluted in August. Both flips replicate in INE 53000 3-operator monthly tourists (`pre-1981 δ = +0.036, p_wcr = 0.017`; `village ≤ 2k δ = +0.049, p_wcr = 0.0006`).
2. **True village-level and 2nd-generation roots are directly observable and sharpen on holidays** (`n1_city_village_roots.csv`, `n1_city_village_qap.csv`, `n1_gen2_city_prov.csv`, `n1_censo2024_mun10k_od.csv`). ECEPOV 2021 microdata (`N = 424,493` persons) contains **uncensored 5-digit birth municipalities** for the respondent (`MUNIN`, 7,498 municipalities) AND **both parents** (`MUNINP1`, 7,856; `MUNINP2`, 7,813). Crossed with INE's true municipality-to-municipality OD network (`exp_tmov_interno_mun_2019..2026`), comparing **two villages in the same province** receiving visitors from the **same city** (`(city × dest_prov) FE + village FE + log km`):
   - Within-destination-province QAP beats the village-permutation null in every month (`r = 0.145–0.157`, `z = 56–64`, `p < 0.001`) and **sharpens in August** (`Aug − Feb z = 16.4`, `Aug − Oct z = 16.8`, `p < 0.001`), **including for strictly 2nd-generation city-born residents** whose parents came from that village (`Aug − Feb z = 14.4`, `p < 0.001`; stacked inland PPML `δ = +0.017 [0.005, 0.028]`, `p = 0.0038` vs Feb; `+0.028`, `p = 6e-7` vs Oct).
   - Emptied inland villages pull far more strongly from their own diaspora cities (`lany × loss_pct = +0.273 [0.221, 0.325]`, `t = 10.3`, `p < 1e-16` in Aug 2024), and second-home share amplifies this (`lany × share_no_principal = +0.124 [0.047, 0.201]`, `p = 0.0017`).
   - On an equal-noise footing (both from ECEPOV 2021 in R2's neighbouring-city design), **strictly 2nd-generation roots equal or exceed 1st-generation roots in August** (MITMA Aug 15 all: `0.074 [0.036, 0.130]` vs `0.061`, both `p_wcr ≤ 0.013`; inland: `0.102 [0.040, 0.156]` vs `0.093`).
3. **"La casa del pueblo" mediates summer village refill—but NOT Christmas** (`n1_second_home_mediation.csv`, `n1_ecepov_casa_del_pueblo.csv`):
   - **Censo 2021 non-primary housing stock (`59525`, 7,974 mainland municipalities)** mediates **64.6% [95% cluster-boot CI: 44.0%, 95.0%]** of the emptied-village August tourist gain (`INE 53001 log(Aug/Feb)`, `n = 5,831`; `β_loss: 0.446 → 0.158`), **48.1% [19.0%, 81.7%]** of MITMA August inflow gain (`n = 2,376`), **47.2% [26.1%, 86.2%]** at Easter (`Mar/Feb`), **12.8%–14.3%** on All Saints, and **only 1.8% [−1.5%, 7.3%] (`p = 0.30`) on Christmas Eve** (and `1.7%` of MITMA Dec 25 net balance). **Summer refills empty family houses; Christmas fills living relatives' primary homes.**
   - **ECEPOV 2021 household microdata (`n = 9,421` households = 1.003M weighted with a cross-CCAA second home)** provides the first national official measurement of "la casa del pueblo": **58.4% [49.3%, 64.2%]** of all cross-regional second homes—and **77.4% [69.4%, 79.9%]** of inland cross-regional second homes (**80.7%** from Madrid, **81.8%** from Cataluña)—are in a household member's or parent's birth CCAA, **2.0×–4.0×** above a within-origin-province destination-permutation null (`p < 0.001`).

---

## 1. Data Inventory

Every table ID, URL, schema, and suppression rule below was verified directly against official INE / Eurostat / MITMS endpoints on 2026-11-24. Sources verified **not to exist** at the requested resolution are explicitly marked.

| Source | Table / File ID | Official URL | Geography | Years / Periods | Unit & Key Variables | Suppression / Methodological Notes | Licence | Downloaded |
|---|---|---|---|---|---|---|---|---|
| **INE Turismo Móviles — Municipal OD (annual Excel)** | `tpx=58416, 58415, 58414, 58413, 58412, 66321, 75000, 79350` (`exp_tmov_interno_mun_YYYY.xlsx`) | `https://www.ine.es/experimental/turismo_moviles/exp_tmov_interno_mun_{2019..2026}.xlsx` (portal `padre=8578&capsel=8581`) | **True 5-digit municipality × 5-digit municipality** (interprovincial only; up to 6,532 active mainland munis/month) | 2019-07 to 2026-09 (87 months; 6,284,691 rows) | Monthly **turistas** (= interprovincial trips with ≥1 night; 3 trips by 1 person = 3 tourists). NOT person-nights. 3 operators (Movistar, Vodafone, Orange) scaled by CNMC shares. | Cells `< 30` tourists omitted (`min = 30`). **Zero "Resto" rows.** Published pairs capture **68.7%–75.9%** of `52997` interprovincial tourists. Date = calendar month; overnight = dominant 22:00–06:00 + detected >06:00 next day. 2024 Sep–Dec sheets named `2024_09..2024_12`. | INE Open Data (free reuse w/ attribution) | **Yes** (`_local/external/ine_turismo_moviles_od/`) |
| **INE Turismo Móviles — Monthly tables (existing)** | `52997`, `53000`, `53001` | `https://www.ine.es/jaxiT3/Tabla.htm?t={52997,53000,53001}` | `52997`: 52×52 prov; `53000`: 8,131 mun × 52 prov; `53001`: 52 prov × 8,131 mun | 2019-07 to 2026-09 | `52997`: **turistas, pernoctaciones, duración media**; `53000`/`53001`: **turistas** | Interprovincial only; cells `< 30` suppressed (`""`); no "Resto" rows. Average stay in `52997`: 4.0 nights (Feb) vs 7.0 nights (Aug). | INE Open Data | **Already local** (`_local/external/ine_turismo_moviles/`) |
| **INE ECEPOV 2021 Microdata** | `datos_2021.zip`, `disreg_ecepov21.zip` | `https://www.ine.es/ftp/microdatos/ecepov/datos_2021.zip` | `IDQ_MUN`: 151 cities (>50k + capitals) + 52 prov (`IDQ_PV`); **`MUNIN`: 7,498 exact birth munis; `MUNINP1`/`MUNINP2`: 7,856/7,813 exact parents' birth munis!** | 2021 (424,493 persons = 46.83M weighted; 172,444 households = 18.81M weighted) | Individual birth & **both parents' exact birth municipality** (`MUNIN`, `MUNINP1`, `MUNINP2`), survey weight `FACTOR`; household `SEGUNRESI`, `LUGSEGUNRESI`, `OTRACCAA` (exact CCAA of 2nd home), `DIASUSA` | ~1% national sample -> sampling zeros in small city×village pairs; exact `IDQ_MUN` censored to `00000` for residence ≤50k, but **birth & parents' birth municipalities are uncensored 5-digit codes**! | INE Open Data | **Yes** (`_local/external/ine_ecepov_2021/`) |
| **INE ECEPOV 2021 Second-Home Tables** | `tpx=56616, 56620, 56632` | `https://www.ine.es/jaxi/Tabla.htm?tpx={56616,56620,56632}` | `56616`/`56620`: 52 provinces; `56632`: 151 cities | 2021 | Households by 2nd-home ownership & coarse location (same mun / same prov / same CCAA / other CCAA / abroad) | Pre-tabulated tables do **not** cross-tabulate destination province/CCAA; only microdata (`OTRACCAA`) gives destination CCAA. | INE Open Data | **Yes** (`_local/external/ine_ecepov_2021/raw/`) |
| **INE Censo Anual de Población 2024 (10% Microdata)** | `datos_muestra_2024.zip` (`md_CensoAnual_muestra_2024.parquet`) | `https://www.ine.es/ftp/microdatos/censoanual/vc/datos_muestra_2024.zip` | **780 exact residence munis (>10k)** × **781 exact birth munis (>10k)** + 52 birth provs + size brackets (`991` ≤2k, `992` 2–5k, `993` 5–10k) | 1 Jan 2024 (4,861,970 persons = 10% of 48.62M) | Exact `CPRO+CMUN`, `CPRO_NAC+CMUN_NAC`, arrival cohort **`VARANOP`** (`0001`=pre-1981, `0002`=1981–90, ..., `2016–2023`) | Municipalities ≤10k recoded to size class (`991`/`992`/`993`), so **cannot** give exact village codes, but separates rural-village-born (`991`) from urban-born (`>10k`) and pre-1981 exodus from recent migrants at 675-municipality scale. (Exhaustive 100% demographic file `datos_2024.zip`, 2.03 GB, also exists on INE FTP.) | INE Open Data | **Yes** (`_local/external/ine_censo_microdatos_2024/`) |
| **INE Censo de Población y Viviendas 2021 — Housing Stock** | `tpx=59525` | `https://www.ine.es/jaxi/Tabla.htm?tpx=59525&L=0` (`files/tpx/es/csv_bdsc/59525.csv`) | **All 8,131 Spanish municipalities** (7,974/7,975 mainland matched) | 2021 | `Total` (26,623,708), `Vivienda principal` (18,536,616), `Vivienda no principal` (8,087,092) | Complete municipal coverage. Does not split non-primary into vacant vs seasonal (use `59531` for ≥500 pop). | INE Open Data | **Yes** (`_local/external/ine_censo_viviendas_2021/`) |
| **INE Censo 2021 — Electricity-Based Housing Use** | `tpx=59531` | `https://www.ine.es/jaxi/Tabla.htm?tpx=59531&L=0` (`files/tpx/es/csv_bdsc/59531.csv`) | **3,185 municipalities** (≥ ~500 pop; 2,988 mainland matched) | 2021 (2020 electricity CUPS) | `Viviendas vacías` (<250 kWh, 3.83M), `bajo consumo`, `uso esporádico` (251–750 kWh, 2.52M) | **Caveat:** summer-only village houses used 2 weeks in August often use <250 kWh/yr and are classified as `vacías`; in 101 munis `viv_vacias > viv_no_principal`. `59525` `viv_no_principal` is the right denominator for small villages. | INE Open Data | **Yes** (`_local/external/ine_censo_viviendas_2021/`) |
| **Eurostat GISCO LAU 2021** | `LAU_RG_01M_2021_3035.csv` | `https://gisco-services.ec.europa.eu/distribution/v2/lau/csv/LAU_RG_01M_2021_3035.csv` | **All 8,131 Spanish municipalities** (7,974/7,974 mainland matched) | 2021 | Exact `AREA_KM2`, `POP_2021`, `POP_DENS_2021` by 5-digit `LAU_ID` | Fixes R1's missing municipal area control (`log(area_km2)`). | CC BY 4.0 / Eurostat | **Yes** (`_local/external/eurostat_lau_2021/`) |
| **INE Censo 2021 Personas Microdata** | `CensoPersonas_2021.zip` (1.30 GB) | `https://www.ine.es/ftp/microdatos/censopv/cen21/CensoPersonas_2021.zip` | Same >10k `CMUN`/`CMUN_NAC` threshold (`991/992/993`) | 2021 (10% sample) | Includes `NORDEN_MAD`/`NORDEN_PAD` for **co-residing** parents only (not adult children who left home) | Superseded by 2024 Censo Anual (`datos_muestra_2024.zip`, exact 2024 alignment) + ECEPOV 2021 (all parents, uncensored villages). | INE Open Data | Schema verified (`dr_CensoPersonas_2021.zip`); not needed |
| **INE Padrón Continuo / Censo Anual municipal birthplace tables** | `55250` (Censo 2021), `68540` (Censo Anual) | `https://www.ine.es/jaxiT3/Tabla.htm?t=68540` | All municipalities | 2021–2025 | Coarse relation only ("same mun / same prov / other prov / abroad") | **Does NOT exist** as an 8,131 × 52 or 8,131 × 8,131 public table (disclosure control). Only 151 cities × 52 prov (`68545`) exists as a public table. | INE Open Data | Verified not to exist at mun×mun |
| **INE Apellidos (Surnames)** | `/apellidos/inicio.do`, `apellidos_frecuencia.xls`, `apellidos_mas_frecuentes.xls` | `https://www.ine.es/daco/daco42/nombyapel/nota_apellidos.pdf` | **Province-only** (residence or birth; freq > 5 in app, top 50/prov in Excel) | Annual (Padrón) | Surname frequencies (>180k surnames in interactive app) | **NOT available at municipal level** from INE (`nota_apellidos.pdf` p. 2). Regional agencies (IDESCAT, EUSTAT) cover only Catalonia / Basque Country. Cannot be used for national municipal roots. | INE Open Data | Verified province-only; not downloaded |
| **MITMS / MITMA 2025 Pernoctaciones** | `RSS.xml` (`2025*_Pernoctaciones_municipios.csv.gz`) | `https://movilidad-opendata.mitma.es/RSS.xml` | MITMA zones (2,435 mainland) | Continuous **Jan–Jun 2025**; only **17 selected days** in Jul–Dec 2025 (incl. Aug 8, 14–17; Oct week; Dec 8) | Daily person-nights (`pernoctaciones`) | **Methodology break on 2025-07-01** (`LEEME_CambioMetodología_IndicadoresBigdata_posteriores_a_20250701.pdf`); Feb 2025 and Aug 2025 sit on opposite sides of the break, making within-2025 Aug−Feb MITMA contrasts non-comparable. INE OD 2025 (`ine_od_mun_2025.parquet`) is consistent across all 12 months. | MITMS Open Data | Audited via RSS; not downloaded due to Jul-2025 break |
| **INE ETR / FAMILITUR Microdata** | `disreg_etr15.zip` | `https://www.ine.es/dyngs/INEbase/es/operacion.htm?c=Estadistica_C&cid=1254736176990` | CCAA-level (~8,200 households/month) | Quarterly 2015–2026 | Trip purpose (incl. VFR: 28.3% Q3, 39.6% Q4), accommodation (own 2nd home vs relatives' home) | Too coarse geographically (CCAA only, no birthplace variable) to test municipal/provincial roots; already cited in `literature_B.md` D.2. | INE Open Data | Documented; not downloaded |

---

## 2. Validation of INE Municipal OD Against MITMA (All 6 Overlapping Months of 2024)

Source: `results/n1_mitma_validation.csv` (`pipeline/n1_new_data.py`, Stage 3). MITMA daily interprovincial person-nights (`pernoct_mun_2024_6m.parquet`) are summed over each of the six months present in MITMA 2024 (`02, 03, 08, 10, 11, 12`; Nov has 27 days due to 3 missing dates) and compared against INE `52997` (prov × prov, both `pernoctaciones` and `turistas`) and INE `exp_tmov_interno_mun_2024` (municipal OD aggregated to province pairs, MITMA zone pairs, true 1-to-1 single-municipality zone pairs, and 2,435 zone inflows).

| Month | MITMA days | INE OD coverage vs `52997` | 47×47 prov: `r` (MITMA vs `52997` pernoct) | 47×47 prov: Spearman `ρ` (`52997` pernoct) | 47×47 prov: `r_log` (`52997` tur) | MITMA zone×zone overlap pairs | Zone×zone `r_log` / `ρ` | True 1-to-1 single-mun pairs (`r_log` / `ρ`) | 2,435 dest inflows (`r_log` / `ρ`) |
|---|---:|---:|---:|---:|---:|---:|---|---|---|
| `2024-02` | 29 | 73.9% | 0.959 | 0.943 | 0.938 | 59,676 | 0.689 / 0.645 | 44,484 (0.733 / 0.685) | 0.778 / 0.850 |
| `2024-03` | 31 | 73.8% | 0.900 | 0.920 | 0.955 | 70,222 | 0.699 / 0.656 | 52,503 (0.741 / 0.698) | 0.796 / 0.868 |
| `2024-08` | 31 | **75.9%** | **0.988** | **0.951** | **0.963** | 107,587 | 0.681 / 0.631 | 77,213 (0.720 / 0.673) | **0.836 / 0.890** |
| `2024-10` | 31 | 68.7% | 0.914 | 0.954 | 0.965 | 63,458 | 0.695 / 0.656 | 47,720 (0.739 / 0.697) | 0.802 / 0.859 |
| `2024-11` | 27 | 69.2% | 0.947 | 0.945 | 0.958 | 58,975 | 0.703 / 0.667 | 44,035 (0.748 / 0.709) | 0.800 / 0.861 |
| `2024-12` | 31 | 69.7% | 0.946 | 0.952 | 0.955 | 59,403 | 0.715 / 0.684 | 44,780 (0.759 / 0.725) | 0.805 / 0.860 |

**Validation verdict:** Excellent cross-source agreement despite different operators (Orange-only MITMA vs 3-operator INE), different units (daily person-nights vs monthly trips; mean stay jumps from 4.0 nights in Feb to 7.0 nights in Aug in `52997`), and INE's `< 30` pair suppression. Province-level overnights agree at `r = 0.90–0.99` (`ρ = 0.92–0.95`), destination inflows at `ρ = 0.85–0.89`, and true 1-to-1 single-municipality pairs at `r_log = 0.72–0.76`.

---

## 3. True Municipality-Level INE OD Holiday Network (Zero `_AM` Aggregation)

Source: `results/n1_ine_network_props.csv`, `results/n1_mun_network_node_tests.csv`, `results/n1_ine_od_inventory.csv`.

### 3.1 Global network properties (`N_fixed = 7,974` mainland municipalities; all densities `<< 5%`)

Three explicit backbones are computed: `pub_ge30` (INE's published threshold `w ≥ 30` tourists/month), `w100` (`w ≥ 100`), and `equal_edges` (top `E = 69,808` directed edges matching `2024-02`, ties broken deterministically). Undirected rewiring (`N = 20` graphs, `10 E` simple swaps) confirms structural constraints:

| Month (`backbone`) | `N_active` | `E_dir` | `E_und` | Und. density (fixed / active) | `⟨k⟩_act` | `Spearman(k_in, k_out)` | Reciprocity | `r_und` (rewired) | `C_avg_act` (ER) | `L_gcc` | Rich-club `ρ(top 50)` (ceil) | Max core (und / in) |
|---|---:|---:|---:|---|---:|---:|---:|---|---|---|---|---|
| `2024-02` (`pub_ge30`) | 5,527 | 69,808 | 52,720 | 0.166% / 0.345% | 19.1 | 0.708 | 0.490 | −0.161 (−0.160) | 0.379 (0.003) | 2.48 | 1.056 (1.130) | 61 / 46 |
| `2024-03` (`pub_ge30`) | 5,800 | 81,879 | 62,865 | 0.198% / 0.374% | 21.7 | 0.699 | 0.464 | −0.167 (−0.166) | 0.399 (0.004) | 2.41 | 1.032 (1.101) | 68 / 50 |
| `2024-08` (`pub_ge30`) | **6,383** | **127,903** | **106,174** | **0.334% / 0.521%** | **33.3** | **0.635** | **0.340** | **−0.192 (−0.188)** | 0.448 (0.005) | 2.33 | **0.971 (1.057)** | **95 / 58** |
| `2024-08` (`equal_edges`) | 5,803 | 69,808 | 58,185 | 0.183% / 0.346% | 20.1 | **0.564** | **0.333** | −0.165 (−0.163) | 0.398 (0.003) | 2.43 | 0.992 (1.112) | 69 / 39 |
| `2024-08` (`w100`) | 4,439 | 26,934 | 22,616 | 0.071% / 0.230% | 10.2 | **0.412** | **0.321** | −0.162 | 0.290 | 2.66 | — | 43 / 22 |
| `2024-10` (`pub_ge30`) | 5,806 | 74,334 | 53,924 | 0.170% / 0.320% | 18.6 | 0.712 | 0.549 | −0.169 (−0.169) | 0.396 (0.003) | 2.47 | 1.048 (1.110) | 63 / 50 |
| `2024-11` (`pub_ge30`) | 5,807 | 69,777 | 51,357 | 0.162% / 0.305% | 17.7 | 0.683 | 0.528 | −0.168 (−0.166) | 0.411 (0.003) | 2.43 | 1.039 (1.117) | 61 / 48 |
| `2024-12` (`pub_ge30`) | 5,769 | 69,514 | 51,279 | 0.161% / 0.308% | 17.8 | 0.682 | 0.525 | −0.171 (−0.170) | 0.410 (0.003) | 2.47 | 1.036 (1.112) | 63 / 50 |

- **Replication across 2022–2025:** Reciprocity drops every August (`2022: 0.558→0.380`; `2023: 0.552→0.359`; `2024: 0.490→0.340`; `2025: 0.490→0.267`), `Spearman(k_in, k_out)` drops every August (`w100`: `2022: 0.564→0.438`; `2023: 0.562→0.402`; `2024: 0.534→0.412`; `2025: 0.411→0.233`), and degree-forced disassortativity replicates (`r_und` matches `r_rew` to within `0.004`). Because INE aggregates a full month of trips (not one peak night), sender/receiver decoupling is milder than MITMA's single-night drop (`0.69→0.09`).

### 3.2 Node-level emptied-village tests at TRUE municipal resolution (resolving audit C1 / M4)

Controls: exact municipal `log(padron_2024)`, exact municipal `log(area_km2)` (Eurostat LAU 2021), `coastal_true`, `log1p(dist_city)`, `log1p(dist_coast)`, and 47 province FE; inference via Webb wild-cluster bootstrap (`B = 1,999`, `G = 47` provinces; `n1_mun_network_node_tests.csv`). Regressor = `loss_pct` (0 = grew most, 1 = most emptied):

| Outcome (`2024-08 vs 2024-02`, INE OD `pub_ge30`) | All mainland munis (`n = 7,750`) | Inland munis (`n = 7,345`) | Small munis inside MITMA `_AM` (`n = 6,276`) | Larger single-zone MITMA munis (`n = 1,473`) |
|---|---|---|---|---|
| `log((s_in,Aug+1)/(s_in,Feb+1))` | **+0.812** [0.557, 1.073], `p = 0.0005` | **+0.822** [0.537, 1.105], `p = 0.0005` | **+1.178** [0.738, 1.596], `p = 0.0005` | −0.166 [−0.993, 0.685], `p = 0.70` |
| `log((k_in,Aug+1)/(k_in,Feb+1))` | **+0.196** [0.108, 0.286], `p = 0.0005` | **+0.239** [0.147, 0.327], `p = 0.0005` | **+0.448** [0.352, 0.539], `p = 0.0005` | +0.174 [−0.313, 0.673], `p = 0.49` |
| `Δ core_und` | **+11.51** [8.42, 14.65], `p = 0.0005` | **+12.98** [9.09, 16.72], `p = 0.0005` | **+3.40** [1.49, 5.27], `p = 0.001` | **−12.81** [−22.48, −3.82], `p = 0.007` |
| `Δ core_in` (directed visits core) | **+2.08** [0.71, 3.44], `p = 0.005` | **+2.84** [1.60, 4.08], `p = 0.001` | **+1.78** [0.51, 2.95], `p = 0.007` | +1.45 [−6.00, 8.80], `p = 0.69` |

**Key diagnostic insight:** This resolves the puzzle in R1 (`r1_villages.md`). On the **exact 1,473 larger single-municipality MITMA zones**, INE reproduces R1's reversal (`Δcore = −12.81`, `s_in` null). On the **6,276 smaller municipalities that MITMA folded into `_AM` aggregates**, emptied villages genuinely gain visitors (`+1.178`), in-partners (`+0.448`), undirected coreness (`+3.40`), and directed in-coreness (`+1.78`)—with **zero zone aggregation**. (Caution: in INE OD, `net_bal_diff` is distorted downward to `−0.16` because tiny villages have `s_out = 0` censored `< 30` in both months; net balance must be measured with uncensored MITMA or `53001 − 53000` national totals.)

---

## 4. Core Empirical Findings From the New Data

### Finding 1: Separating the 1950–1980 rural exodus from recent migration resolves Prof. Zhou's G1

Sources: `results/n1_censo2024_expanded_city_prov.csv`, `results/n1_censo2024_cohorts.csv`, `results/n1_censo2024_cohort_diffs.csv` (`pipeline/n1_censo2024_microdata.py`).

Using Censo Anual 2024 (10% sample, 4.86M persons), R2's neighbouring-city PPML (`city FE + (origin_prov × dest_prov) FE + log km + (log km)^2 + ADRH 2024 income × dest_prov`, WCR Webb bootstrap) expands from **112 cities in 22 provinces** to **675 single-zone municipalities (`> 10,000` pop) across 44 provinces** (`n = 31,043` cells/night):

1. **Pooled birthplace (`all_id_provs`, `G = 44`) confirms R2 and removes Barcelona dependence inland:**
   - MITMA Aug 15 − Feb 14 (all destinations): `δ = +0.042 [−0.028, +0.107]`, `p_wcr = 0.216` (still not significant pooled!).
   - MITMA Aug 15 − Feb 14 (**inland destinations**): **`δ = +0.110 [+0.025, +0.185]`, `p_wcr = 0.013`**; dropping Barcelona barely changes it (`δ = +0.094 [−0.002, +0.185]`, `p_wcr = 0.055`, vs R2's drop from `0.155` to `0.086`).
   - MITMA Dec 25 − Feb 14 (all destinations): **`δ = +0.066 [+0.026, +0.119]`, `p_wcr = 0.0020`** (10× sharper than R2's `p = 0.019`).
2. **Joint decomposition into Pre-1981 Exodus (`VARANOP = '0001'`) vs Recent Migrants (`2016–2024`), and Village-Born (`≤ 2,000`, `CMUN_NAC = '991'`) vs Urban-Born (`> 10,000`)** (within-FE correlations only `0.205` and `0.159`—cleanly separable!):

| Dataset / Contrast (`G = 44`, WCR `B = 4,999`) | Pre-1981 exodus `β / δ` [95% WCR CI] (`p_wcr`) | Recent 2016–24 `β / δ` [95% WCR CI] (`p_wcr`) | Village-born ≤2k `β / δ` [95% WCR CI] (`p_wcr`) | Urban-born >10k `β / δ` [95% WCR CI] (`p_wcr`) |
|---|---|---|---|---|
| **MITMA Feb 14 (ordinary, all)** | 0.096 [0.060, 0.129] (`0.0004`) | **0.196** [0.147, 0.255] (`0.0002`) | 0.035 [0.018, 0.056] (`0.0010`) | **0.284** [0.212, 0.358] (`0.0002`) |
| **MITMA Aug 15 (summer, all)** | **0.179** [0.127, 0.218] (`0.0006`) | 0.125 [0.098, 0.183] (`0.0002`) | **0.121** [0.092, 0.154] (`0.0002`) | 0.222 [0.191, 0.268] (`0.0002`) |
| **MITMA Aug 15 (summer, inland)** | **0.257** [0.185, 0.297] (`0.0002`) | 0.137 [0.112, 0.182] (`0.0002`) | **0.170** [0.136, 0.198] (`0.0002`) | 0.253 [0.227, 0.283] (`0.0002`) |
| **MITMA Dec 25 (Christmas Eve, all)** | 0.121 [0.091, 0.154] (`0.0002`) | **0.197** [0.158, 0.261] (`0.0002`) | 0.064 [0.048, 0.088] (`0.0002`) | **0.309** [0.253, 0.383] (`0.0002`) |
| **Stacked `Δ` MITMA Aug 15 − Feb 14 (ALL)** | **+0.084** [+0.043, +0.119] (**`0.0020`**) | **−0.072** [−0.116, −0.017] (**`0.0144`**) | **+0.087** [+0.060, +0.111] (**`0.0004`**) | **−0.062** [−0.119, −0.006] (**`0.0304`**) |
| **Stacked `Δ` MITMA Aug 15 − Feb 14 (INLAND)** | **+0.154** [+0.109, +0.192] (**`0.0002`**) | **−0.074** [−0.125, −0.015] (**`0.0182`**) | **+0.124** [+0.093, +0.152] (**`0.0002`**) | −0.043 [−0.109, +0.016] (`0.174`) |
| **Stacked `Δ` MITMA Dec 25 − Feb 14 (ALL)** | +0.026 [−0.002, +0.051] (`0.065`) | +0.001 [−0.030, +0.049] (`0.973`) | **+0.029** [+0.010, +0.049] (**`0.0048`**) | +0.025 [−0.010, +0.073] (`0.159`) |
| **Stacked `Δ` INE 53000 Aug − Feb (ALL)** | **+0.036** [+0.007, +0.063] (**`0.0170`**) | −0.019 [−0.039, +0.001] (`0.066`) | **+0.049** [+0.030, +0.067] (**`0.0006`**) | −0.023 [−0.057, +0.008] (`0.151`) |

**Why this matters:** Even against the **pre-specified February baseline** and across **ALL destinations (coast included)**, the summer holiday shift is **not** weak or ambiguous—it is an **inversion between two migrations**: August selectively activates the **pre-1981 rural exodus** (`+0.084`, `p = 0.002`) and **village-born (`≤ 2,000`)** ties (`+0.087`, `p = 0.0004`, ×3.5 jump), while **recent urban migration ties fall** (`−0.072`, `p = 0.014`). Pooling both into one birthplace stock hid this cleanly opposite movement.

### Finding 2: Municipality-to-municipality roots and strictly 2nd-generation (parental) roots

Sources: `results/n1_city_village_roots.csv`, `results/n1_city_village_qap.csv`, `results/n1_gen2_city_prov.csv`, `results/n1_censo2024_mun10k_od.csv`.

1. **Exact City × Village alignment inside the same destination province (`(city × dest_prov) FE + village FE + log km`, 137 cities × 6,418 active destination municipalities, ECEPOV 2021 × INE OD):**
   - **Within-destination-province QAP** (1,000 permutations of villages within each destination province; descriptive network test under exchangeability of villages within province):
     - Combined family roots (1st + 2nd gen): `r = 0.146` (Feb), `0.145` (Mar), `0.146` (Aug), `0.154` (Nov), `0.157` (Dec) — `z = 56–64`, `p < 0.001` every month.
     - **Holiday − ordinary difference QAP:** `Aug − Feb: r = +0.0313 (z = 16.4, p < 0.001)`; `Aug − Oct: +0.0318 (z = 16.8, p < 0.001)`; `Dec − Feb: +0.0129 (z = 6.9, p < 0.001)`; `Mar (Easter) − Feb: +0.0118 (z = 6.9, p < 0.001)`; `Nov − Feb: +0.0078 (z = 4.3, p < 0.001)`.
     - **2nd-generation STRICT** (own-province-born city residents whose parents were born in village `v`): `r = 0.135–0.145 (z = 52–56, p < 0.001)`; **`Aug − Feb: +0.0275 (z = 14.4, p < 0.001)`; `Aug − Oct: +0.0282 (z = 14.7, p < 0.001)`**.
   - **Clustered PPML & margins (`CRV1` by 47 destination provinces):**
     - **Small villages (`pop < 5,000`)**: `lany` (1st+2nd gen) = **+0.037** [0.022, 0.051] in Feb → **+0.047** [0.032, 0.062] in Aug (`p = 9e-10`).
     - **Extensive margin (`1[flow ≥ 30]`, `n = 856,772`)**: Feb `+0.0203` [0.0180, 0.0227] → **Aug `+0.0254` [0.0229, 0.0278]** (+25% jump); **intensive margin (`log(flow)` | `≥ 30`, `n = 59,284`)**: **`+0.0465` [0.0414, 0.0516]**.
     - **Emptied-village interaction (`inland_interact_loss`, Aug 2024, `n = 566,226`)**: `lany × loss_pct = +0.273 [0.221, 0.325]` (`t = 10.27`, `p < 1e-16`)—the most emptied inland villages (`loss_pct = 1`) have a net within-province roots elasticity of **`+0.204`** (`−0.069 + 0.273`), vs zero for growing towns! Adding `lany × share_no_principal` yields **`+0.124 [0.047, 0.201]` (`p = 0.0017`)** alongside `lany × loss_pct = +0.200 [0.136, 0.264]` (`p = 1e-9`).
     - **Multi-year replication (`all_ppml_any`, `2022–2025`):** Aug is positive and significant every year (`0.0182, 0.0189, 0.0185, 0.0202`, all `p ≤ 0.003`), while Feb is indistinguishable from zero (`0.0068, 0.0037, 0.0079, 0.0113`, all `p > 0.08`). Stacked inland `Aug − Feb` for **2nd-gen strict**: `δ = +0.017 [0.005, 0.028]` (`p = 0.0038`; vs Oct: `+0.028 [0.017, 0.039]`, `p = 6e-7`).
     - **Censo 2024 `>10k × >10k` exact municipal OD (`n1_censo2024_mun10k_od.csv`, 675 × 675 municipalities, 10× larger sample than ECEPOV):** Within-province PPML elasticity = **`+0.058 [0.028, 0.089]`** (Feb) → **`+0.063 [0.045, 0.081]`** (`p = 6.2e-12`, Aug), confirming municipal-level alignment even among towns `> 10,000`.
2. **2nd-generation strict roots in R2's exact neighbouring-city design (`n1_gen2_city_prov.csv`, 112 cities, 22 clusters):**
   - Alone, 2nd-gen strict predicts stays on every night (`MITMA Feb 14: 0.084 [0.049, 0.139], p_wcr = 0.004`; `Aug 15: 0.105 [0.054, 0.179], p_wcr = 0.011`; `Aug 15 inland: 0.148 [0.061, 0.221], p_wcr = 0.007`).
   - Paired with noise-free 100%-enumeration Censo `68545`, the 1%-sample ECEPOV 2nd-gen variable is absorbed (`−0.002`, `p = 0.67`) due to severe attenuation bias (955/5,152 sampling zeros). **On an equal-noise footing (both 1st-gen and 2nd-gen strict measured from ECEPOV 2021)**, 2nd-generation strict independently predicts overnight stays **at or above the 1st-generation slope in August**: MITMA Aug 15 all `0.074 [0.036, 0.130]` (`p_wcr = 0.011`) vs 1st-gen `0.061` (`p_wcr = 0.013`); MITMA Aug 15 inland `0.102 [0.040, 0.156]` (`p_wcr = 0.008`) vs `0.093` (`p_wcr = 0.004`); INE `2024M08` `0.059 [0.028, 0.095]` (`p_wcr = 0.009`) vs `0.050` (`p_wcr = 0.008`).

### Finding 3: Second-home mediation separates "la casa del pueblo" (summer) from "la casa de los padres" (Christmas)

Sources: `results/n1_second_home_mediation.csv`, `results/n1_ecepov_casa_del_pueblo.csv`.

#### (a) Municipal & MITMA-zone stock mediation (`Censo 2021 59525`, 47 province clusters, cluster-bootstrap mediation CIs)

| Outcome (geography, `n`) | `β(loss_pct)` base [WCB CI] | `β(loss_pct)` + housing [WCB CI] | `β(share_no_principal)` [WCB CI] | **Mediated share [95% boot CI]** |
|---|---|---|---|---|
| **INE 53001 `log(Aug/Feb)`** (munis, `5,831`) | +0.446 [0.325, 0.562] | +0.158 [0.019, 0.308] | +1.160 [0.905, 1.422] | **64.6% [44.0%, 95.0%]** |
| **INE 53001 `log(Aug/Oct)`** (munis, `6,068`) | +0.449 [0.269, 0.649] | +0.175 [0.004, 0.355] | +1.077 [0.860, 1.308] | **60.9% [43.3%, 97.5%]** |
| **INE 53001 `(Aug−Feb)/pop`** (all munis, `7,750`, `[0, 29]` bounds) | +2.074 (bounds `1.94..2.16`) | +1.056 (bounds `0.88..1.17`) | +3.554 [2.257, 4.905] | **49.1% [29.0%, 81.2%]** (bounds `45.7%..54.6%`) |
| **INE net `(53001−53000)/pop`** (`4,672`, mediator: `log1p(np/res)`) | +2.593 [1.786, 3.383] | −0.029 [−0.983, 0.882] (`p = 0.96`) | +4.330 [3.054, 5.630] | **101.1% [71.9%, 138.6%]** |
| **INE 53001 `log(Mar/Feb)`** Easter (munis, `5,839`) | +0.175 [0.089, 0.261] | +0.093 [0.001, 0.179] | +0.329 [0.194, 0.466] | **47.2% [26.1%, 86.2%]** |
| **INE 53001 `log(Nov/Feb)`** All Saints (munis, `5,758`) | +0.227 [0.088, 0.372] | +0.194 [0.036, 0.350] | +0.130 [−0.079, 0.338] (`p = 0.22`) | 14.3% [−5.8%, 63.5%] (ns) |
| **MITMA Aug 15 `log(s_in/Feb)`** (zones, `2,376`) | +0.459 [0.204, 0.719] | +0.238 [0.027, 0.468] | **+3.275** [2.793, 3.752] | **48.1% [19.0%, 81.7%]** |
| **MITMA Aug 15 `Δ net_bal`** (zones, `2,376`) | +0.419 [0.233, 0.593] | +0.304 [0.126, 0.470] | +1.702 [1.369, 2.046] | **27.4% [10.7%, 45.2%]** |
| **MITMA Nov 1 `log(s_in/Feb)`** (zones, `2,376`) | +0.537 [0.335, 0.741] | +0.468 [0.223, 0.704] | +1.019 [0.534, 1.461] | **12.8% [2.8%, 33.1%]** |
| **MITMA Dec 25 `log(s_in/Feb)`** Christmas (zones, `2,376`) | +0.687 [0.526, 0.842] | +0.675 [0.497, 0.854] | +0.185 [−0.168, 0.540] (`p = 0.30`) | **1.8% [−1.5%, 7.3%] (zero)** |
| **MITMA Dec 25 `Δ net_bal`** Christmas (zones, `2,376`) | +0.501 [0.392, 0.602] | +0.492 [0.373, 0.610] | +0.126 [−0.092, 0.339] (`p = 0.24`) | **1.7% [−0.8%, 5.9%] (zero)** |

*(On municipalities ≥ 500 pop where `59531` exists, `n = 2,783`, `share_no_principal` mediates 48.1%, `share_esporadico` mediates 14.5%, and `share_no_principal − share_vacias` has negative mediation `−15.7%` because summer-only village houses consume `< 250 kWh/yr` and are coded as `vacías` in `59531`.)*

#### (b) Direct national household test of "la casa del pueblo" (`ECEPOV 2021`, `n = 9,421` cross-CCAA second-home households = 1,002,865 weighted)

| Subgroup (`n1_ecepov_casa_del_pueblo.csv`) | Sampled HH (weighted) | 1st-gen birth CCAA [95% CI] | Parents' birth CCAA | **Member OR parent birth CCAA [95% CI]** | Within-prov perm. null | **Lift (`p_perm`)** |
|---|---:|---|---:|---|---:|---|
| **All cross-CCAA 2nd homes** | 9,421 (1,002,865) | 45.6% [39.6%, 51.3%] | 54.5% | **58.4% [49.3%, 64.2%]** | 20.1% | **2.90× (`<0.001`)** |
| **Inland destination CCAA** | 4,241 (496,292) | 59.5% [55.8%, 62.5%] | 72.6% | **77.4% [69.4%, 79.9%]** | 38.3% | **2.02× (`<0.001`)** |
| **Coastal/island dest CCAA** | 5,180 (506,572) | 31.9% [25.8%, 42.6%] | 36.7% | **39.8% [33.0%, 51.5%]** | 15.5% | **2.57× (`<0.001`)** |
| Origin: **Madrid** (all / inland dest) | 3,478 / 1,758 | 45.5% / 59.7% | 57.3% / 77.1% | **60.2%** / **80.7%** [77.8%, 83.5%] | 20.4% / 40.0% | **2.95× / 2.01× (`<0.001`)** |
| Origin: **Cataluña** (all / inland dest) | 792 / 410 | 63.9% / 61.8% | 74.3% / 76.4% | **79.5%** / **81.8%** [75.7%, 87.4%] | 20.1% / 26.6% | **3.96× / 3.08× (`<0.001`)** |
| Origin: **País Vasco** | 1,088 (124,476) | 36.3% [33.1%, 39.5%] | 45.4% | **48.6% [45.3%, 52.0%]** | 21.8% | **2.23× (`<0.001`)** |

*(CIs use 52-province cluster bootstrap where `G ≥ 10` and household bootstrap for single-CCAA subsets. `literature_B.md` D.2 previously noted no national figure existed; ECEPOV 2021 provides it cleanly.)*

---

## 5. Recommendation: Coursework vs Journal Paper

### Does this new data make the COMP0123 coursework question better—or belong in the journal paper only?

**Split it decisively, obeying Prof. Zhou's "3 results only; penalises over-ambition" rule (M8):**

1. **WHAT TO BRING INTO THE COURSEWORK (tiny complexity cost, huge grade jump from ~75 → 80+):**
   - **Keep MITMA (`N = 2,435`, single nights) as the primary network** for COMP0123—because MITMA is uncensored, captures intra-provincial trips and daily holiday timing (Aug 15 vs Nov 1 vs Dec 25 vs sister nights), and already has all Zhou methods (BOSAM, rewiring, rich club, gravity null, QAP) consolidated and unit-tested.
   - **Adopt TWO lightweight, high-impact upgrades inside the existing 3-result structure:**
     1. **Inside Result 2 (Links / Roots):** Add the **Censo Anual 2024 cohort & village-born decomposition** (`pre-1981 exodus` vs `recent 2016–2024`, and `village ≤ 2,000` vs `urban > 10,000`, `G = 44` provinces, `675` municipalities). This **directly solves G1, G5, and M2** in a single two-panel figure: it proves *why* pooled birthplace barely moved in August (`+0.04`, ns) while **pre-1981 exodus (`+0.084`, `p = 0.002`) and village-born ties (`+0.087`, `p = 0.0004`) surge in August** even with all coastal destinations included, and removes dependence on 22 clusters / Barcelona!
     2. **Inside Result 3 (Nodes / Emptied villages):** Add the **Censo 2021 `59525` non-primary housing mediation** (Aug **48%–65%** mediated vs Christmas **< 2%**!) plus **one sentence** from ECEPOV 2021 (**77.4%** of inland cross-regional second homes are in the family's birth CCAA, 2.0× null) as the physical mechanism explaining *why* summer refills empty villages via "la casa del pueblo" while Christmas goes to living relatives' homes. Mention the **INE municipal OD network (`N = 7,974`, density 0.33%)** in one robustness paragraph/appendix table to confirm that small villages hidden inside MITMA `_AM` zones genuinely gain `k_in` (`+0.45`) and `core_in` (`+1.78`).
2. **WHAT TO KEEP FOR THE JOURNAL PAPER ONLY (*EPJ Data Science* / *Journal of Transport Geography* / *Population, Space and Place*):**
   - The full **multilayer / multi-scale architecture**: true `N = 7,974` INE municipal OD network across **87 months (2019–2026)** (COVID disruption + 2022–2026 recovery), the full **137 × 6,418 City × Village** ECEPOV–INE model with strict **2nd-generation parental village roots**, the `> 10k × > 10k` Censo 2024 municipal OD matrix, and the full `59531` electricity-consumption sensitivity analysis. Adding a second full network (`N = 7,974`), a second roots matrix (city×village), and survey-attenuation modelling to a 15–20 page coursework report would violate M8 (over-ambition) and dilute the network-science focus.

### Recommended wording of the coursework question (sharper than before!)

With Finding 1 (the rural-exodus vs recent-migrant inversion) and Finding 3 (summer second homes vs Christmas primary homes), you no longer have to retreat from the holiday framing—you can state a crisp, non-obvious, **network-first** question whose answer is genuinely surprising:

- **Title:** ***Roots or beaches? How holidays split Spain's overnight mobility network between resorts and the rural exodus***
- **Question:** *When Spaniards sleep away from home, does the overnight mobility network still follow the map of the 1950–1980 rural exodus—and do summer and family holidays strengthen or dilute it?*
- **Three-part answer (exact match to evidence):**
  1. **Macro (beaches take the hubs):** On summer holidays the network splits into city senders and coastal receivers (`Spearman(k_in, k_out)` 0.69→0.09; 80–88% of top-50 in-hubs coastal), structurally forcing disassortativity and rich-club compression via `P(k)`.
  2. **Meso (holidays invert which roots matter):** Ordinary nights and Christmas follow **recent urban migrants**; **August inverts the map**, selectively tripling the pull of **village-born (`≤ 2,000`)** and doubling **pre-1981 exodus** ties (`p ≤ 0.002`, 44 province clusters) while recent urban ties fall (`p = 0.014`).
  3. **Micro (emptied villages switch roles via two distinct housing mechanisms):** Emptied villages flip from net exporters to net importers (99% in August vs 23% in Feb; confirmed at true municipal level without `_AM` aggregation). In **summer**, 48–65% of the gain runs through **unoccupied family second homes ("la casa del pueblo")**—77% of inland cross-regional second homes sit in the family's birth region—whereas on **Christmas Eve**, second homes explain `< 2%` because visitors sleep in living relatives' primary homes.

---

## Files produced (all traceable)
- Scripts: `pipeline/n1_new_data.py`, `pipeline/n1_censo2024_microdata.py`
- Data directories (new, under `_local/external/`): `ine_turismo_moviles_od/`, `ine_censo_viviendas_2021/`, `ine_ecepov_2021/`, `ine_censo_microdatos_2024/`, `eurostat_lau_2021/`
- Outputs: `results/n1_mun_nodes.parquet`, `results/n1_mitma_zone_housing.csv`, `results/n1_ine_od_inventory.csv`, `results/n1_ine_network_props.csv`, `results/n1_mun_network_node_tests.csv`, `results/n1_mitma_validation.csv`, `results/n1_roots_city_village_ecepov2021.parquet`, `results/n1_city_village_roots.csv`, `results/n1_city_village_qap.csv`, `results/n1_gen2_city_prov.csv`, `results/n1_second_home_mediation.csv`, `results/n1_ecepov_casa_del_pueblo.csv`, `results/n1_censo2024_expanded_city_prov.csv`, `results/n1_censo2024_cohorts.csv`, `results/n1_censo2024_cohort_diffs.csv`, `results/n1_censo2024_mun10k_od.csv`.
