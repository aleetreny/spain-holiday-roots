# Replication package: Roots or beaches? Historical internal migration and holiday mobility in Spain

Alejandro Treny Ortega, University College London. DOI: [10.5281/zenodo.24624569](https://doi.org/10.5281/zenodo.24624569).

This archive holds the code, the pre-registration and the derived data behind the article *Roots or beaches? Historical internal migration and holiday mobility in Spain*. It is a snapshot of the `journal` branch of <https://github.com/aleetreny/spain-holiday-roots>. `VERSION.txt` gives the exact commit.

## What is in the archive

| Path | Contents |
|---|---|
| `journal/analysis/PREREGISTRATION.md`, `design.yml` | Pre-registration, frozen before the new MITMA nights were opened. The Confirmation II addendum and the two deviations are at the end of the file. The public commit history on GitHub dates them. |
| `journal/analysis/A_calendar/` | Analysis A: the night calendar, the nightly roots elasticities for 2022–2024 and the pre-registered tests. |
| `journal/analysis/A_calendar_noart/` | The pre-registered 2022–2024 tests run without the three artefact zones (Deviation 2). The article reports these values. |
| `journal/analysis/A2_confirmation2025/` | Confirmation II: the method bridge and the tests on unseen 2025 nights. |
| `journal/analysis/B_cohorts/`, `D_villages/` | Cohort models (arrival period and birthplace size) and the receiving end in the most emptied zones. |
| `journal/analysis/C_municipal/` | Parents' villages (ECEPOV 2021). This analysis failed its pre-registered gate and is reported in the Supplementary Information only. |
| `journal/analysis/data_build/` | Downloaders for the MITMA files of 2022–2025. |
| `journal/analysis/RESULTS_SUMMARY.md`, `REPORT_*.md` | One-page results summary and one report per analysis. Each number in the reports traces to a CSV in an `out/` folder. |
| `journal/figures/` | Figure scripts. `out/` holds every figure as PDF and PNG, and `figN_data.csv` holds the values plotted in figure N. |
| `journal/pilots/design/` | Feasibility pilots run before the pre-registration. |
| `pipeline/`, `explore/`, `results/`, `tests/` | The earlier analysis, first written as coursework for UCL COMP0123 *Complex Networks*. The journal analyses import its modules and read its outputs. `results/` holds its tables and diagnostic plots. |

## Where each result in the article comes from

| Article item | File |
|---|---|
| Table 2, 2022–2024 columns | `journal/analysis/A_calendar_noart/out/preregistered_tests_noart.csv` |
| Table 2, 2025 columns and decisions | `journal/analysis/A2_confirmation2025/out/preregistered_tests.csv`, `combined_decision.csv` |
| Figures 1–4 and S1–S5 | `journal/figures/out/fig*.pdf` and `fig*_data.csv` |
| Supplementary Table S1 (robustness) | `journal/analysis/A_calendar/out/robustness.csv` |
| Supplementary Table S2 (receiving zones) | `journal/analysis/A_calendar/out/descriptive_q5_villages_daily.csv` |
| Supplementary Table S3 (discovery and all nights) | `journal/analysis/A_calendar_noart/out/preregistered_tests_noart.csv` |
| Deviation 2, old against new values | `journal/analysis/A_calendar_noart/out/comparison_noart.csv` |
| MITMA method bridge | `journal/analysis/A2_confirmation2025/out/method_bridge_dates.csv`, `method_bridge_summary.csv` |
| Placebo windows | `placebo_summary*.csv` and `placebo_windows*.csv` in the A, A2 and A_calendar_noart `out/` folders |
| Nightly elasticities, all nights | `nightly_pooled.csv` and `nightly_cohorts.csv` in the A and A2 `out/` folders |

## Data that are not in the archive

All inputs are public. They are not redistributed here because of their size. Download them from their providers:

| Source | Used for | Folder expected by the code |
|---|---|---|
| MITMA, *Estudios de movilidad con Big Data*, daily overnight stays by municipality (<https://movilidad-opendata.mitma.es/RSS.xml>) | All flows | `_local/external/mitma_pernoct/` (nights used in the coursework), `_local/external_journal/mitma_pernoct/` (other 2022–2024 nights), `_local/external_journal/mitma_pernoct_2025/` |
| MITMA zoning files (same portal) | Zones and municipality links | `raw_data_cache/v2/zonificacion/`, `clean_data/v2/zones/` |
| INE, Censo Anual table 68545, residents by province of birth (<https://www.ine.es/jaxiT3/Tabla.htm?t=68545>) | Roots shares, 137 cities | `_local/external/ine_birth_residence/` |
| INE, Censo Anual 2024, 10 % public-use microdata (<https://www.ine.es/ftp/microdatos/censoanual/vc/datos_muestra_2024.zip>) | Cohort models, 675 municipalities | `_local/external/ine_censo_microdatos_2024/` |
| INE, intercensal municipal populations (<https://www.ine.es/intercensal/>) | Depopulation ranking | `_local/external/ine_historical_pop/` |
| INE, household income atlas ADRH (<https://www.ine.es/jaxiT3/Tabla.htm?t=30656>) | Income controls | `_local/external/ine_adrh/` |
| INE, Census 2021 dwellings, table 59525 | Second homes | `_local/external/ine_censo_viviendas_2021/` |
| INE, tourism measured with mobile phones, municipal origin–destination matrices | Parents' villages (SI) | `_local/external/ine_turismo_moviles_od/` |
| INE, ECEPOV 2021 microdata (<https://www.ine.es/ftp/microdatos/ecepov/datos_2021.zip>) | Parents' villages (SI) | `_local/external/ine_ecepov_2021/` |
| Eurostat GISCO, Local Administrative Units 2021 | Parents' villages (SI) | `_local/external_journal/gisco_lau_2021/` |
| BOE resolutions on the labour calendar, 2022–2025 | Holiday calendar | Downloaded by `journal/analysis/A_calendar/a01_calendar.py`. The parsed result is in `out/boe_holidays_*.csv`. |

Some downloaders are in `explore/` (`p80_download_mitma_22_23.py`, `p100_download_summer.py`, `p81_download_adrh.py`, `p83_download_ine_turismo.py`) and in `journal/analysis/data_build/`. The night manifests (`night_manifest.csv`) list every MITMA night and the reason for each exclusion. The 2025 manifest also gives the SHA-256 checksum of each file.

## How to run the code

1. Create the environment with Python 3.12: `uv venv --python 3.12 && uv pip install -r requirements.txt`. Run every script as `.venv/bin/python -I <script>`.
2. Put the data in the folders listed above, under the archive root. If the earlier-analysis inputs (`data/`, `raw_data_cache/`, `clean_data/`, `_local/external/`, `results/`) are in another copy of the repository, set `SHR_MAIN_ROOT` to the path of that copy.
3. Build the earlier-analysis inputs. `explore/p60_pernoct_build.py`, `p61_pernoct_mun.py` and `p88_replication_years.py` write `data/derived/`. The `pipeline/` scripts then write `results/`. Start with `b01_build_networks.py`, then run `p_patterns.py`, `r1_villages.py`, `r2_roots.py`, `n1_new_data.py`, `n1_censo2024_microdata.py` and `v1_verify_n1.py`.
4. Run the journal analyses in this order:
   - A: `a01_calendar.py` to `a06_descriptive.py`.
   - A2: `a2_data.py`, then `a2_models.py`.
   - A_calendar_noart: `run_noart.py`, then `compare_noart.py`.
   - B: `run_b_cohorts.py`.
   - D: `run_d_villages.py`, then `d_inland_decomp.py`.
   - C: `run.py`.
   - Figures: `journal/figures/make_all.py`.
5. Large intermediates are written to `_local/journal_data/`. The wild cluster bootstrap steps take several hours on a laptop.

You do not need to run anything to check the article. Every reported estimate and every plotted value is already in the CSV files of this archive.

## Licence

The code is under the MIT licence. Derived data, tables, figures and documentation are under CC BY 4.0. See `LICENSE.md`. The source data stay under the terms of their providers: MITMA, INE, Eurostat GISCO and the BOE.

## Citation

Please cite the article and this archive. `CITATION.cff` has the details.
