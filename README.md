# Spain holiday roots

**Question.** On holiday, does Spain go back to the villages its families left? Does Spain's holiday mobility network still carry the imprint of the great rural-to-urban migration of 1950–1975 ("rural exodus")?

Coursework for UCL COMP0123 *Complex Networks*, to be followed by a journal paper.

**Report:** [`reports/manuscript/output/main.pdf`](reports/manuscript/output/main.pdf) (LaTeX sources in `reports/manuscript/`). **Agents and collaborators: start with [`AGENTS.md`](AGENTS.md).**

## Three networks over the same ~2,500 Spanish municipal zones
1. **Everyday network**: overnight stays and trips on ordinary days (MITMA Big Data mobility study).
2. **Holiday network**: where residents sleep in mid-August, at Easter, on All Saints' Day and at Christmas (MITMA *pernoctaciones*, 2022–2024).
3. **Roots network**: where people were born versus where they live (INE census and register).

## Main results so far (pilots plus robustness checks)
1. **The network changes shape on holidays.**
   - Cities connect mostly to villages: degree assortativity falls from −0.13 to −0.21.
   - The city rich club disappears.
   - Emptied villages jump into the network core (k-core ×4.5).
2. **The holiday network mirrors the roots network.** Comparing neighbouring cities: roots elasticity ≈ 0.72; partial QAP r = 0.49. This survives controls for income, distance, road time and beaches, replicates in 2022–2024, and is confirmed by an independent source (INE three-operator tourism data).
3. **Summer territories.** Bipartite communities under a spatial null; 12 of 16 are non-contiguous "diaspora archipelagos".

## Layout
| Path | Contents |
|---|---|
| `explore/` | Analysis scripts (`p60`–`p94`, `p100`); `common.py` holds helpers |
| `reports/` | Write-ups and figures |
| `research_notes/` | Agent reports: robustness checks, network layer, data scout, literature |
| `.planning/` | Planning and findings log |
| `data/`, `_local/`, `raw_data_cache/`, `clean_data/`, `metadata_cache/` | Data, not in git. Rebuild from the download scripts (`explore/p80_*`, `p81_*`, `p83_*`, `p100_download_summer.py`) and the MITMA/INE URLs |

Environment: `uv venv && uv pip install -r requirements.txt`; run scripts with `.venv/bin/python -I`.
