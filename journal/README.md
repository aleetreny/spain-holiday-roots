# Journal article

*Roots or beaches? Historical internal migration and holiday mobility in Spain* (in preparation, EPJ Data Science).

The article extends the coursework in the rest of this repository. It covers 1,275 nights from 2022 to 2025, pre-registers its tests before the new nights were opened, and checks the August result on unseen 2025 nights.

| Folder | Contents |
|---|---|
| [`analysis/`](analysis/) | [Pre-registration](analysis/PREREGISTRATION.md), analyses A–D with their reports and output tables, and the [one-page results summary](analysis/RESULTS_SUMMARY.md) |
| [`figures/`](figures/) | Figure scripts. `out/` holds each figure as PDF and PNG, with the plotted values in `figN_data.csv` |
| [`manuscript/`](manuscript/) | LaTeX sources of the article and its supplementary information (Springer Nature template) |
| [`replication/`](replication/) | README, licence and citation file of the Zenodo archive ([10.5281/zenodo.24624569](https://doi.org/10.5281/zenodo.24624569)), and the script that builds it |
| [`arxiv/`](arxiv/) | Script that builds the arXiv source bundle |
| [`pilots/`](pilots/) | Feasibility pilots run before the pre-registration |

Scripts read code and data from the repository root. Set `SHR_MAIN_ROOT` if the earlier-analysis inputs live in another copy of the repository.
