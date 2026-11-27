"""Shared read-only paths and coursework module imports for Analysis A (Roots Calendar).

Reuses coursework code from the main checkout by import (never modifies coursework files).
Heavy outputs live under <journal_repo>/_local/journal_data/A_calendar/ (outside git).
Small tidy CSVs and code live under <worktree>/journal/analysis/A_calendar/.
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

WORKTREE = Path(__file__).resolve().parents[3]
A_DIR = Path(__file__).resolve().parent
OUT_DIR = A_DIR / "out"

MAIN_REPO = Path("/Users/alejandrotreny/Proyectos/spain-holiday-roots")
JOURNAL_REPO = Path("/Users/alejandrotreny/Proyectos/spain-holiday-roots-journal")

# Raw MITMA daily files
DISCOVERY_RAW_DIR = MAIN_REPO / "_local" / "external" / "mitma_pernoct"
CONFIRMATION_RAW_DIR = JOURNAL_REPO / "_local" / "external_journal" / "mitma_pernoct"

# Coursework derived parquets (for exact validation on discovery dates)
CW_PERNOCT_2024 = MAIN_REPO / "data" / "derived" / "pernoct_mun_2024_6m.parquet"
CW_PERNOCT_2223 = MAIN_REPO / "data" / "derived" / "pernoct_mun_2022_2023.parquet"
CW_RESULTS_DIR = MAIN_REPO / "results"
CW_EXT_DIR = MAIN_REPO / "_local" / "external"
CW_REL_CSV = MAIN_REPO / "raw_data_cache" / "v2" / "zonificacion" / "relacion_ine_zonificacionMitma.csv"
CW_TT_P84 = MAIN_REPO / "explore" / "out_p84_city_zone_tt.parquet"

# Heavy journal outputs (never committed)
HEAVY_DIR = JOURNAL_REPO / "_local" / "journal_data" / "A_calendar"
PARQUET_DIR = HEAVY_DIR / "pernoct_mun_2022_2024"
PANEL_DIR = HEAVY_DIR / "panel"

# Import coursework modules from MAIN_REPO so C.ROOT resolves to main checkout data
_main_pipe = str(MAIN_REPO / "pipeline")
if _main_pipe not in sys.path:
    sys.path.insert(0, _main_pipe)

import config as C  # noqa: E402
import b00_common as B0  # noqa: E402
import v1_verify_n1 as V1  # noqa: E402

assert C.ROOT == MAIN_REPO, f"Expected C.ROOT={MAIN_REPO}, got {C.ROOT}"


def ensure_dirs() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    PARQUET_DIR.mkdir(parents=True, exist_ok=True)
    PANEL_DIR.mkdir(parents=True, exist_ok=True)


def make_logger(t0: float | None = None):
    start = time.time() if t0 is None else t0

    def log(*args) -> None:
        print(f"[{time.time() - start:7.1f}s]", *args, flush=True)

    return log
