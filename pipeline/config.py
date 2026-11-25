"""Frozen specification for question B ("Roots or beaches?").

Single source of truth for every pipeline step. Change a value here only with a
note in .planning/2026-11-24-consolidation/decisions.md explaining why.
"""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "derived"
RESULTS = ROOT / "results"
FIGS = ROOT / "reports" / "figs_final"

# ---- Inputs ---------------------------------------------------------------
PERNOCT_2024 = DATA / "pernoct_mun_2024_6m.parquet"       # cols: date, r (residence), s (sleep zone), p (persons)
PERNOCT_2022_23 = DATA / "pernoct_mun_2022_2023.parquet"
ZONES_GPKG = ROOT / "clean_data" / "v2" / "zones" / "municipios_mitma.gpkg"
COASTAL_FLAGS = ROOT / "explore" / "out_p86_coastal_flags.csv"  # use column COASTAL_COL; never recompute with a buffer
COASTAL_COL = "coastal_true"

# ---- Node set ---------------------------------------------------------------
# Peninsular mainland MITMA zones only: drop foreign zones, Balearics (07), Canaries (35, 38), Ceuta (51), Melilla (52).
EXCLUDED_PROVINCES = ("07", "35", "38", "51", "52")
# Also drop the 49 uninhabited shared-territory zones with 11-digit ids (no population, no flows; audit M3).
DROP_LONG_IDS = True          # drop zone ids longer than 8 characters
EXPECTED_N_NODES = 2435       # assert this after filtering (was 2484 before audit M3)
# Zone aggregation (audit C1/M4): MITMA '_AM' zones group several small municipalities.
# Every node-level comparison must control log(n_municipalities) and log(area), and be repeated on single-municipality zones.

# ---- Networks ---------------------------------------------------------------
# Main nights (2024). Directed, weighted: edge i->j = residents of i sleeping in j, i != j.
MAIN_NIGHTS = {
    "ordinary": "2024-02-14",   # Wednesday
    "summer":   "2024-08-15",   # holiday peak
    "allsaints": "2024-11-01",  # family holiday, no beach (pattern-mining #2: behaves like a big weekend)
}
# Sister nights: every result must be reported as median and range over these.
SISTER_NIGHTS = {
    "ordinary": "all Tue/Wed/Thu of Feb 2024 present in data",
    "summer":   [f"2024-08-{d:02d}" for d in range(10, 19)],
    "allsaints": ["2024-11-01", "2024-11-02", "2024-11-03"],  # long weekend; was 10-31 (an ordinary Wed->Thu night)
}
# DATE CONVENTION (pattern-mining #1, empirically robust; the MITMA methodology v8 does not state it):
# the file dated d holds the night ENDING on d (d-1 -> d). So "2024-08-15" = night 14->15 Aug (eve of the
# Assumption holiday), "2024-11-01" = night 31 Oct -> 1 Nov, "2024-12-25" = Christmas Eve night.
# Every main night is therefore "the eve of the holiday". Say this in the methods.

# Appendix-only comparison nights (same "eve" convention): family holidays.
APPENDIX_NIGHTS = {
    "christmas": "2024-12-25",   # Christmas Eve night
    "goodfriday": "2024-03-29",  # Holy Thursday -> Good Friday night
}

# Zones with data drop-outs in 2024 (pattern-mining #13): exclude from node-level results.
ARTEFACT_ZONES = ("31905_AM", "19103_AM", "40149_AM")  # Beriáin, Checa, Olombrada
# Cross-year replication (Nov 1 2023 is missing: MITMA outage 30 Oct - 3 Nov 2023).
YEAR_REPLICATION = {
    "ordinary": ["2022-02-16", "2023-02-15"],   # Wednesdays
    "summer":   ["2022-08-15", "2023-08-15"],
    "allsaints": ["2022-11-01"],
}

# Binary backbones. Both are computed; MAIN is used in the text, the other is a robustness check.
BACKBONE_MAIN = {"kind": "weight_threshold", "w_min": 10}
BACKBONE_CHECK = {"kind": "equal_edges", "reference": "ordinary"}  # keep top-E edges by weight, E = E(ordinary main night)
MAX_DENSITY = 0.05

# ---- Window means (audit C3) ------------------------------------------------
# Never use sum(p)/count(distinct date) per pair. Mean over a window = sum(p) / number of nights in the window
# (missing nights count as zero). Preferred: one network per night, then median and range.

# ---- Rich club (audit M1) ---------------------------------------------------
# Compare networks by RANK (top-n nodes, n in RICH_CLUB_TOP_N), not by fixed degree; report the rewired ceiling.
RICH_CLUB_TOP_N = (25, 50, 100, 200)

# ---- Null models ------------------------------------------------------------
SEED = 20261124
N_REWIRE = 50            # degree-preserving rewirings (directed, preserve k_in and k_out)
REWIRE_SWAPS_PER_EDGE = 10
N_QAP_PERM = 2000
N_SPATIAL_DRAWS = 100    # draws from the doubly-constrained gravity null

# ---- Stratification ---------------------------------------------------------
SMALL_VILLAGE_POP = 6000   # size-matched comparison for the k-core promotion
