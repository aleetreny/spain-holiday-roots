"""b00_common: node table, edge builder and backbones for question B.

Everything reads its inputs and parameters from ``config.py``.

Residence-code mapping (documented, not re-done here)
-----------------------------------------------------
The raw MITMA *pernoctaciones* files give the residence zone at MITMA
*district* level and the overnight zone at MITMA *municipal* level. The legacy
scripts ``explore/p61_pernoct_mun.py`` (2024) and
``explore/p88_replication_years.py`` (2022-23) mapped residence district ->
municipal zone with the DISTINCT pairs (distrito_mitma, municipio_mitma) of
``raw_data_cache/v2/zonificacion/relacion_ine_zonificacionMitma.csv`` (inner
join, so unmapped residence codes are dropped) and summed persons by
(date, r, s). The parquet files named in config (PERNOCT_2024,
PERNOCT_2022_23) are the output of that step, so ``r`` and ``s`` are already in
the same MITMA municipal-zone id space as the zones GeoPackage. This module
therefore applies no further mapping: exactly as ``p92`` does, it keeps rows with
both ``r`` and ``s`` in the mainland node set and ``r != s``, and sums ``p``
per (r, s) (a no-op on these files, where (date, r, s) is unique).
"""
from __future__ import annotations

import sys
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C  # noqa: E402

# Province -> autonomous community (copied from explore/common.py; data only).
PROV2CCAA = {
    **{p: "AND" for p in ["04", "11", "14", "18", "21", "23", "29", "41"]},
    **{p: "ARA" for p in ["22", "44", "50"]},
    "33": "AST", "07": "BAL", "35": "CAN", "38": "CAN", "39": "CANT",
    **{p: "CYL" for p in ["05", "09", "24", "34", "37", "40", "42", "47", "49"]},
    **{p: "CLM" for p in ["02", "13", "16", "19", "45"]},
    **{p: "CAT" for p in ["08", "17", "25", "43"]},
    **{p: "VAL" for p in ["03", "12", "46"]},
    "06": "EXT", "10": "EXT",
    **{p: "GAL" for p in ["15", "27", "32", "36"]},
    "28": "MAD", "30": "MUR", "31": "NAV",
    **{p: "PVA" for p in ["01", "20", "48"]},
    "26": "RIO", "51": "CEU", "52": "MEL",
}

HIST_POP_CSV = C.ROOT / "_local" / "external" / "ine_historical_pop" / "municipal_population_1900_2025.csv"
INE_MITMA_REL = C.ROOT / "raw_data_cache" / "v2" / "zonificacion" / "relacion_ine_zonificacionMitma.csv"
# Exactly the columns p92 sums per zone (pop_1991 and pop_2011 are not used there).
_HCOLS = ["pop_1900", "pop_1910", "pop_1920", "pop_1930", "pop_1940", "pop_1950",
          "pop_1960", "pop_1970", "pop_1981", "pop_2001", "pop_2021"]
_PEAK_COLS = ["pop_1900", "pop_1910", "pop_1920", "pop_1930", "pop_1940", "pop_1950", "pop_1960"]


# ---------------------------------------------------------------------------
# Nodes
# ---------------------------------------------------------------------------
def historical_loss() -> pd.Series:
    """loss_peak_2021 = 1 - pop_2021 / max(pop_1900..pop_1960), per MITMA municipal zone.

    Identical to p92: INE municipalities are inner-joined to MITMA zones and each
    census column is summed per zone with ``min_count=1``.
    """
    H = pd.read_csv(HIST_POP_CSV, dtype={"ine_code": str, "province_code": str})
    rel = pd.read_csv(INE_MITMA_REL, sep="|", dtype=str)[["municipio_ine", "municipio_mitma"]].drop_duplicates()
    H = H.merge(rel, left_on="ine_code", right_on="municipio_ine", how="inner")
    Z = H.groupby("municipio_mitma")[_HCOLS].sum(min_count=1)
    peak = Z[_PEAK_COLS].max(axis=1)
    return (1.0 - Z["pop_2021"] / peak).rename("loss_peak_2021")


def load_nodes() -> pd.DataFrame:
    """Mainland node table, sorted by zone id; row position == node index.

    Columns: idx, id, name, population, prov, ccaa, coastal_true, loss_peak_2021, x, y,
    n_municipalities, area_km2, is_am.
    x, y: representative point of the zone polygon, EPSG:3042 metres (as explore/common.py).
    n_municipalities: number of INE municipalities listed in the gpkg ``municipalities`` column
    ('; '-separated); '_AM' zones aggregate several. area_km2: polygon area in EPSG:3042.
    config.DROP_LONG_IDS drops the 49 uninhabited shared-territory zones with 11-digit ids (audit M3).
    """
    import geopandas as gpd

    g = gpd.read_file(C.ZONES_GPKG)
    g = g[g["id"].str.match(r"^\d")]                       # drop foreign zones
    g = g[~g["id"].str[:2].isin(C.EXCLUDED_PROVINCES)]      # drop islands, Ceuta, Melilla
    if C.DROP_LONG_IDS:
        g = g[g["id"].str.len() <= 8]                       # drop 11-digit shared territories
    pt = g.geometry.representative_point()
    nd = pd.DataFrame({
        "id": g["id"].values,
        "name": g["name"].values,
        "population": pd.to_numeric(g["population"], errors="coerce").values,
        "x": pt.x.values,
        "y": pt.y.values,
        "area_km2": (g.geometry.area / 1e6).values,
        "municipalities": g["municipalities"].values,
    })
    nd["n_municipalities"] = nd["municipalities"].map(
        lambda v: len([m for m in str(v).split(";") if m.strip()]) if isinstance(v, str) and v.strip() not in ("", "NA") else np.nan)
    nd["is_am"] = nd["id"].str.endswith("_AM")
    nd = nd.drop(columns="municipalities")
    nd["prov"] = nd["id"].str[:2]
    nd["ccaa"] = nd["prov"].map(PROV2CCAA)
    coast = pd.read_csv(C.COASTAL_FLAGS, dtype={"id": str}).set_index("id")[C.COASTAL_COL]
    nd["coastal_true"] = nd["id"].map(coast)
    nd["loss_peak_2021"] = nd["id"].map(historical_loss())
    nd = nd.sort_values("id").reset_index(drop=True)
    nd.insert(0, "idx", np.arange(len(nd)))

    assert len(nd) == C.EXPECTED_N_NODES, f"N={len(nd)} != {C.EXPECTED_N_NODES}"
    assert nd["id"].is_unique
    assert nd["ccaa"].notna().all()
    assert nd["coastal_true"].notna().all(), "zones without a coastal_true flag"
    assert nd["n_municipalities"].notna().all(), "zones without a municipalities list"
    nd["coastal_true"] = nd["coastal_true"].astype(bool)
    return nd


# ---------------------------------------------------------------------------
# Edges
# ---------------------------------------------------------------------------
def pernoct_file(date: str) -> Path:
    year = int(str(date)[:4])
    if year == 2024:
        return C.PERNOCT_2024
    if year in (2022, 2023):
        return C.PERNOCT_2022_23
    raise ValueError(f"no pernoct file for {date}")


def available_dates(path: Path) -> list[str]:
    con = duckdb.connect()
    d = con.execute("SELECT DISTINCT date FROM read_parquet(?) ORDER BY 1", [str(path)]).fetchall()
    return [str(x[0]) for x in d]


def edges_from_frame(df: pd.DataFrame, node_ids, expected_n: int | None = None) -> pd.DataFrame:
    """Pure edge builder (tested). ``df`` has columns r, s, p (one date).

    Keeps r, s in ``node_ids``; drops self-loops r == s; sums p per (r, s).
    Returns columns i, j (node indices), w, sorted by (i, j).
    """
    node_ids = list(node_ids)
    if expected_n is not None:
        assert len(node_ids) == expected_n, f"N={len(node_ids)} != {expected_n}"
    pos = {z: k for k, z in enumerate(node_ids)}
    d = df[df["r"].isin(pos) & df["s"].isin(pos) & (df["r"] != df["s"])]
    d = d.groupby(["r", "s"], as_index=False)["p"].sum()
    out = pd.DataFrame({
        "i": d["r"].map(pos).astype(np.int32).values,
        "j": d["s"].map(pos).astype(np.int32).values,
        "w": d["p"].astype(float).values,
    })
    out = out[out["w"] > 0].sort_values(["i", "j"]).reset_index(drop=True)
    assert (out["i"] != out["j"]).all()
    return out


def build_edges(date: str, nodes: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Directed weighted edge list r -> s for one date on the mainland node set.

    Returns (edges[i, j, w], info) where info holds totals used for bookkeeping.
    """
    path = pernoct_file(date)
    con = duckdb.connect()
    raw = con.execute(
        "SELECT r, s, sum(p) AS p FROM read_parquet(?) WHERE date = CAST(? AS DATE) GROUP BY 1, 2",
        [str(path), str(date)],
    ).df()
    if raw.empty:
        raise ValueError(f"date {date} not present in {path.name}")
    ids = set(nodes["id"])
    in_ml = raw["r"].isin(ids) & raw["s"].isin(ids)
    info = {
        "date": str(date),
        "total_p_all": float(raw["p"].sum()),
        "total_p_mainland": float(raw.loc[in_ml, "p"].sum()),
        "total_p_mainland_away": float(raw.loc[in_ml & (raw["r"] != raw["s"]), "p"].sum()),
    }
    edges = edges_from_frame(raw, nodes["id"], expected_n=C.EXPECTED_N_NODES)
    return edges, info


def window_mean_from_frame(df: pd.DataFrame, n_nights: int) -> pd.DataFrame:
    """Mean weight per (r, s) over a window of ``n_nights`` nights (audit C3).

    ``df`` has columns date, r, s, p for the nights of the window. A pair missing on a
    night counts as 0 on that night: mean = sum(p) / n_nights, NEVER sum(p) / count(distinct date).
    The pipeline itself works on single nights; this helper exists so that any window mean
    computed later uses the correct denominator.
    """
    present = df["date"].nunique()
    assert n_nights >= present, f"window has {present} nights present but n_nights={n_nights}"
    out = df.groupby(["r", "s"], as_index=False)["p"].sum()
    out["p"] = out["p"] / float(n_nights)
    return out


# ---------------------------------------------------------------------------
# Backbones
# ---------------------------------------------------------------------------
def backbone_threshold(edges: pd.DataFrame, w_min: float) -> pd.DataFrame:
    """Keep directed edges with w >= w_min."""
    return edges[edges["w"] >= w_min].reset_index(drop=True)


def backbone_equal_edges(edges: pd.DataFrame, E: int) -> pd.DataFrame:
    """Keep exactly the E heaviest directed edges.

    Ties in w are broken deterministically by (i, j) ascending, so the result has
    exactly min(E, len(edges)) edges.
    """
    d = edges.sort_values(["w", "i", "j"], ascending=[False, True, True], kind="mergesort")
    return d.head(int(E)).sort_values(["i", "j"]).reset_index(drop=True)


def to_igraph(edges: pd.DataFrame, n: int, directed: bool = True):
    """Binary igraph graph. Undirected: an edge {a, b} exists if a->b or b->a is in ``edges``."""
    import igraph as ig

    if directed:
        el = list(zip(edges["i"].tolist(), edges["j"].tolist()))
        return ig.Graph(n=n, edges=el, directed=True)
    a = np.minimum(edges["i"].values, edges["j"].values)
    b = np.maximum(edges["i"].values, edges["j"].values)
    und = np.unique(np.stack([a, b], axis=1), axis=0) if len(a) else np.empty((0, 2), int)
    return ig.Graph(n=n, edges=[tuple(x) for x in und.tolist()], directed=False)
