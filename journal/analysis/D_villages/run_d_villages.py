"""Analysis D: The receiving end (emptied villages) for EPJ Data Science.

Implements items D.1 through D.3 strictly on the coursework MITMA cache (2022-2024):
  D.1 Fixed sample (R1 loss rank, 3 artefact zones & 56 missing-loss zones excluded:
      n=2,376 all zones, n=1,473 single-municipality zones, 47 provinces).
      For main nights, every cached 2024 night (180), and every cached 2022-2023 night (126):
      share of zones with inflow > outflow by loss quintile, and arrivals per resident and
      departures per resident separately, to decompose August vs Christmas Eve.
  D.2 Uncertainty via stratified bootstrap of zones within province (B=4,999; primary)
      and 47-province cluster bootstrap (sensitivity), for both all zones and single-municipality
      zones only, covering quintile levels, paired holiday-minus-ordinary changes, Q5-Q1 gaps,
      and controlled FE-OLS gradients (exact linear decomposition net = arrivals - departures).
  D.3 Housing slope attenuation (Census 2021 table 59525 non-primary dwellings):
      August vs Christmas Eve on a common sample with raw slopes, housing-adjusted slopes,
      absolute attenuation, relative attenuation share, and formal paired August-minus-Christmas
      differences with bootstrap intervals. Labelled strictly as attenuation, never mediation.

Run:
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
    .venv/bin/python -I \
    journal/analysis/D_villages/run_d_villages.py
"""
from __future__ import annotations

import re
import os
import sys
import time
import zlib
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
import scipy.sparse as sps
from scipy import stats
from scipy.linalg import qr

# One repository root holds code and data. SHR_MAIN_ROOT may point the coursework inputs
# (data/, raw_data_cache/, clean_data/, _local/external/, results/) at another checkout.
REPO = Path(__file__).resolve().parents[3]
MAIN_ROOT = Path(os.environ.get("SHR_MAIN_ROOT", REPO)).resolve()
JOURNAL_DATA_BD = REPO / "_local" / "journal_data" / "BD"
HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
OUT.mkdir(parents=True, exist_ok=True)
JOURNAL_DATA_BD.mkdir(parents=True, exist_ok=True)

# Reuse coursework pipeline by import only (never edit pipeline/)
sys.path.insert(0, str(MAIN_ROOT / "pipeline"))
import config as C  # noqa: E402
import b00_common as B0  # noqa: E402
import r1_villages as R1  # noqa: E402
import v1_verify_n1 as V1  # noqa: E402

PERNOCT_2024 = MAIN_ROOT / "data" / "derived" / "pernoct_mun_2024_6m.parquet"
PERNOCT_2022_23 = MAIN_ROOT / "data" / "derived" / "pernoct_mun_2022_2023.parquet"
NODES_PQ = MAIN_ROOT / "results" / "b01_networks" / "nodes.parquet"
NET_DIR = MAIN_ROOT / "results" / "b01_networks"
P_ZONE_DAILY = MAIN_ROOT / "results" / "p_zone_daily.parquet"
P_VIV59525 = MAIN_ROOT / "_local" / "external" / "ine_censo_viviendas_2021" / "raw" / "censo2021_59525.csv"
INE_MITMA_REL = MAIN_ROOT / "raw_data_cache" / "v2" / "zonificacion" / "relacion_ine_zonificacionMitma.csv"

CONTROLS = ["lpop", "lnmun", "larea", "coastal", "ldcity", "ldcoast"]
DANA_DATES = {"2024-10-29", "2024-10-30", "2024-10-31"}
N_BOOT_ZONE = 4999
N_BOOT_PROV = 4999
N_WCB_OLS = 9999

T0 = time.time()


def log(*a):
    print(f"[{time.time() - T0:7.1f}s]", *a, flush=True)


def rng_for(*parts) -> np.random.Generator:
    ints = [C.SEED] + [
        zlib.crc32(str(p).encode()) if not isinstance(p, (int, np.integer)) else int(p)
        for p in parts
    ]
    return np.random.default_rng(np.random.SeedSequence(ints))


def load_fixed_nodes_and_housing() -> pd.DataFrame:
    """Load the exact R1 node table (2,435 zones) with frozen R1 loss rank/quintiles and 59525 housing stock."""
    nd = pd.read_parquet(NODES_PQ).sort_values("idx").reset_index(drop=True)
    assert len(nd) == C.EXPECTED_N_NODES and (nd.idx.to_numpy() == np.arange(len(nd))).all()
    xy = nd[["x", "y"]].to_numpy(dtype=float) / 1000.0
    D = np.hypot(xy[:, None, 0] - xy[None, :, 0], xy[:, None, 1] - xy[None, :, 1])
    big = nd["population"].to_numpy(dtype=float) >= R1.BIG_CITY_POP
    coast = nd["coastal_true"].to_numpy(dtype=bool)
    nd["dist_city_km"] = D[:, big].min(axis=1)
    nd["dist_coast_km"] = D[:, coast].min(axis=1)
    nd["lpop"] = np.log(nd["population"].to_numpy(dtype=float))
    nd["lnmun"] = np.log(nd["n_municipalities"].to_numpy(dtype=float))
    nd["larea"] = np.log(nd["area_km2"].to_numpy(dtype=float))
    nd["coastal"] = coast.astype(float)
    nd["ldcity"] = np.log1p(nd["dist_city_km"].to_numpy(dtype=float))
    nd["ldcoast"] = np.log1p(nd["dist_coast_km"].to_numpy(dtype=float))
    nd["artefact"] = nd["id"].isin(C.ARTEFACT_ZONES)
    assert int(nd["artefact"].sum()) == 3

    ok = nd["loss_peak_2021"].notna() & ~nd["artefact"]
    nd["eligible"] = ok
    nd["loss_pct"] = np.nan
    nd.loc[ok, "loss_pct"] = nd.loc[ok, "loss_peak_2021"].rank(pct=True, method="average")
    nd["loss_q"] = np.nan
    nd.loc[ok, "loss_q"] = pd.qcut(nd.loc[ok, "loss_peak_2021"], 5, labels=False) + 1
    nd["emptied_inland_50"] = ok & (nd["loss_peak_2021"] >= 0.5) & (~nd["coastal_true"])
    assert int(ok.sum()) == 2376
    assert int((ok & (nd["n_municipalities"] == 1)).sum()) == 1473
    assert int(nd["emptied_inland_50"].sum()) == 455

    # Parse raw Census 2021 table 59525 (exact V1 C3 mapping)
    rows_59 = []
    with open(P_VIV59525, encoding="utf-8-sig") as f:
        f.readline()
        for line in f:
            p = line.rstrip("\n").split(";")
            if len(p) < 4 or not p[1].strip():
                continue
            code = p[1].strip()[:5]
            if not re.match(r"^\d{5}$", code):
                continue
            val = V1.parse_es_num(p[3])
            if val is not None:
                rows_59.append((code, p[2].strip(), val))
    viv59 = pd.DataFrame(rows_59, columns=["ine_code", "tipo", "val"]).pivot(index="ine_code", columns="tipo", values="val")
    assert len(viv59) == 8131

    rel_full = pd.read_csv(INE_MITMA_REL, sep="|", dtype=str)
    zone_of = (
        rel_full.groupby(["municipio_ine", "municipio_mitma"]).size().reset_index(name="n")
        .sort_values(["municipio_ine", "n", "municipio_mitma"], ascending=[True, False, True])
        .drop_duplicates("municipio_ine").set_index("municipio_ine")["municipio_mitma"]
    )
    zv = pd.DataFrame({"municipio_ine": zone_of.index, "municipio_mitma": zone_of.values}).merge(
        viv59.reset_index(), left_on="municipio_ine", right_on="ine_code", how="inner"
    )
    zg = zv.groupby("municipio_mitma")[["Vivienda no principal", "Total"]].sum()
    nd["viv_no_principal"] = nd["id"].map(zg["Vivienda no principal"])
    nd["viv_total"] = nd["id"].map(zg["Total"])
    nd["share_np"] = nd["viv_no_principal"] / nd["viv_total"]
    nd["log1p_np_per_res"] = np.log1p(nd["viv_no_principal"] / nd["population"])
    assert nd.loc[ok, "share_np"].notna().all(), "Missing housing share_np on eligible zones"
    assert nd.loc[ok, "log1p_np_per_res"].notna().all(), "Missing log1p_np_per_res on eligible zones"
    return nd


def build_all_cached_zone_panel(nd: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compute s_in and s_out directly from coursework cache parquets for all 306 nights & verify against p_zone_daily."""
    con = duckdb.connect()
    valid_df = nd[["id"]].copy()
    con.register("valid_zones", valid_df)

    parts = []
    tot_parts = []
    for path in (PERNOCT_2024, PERNOCT_2022_23):
        # s_in and s_out on full mainland node set (r != s), matching R1 / p_zone_daily
        q_in = con.execute(f"""
            SELECT CAST(t.date AS VARCHAR) AS date, t.s AS id, sum(t.p) AS s_in
            FROM read_parquet('{path}') t
            JOIN valid_zones vr ON t.r = vr.id
            JOIN valid_zones vs ON t.s = vs.id
            WHERE t.r != t.s
            GROUP BY 1, 2
        """).df()
        q_out = con.execute(f"""
            SELECT CAST(t.date AS VARCHAR) AS date, t.r AS id, sum(t.p) AS s_out
            FROM read_parquet('{path}') t
            JOIN valid_zones vr ON t.r = vr.id
            JOIN valid_zones vs ON t.s = vs.id
            WHERE t.r != t.s
            GROUP BY 1, 2
        """).df()
        q_tot = con.execute(f"""
            SELECT CAST(t.date AS VARCHAR) AS date,
                   count(DISTINCT t.r) AS n_res_zones,
                   sum(t.p) AS total_mainland_p,
                   sum(CASE WHEN t.r != t.s THEN t.p ELSE 0 END) AS total_away_p
            FROM read_parquet('{path}') t
            JOIN valid_zones vr ON t.r = vr.id
            JOIN valid_zones vs ON t.s = vs.id
            GROUP BY 1 ORDER BY 1
        """).df()
        m = q_in.merge(q_out, on=["date", "id"], how="outer").fillna(0.0)
        parts.append(m)
        tot_parts.append(q_tot)

    raw_panel = pd.concat(parts, ignore_index=True)
    night_meta = pd.concat(tot_parts, ignore_index=True).sort_values("date").reset_index(drop=True)
    dates = sorted(night_meta["date"].unique())
    assert len(dates) == 306, f"Expected 306 cached nights, got {len(dates)}"

    # Complete cartesian grid (306 dates x 2,435 zones) so zero-flow zones keep explicit 0.0
    grid = pd.MultiIndex.from_product([dates, nd["id"].tolist()], names=["date", "id"]).to_frame(index=False)
    panel = grid.merge(raw_panel, on=["date", "id"], how="left").fillna({"s_in": 0.0, "s_out": 0.0})

    # Independent cross-check against results/p_zone_daily.parquet on all 306 nights!
    pz = pd.read_parquet(P_ZONE_DAILY)
    pz["date"] = pz["date"].astype(str)
    chk = panel.merge(pz[["date", "id", "in_away_ml", "out_away_ml"]], on=["date", "id"])
    assert len(chk) == len(panel)
    assert np.allclose(chk["s_in"], chk["in_away_ml"], rtol=1e-9, atol=1e-6), "s_in mismatch vs p_zone_daily"
    assert np.allclose(chk["s_out"], chk["out_away_ml"], rtol=1e-9, atol=1e-6), "s_out mismatch vs p_zone_daily"
    log("Verified independently computed s_in and s_out across all 306 nights x 2,435 zones against p_zone_daily (diff < 1e-6)")

    # Cross-check Q5 net-importer counts on Feb 14, Aug 15, Dec 25 against results/r1_net_balance_by_quintile.csv
    r1_q = pd.read_csv(MAIN_ROOT / "results" / "r1_net_balance_by_quintile.csv")
    q5_ids = set(nd.loc[nd["loss_q"] == 5, "id"])
    for d in ("2024-02-14", "2024-08-15", "2024-12-25"):
        sub = panel[(panel["date"] == d) & (panel["id"].isin(q5_ids))]
        sh = float((sub["s_in"] > sub["s_out"]).mean())
        exp = float(r1_q[(r1_q["date"] == d) & (r1_q["sample"] == "all") & (r1_q["cls"] == "Q5") & (r1_q["role"] == "main")]["share_net_importer"].iloc[0])
        assert abs(sh - exp) < 1e-12, (d, sh, exp)
    log("Verified Q5 net-importer shares on Feb14 (24.00%), Aug15 (98.95%), Dec25 (99.58%) match R1 exactly")

    # Add calendar & quality metadata to night_meta
    med_res_zones = night_meta["n_res_zones"].median()
    night_meta["year"] = night_meta["date"].str[:4].astype(int)
    night_meta["month"] = night_meta["date"].str[5:7].astype(int)
    night_meta["weekday_file"] = pd.to_datetime(night_meta["date"]).dt.day_name().str[:3]
    night_meta["is_dana"] = night_meta["date"].isin(DANA_DATES)
    night_meta["coverage_ratio"] = night_meta["n_res_zones"] / med_res_zones
    night_meta["fails_coverage_99"] = night_meta["coverage_ratio"] < 0.99
    # Year-specific median total mainland persons for coarse deviation flag
    yr_med_p = night_meta.groupby("year")["total_mainland_p"].transform("median")
    night_meta["total_dev_vs_yrmed"] = np.abs(night_meta["total_mainland_p"] - yr_med_p) / yr_med_p
    night_meta["fails_total_dev_10pct"] = night_meta["total_dev_vs_yrmed"] > 0.10

    NW = pd.read_csv(NET_DIR / "networks.csv", dtype=str)
    feb24_sis = set(NW.loc[(NW.group == "ordinary") & (NW.role == "sister"), "date"])
    aug24_sis = set(C.SISTER_NIGHTS["summer"])
    nov24_sis = set(C.SISTER_NIGHTS["allsaints"])

    def classify_date(d: str) -> str:
        if d == "2024-02-14":
            return "main_ordinary_2024"
        if d == "2024-08-15":
            return "main_summer_2024"
        if d == "2024-11-01":
            return "main_allsaints_2024"
        if d == "2024-12-25":
            return "main_christmas_2024"
        if d == "2024-03-29":
            return "main_goodfriday_2024"
        if d in DANA_DATES:
            return "dana_excluded_oct2024"
        if d in feb24_sis:
            return "sister_feb_2024"
        if d in aug24_sis:
            return "sister_aug_2024"
        if d in nov24_sis:
            return "sister_allsaints_2024"
        if d.startswith("2024-10") and pd.Timestamp(d).day_of_week in (1, 2, 3):
            return "ordinary_oct_clean_2024"
        if d in ("2022-02-16", "2023-02-15"):
            return "repl_ordinary"
        if d in ("2022-08-15", "2023-08-15"):
            return "repl_summer"
        if d == "2022-11-01":
            return "repl_allsaints"
        return f"cached_{d[:4]}_{d[5:7]}"

    night_meta["night_class"] = [classify_date(d) for d in night_meta["date"]]
    return panel, night_meta


def make_zone_bootstrap_indices(provs: np.ndarray, n_boot: int, rng: np.random.Generator) -> np.ndarray:
    """Stratified within-province zone bootstrap: resample n_p zones with replacement within each province p."""
    n = len(provs)
    out = np.empty((n_boot, n), dtype=np.int32)
    pos = 0
    for p in sorted(np.unique(provs)):
        idx = np.where(provs == p)[0]
        m = len(idx)
        draws = rng.integers(0, m, size=(n_boot, m))
        out[:, pos:pos + m] = idx[draws]
        pos += m
    assert pos == n
    return out


def make_prov_bootstrap_draws(provs: np.ndarray, n_boot: int, rng: np.random.Generator) -> tuple[list[np.ndarray], np.ndarray]:
    """Cluster bootstrap resampling whole provinces with replacement (G provinces)."""
    u_provs = sorted(np.unique(provs))
    idx_by_p = [np.where(provs == p)[0] for p in u_provs]
    G = len(u_provs)
    picks = rng.integers(0, G, size=(n_boot, G))
    return idx_by_p, picks


def summarize_boot(est: float, bvals: np.ndarray) -> dict:
    bvals = bvals[np.isfinite(bvals)]
    lo = float(np.quantile(bvals, 0.025))
    hi = float(np.quantile(bvals, 0.975))
    se = float(np.std(bvals, ddof=1))
    # Two-sided bootstrap p-value against 0 (sign-tail doubling)
    p_lo = (1 + np.sum(bvals <= 0.0)) / (len(bvals) + 1)
    p_hi = (1 + np.sum(bvals >= 0.0)) / (len(bvals) + 1)
    p2 = float(min(1.0, 2.0 * min(p_lo, p_hi)))
    return {"est": float(est), "se_boot": se, "ci_lo": lo, "ci_hi": hi, "p_boot": p2}


def fast_ols_fe_zone_boot(
    Y_mat: np.ndarray,
    x: np.ndarray,
    Z: np.ndarray,
    provs: np.ndarray,
    boot_idx: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Vectorized FE-OLS slopes of multiple outcomes Y_mat (n x M) on x + Z | prov FE
    for both full sample and every stratified within-province zone bootstrap draw.

    Because within-province resampling preserves exact province sizes and province blocks,
    each draw's province FE demeaning is simply demeaning within fixed slice [pos:pos+m]!
    """
    n, M = Y_mat.shape
    u_provs = sorted(np.unique(provs))
    # Reorder data into contiguous province slices matching make_zone_bootstrap_indices order!
    order = np.concatenate([np.where(provs == p)[0] for p in u_provs])
    inv_order = np.argsort(order)
    # Map boot_idx (which references original row indices) into reordered indices
    # Actually even simpler: directly index original arrays using boot_idx, and demean within
    # the known contiguous slices [s:e] because make_zone_bootstrap_indices places province p in [s:e]!
    slices = []
    pos = 0
    for p in u_provs:
        m = int(np.sum(provs == p))
        slices.append((pos, pos + m))
        pos += m

    def solve_one(Y_in: np.ndarray, x_in: np.ndarray, Z_in: np.ndarray) -> np.ndarray:
        # Demean [Y, x, Z] within each province slice
        A = np.column_stack([Y_in, x_in[:, None], Z_in])
        A_dm = A.copy()
        for s, e in slices:
            A_dm[s:e] -= A_dm[s:e].mean(axis=0, keepdims=True)
        Yd = A_dm[:, :M]
        xd = A_dm[:, M]
        Zd = A_dm[:, M + 1:]
        # Keep non-constant columns in Zd
        keep_z = np.std(Zd, axis=0) > 1e-12
        Zd = Zd[:, keep_z]
        if Zd.shape[1]:
            # Project Yd and xd out of Zd via QR
            Q, R = np.linalg.qr(Zd, mode="reduced")
            diag_ok = np.abs(np.diag(R)) > 1e-10 * max(1e-12, float(np.abs(np.diag(R)).max()))
            Q = Q[:, diag_ok]
            xd = xd - Q @ (Q.T @ xd)
            Yd = Yd - Q @ (Q.T @ Yd)
        return (xd @ Yd) / float(xd @ xd)

    # Full-sample estimate (on province-ordered arrays)
    b_full = solve_one(Y_mat[order], x[order], Z[order])
    B = boot_idx.shape[0]
    b_boot = np.empty((B, M), dtype=float)
    for b in range(B):
        ib = boot_idx[b]
        b_boot[b] = solve_one(Y_mat[ib], x[ib], Z[ib])
    return b_full, b_boot


def fast_ols_fe_prov_boot(
    Y_mat: np.ndarray,
    x: np.ndarray,
    Z: np.ndarray,
    idx_by_p: list[np.ndarray],
    picks: np.ndarray,
) -> np.ndarray:
    """Province-cluster bootstrap for FE-OLS (each drawn province gets its own FE block, matching V1 C3)."""
    B, G = picks.shape
    M = Y_mat.shape[1]
    b_boot = np.empty((B, M), dtype=float)
    for b in range(B):
        pick = picks[b]
        A_parts = []
        for g_idx in pick:
            ii = idx_by_p[g_idx]
            Ab = np.column_stack([Y_mat[ii], x[ii, None], Z[ii]])
            Ab -= Ab.mean(axis=0, keepdims=True)
            A_parts.append(Ab)
        A_dm = np.vstack(A_parts)
        Yd = A_dm[:, :M]
        xd = A_dm[:, M]
        Zd = A_dm[:, M + 1:]
        keep_z = np.std(Zd, axis=0) > 1e-12
        Zd = Zd[:, keep_z]
        if Zd.shape[1]:
            Q, R = np.linalg.qr(Zd, mode="reduced")
            diag_ok = np.abs(np.diag(R)) > 1e-10 * max(1e-12, float(np.abs(np.diag(R)).max()))
            Q = Q[:, diag_ok]
            xd = xd - Q @ (Q.T @ xd)
            Yd = Yd - Q @ (Q.T @ Yd)
        b_boot[b] = (xd @ Yd) / float(xd @ xd)
    return b_boot


def main():
    log("=== Analysis D: The Receiving End (Emptied Villages) ===")
    nd = load_fixed_nodes_and_housing()
    panel_all2435, night_meta = build_all_cached_zone_panel(nd)
    night_meta.to_csv(OUT / "d_night_manifest.csv", index=False)

    # Filter panel to fixed R1 eligible sample (n = 2,376 zones, 47 provinces)
    nd_el = nd[nd["eligible"]].copy().reset_index(drop=True)
    el_ids = nd_el["id"].tolist()
    id2pos = {z: k for k, z in enumerate(el_ids)}
    pop_el = nd_el["population"].to_numpy(dtype=float)

    panel = panel_all2435[panel_all2435["id"].isin(id2pos)].copy()
    panel["pop"] = panel["id"].map(nd_el.set_index("id")["population"])
    panel["loss_q"] = panel["id"].map(nd_el.set_index("id")["loss_q"]).astype(int)
    panel["is_single_mun"] = panel["id"].map(nd_el.set_index("id")["n_municipalities"] == 1)
    panel["emptied_inland_50"] = panel["id"].map(nd_el.set_index("id")["emptied_inland_50"])
    panel["arrivals_per_res"] = panel["s_in"] / panel["pop"]
    panel["departures_per_res"] = panel["s_out"] / panel["pop"]
    panel["net_per_res"] = (panel["s_in"] - panel["s_out"]) / panel["pop"]
    tot_io = panel["s_in"] + panel["s_out"]
    panel["net_bal"] = np.where(tot_io > 0, (panel["s_in"] - panel["s_out"]) / tot_io, np.nan)
    panel["net_importer"] = (panel["s_in"] > panel["s_out"]).astype(float)

    # Save heavy 306-night x 2,376-zone panel to _local/journal_data/BD/
    panel.to_parquet(JOURNAL_DATA_BD / "d_zone_daily_panel.parquet", index=False)
    log(f"Saved heavy per-zone daily panel ({len(panel):,} rows) -> {JOURNAL_DATA_BD / 'd_zone_daily_panel.parquet'}")

    # Build fast (n_dates, 2376) matrices for vectorized summaries and window means
    all_dates = sorted(night_meta["date"].unique())
    d2row = {d: i for i, d in enumerate(all_dates)}
    panel_sorted = panel.sort_values(["date", "id"])
    # Verify exact alignment with el_ids (sorted by id!)
    assert el_ids == sorted(el_ids)
    N_D, N_Z = len(all_dates), len(el_ids)
    SIN = panel_sorted["s_in"].to_numpy(dtype=float).reshape(N_D, N_Z)
    SOUT = panel_sorted["s_out"].to_numpy(dtype=float).reshape(N_D, N_Z)
    ARR = SIN / pop_el[None, :]
    DEP = SOUT / pop_el[None, :]
    NET = ARR - DEP
    TOT_IO = SIN + SOUT
    NBAL = np.where(TOT_IO > 0, (SIN - SOUT) / TOT_IO, np.nan)
    IMP = (SIN > SOUT).astype(float)

    # =========================================================================
    # D.1: Quintile summary across ALL 306 cached nights (all & single_mun)
    # =========================================================================
    log("=== D.1: Quintile receiver share & arrivals/departures across all 306 cached nights ===")
    loss_q = nd_el["loss_q"].to_numpy(dtype=int)
    single_mask = (nd_el["n_municipalities"] == 1).to_numpy(dtype=bool)
    emp50_mask = nd_el["emptied_inland_50"].to_numpy(dtype=bool)
    meta_map = night_meta.set_index("date").to_dict("index")

    q_all_rows = []
    for sname, smask in (("all", np.ones(N_Z, bool)), ("single_mun", single_mask)):
        cls_masks = [(f"Q{q}", smask & (loss_q == q)) for q in (1, 2, 3, 4, 5)]
        cls_masks.append(("emptied_inland_loss_ge_0.5", smask & emp50_mask))
        for cls_name, cmask in cls_masks:
            n_c = int(cmask.sum())
            imp_c = IMP[:, cmask]
            arr_c = ARR[:, cmask]
            dep_c = DEP[:, cmask]
            net_c = NET[:, cmask]
            nbal_c = NBAL[:, cmask]
            sh_imp = imp_c.mean(axis=1)
            arr_mean, arr_med = arr_c.mean(axis=1), np.median(arr_c, axis=1)
            arr_q25, arr_q75 = np.quantile(arr_c, 0.25, axis=1), np.quantile(arr_c, 0.75, axis=1)
            dep_mean, dep_med = dep_c.mean(axis=1), np.median(dep_c, axis=1)
            dep_q25, dep_q75 = np.quantile(dep_c, 0.25, axis=1), np.quantile(dep_c, 0.75, axis=1)
            net_mean, net_med = net_c.mean(axis=1), np.median(net_c, axis=1)
            nbal_med = np.nanmedian(nbal_c, axis=1)
            for ti, d in enumerate(all_dates):
                mm = meta_map[d]
                q_all_rows.append({
                    "date": d, "year": mm["year"], "month": mm["month"],
                    "weekday_file": mm["weekday_file"], "night_class": mm["night_class"],
                    "is_dana": mm["is_dana"],
                    "sample": sname, "cls": cls_name, "n": n_c,
                    "share_net_importer": float(sh_imp[ti]),
                    "mean_arrivals_per_res": float(arr_mean[ti]),
                    "median_arrivals_per_res": float(arr_med[ti]),
                    "q25_arrivals_per_res": float(arr_q25[ti]),
                    "q75_arrivals_per_res": float(arr_q75[ti]),
                    "mean_departures_per_res": float(dep_mean[ti]),
                    "median_departures_per_res": float(dep_med[ti]),
                    "q25_departures_per_res": float(dep_q25[ti]),
                    "q75_departures_per_res": float(dep_q75[ti]),
                    "mean_net_per_res": float(net_mean[ti]),
                    "median_net_per_res": float(net_med[ti]),
                    "median_net_bal": float(nbal_med[ti]),
                })
    df_q_all = pd.DataFrame(q_all_rows)
    df_q_all.to_csv(OUT / "d_quintiles_all_nights.csv", index=False)
    log(f"Saved D.1 all-nights quintile summary ({len(df_q_all)} rows) -> out/d_quintiles_all_nights.csv")

    # Build named night/window profiles (s_in, s_out, arr, dep, net, nbal, imp) for key comparisons
    NW = pd.read_csv(NET_DIR / "networks.csv", dtype=str)
    feb24_sis = sorted(NW.loc[(NW.group == "ordinary") & (NW.role == "sister"), "date"])
    aug24_sis = list(C.SISTER_NIGHTS["summer"])
    oct24_clean12 = [d for d in all_dates if d.startswith("2024-10") and pd.Timestamp(d).day_of_week in (1, 2, 3) and d not in DANA_DATES]
    feb22_sis = [d for d in all_dates if d.startswith("2022-02") and pd.Timestamp(d).day_of_week in (1, 2, 3)]
    aug22_sis = [d for d in all_dates if d.startswith("2022-08") and 10 <= int(d[8:10]) <= 18]
    feb23_sis = [d for d in all_dates if d.startswith("2023-02") and pd.Timestamp(d).day_of_week in (1, 2, 3)]
    aug23_sis = [d for d in all_dates if d.startswith("2023-08") and 10 <= int(d[8:10]) <= 18]
    assert len(feb24_sis) == 13 and len(oct24_clean12) == 12 and len(aug24_sis) == 9

    KEY_DEFS = {
        "2024-02-14": ["2024-02-14"],
        "2024-02-mean13": feb24_sis,
        "2024-10-mean12": oct24_clean12,
        "2024-08-15": ["2024-08-15"],
        "2024-08-mean9": aug24_sis,
        "2024-11-01": ["2024-11-01"],
        "2024-12-25": ["2024-12-25"],
        "2024-03-29": ["2024-03-29"],
        "2022-02-16": ["2022-02-16"],
        "2022-08-15": ["2022-08-15"],
        "2022-11-01": ["2022-11-01"],
        "2022-02-mean12": feb22_sis,
        "2022-08-mean9": aug22_sis,
        "2023-02-15": ["2023-02-15"],
        "2023-08-15": ["2023-08-15"],
        "2023-02-mean12": feb23_sis,
        "2023-08-mean9": aug23_sis,
    }
    PROF = {}
    for kname, dlist in KEY_DEFS.items():
        idxs = [d2row[d] for d in dlist]
        s_in = SIN[idxs].mean(axis=0)
        s_out = SOUT[idxs].mean(axis=0)
        arr = s_in / pop_el
        dep = s_out / pop_el
        net = arr - dep
        tot = s_in + s_out
        nbal = np.where(tot > 0, (s_in - s_out) / tot, np.nan)
        imp = (s_in > s_out).astype(float)
        PROF[kname] = {
            "s_in": s_in, "s_out": s_out, "arr": arr, "dep": dep, "net": net, "nbal": nbal, "imp": imp,
        }

    # Save compact per-zone CSV for the main nights (for figure plotting without heavy parquet)
    z_out = nd_el[[
        "id", "name", "prov", "ccaa", "population", "n_municipalities", "is_am", "area_km2",
        "coastal_true", "loss_peak_2021", "loss_pct", "loss_q", "emptied_inland_50", "share_np", "log1p_np_per_res",
    ]].copy()
    for tag, kname in [
        ("feb14", "2024-02-14"), ("oct12", "2024-10-mean12"), ("aug15", "2024-08-15"),
        ("nov01", "2024-11-01"), ("dec25", "2024-12-25"), ("mar29", "2024-03-29"),
    ]:
        p = PROF[kname]
        z_out[f"arr_res_{tag}"] = p["arr"]
        z_out[f"dep_res_{tag}"] = p["dep"]
        z_out[f"net_res_{tag}"] = p["net"]
        z_out[f"imp_{tag}"] = p["imp"].astype(int)
    z_out.to_csv(OUT / "d_zones_main_nights.csv", index=False)
    log("Saved compact main-nights per-zone CSV -> out/d_zones_main_nights.csv")

    # =========================================================================
    # D.2: Within-province zone bootstrap (B=4,999) for quintiles & paired changes
    # =========================================================================
    log("=== D.2: Stratified within-province zone bootstrap (all & single_mun) ===")
    BOOT_NIGHTS = [
        "2024-02-14", "2024-02-mean13", "2024-10-mean12",
        "2024-08-15", "2024-08-mean9", "2024-11-01", "2024-12-25", "2024-03-29",
        "2022-02-16", "2022-08-15", "2022-11-01", "2023-02-15", "2023-08-15",
    ]
    BOOT_PAIRS = [
        ("Aug15_minus_Feb14", "2024-08-15", "2024-02-14"),
        ("Aug15_minus_OctMean12", "2024-08-15", "2024-10-mean12"),
        ("AugMean9_minus_FebMean13", "2024-08-mean9", "2024-02-mean13"),
        ("AugMean9_minus_OctMean12", "2024-08-mean9", "2024-10-mean12"),
        ("Dec25_minus_Feb14", "2024-12-25", "2024-02-14"),
        ("Dec25_minus_OctMean12", "2024-12-25", "2024-10-mean12"),
        ("Nov01_minus_Feb14", "2024-11-01", "2024-02-14"),
        ("Mar29_minus_Feb14", "2024-03-29", "2024-02-14"),
        ("Y2022_Aug15_minus_Feb16", "2022-08-15", "2022-02-16"),
        ("Y2023_Aug15_minus_Feb15", "2023-08-15", "2023-02-15"),
    ]

    q_boot_rows = []
    # Pre-generate bootstrap index arrays for both samples (reused consistently across D.2 & D.3!)
    SAMPLE_MASKS = {"all": np.ones(N_Z, bool), "single_mun": single_mask}
    ZONE_BOOTS = {}
    PROV_BOOTS = {}
    for sname, smask in SAMPLE_MASKS.items():
        provs_s = nd_el.loc[smask, "prov"].to_numpy(dtype=str)
        ZONE_BOOTS[sname] = make_zone_bootstrap_indices(provs_s, N_BOOT_ZONE, rng_for("d_zone_boot", sname))
        PROV_BOOTS[sname] = make_prov_bootstrap_draws(provs_s, N_BOOT_PROV, rng_for("d_prov_boot", sname))

    def boot_quintile_stats(vals_dict: dict[str, np.ndarray], lq_s: np.ndarray, emp_s: np.ndarray,
                            b_idx: np.ndarray, mode_label: str, sample: str, period: str) -> list[dict]:
        """Compute point estimates and stratified within-province zone bootstrap 95% CIs for Q1..Q5, Q5-Q1, emptied_inland."""
        B = b_idx.shape[0]
        lq_b = lq_s[b_idx]  # (B, n_s)
        emp_b = emp_s[b_idx]
        # Precompute class boolean masks in bootstrap draws
        cls_list = [(f"Q{q}", lq_s == q, lq_b == q) for q in (1, 2, 3, 4, 5)]
        cls_list.append(("emptied_inland_loss_ge_0.5", emp_s, emp_b))

        rows = []
        # Store per-class bootstrap arrays for Q5 - Q1 difference
        store_mean = {}
        store_med = {}
        for mname, vec in vals_dict.items():
            vb = vec[b_idx]  # (B, n_s)
            for cname, cm, cmb in cls_list:
                n_c = int(cm.sum())
                cnt_b = cmb.sum(axis=1)
                # Mean
                est_mean = float(np.nanmean(vec[cm]))
                # Vectorized masked mean across B draws
                valid_b = cmb & np.isfinite(vb)
                vcnt = valid_b.sum(axis=1)
                mean_b = np.where(vcnt > 0, np.where(valid_b, vb, 0.0).sum(axis=1) / np.maximum(vcnt, 1), np.nan)
                sm = summarize_boot(est_mean, mean_b)
                store_mean[(mname, cname)] = (est_mean, mean_b)

                # Median (for non-binary metrics)
                if mname != "share_net_importer":
                    est_med = float(np.nanmedian(vec[cm]))
                    # Fast median per row via masked nan
                    vb_nan = np.where(cmb, vb, np.nan)
                    med_b = np.nanmedian(vb_nan, axis=1)
                    smed = summarize_boot(est_med, med_b)
                    store_med[(mname, cname)] = (est_med, med_b)
                else:
                    smed = {"est": np.nan, "se_boot": np.nan, "ci_lo": np.nan, "ci_hi": np.nan, "p_boot": np.nan}

                rows.append({
                    "mode": mode_label, "period": period, "sample": sample, "cls": cname, "n": n_c,
                    "metric": mname,
                    "mean": sm["est"], "mean_se_boot": sm["se_boot"],
                    "mean_ci_lo": sm["ci_lo"], "mean_ci_hi": sm["ci_hi"], "mean_p_boot": sm["p_boot"],
                    "median": smed["est"], "median_se_boot": smed["se_boot"],
                    "median_ci_lo": smed["ci_lo"], "median_ci_hi": smed["ci_hi"], "median_p_boot": smed["p_boot"],
                })

            # Q5 - Q1 difference
            e5, b5 = store_mean[(mname, "Q5")]
            e1, b1 = store_mean[(mname, "Q1")]
            sm_d = summarize_boot(e5 - e1, b5 - b1)
            if mname != "share_net_importer":
                em5, bm5 = store_med[(mname, "Q5")]
                em1, bm1 = store_med[(mname, "Q1")]
                smed_d = summarize_boot(em5 - em1, bm5 - bm1)
            else:
                smed_d = {"est": np.nan, "se_boot": np.nan, "ci_lo": np.nan, "ci_hi": np.nan, "p_boot": np.nan}
            rows.append({
                "mode": mode_label, "period": period, "sample": sample, "cls": "Q5_minus_Q1",
                "n": int((lq_s == 5).sum() + (lq_s == 1).sum()),
                "metric": mname,
                "mean": sm_d["est"], "mean_se_boot": sm_d["se_boot"],
                "mean_ci_lo": sm_d["ci_lo"], "mean_ci_hi": sm_d["ci_hi"], "mean_p_boot": sm_d["p_boot"],
                "median": smed_d["est"], "median_se_boot": smed_d["se_boot"],
                "median_ci_lo": smed_d["ci_lo"], "median_ci_hi": smed_d["ci_hi"], "median_p_boot": smed_d["p_boot"],
            })
        return rows

    for sname, smask in SAMPLE_MASKS.items():
        lq_s = loss_q[smask]
        emp_s = emp50_mask[smask]
        b_idx = ZONE_BOOTS[sname]
        for kname in BOOT_NIGHTS:
            p = PROF[kname]
            vdict = {
                "share_net_importer": p["imp"][smask],
                "arrivals_per_res": p["arr"][smask],
                "departures_per_res": p["dep"][smask],
                "net_per_res": p["net"][smask],
                "net_bal": p["nbal"][smask],
            }
            q_boot_rows.extend(boot_quintile_stats(vdict, lq_s, emp_s, b_idx, "level", sname, kname))
        for pname, h, o in BOOT_PAIRS:
            ph, po = PROF[h], PROF[o]
            vdict = {
                "share_net_importer": (ph["imp"] - po["imp"])[smask],
                "arrivals_per_res": (ph["arr"] - po["arr"])[smask],
                "departures_per_res": (ph["dep"] - po["dep"])[smask],
                "net_per_res": (ph["net"] - po["net"])[smask],
                "net_bal": (ph["nbal"] - po["nbal"])[smask],
            }
            q_boot_rows.extend(boot_quintile_stats(vdict, lq_s, emp_s, b_idx, "change", sname, pname))
        log(f"  D.2 quintile bootstrap done: {sname}")
    df_q_boot = pd.DataFrame(q_boot_rows)
    df_q_boot.to_csv(OUT / "d_quintiles_main_boot.csv", index=False)

    # =========================================================================
    # D.2b: Controlled FE-OLS decomposition of arrivals vs departures per resident
    # =========================================================================
    log("=== D.2b: Controlled FE-OLS decomposition (arrivals vs departures per resident) ===")
    ctrl_rows = []
    for sname, smask in SAMPLE_MASKS.items():
        d_sub = nd_el.loc[smask].reset_index(drop=True)
        x_loss = d_sub["loss_pct"].to_numpy(dtype=float)
        Z_ctrl = d_sub[CONTROLS].to_numpy(dtype=float)
        provs_s = d_sub["prov"].to_numpy(dtype=str)
        b_idx = ZONE_BOOTS[sname]
        idx_by_p, picks = PROV_BOOTS[sname]
        pop_s = d_sub["population"].to_numpy(dtype=float)

        # Build matrix of all level & change outcomes to solve jointly in one vectorized pass!
        col_meta = []
        Y_cols = []
        for kname in ("2024-02-14", "2024-10-mean12", "2024-08-15", "2024-08-mean9",
                      "2024-11-01", "2024-12-25", "2024-03-29", "2022-08-15", "2023-08-15"):
            p = PROF[kname]
            for oname, yv in [
                ("arrivals_per_res", p["arr"][smask]),
                ("departures_per_res", p["dep"][smask]),
                ("net_per_res", p["net"][smask]),
                ("log_arrivals_per_res", np.log((p["s_in"][smask] + 1.0) / pop_s)),
                ("log_departures_per_res", np.log((p["s_out"][smask] + 1.0) / pop_s)),
                ("net_importer", p["imp"][smask]),
                ("net_bal", p["nbal"][smask]),
            ]:
                col_meta.append(("level", kname, oname))
                Y_cols.append(yv)

        for pname, h, o in BOOT_PAIRS:
            ph, po = PROF[h], PROF[o]
            for oname, yv in [
                ("arrivals_per_res_diff", (ph["arr"] - po["arr"])[smask]),
                ("departures_per_res_diff", (ph["dep"] - po["dep"])[smask]),
                ("net_per_res_diff", (ph["net"] - po["net"])[smask]),
                ("log_s_in_ratio", np.log((ph["s_in"][smask] + 1.0) / (po["s_in"][smask] + 1.0))),
                ("log_s_out_ratio", np.log((ph["s_out"][smask] + 1.0) / (po["s_out"][smask] + 1.0))),
                ("net_importer_diff", (ph["imp"] - po["imp"])[smask]),
                ("net_bal_diff", (ph["nbal"] - po["nbal"])[smask]),
            ]:
                col_meta.append(("change", pname, oname))
                Y_cols.append(yv)

        Y_mat = np.column_stack(Y_cols)
        assert np.isfinite(Y_mat).all(), "Non-finite value in controlled decomposition outcomes"
        b_full, b_zone = fast_ols_fe_zone_boot(Y_mat, x_loss, Z_ctrl, provs_s, b_idx)
        b_prov = fast_ols_fe_prov_boot(Y_mat, x_loss, Z_ctrl, idx_by_p, picks)

        # Also run R1's exact CRV1 + Webb WCR on key rows for direct reconciliation with R1!
        for j, (mode_l, per_l, oname) in enumerate(col_meta):
            sz = summarize_boot(b_full[j], b_zone[:, j])
            sp = summarize_boot(b_full[j], b_prov[:, j])
            r_wcb = R1.fe_ols(Y_mat[:, j], x_loss, Z_ctrl, provs_s, provs_s,
                              n_boot=N_WCB_OLS if per_l in ("2024-02-14", "2024-08-15", "2024-12-25", "Aug15_minus_Feb14", "Aug15_minus_OctMean12", "Dec25_minus_Feb14", "Dec25_minus_OctMean12") else 0,
                              seed=int(rng_for("d_wcb", sname, mode_l, per_l, oname).integers(1, 2**31)))
            assert abs(r_wcb["coef"] - b_full[j]) < 1e-10
            ctrl_rows.append({
                "sample": sname, "mode": mode_l, "period": per_l, "outcome": oname,
                "n": len(d_sub), "G": r_wcb["G"],
                "coef": float(b_full[j]),
                "zone_boot_se": sz["se_boot"], "zone_boot_ci_lo": sz["ci_lo"], "zone_boot_ci_hi": sz["ci_hi"], "zone_boot_p": sz["p_boot"],
                "prov_boot_se": sp["se_boot"], "prov_boot_ci_lo": sp["ci_lo"], "prov_boot_ci_hi": sp["ci_hi"], "prov_boot_p": sp["p_boot"],
                "se_crv1": float(r_wcb["se"]), "p_t_crv1": float(r_wcb["p_t"]),
                "p_wcr_prov": float(r_wcb.get("p_wcr", np.nan)),
                "ci_wcu_prov_lo": float(r_wcb.get("ci_wcu_lo", np.nan)),
                "ci_wcu_prov_hi": float(r_wcb.get("ci_wcu_hi", np.nan)),
            })
        log(f"  D.2b controlled decomposition done: {sname}")
    df_ctrl = pd.DataFrame(ctrl_rows)
    # Assert exact linear identity: coef(net_per_res) == coef(arrivals_per_res) - coef(departures_per_res)
    for (sname, per_l), g in df_ctrl[df_ctrl["mode"] == "level"].groupby(["sample", "period"]):
        gm = g.set_index("outcome")["coef"]
        assert abs(gm["net_per_res"] - (gm["arrivals_per_res"] - gm["departures_per_res"])) < 1e-12
    for (sname, per_l), g in df_ctrl[df_ctrl["mode"] == "change"].groupby(["sample", "period"]):
        gm = g.set_index("outcome")["coef"]
        assert abs(gm["net_per_res_diff"] - (gm["arrivals_per_res_diff"] - gm["departures_per_res_diff"])) < 1e-12
    log("Verified exact linear decomposition identity beta(net_per_res) == beta(arrivals_per_res) - beta(departures_per_res)")
    df_ctrl.to_csv(OUT / "d_decomp_controlled.csv", index=False)

    # =========================================================================
    # D.3: Housing slope attenuation (Census 2021 table 59525) — August vs Christmas Eve
    # =========================================================================
    log("=== D.3: Common-sample housing slope attenuation (August vs Christmas Eve) ===")
    HOUSING_SAMPLES = {
        "all": np.ones(N_Z, bool),
        "single_mun": single_mask,
        "inland": ~nd_el["coastal_true"].to_numpy(dtype=bool),
    }
    # Add inland bootstrap indices
    provs_inl = nd_el.loc[HOUSING_SAMPLES["inland"], "prov"].to_numpy(dtype=str)
    ZONE_BOOTS["inland"] = make_zone_bootstrap_indices(provs_inl, N_BOOT_ZONE, rng_for("d_zone_boot", "inland"))
    PROV_BOOTS["inland"] = make_prov_bootstrap_draws(provs_inl, N_BOOT_PROV, rng_for("d_prov_boot", "inland"))

    def fit_attenuation_multi(
        Y_mat: np.ndarray,
        x: np.ndarray,
        h_var: np.ndarray,
        Z: np.ndarray,
        provs: np.ndarray,
        boot_idx: np.ndarray,
        idx_by_p: list[np.ndarray],
        picks: np.ndarray,
    ) -> dict:
        """Compute raw slope, housing-adjusted slope, housing slope, absolute attenuation,
        and relative attenuation share across M outcomes simultaneously on the EXACT same bootstrap draws."""
        n, M = Y_mat.shape
        u_provs = sorted(np.unique(provs))
        order = np.concatenate([np.where(provs == p)[0] for p in u_provs])
        slices = []
        pos = 0
        for p in u_provs:
            m = int(np.sum(provs == p))
            slices.append((pos, pos + m))
            pos += m

        def solve_from_A(A_dm: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
            Yd = A_dm[:, :M]
            xd = A_dm[:, M]
            hd = A_dm[:, M + 1]
            Zd = A_dm[:, M + 2:]
            keep_z = np.std(Zd, axis=0) > 1e-12
            Zd = Zd[:, keep_z]
            if Zd.shape[1]:
                Q, R = np.linalg.qr(Zd, mode="reduced")
                diag_ok = np.abs(np.diag(R)) > 1e-10 * max(1e-12, float(np.abs(np.diag(R)).max()))
                Q = Q[:, diag_ok]
                xd = xd - Q @ (Q.T @ xd)
                hd = hd - Q @ (Q.T @ hd)
                Yd = Yd - Q @ (Q.T @ Yd)
            b_raw = (xd @ Yd) / float(xd @ xd)
            X2 = np.column_stack([xd, hd])
            B2 = np.linalg.solve(X2.T @ X2, X2.T @ Yd)  # (2, M)
            return b_raw, B2[0], B2[1]

        def solve_zone_order(Y_in, x_in, h_in, Z_in):
            A = np.column_stack([Y_in, x_in[:, None], h_in[:, None], Z_in])
            A_dm = A.copy()
            for s, e in slices:
                A_dm[s:e] -= A_dm[s:e].mean(axis=0, keepdims=True)
            return solve_from_A(A_dm)

        b_raw0, b_adj0, b_h0 = solve_zone_order(Y_mat[order], x[order], h_var[order], Z[order])

        Bz = boot_idx.shape[0]
        raw_z, adj_z, h_z = np.empty((Bz, M)), np.empty((Bz, M)), np.empty((Bz, M))
        for b in range(Bz):
            ib = boot_idx[b]
            raw_z[b], adj_z[b], h_z[b] = solve_zone_order(Y_mat[ib], x[ib], h_var[ib], Z[ib])

        Bp = picks.shape[0]
        raw_p, adj_p, h_p = np.empty((Bp, M)), np.empty((Bp, M)), np.empty((Bp, M))
        for b in range(Bp):
            parts = []
            for g_idx in picks[b]:
                ii = idx_by_p[g_idx]
                Ab = np.column_stack([Y_mat[ii], x[ii, None], h_var[ii, None], Z[ii]])
                Ab -= Ab.mean(axis=0, keepdims=True)
                parts.append(Ab)
            raw_p[b], adj_p[b], h_p[b] = solve_from_A(np.vstack(parts))

        return {
            "raw0": b_raw0, "adj0": b_adj0, "h0": b_h0,
            "raw_z": raw_z, "adj_z": adj_z, "h_z": h_z,
            "raw_p": raw_p, "adj_p": adj_p, "h_p": h_p,
        }

    h_rows = []
    h_diff_rows = []
    HOLIDAYS_H = [
        ("Aug15", "2024-08-15"),
        ("Dec25", "2024-12-25"),
        ("Nov01", "2024-11-01"),
        ("Mar29", "2024-03-29"),
    ]
    OUTCOME_NAMES_H = [
        "log_s_in_ratio",
        "arrivals_per_res_diff",
        "departures_per_res_diff",
        "net_per_res_diff",
        "net_bal_diff",
        "net_importer_level",
        "net_importer_diff",
    ]
    for sname, smask in HOUSING_SAMPLES.items():
        d_sub = nd_el.loc[smask].reset_index(drop=True)
        x_loss = d_sub["loss_pct"].to_numpy(dtype=float)
        ctrls_s = [c for c in CONTROLS if not (sname == "inland" and c == "coastal")]
        Z_ctrl = d_sub[ctrls_s].to_numpy(dtype=float)
        provs_s = d_sub["prov"].to_numpy(dtype=str)
        b_idx = ZONE_BOOTS[sname]
        idx_by_p, picks = PROV_BOOTS[sname]

        for bname, o_date in (("vs_Feb14", "2024-02-14"), ("vs_OctMean12", "2024-10-mean12")):
            po = PROF[o_date]
            # Stack all 4 holidays x 7 outcomes = 28 columns on the EXACT same sample!
            col_keys = []
            Y_cols = []
            for hlabel, h_date in HOLIDAYS_H:
                ph = PROF[h_date]
                ymap = {
                    "log_s_in_ratio": np.log((ph["s_in"][smask] + 1.0) / (po["s_in"][smask] + 1.0)),
                    "arrivals_per_res_diff": (ph["arr"] - po["arr"])[smask],
                    "departures_per_res_diff": (ph["dep"] - po["dep"])[smask],
                    "net_per_res_diff": (ph["net"] - po["net"])[smask],
                    "net_bal_diff": (ph["nbal"] - po["nbal"])[smask],
                    "net_importer_level": ph["imp"][smask],
                    "net_importer_diff": (ph["imp"] - po["imp"])[smask],
                }
                for oname in OUTCOME_NAMES_H:
                    col_keys.append((hlabel, oname))
                    Y_cols.append(ymap[oname])
            Y_mat = np.column_stack(Y_cols)
            key2j = {k: j for j, k in enumerate(col_keys)}

            for hproxy in ("share_np", "log1p_np_per_res"):
                if hproxy == "log1p_np_per_res" and (sname != "all" or bname != "vs_Feb14"):
                    continue
                h_vec = d_sub[hproxy].to_numpy(dtype=float)
                res = fit_attenuation_multi(Y_mat, x_loss, h_vec, Z_ctrl, provs_s, b_idx, idx_by_p, picks)

                # Per-holiday attenuation rows
                for j, (hlabel, oname) in enumerate(col_keys):
                    r0, a0, hh0 = float(res["raw0"][j]), float(res["adj0"][j]), float(res["h0"][j])
                    att0 = r0 - a0
                    sh0 = att0 / r0 if abs(r0) > 1e-8 else np.nan
                    att_z = res["raw_z"][:, j] - res["adj_z"][:, j]
                    sh_z = np.where(np.abs(res["raw_z"][:, j]) > 1e-8, att_z / res["raw_z"][:, j], np.nan)
                    att_p = res["raw_p"][:, j] - res["adj_p"][:, j]
                    sh_p = np.where(np.abs(res["raw_p"][:, j]) > 1e-8, att_p / res["raw_p"][:, j], np.nan)

                    s_raw_z = summarize_boot(r0, res["raw_z"][:, j])
                    s_adj_z = summarize_boot(a0, res["adj_z"][:, j])
                    s_h_z = summarize_boot(hh0, res["h_z"][:, j])
                    s_att_z = summarize_boot(att0, att_z)
                    s_sh_z = summarize_boot(sh0, sh_z)
                    s_raw_p = summarize_boot(r0, res["raw_p"][:, j])
                    s_adj_p = summarize_boot(a0, res["adj_p"][:, j])
                    s_att_p = summarize_boot(att0, att_p)
                    s_sh_p = summarize_boot(sh0, sh_p)

                    h_rows.append({
                        "sample": sname, "baseline": bname, "holiday": hlabel,
                        "outcome": oname, "housing_proxy": hproxy,
                        "n": len(d_sub), "G": int(d_sub["prov"].nunique()),
                        "beta_raw": r0,
                        "raw_z_ci_lo": s_raw_z["ci_lo"], "raw_z_ci_hi": s_raw_z["ci_hi"],
                        "raw_p_ci_lo": s_raw_p["ci_lo"], "raw_p_ci_hi": s_raw_p["ci_hi"],
                        "beta_adj_housing": a0,
                        "adj_z_ci_lo": s_adj_z["ci_lo"], "adj_z_ci_hi": s_adj_z["ci_hi"],
                        "adj_p_ci_lo": s_adj_p["ci_lo"], "adj_p_ci_hi": s_adj_p["ci_hi"],
                        "beta_housing": hh0,
                        "housing_z_ci_lo": s_h_z["ci_lo"], "housing_z_ci_hi": s_h_z["ci_hi"],
                        "attenuation_abs": att0,
                        "atten_abs_z_ci_lo": s_att_z["ci_lo"], "atten_abs_z_ci_hi": s_att_z["ci_hi"], "atten_abs_z_p": s_att_z["p_boot"],
                        "atten_abs_p_ci_lo": s_att_p["ci_lo"], "atten_abs_p_ci_hi": s_att_p["ci_hi"], "atten_abs_p_p": s_att_p["p_boot"],
                        "attenuation_share": sh0,
                        "atten_share_z_ci_lo": s_sh_z["ci_lo"], "atten_share_z_ci_hi": s_sh_z["ci_hi"],
                        "atten_share_p_ci_lo": s_sh_p["ci_lo"], "atten_share_p_ci_hi": s_sh_p["ci_hi"],
                    })

                # Formal paired August 15 minus Christmas Eve differences on the EXACT same common sample & draws!
                for oname in OUTCOME_NAMES_H:
                    ja = key2j[("Aug15", oname)]
                    jx = key2j[("Dec25", oname)]
                    ra0, aa0, ha0 = float(res["raw0"][ja]), float(res["adj0"][ja]), float(res["h0"][ja])
                    rx0, ax0, hx0 = float(res["raw0"][jx]), float(res["adj0"][jx]), float(res["h0"][jx])
                    atta0, attx0 = ra0 - aa0, rx0 - ax0
                    sha0 = atta0 / ra0 if abs(ra0) > 1e-8 else np.nan
                    shx0 = attx0 / rx0 if abs(rx0) > 1e-8 else np.nan

                    atta_z = res["raw_z"][:, ja] - res["adj_z"][:, ja]
                    attx_z = res["raw_z"][:, jx] - res["adj_z"][:, jx]
                    atta_p = res["raw_p"][:, ja] - res["adj_p"][:, ja]
                    attx_p = res["raw_p"][:, jx] - res["adj_p"][:, jx]

                    sha_z = np.where(np.abs(res["raw_z"][:, ja]) > 1e-8, atta_z / res["raw_z"][:, ja], np.nan)
                    shx_z = np.where(np.abs(res["raw_z"][:, jx]) > 1e-8, attx_z / res["raw_z"][:, jx], np.nan)
                    sha_p = np.where(np.abs(res["raw_p"][:, ja]) > 1e-8, atta_p / res["raw_p"][:, ja], np.nan)
                    shx_p = np.where(np.abs(res["raw_p"][:, jx]) > 1e-8, attx_p / res["raw_p"][:, jx], np.nan)

                    sd_abs_z = summarize_boot(atta0 - attx0, atta_z - attx_z)
                    sd_abs_p = summarize_boot(atta0 - attx0, atta_p - attx_p)
                    sd_sh_z = summarize_boot(sha0 - shx0, sha_z - shx_z)
                    sd_sh_p = summarize_boot(sha0 - shx0, sha_p - shx_p)
                    sd_raw_z = summarize_boot(ra0 - rx0, res["raw_z"][:, ja] - res["raw_z"][:, jx])
                    sd_raw_p = summarize_boot(ra0 - rx0, res["raw_p"][:, ja] - res["raw_p"][:, jx])
                    sd_adj_z = summarize_boot(aa0 - ax0, res["adj_z"][:, ja] - res["adj_z"][:, jx])
                    sd_h_z = summarize_boot(ha0 - hx0, res["h_z"][:, ja] - res["h_z"][:, jx])

                    h_diff_rows.append({
                        "sample": sname, "baseline": bname, "outcome": oname, "housing_proxy": hproxy,
                        "n_common": len(d_sub), "G": int(d_sub["prov"].nunique()),
                        "beta_raw_aug": ra0, "beta_raw_xmas": rx0,
                        "raw_diff_aug_minus_xmas": sd_raw_z["est"],
                        "raw_diff_z_ci_lo": sd_raw_z["ci_lo"], "raw_diff_z_ci_hi": sd_raw_z["ci_hi"],
                        "raw_diff_p_ci_lo": sd_raw_p["ci_lo"], "raw_diff_p_ci_hi": sd_raw_p["ci_hi"],
                        "beta_adj_aug": aa0, "beta_adj_xmas": ax0,
                        "adj_diff_aug_minus_xmas": sd_adj_z["est"],
                        "adj_diff_z_ci_lo": sd_adj_z["ci_lo"], "adj_diff_z_ci_hi": sd_adj_z["ci_hi"],
                        "beta_housing_aug": ha0, "beta_housing_xmas": hx0,
                        "housing_diff_aug_minus_xmas": sd_h_z["est"],
                        "housing_diff_z_ci_lo": sd_h_z["ci_lo"], "housing_diff_z_ci_hi": sd_h_z["ci_hi"],
                        "atten_abs_aug": atta0, "atten_abs_xmas": attx0,
                        "atten_abs_diff_aug_minus_xmas": sd_abs_z["est"],
                        "atten_abs_diff_z_ci_lo": sd_abs_z["ci_lo"], "atten_abs_diff_z_ci_hi": sd_abs_z["ci_hi"], "atten_abs_diff_z_p": sd_abs_z["p_boot"],
                        "atten_abs_diff_p_ci_lo": sd_abs_p["ci_lo"], "atten_abs_diff_p_ci_hi": sd_abs_p["ci_hi"], "atten_abs_diff_p_p": sd_abs_p["p_boot"],
                        "atten_share_aug": sha0, "atten_share_xmas": shx0,
                        "atten_share_diff_aug_minus_xmas": sd_sh_z["est"],
                        "atten_share_diff_z_ci_lo": sd_sh_z["ci_lo"], "atten_share_diff_z_ci_hi": sd_sh_z["ci_hi"], "atten_share_diff_z_p": sd_sh_z["p_boot"],
                        "atten_share_diff_p_ci_lo": sd_sh_p["ci_lo"], "atten_share_diff_p_ci_hi": sd_sh_p["ci_hi"], "atten_share_diff_p_p": sd_sh_p["p_boot"],
                    })
        log(f"  D.3 housing attenuation done: {sname}")

    df_h = pd.DataFrame(h_rows)
    df_h_diff = pd.DataFrame(h_diff_rows)
    df_h.to_csv(OUT / "d_housing_attenuation.csv", index=False)
    df_h_diff.to_csv(OUT / "d_housing_aug_vs_xmas_diff.csv", index=False)

    # Cross-check exact reproduction of V1 C3 MITMA log_s_in_ratio numbers
    v1_c3 = pd.read_csv(MAIN_ROOT / "results" / "v1_c3_mediation.csv").set_index("spec")
    for hlabel, v1_key in [
        ("Aug15", "MITMA_Aug15_vs_Feb14_log_sin"),
        ("Nov01", "MITMA_Nov01_vs_Feb14_log_sin"),
        ("Dec25", "MITMA_Dec25_vs_Feb14_log_sin"),
        ("Mar29", "MITMA_Mar29_vs_Feb14_log_sin"),
    ]:
        row = df_h[(df_h["sample"] == "all") & (df_h["baseline"] == "vs_Feb14") & (df_h["holiday"] == hlabel)
                   & (df_h["outcome"] == "log_s_in_ratio") & (df_h["housing_proxy"] == "share_np")].iloc[0]
        assert abs(row["beta_raw"] - float(v1_c3.loc[v1_key, "beta_base"])) < 1e-10
        assert abs(row["beta_adj_housing"] - float(v1_c3.loc[v1_key, "beta_after_med"])) < 1e-10
        assert abs(row["attenuation_share"] - float(v1_c3.loc[v1_key, "decomp_share"])) < 1e-10
    log("Verified D.3 housing attenuation slopes & shares reproduce V1 C3 to < 1e-10 on all 4 holidays")

    # =========================================================================
    # Tidy plotting table for Figure "The receiving end"
    # =========================================================================
    fig_d = []
    # Panel B: Receiver share by loss quintile across key nights (all & single_mun)
    for _, r in df_q_boot[
        (df_q_boot["mode"] == "level") & (df_q_boot["metric"] == "share_net_importer")
        & (df_q_boot["cls"].isin(["Q1", "Q2", "Q3", "Q4", "Q5", "Q5_minus_Q1"]))
        & (df_q_boot["period"].isin(["2024-02-14", "2024-10-mean12", "2024-08-15", "2024-11-01", "2024-12-25", "2024-03-29", "2022-08-15", "2023-08-15"]))
    ].iterrows():
        fig_d.append({
            "panel": "B_receiver_share_by_quintile", "sample": r["sample"], "period": r["period"],
            "group": r["cls"], "metric": "share_net_importer", "stat": "mean",
            "point": r["mean"], "ci_lo": r["mean_ci_lo"], "ci_hi": r["mean_ci_hi"], "n": r["n"],
        })
    # Panel C: Arrivals vs departures vs net per resident (quintile levels & changes, means + medians)
    for _, r in df_q_boot[
        (df_q_boot["metric"].isin(["arrivals_per_res", "departures_per_res", "net_per_res"]))
        & (df_q_boot["cls"].isin(["Q1", "Q2", "Q3", "Q4", "Q5", "Q5_minus_Q1"]))
        & (df_q_boot["period"].isin(["2024-02-14", "2024-10-mean12", "2024-08-15", "2024-12-25", "Aug15_minus_Feb14", "Aug15_minus_OctMean12", "Dec25_minus_Feb14", "Dec25_minus_OctMean12"]))
    ].iterrows():
        fig_d.append({
            "panel": f"C_decomp_quintile_{r['mode']}", "sample": r["sample"], "period": r["period"],
            "group": r["cls"], "metric": r["metric"], "stat": "mean",
            "point": r["mean"], "ci_lo": r["mean_ci_lo"], "ci_hi": r["mean_ci_hi"], "n": r["n"],
        })
        fig_d.append({
            "panel": f"C_decomp_quintile_{r['mode']}", "sample": r["sample"], "period": r["period"],
            "group": r["cls"], "metric": r["metric"], "stat": "median",
            "point": r["median"], "ci_lo": r["median_ci_lo"], "ci_hi": r["median_ci_hi"], "n": r["n"],
        })
    # Panel D: Controlled FE-OLS gradients (loss_pct slope for arrivals/departures/net per resident)
    for _, r in df_ctrl[
        (df_ctrl["outcome"].isin(["arrivals_per_res", "departures_per_res", "net_per_res",
                                  "arrivals_per_res_diff", "departures_per_res_diff", "net_per_res_diff",
                                  "net_importer", "log_s_in_ratio"]))
        & (df_ctrl["period"].isin(["2024-02-14", "2024-10-mean12", "2024-08-15", "2024-12-25",
                                   "Aug15_minus_Feb14", "Aug15_minus_OctMean12", "Dec25_minus_Feb14", "Dec25_minus_OctMean12"]))
    ].iterrows():
        fig_d.append({
            "panel": f"D_controlled_gradient_{r['mode']}", "sample": r["sample"], "period": r["period"],
            "group": "loss_pct_0_to_1", "metric": r["outcome"], "stat": "fe_ols_zone_boot",
            "point": r["coef"], "ci_lo": r["zone_boot_ci_lo"], "ci_hi": r["zone_boot_ci_hi"], "n": r["n"],
        })
    # Panel SI: Housing attenuation Aug vs Christmas
    for _, r in df_h[
        (df_h["housing_proxy"] == "share_np") & (df_h["holiday"].isin(["Aug15", "Dec25", "Nov01", "Mar29"]))
        & (df_h["outcome"].isin(["log_s_in_ratio", "arrivals_per_res_diff", "net_per_res_diff", "net_bal_diff"]))
    ].iterrows():
        fig_d.append({
            "panel": "SI_housing_attenuation", "sample": r["sample"], "period": f"{r['holiday']}_{r['baseline']}",
            "group": r["outcome"], "metric": "attenuation_abs", "stat": "zone_boot",
            "point": r["attenuation_abs"], "ci_lo": r["atten_abs_z_ci_lo"], "ci_hi": r["atten_abs_z_ci_hi"], "n": r["n"],
        })
    pd.DataFrame(fig_d).to_csv(OUT / "fig_receiving_end_tidy.csv", index=False)
    log("ALL ANALYSIS D OUTPUTS COMPLETE.")


if __name__ == "__main__":
    main()
