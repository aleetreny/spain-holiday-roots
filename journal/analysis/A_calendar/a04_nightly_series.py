"""a04_nightly_series: Design builders (112-city pooled R2 + 675-city cohort/birthsize V1),
strict reproduction check against coursework results (r2_ppml_nights.csv, r2_diff_tests.csv,
v1_c2_main_contrasts.csv), and nightly PPML series across all valid 2022-2024 nights.

Run:
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
    .venv/bin/python -I \
    journal/analysis/A_calendar/a04_nightly_series.py
"""
from __future__ import annotations

import argparse
import os
import sys
import time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")

import duckdb
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from a00_paths import (  # noqa: E402
    B0,
    C,
    CW_EXT_DIR,
    CW_RESULTS_DIR,
    CW_TT_P84,
    OUT_DIR,
    PANEL_DIR,
    V1,
    ensure_dirs,
    make_logger,
)
from a05_solver import (  # noqa: E402
    FastStackedMultiNight,
    build_equal_year_loadings,
    fit_single_block,
    make_block_design,
)


# =============================================================================
# 1. Static Design Builders (reused by a05_preregistered_tests.py)
# =============================================================================
def build_pooled_meta(
    nodes: pd.DataFrame,
    roots_year: int = 2024,
    income_year: int = 2024,
    use_road_time: bool = False,
    drop_orig_provs: set[str] | None = None,
    drop_dest_provs: set[str] | None = None,
    only_dest_provs: set[str] | None = None,
) -> dict:
    """Exact R2 137-city / 112-identifying-city pooled PPML static design."""
    PROVS = np.array(sorted(nodes["prov"].unique()), dtype=str)
    assert len(PROVS) == 47
    P2I = {p: k for k, p in enumerate(PROVS)}
    NID = nodes["id"].to_numpy(dtype=str)
    NPI = np.array([P2I[p] for p in nodes["prov"].to_numpy(dtype=str)])
    COASTAL_PROVS = sorted(nodes.loc[nodes["coastal_true"], "prov"].unique())
    INLAND_PROVS = [p for p in PROVS if p not in COASTAL_PROVS]

    # Census 68545 roots
    bm_all = pd.read_csv(
        CW_EXT_DIR / "ine_birth_residence" / "mun_residence_x_prov_birth_censo_2021_2025.csv",
        dtype={"ine_code": str, "prov_birth": str},
    )
    bm = bm_all[bm_all["year"] == roots_year]
    assert bm["persons"].notna().all()
    cities_137 = np.array(sorted(c for c in bm["ine_code"].unique() if c[:2] not in C.EXCLUDED_PROVINCES))
    assert len(cities_137) == 137
    cprov = np.array([c[:2] for c in cities_137])
    pcount = Counter(cprov)
    id_provs = sorted(p for p, k in pcount.items() if k >= 2)
    idrow = np.array([pcount[p] >= 2 for p in cprov])
    assert len(id_provs) == 22 and int(idrow.sum()) == 112

    tot_res = bm.groupby("ine_code")["persons"].sum().reindex(cities_137).to_numpy(dtype=float)
    bm_s = bm[bm["prov_birth"].isin(PROVS)].set_index(["ine_code", "prov_birth"])["persons"]
    BORN = np.array([[bm_s.get((c, p), np.nan) for p in PROVS] for c in cities_137], dtype=float)
    assert np.isfinite(BORN).all()
    OWN = cprov[:, None] == PROVS[None, :]
    MASK_C = ~OWN
    assert (BORN[MASK_C] > 0).all()
    SHARE = BORN / tot_res[:, None]

    # Distances (city to pop-weighted province centroid, EPSG:3042 km)
    pw = nodes["population"].to_numpy(dtype=float)
    nx, ny = nodes["x"].to_numpy(dtype=float), nodes["y"].to_numpy(dtype=float)
    PCX = np.array([np.average(nx[NPI == k], weights=pw[NPI == k]) for k in range(47)])
    PCY = np.array([np.average(ny[NPI == k], weights=pw[NPI == k]) for k in range(47)])
    nid2pos = {z: k for k, z in enumerate(NID)}
    city_node = np.array([nid2pos[c] for c in cities_137])
    KM_C = np.hypot(nx[city_node][:, None] - PCX[None, :], ny[city_node][:, None] - PCY[None, :]) / 1000.0

    # Road time (pop-weighted mean of legacy p84 tt_adj)
    TT_C = None
    if use_road_time:
        tt = pd.read_parquet(CW_TT_P84, columns=["c", "z", "tt_adj"])
        tt = tt[tt["z"].isin(set(NID)) & tt["c"].isin(set(cities_137))].copy()
        c2r = {c: k for k, c in enumerate(cities_137)}
        zpos = pd.Series(np.arange(len(NID)), index=NID)
        tt["pi"] = NPI[zpos.loc[tt["z"]].to_numpy()]
        tt["pop"] = pw[zpos.loc[tt["z"]].to_numpy()]
        tt["ri"] = tt["c"].map(c2r)
        tt["wtt"] = tt["tt_adj"] * tt["pop"]
        g = tt.groupby(["ri", "pi"])[["wtt", "pop"]].sum()
        TT_C = np.full((len(cities_137), 47), np.nan)
        TT_C[g.index.get_level_values(0), g.index.get_level_values(1)] = (g["wtt"] / g["pop"]).to_numpy()
        assert np.isfinite(TT_C[MASK_C]).all()

    # ADRH income
    adrh = pd.read_csv(CW_EXT_DIR / "ine_adrh" / "adrh_municipios.csv", dtype={"ine_code": str})
    inc = adrh[(adrh["year"] == income_year) & (adrh["indicator"] == "Renta neta media por persona")].set_index("ine_code")["value"]
    INC = inc.reindex(cities_137).to_numpy(dtype=float)
    assert np.isfinite(INC).all(), f"Missing ADRH {income_year} income for 137 cities"
    LINC = np.log(INC) - np.log(INC).mean()

    # Dyad mask (identifying cities, off-diagonal, optional origin/dest exclusions)
    mask = MASK_C & idrow[:, None]
    if drop_orig_provs:
        mask &= ~np.isin(cprov, sorted(drop_orig_provs))[:, None]
    if drop_dest_provs:
        mask &= ~np.isin(PROVS, sorted(drop_dest_provs))[None, :]
    if only_dest_provs is not None:
        mask &= np.isin(PROVS, sorted(only_dest_provs))[None, :]

    ri_all, pj_all = np.where(mask)
    cities_dyad = cities_137[ri_all]
    oprov_dyad = cprov[ri_all]
    dprov_dyad = PROVS[pj_all]

    x_lshare = np.log(SHARE[ri_all, pj_all])[:, None]
    dist_mat = TT_C if use_road_time else KM_C
    ld = np.log(dist_mat[ri_all, pj_all])
    ld_c = ld - ld.mean()
    cov_dist = np.column_stack([ld_c, ld_c ** 2])

    c_code = pd.factorize(cities_dyad)[0]
    p_code = pd.factorize(np.char.add(np.char.add(oprov_dyad, "_"), dprov_dyad))[0]
    cl_levels = np.array(sorted(set(oprov_dyad)), dtype=str)
    cl_code = pd.Categorical(oprov_dyad, categories=cl_levels).codes.astype(np.int32)
    slope_raw = [(pj_all.astype(np.int32), LINC[ri_all])]

    c2i = {c: k for k, c in enumerate(cities_137)}
    return {
        "cities_all": cities_137,
        "c2i": c2i,
        "PROVS": PROVS,
        "P2I": P2I,
        "COASTAL_PROVS": COASTAL_PROVS,
        "INLAND_PROVS": INLAND_PROVS,
        "ri_all": ri_all,
        "pj_all": pj_all,
        "cities_dyad": cities_dyad,
        "oprov_dyad": oprov_dyad,
        "dprov_dyad": dprov_dyad,
        "Xr": x_lshare,
        "cov_dist": cov_dist,
        "c_code": c_code,
        "p_code": p_code,
        "cl_code": cl_code,
        "cl_levels": cl_levels,
        "G_tot": len(cl_levels),
        "slope_raw": slope_raw,
        "n_cells": len(ri_all),
        "n_cities_all": len(cities_137),
    }


def build_cohort_meta(
    nodes: pd.DataFrame,
    income_year: int = 2024,
    zero_rule: str = "haldane_0.5",
    restrict_to_112: bool = False,
    drop_orig_provs: set[str] | None = None,
    drop_dest_provs: set[str] | None = None,
    only_dest_provs: set[str] | None = None,
) -> dict:
    """Exact V1 675-city cohort_joint & birthplace_joint static design.

    Arrival codes: ONLY explicit codes enter each group (VARANOP='0001' for pre1981,
    VARANOP IN '2016'..'2023' for recent; special codes '7777','8888' dropped from both).
    """
    PROVS = np.array(sorted(nodes["prov"].unique()), dtype=str)
    P2I = {p: k for k, p in enumerate(PROVS)}
    NID = nodes["id"].to_numpy(dtype=str)
    NPI = np.array([P2I[p] for p in nodes["prov"].to_numpy(dtype=str)])
    COASTAL_PROVS = sorted(nodes.loc[nodes["coastal_true"], "prov"].unique())
    INLAND_PROVS = [p for p in PROVS if p not in COASTAL_PROVS]

    pq = CW_EXT_DIR / "ine_censo_microdatos_2024" / "md_CensoAnual_muestra_2024.parquet"
    con = duckdb.connect()
    spec_muni = "('991','992','993','666','966','777','888','999')"
    m_tot = con.execute(f"""
        SELECT CPRO || CMUN AS mun, CPRO AS prov, count(*) AS n_tot
        FROM read_parquet('{pq}')
        WHERE CMUN NOT IN {spec_muni}
          AND CPRO NOT IN ('07','35','38','51','52')
        GROUP BY 1, 2 ORDER BY 1
    """).df()
    single_set = set(nodes.loc[nodes["n_municipalities"] == 1, "id"])
    m_single = m_tot[m_tot["mun"].isin(single_set)].copy()
    pcounts = m_single.groupby("prov").size()
    id_provs = sorted(pcounts[pcounts >= 2].index)
    m675 = m_single[m_single["prov"].isin(id_provs)].sort_values("mun").reset_index(drop=True)
    assert len(m675) == 675 and len(id_provs) == 44

    cp = con.execute(f"""
        SELECT CPRO || CMUN AS mun, CPRO_NAC AS dprov,
               count(*) AS n_born,
               sum(CASE WHEN VARANOP = '0001' THEN 1 ELSE 0 END) AS n_pre1981,
               sum(CASE WHEN VARANOP IN ('2016','2017','2018','2019','2020','2021','2022','2023') THEN 1 ELSE 0 END) AS n_rec2016,
               sum(CASE WHEN CMUN_NAC = '991' THEN 1 ELSE 0 END) AS n_vil2k,
               sum(CASE WHEN CMUN_NAC NOT IN {spec_muni} AND CMUN_NAC IS NOT NULL AND CMUN_NAC != '' THEN 1 ELSE 0 END) AS n_urb10k
        FROM read_parquet('{pq}')
        WHERE CMUN NOT IN {spec_muni}
          AND CPRO NOT IN ('07','35','38','51','52')
          AND CPRO_NAC BETWEEN '01' AND '52'
          AND CPRO_NAC NOT IN ('07','35','38','51','52')
          AND CPRO != CPRO_NAC
        GROUP BY 1, 2
    """).df()

    MMUN = m675["mun"].to_numpy(dtype=str)
    MPROV = m675["prov"].to_numpy(dtype=str)
    NTOT = m675["n_tot"].to_numpy(dtype=float)[:, None]
    M2R = {m: r for r, m in enumerate(MMUN)}
    C_M = len(MMUN)

    cp675 = cp[cp["mun"].isin(M2R)].copy()
    ri_cp = cp675["mun"].map(M2R).to_numpy(dtype=int)
    pj_cp = cp675["dprov"].map(P2I).to_numpy(dtype=int)

    MATS = {}
    for col in ("n_born", "n_pre1981", "n_rec2016", "n_vil2k", "n_urb10k"):
        M = np.zeros((C_M, 47))
        M[ri_cp, pj_cp] = cp675[col].to_numpy(dtype=float)
        MATS[col] = M

    # Optional restriction to R2's 112 large cities (Table A2 robustness)
    row_keep = np.ones(C_M, dtype=bool)
    if restrict_to_112:
        bm = pd.read_csv(
            CW_EXT_DIR / "ine_birth_residence" / "mun_residence_x_prov_birth_censo_2021_2025.csv",
            dtype={"ine_code": str},
        )
        c137 = sorted(c for c in bm.loc[bm["year"] == 2024, "ine_code"].unique() if c[:2] not in C.EXCLUDED_PROVINCES)
        c137_cnt = Counter(c[:2] for c in c137)
        c112 = {c for c in c137 if c137_cnt[c[:2]] >= 2}
        row_keep = np.isin(MMUN, sorted(c112))
        assert int(row_keep.sum()) == 112

    OWN = MPROV[:, None] == PROVS[None, :]
    mask = (~OWN) & row_keep[:, None]
    if drop_orig_provs:
        mask &= ~np.isin(MPROV, sorted(drop_orig_provs))[:, None]
    if drop_dest_provs:
        mask &= ~np.isin(PROVS, sorted(drop_dest_provs))[None, :]
    if only_dest_provs is not None:
        mask &= np.isin(PROVS, sorted(only_dest_provs))[None, :]

    ri_all, pj_all = np.where(mask)
    cities_dyad = MMUN[ri_all]
    oprov_dyad = MPROV[ri_all]
    dprov_dyad = PROVS[pj_all]

    # Distances & ADRH income
    pw = nodes["population"].to_numpy(dtype=float)
    nx, ny = nodes["x"].to_numpy(dtype=float), nodes["y"].to_numpy(dtype=float)
    PCX = np.array([np.average(nx[NPI == k], weights=pw[NPI == k]) for k in range(47)])
    PCY = np.array([np.average(ny[NPI == k], weights=pw[NPI == k]) for k in range(47)])
    nid2pos = {z: k for k, z in enumerate(NID)}
    MNODE = np.array([nid2pos[m] for m in MMUN])
    KM = np.hypot(nx[MNODE][:, None] - PCX[None, :], ny[MNODE][:, None] - PCY[None, :]) / 1000.0
    LKM = np.log(KM[ri_all, pj_all])
    lkm_c = LKM - LKM.mean()
    cov_dist = np.column_stack([lkm_c, lkm_c ** 2])

    adrh = pd.read_csv(CW_EXT_DIR / "ine_adrh" / "adrh_municipios.csv", dtype={"ine_code": str})
    inc_s = adrh[(adrh["year"] == income_year) & (adrh["indicator"] == "Renta neta media por persona")].set_index("ine_code")["value"]
    INC = inc_s.reindex(MMUN).to_numpy(dtype=float)
    assert np.isfinite(INC).all(), f"Missing ADRH {income_year} income for 675 cities"
    LINC = np.log(INC) - np.log(INC).mean()

    # Zero-cell transforms (Table A2)
    def transform(M: np.ndarray) -> np.ndarray:
        if zero_rule == "haldane_0.5":
            return np.log((M + 0.5) / NTOT)[ri_all, pj_all]
        if zero_rule == "plus_1.0":
            return np.log((M + 1.0) / NTOT)[ri_all, pj_all]
        if zero_rule == "plus_0.1":
            return np.log((10.0 * M + 1.0) / (10.0 * NTOT))[ri_all, pj_all]
        if zero_rule == "share_plus_1e4":
            return np.log(M / NTOT + 1e-4)[ri_all, pj_all]
        if zero_rule == "share_plus_1e5":
            return np.log(M / NTOT + 1e-5)[ri_all, pj_all]
        if zero_rule == "asinh_1000_share":
            return np.arcsinh(1000.0 * (M / NTOT))[ri_all, pj_all]
        raise ValueError(f"Unknown zero_rule: {zero_rule}")

    X_COH = np.column_stack([transform(MATS["n_pre1981"]), transform(MATS["n_rec2016"])])
    X_BS = np.column_stack([transform(MATS["n_vil2k"]), transform(MATS["n_urb10k"])])

    c_code = pd.factorize(cities_dyad)[0]
    p_code = pd.factorize(np.char.add(np.char.add(oprov_dyad, "_"), dprov_dyad))[0]
    cl_levels = np.array(sorted(set(oprov_dyad)), dtype=str)
    cl_code = pd.Categorical(oprov_dyad, categories=cl_levels).codes.astype(np.int32)
    slope_raw = [(pj_all.astype(np.int32), LINC[ri_all])]

    return {
        "cities_all": MMUN,
        "c2i": M2R,
        "PROVS": PROVS,
        "P2I": P2I,
        "COASTAL_PROVS": COASTAL_PROVS,
        "INLAND_PROVS": INLAND_PROVS,
        "ri_all": ri_all,
        "pj_all": pj_all,
        "cities_dyad": cities_dyad,
        "oprov_dyad": oprov_dyad,
        "dprov_dyad": dprov_dyad,
        "X_COH": X_COH,
        "X_BS": X_BS,
        "cov_dist": cov_dist,
        "c_code": c_code,
        "p_code": p_code,
        "cl_code": cl_code,
        "cl_levels": cl_levels,
        "G_tot": len(cl_levels),
        "slope_raw": slope_raw,
        "n_cells": len(ri_all),
        "n_cities_all": C_M,
    }


def load_panel_flows(
    panel_parquet: Path,
    meta: dict,
    dates: list[str] | None = None,
    use_noart: bool = False,
) -> dict[str, dict[str, np.ndarray]]:
    """Load per-date flow vectors aligned to meta['ri_all'], meta['pj_all']."""
    con = duckdb.connect()
    if dates is not None:
        dlist = "[" + ",".join(f"'{d}'" for d in sorted(dates)) + "]"
        where = f"WHERE date IN {dlist}"
    else:
        where = ""
    inl_col = "flow_inland_noart" if use_noart else "flow_inland"
    all_expr = "(flow_inland_noart + flow_coastal)" if use_noart else "flow_all"
    df = con.execute(f"""
        SELECT date, city, dprov,
               {inl_col} AS f_inland,
               flow_coastal AS f_coastal,
               {all_expr} AS f_all
        FROM read_parquet('{panel_parquet}')
        {where}
    """).df()

    c2i = meta["c2i"]
    P2I = meta["P2I"]
    C_n = meta["n_cities_all"]
    ri_all = meta["ri_all"]
    pj_all = meta["pj_all"]
    coast_idx = [P2I[p] for p in meta["COASTAL_PROVS"]]

    df = df[df["city"].isin(c2i) & df["dprov"].isin(P2I)].copy()
    ri = df["city"].map(c2i).to_numpy(dtype=np.int32)
    pj = df["dprov"].map(P2I).to_numpy(dtype=np.int32)
    df["ri"] = ri
    df["pj"] = pj

    flows: dict[str, dict[str, np.ndarray]] = {}
    for d, grp in df.groupby("date", sort=True):
        r_g = grp["ri"].to_numpy()
        p_g = grp["pj"].to_numpy()
        M_all = np.zeros((C_n, 47))
        M_inl = np.zeros((C_n, 47))
        M_cst = np.zeros((C_n, 47))
        M_all[r_g, p_g] = grp["f_all"].to_numpy(dtype=float)
        M_inl[r_g, p_g] = grp["f_inland"].to_numpy(dtype=float)
        M_cst[r_g, p_g] = grp["f_coastal"].to_numpy(dtype=float)
        M_inlp = M_all.copy()
        M_inlp[:, coast_idx] = 0.0
        flows[str(d)] = {
            "all": M_all[ri_all, pj_all],
            "inland": M_inl[ri_all, pj_all],
            "coastal": M_cst[ri_all, pj_all],
            "inlandprov": M_inlp[ri_all, pj_all],
        }
    return flows


# =============================================================================
# 2. Strict Coursework Reproduction Check
# =============================================================================
def run_reproduction_check(nodes: pd.DataFrame) -> pd.DataFrame:
    """Reproduce exact coursework values on 2024-02-14 and 2024-08-15 (plus key nights)
    and cross-check FastStackedMultiNight against r2_diff_tests.csv and v1_c2_main_contrasts.csv.
    Stops with AssertionError if any check fails!
    """
    log = make_logger()
    log("=== STRICT COURSEWORK REPRODUCTION CHECK ===")
    p_meta = build_pooled_meta(nodes)
    c_meta = build_cohort_meta(nodes)

    chk_dates = ["2024-02-14", "2024-08-15", "2024-11-01", "2024-12-25", "2024-03-29", "2022-02-16", "2022-08-15", "2023-02-15", "2023-08-15"]
    p_flows = load_panel_flows(PANEL_DIR / "panel_cities137_prov.parquet", p_meta, dates=chk_dates)
    p_flows_noart = load_panel_flows(PANEL_DIR / "panel_cities137_prov.parquet", p_meta, dates=["2024-02-14", "2024-08-15"], use_noart=True)
    c_flows = load_panel_flows(PANEL_DIR / "panel_cities675_prov.parquet", c_meta, dates=["2024-02-14", "2024-08-15"])

    r2_ref = pd.read_csv(CW_RESULTS_DIR / "r2_ppml_nights.csv")
    r2_diff = pd.read_csv(CW_RESULTS_DIR / "r2_diff_tests.csv")
    v1_ref = pd.read_csv(CW_RESULTS_DIR / "v1_c2_main_contrasts.csv")

    rows = []
    # 1. Pooled R2 single nights
    for d in chk_dates:
        for sample in ("all", "inland", "coastal", "inlandprov"):
            ref_row = r2_ref[(r2_ref["date"] == d) & (r2_ref["sample"] == sample) & (r2_ref["spec"] == "main")]
            if ref_row.empty:
                continue
            blk = make_block_design(
                d, int(d[:4]), p_flows[d][sample], p_meta["Xr"], p_meta["cov_dist"],
                p_meta["c_code"], p_meta["p_code"], p_meta["cl_code"], slope_raw=p_meta["slope_raw"],
            )
            fit = fit_single_block(blk, p_meta["G_tot"])
            b_ref = float(ref_row["beta"].iloc[0])
            se_ref = float(ref_row["se_crv1"].iloc[0])
            b_new = float(fit["beta"][0])
            se_new = float(fit["se"][0])
            assert abs(b_new - b_ref) < 1e-6, f"Pooled mismatch {d} {sample}: {b_new} vs {b_ref}"
            assert abs(se_new - se_ref) < 1e-5, f"Pooled SE mismatch {d} {sample}: {se_new} vs {se_ref}"
            rows.append({
                "check": f"pooled_{sample}_{d}", "ref_file": "r2_ppml_nights.csv",
                "expected": b_ref, "reproduced": b_new, "abs_diff": abs(b_new - b_ref),
                "se_expected": se_ref, "se_reproduced": se_new,
            })

    # Artefact-zone sensitivity check on 2024-02-14 & 2024-08-15
    for d in ("2024-02-14", "2024-08-15"):
        blk_na = make_block_design(
            d, 2024, p_flows_noart[d]["all"], p_meta["Xr"], p_meta["cov_dist"],
            p_meta["c_code"], p_meta["p_code"], p_meta["cl_code"], slope_raw=p_meta["slope_raw"],
        )
        fit_na = fit_single_block(blk_na, p_meta["G_tot"])
        b_ref = float(r2_ref[(r2_ref["date"] == d) & (r2_ref["sample"] == "all") & (r2_ref["spec"] == "main")]["beta"].iloc[0])
        log(f"  Artefact-zone dest check {d}: all2435={b_ref:.6f}, drop3art={fit_na['beta'][0]:.6f}, diff={fit_na['beta'][0] - b_ref:+.2e}")

    # 2. Pooled stacked diff + WCR check on 2024-08-15 vs 2024-02-14
    blk_o = make_block_design("2024-02-14", 2024, p_flows["2024-02-14"]["all"], p_meta["Xr"], p_meta["cov_dist"],
                              p_meta["c_code"], p_meta["p_code"], p_meta["cl_code"], slope_raw=p_meta["slope_raw"])
    blk_h = make_block_design("2024-08-15", 2024, p_flows["2024-08-15"]["all"], p_meta["Xr"], p_meta["cov_dist"],
                              p_meta["c_code"], p_meta["p_code"], p_meta["cl_code"], slope_raw=p_meta["slope_raw"])
    blk_o.mu_single = fit_single_block(blk_o, p_meta["G_tot"])["mu"]
    blk_h.mu_single = fit_single_block(blk_h, p_meta["G_tot"])["mu"]
    l_tgt, L_oth = build_equal_year_loadings([blk_o, blk_h], [False, True], K=1, target_vec=np.array([1.0]))
    st = FastStackedMultiNight([blk_o, blk_h], l_tgt, L_oth, p_meta["G_tot"])
    u = st.fit_unrestricted()
    # Exact R2 seed for VB check!
    vb_r2 = V1.rng_for("r2_webb").choice(V1.WEBB, size=(9999, 22))
    p_wcr, _, _ = st.p_wcr(0.0, vb_r2, mu_start=u["mus"])
    drow = r2_diff[(r2_diff["sample"] == "all") & (r2_diff["spec"] == "main") & (r2_diff["holiday"] == "2024-08-15") & (r2_diff["ordinary"] == "2024-02-14")].iloc[0]
    assert abs(u["est"] - float(drow["delta"])) < 1e-6
    assert abs(u["se_crv1"] - float(drow["se_crv1"])) < 1e-5
    log(f"  Stacked pooled Aug15-Feb14: est={u['est']:.6f} (ref {drow['delta']:.6f}), "
        f"se={u['se_crv1']:.6f} (ref {drow['se_crv1']:.6f}), p_wcr={p_wcr:.4f} (ref {drow['p_wcr_delta0']:.4f})")
    rows.append({
        "check": "stacked_pooled_Aug15_vs_Feb14", "ref_file": "r2_diff_tests.csv",
        "expected": float(drow["delta"]), "reproduced": u["est"], "abs_diff": abs(u["est"] - float(drow["delta"])),
        "se_expected": float(drow["se_crv1"]), "se_reproduced": u["se_crv1"],
    })

    # 3. Cohort & Birthsize reproduction on 2024-02-14 & 2024-08-15
    vb_v1 = V1.rng_for("v1_webb_9999").choice(V1.WEBB, size=(9999, 44))
    for model_name, X_mat, terms in [
        ("cohort_joint", c_meta["X_COH"], ["pre1981_exodus", "recent_2016_2024"]),
        ("birthsize_joint", c_meta["X_BS"], ["born_village_le2k", "born_urban_gt10k"]),
    ]:
        bo = make_block_design("2024-02-14", 2024, c_flows["2024-02-14"]["all"], X_mat, c_meta["cov_dist"],
                               c_meta["c_code"], c_meta["p_code"], c_meta["cl_code"], slope_raw=c_meta["slope_raw"])
        bh = make_block_design("2024-08-15", 2024, c_flows["2024-08-15"]["all"], X_mat, c_meta["cov_dist"],
                               c_meta["c_code"], c_meta["p_code"], c_meta["cl_code"], slope_raw=c_meta["slope_raw"])
        fo = fit_single_block(bo, c_meta["G_tot"])
        fh = fit_single_block(bh, c_meta["G_tot"])
        bo.mu_single, bh.mu_single = fo["mu"], fh["mu"]
        for k_i, tnm in enumerate(terms):
            vrow = v1_ref[(v1_ref["contrast"] == "Aug15_vs_Feb14") & (v1_ref["sample"] == "all") & (v1_ref["model"] == model_name) & (v1_ref["term"] == tnm)].iloc[0]
            assert abs(fo["beta"][k_i] - float(vrow["beta_ordinary"])) < 1e-6
            assert abs(fh["beta"][k_i] - float(vrow["beta_holiday"])) < 1e-6
            tvec = np.eye(2)[k_i]
            ovec = [np.eye(2)[1 - k_i]]
            lt, Lo = build_equal_year_loadings([bo, bh], [False, True], K=2, target_vec=tvec, other_vecs=ovec)
            stc = FastStackedMultiNight([bo, bh], lt, Lo, c_meta["G_tot"])
            uc = stc.fit_unrestricted()
            pw, _, _ = stc.p_wcr(0.0, vb_v1, mu_start=uc["mus"])
            assert abs(uc["est"] - float(vrow["delta"])) < 1e-6
            assert abs(uc["se_crv1"] - float(vrow["se_crv1"])) < 1e-5
            assert abs(pw - float(vrow["p_wcr_9999"])) < 1e-9, f"Exact WCR p mismatch {tnm}: {pw} vs {vrow['p_wcr_9999']}"
            log(f"  Cohort {model_name}/{tnm}: Feb={fo['beta'][k_i]:.6f} (ref {vrow['beta_ordinary']:.6f}), "
                f"Aug={fh['beta'][k_i]:.6f} (ref {vrow['beta_holiday']:.6f}), "
                f"delta={uc['est']:+.6f}, p_wcr={pw:.4f} (EXACT ref {vrow['p_wcr_9999']:.4f})")
            rows.append({
                "check": f"{model_name}_{tnm}_Feb14", "ref_file": "v1_c2_main_contrasts.csv",
                "expected": float(vrow["beta_ordinary"]), "reproduced": float(fo["beta"][k_i]),
                "abs_diff": abs(float(fo["beta"][k_i]) - float(vrow["beta_ordinary"])),
                "se_expected": np.nan, "se_reproduced": float(fo["se"][k_i]),
            })
            rows.append({
                "check": f"{model_name}_{tnm}_Aug15", "ref_file": "v1_c2_main_contrasts.csv",
                "expected": float(vrow["beta_holiday"]), "reproduced": float(fh["beta"][k_i]),
                "abs_diff": abs(float(fh["beta"][k_i]) - float(vrow["beta_holiday"])),
                "se_expected": np.nan, "se_reproduced": float(fh["se"][k_i]),
            })
            rows.append({
                "check": f"stacked_{model_name}_{tnm}_Aug15_vs_Feb14", "ref_file": "v1_c2_main_contrasts.csv",
                "expected": float(vrow["delta"]), "reproduced": uc["est"],
                "abs_diff": abs(uc["est"] - float(vrow["delta"])),
                "se_expected": float(vrow["se_crv1"]), "se_reproduced": uc["se_crv1"],
            })

    out_df = pd.DataFrame(rows)
    out_csv = OUT_DIR / "reproduction_check.csv"
    out_df.to_csv(out_csv, index=False)
    log(f"ALL REPRODUCTION CHECKS PASSED (max |diff| = {out_df['abs_diff'].max():.2e})! Saved {out_csv}")
    return out_df


# =============================================================================
# 3. Parallel Workers for Nightly Series
# =============================================================================
_G_META = None


def _init_worker(meta: dict):
    global _G_META
    _G_META = meta


def _fit_pooled_chunk(items: list[tuple[str, dict[str, np.ndarray]]]) -> list[dict]:
    meta = _G_META
    G_tot = meta["G_tot"]
    rows = []
    for d, fdict in items:
        yr = int(d[:4])
        for sample in ("all", "inland", "coastal", "inlandprov"):
            blk = make_block_design(
                d, yr, fdict[sample], meta["Xr"], meta["cov_dist"],
                meta["c_code"], meta["p_code"], meta["cl_code"], slope_raw=meta["slope_raw"],
            )
            fit = fit_single_block(blk, G_tot)
            b = float(fit["beta"][0])
            se = float(fit["se"][0])
            G = int(fit["G"])
            tc = float(stats.t.ppf(0.975, G - 1))
            p_t = float(2.0 * stats.t.sf(abs(b / se), G - 1))
            rows.append({
                "date": d, "sample": sample, "beta": b, "se_crv1": se, "G": G,
                "ci_lo": b - tc * se, "ci_hi": b + tc * se, "p_crv1": p_t,
                "n_obs": fit["n"], "flow_total": fit["flow_total"],
            })
    return rows


def _fit_cohort_chunk(items: list[tuple[str, dict[str, np.ndarray]]]) -> list[dict]:
    meta = _G_META
    G_tot = meta["G_tot"]
    rows = []
    for d, fdict in items:
        yr = int(d[:4])
        for sample in ("all", "inland", "coastal"):
            for mname, X_mat, tnames in (
                ("cohort_joint", meta["X_COH"], ("pre1981", "recent_2016_2023")),
                ("birthplace_joint", meta["X_BS"], ("born_le_2000", "born_gt_10000")),
            ):
                blk = make_block_design(
                    d, yr, fdict[sample], X_mat, meta["cov_dist"],
                    meta["c_code"], meta["p_code"], meta["cl_code"], slope_raw=meta["slope_raw"],
                )
                fit = fit_single_block(blk, G_tot)
                G = int(fit["G"])
                tc = float(stats.t.ppf(0.975, G - 1))
                for k_i, tnm in enumerate(tnames):
                    b = float(fit["beta"][k_i])
                    se = float(fit["se"][k_i])
                    p_t = float(2.0 * stats.t.sf(abs(b / se), G - 1))
                    rows.append({
                        "date": d, "sample": sample, "model": mname, "term": tnm,
                        "beta": b, "se_crv1": se, "G": G,
                        "ci_lo": b - tc * se, "ci_hi": b + tc * se, "p_crv1": p_t,
                        "n_obs": fit["n"], "flow_total": fit["flow_total"],
                    })
    return rows


def run_nightly_series(nodes: pd.DataFrame, workers: int = 4) -> None:
    log = make_logger()
    man = pd.read_csv(OUT_DIR / "night_manifest.csv")
    valid_dates = sorted(man.loc[man["is_valid"], "file_date"].tolist())
    meta_cols = man.loc[man["is_valid"], ["file_date", "evening_date", "year", "evening_weekday", "split", "night_class"]].rename(
        columns={"file_date": "date"}
    )

    # --- 1. Pooled R2 nightly series ---
    t1 = time.time()
    p_meta = build_pooled_meta(nodes)
    p_flows = load_panel_flows(PANEL_DIR / "panel_cities137_prov.parquet", p_meta, dates=valid_dates)
    items_p = [(d, p_flows[d]) for d in valid_dates]
    n_chunks = max(1, min(workers * 4, len(items_p)))
    chunks_p = [items_p[i::n_chunks] for i in range(n_chunks)]
    log(f"Estimating pooled R2 nightly series on {len(valid_dates)} valid nights ({workers} workers)...")
    p_rows = []
    if workers <= 1:
        _init_worker(p_meta)
        for ch in chunks_p:
            p_rows.extend(_fit_pooled_chunk(ch))
    else:
        with ProcessPoolExecutor(max_workers=workers, initializer=_init_worker, initargs=(p_meta,)) as ex:
            for res in ex.map(_fit_pooled_chunk, chunks_p):
                p_rows.extend(res)
    df_p = pd.DataFrame(p_rows).merge(meta_cols, on="date", how="left").sort_values(["sample", "date"]).reset_index(drop=True)
    df_p.to_csv(OUT_DIR / "nightly_pooled.csv", index=False)
    log(f"Saved nightly_pooled.csv ({len(df_p):,} rows, {time.time() - t1:.1f}s)")

    # --- 2. Cohort + Birthsize nightly series ---
    t2 = time.time()
    c_meta = build_cohort_meta(nodes)
    c_flows = load_panel_flows(PANEL_DIR / "panel_cities675_prov.parquet", c_meta, dates=valid_dates)
    items_c = [(d, c_flows[d]) for d in valid_dates]
    chunks_c = [items_c[i::n_chunks] for i in range(n_chunks)]
    log(f"Estimating 675-city cohort & birthplace nightly series on {len(valid_dates)} valid nights ({workers} workers)...")
    c_rows = []
    if workers <= 1:
        _init_worker(c_meta)
        for ch in chunks_c:
            c_rows.extend(_fit_cohort_chunk(ch))
    else:
        with ProcessPoolExecutor(max_workers=workers, initializer=_init_worker, initargs=(c_meta,)) as ex:
            for res in ex.map(_fit_cohort_chunk, chunks_c):
                c_rows.extend(res)
    df_c = pd.DataFrame(c_rows).merge(meta_cols, on="date", how="left").sort_values(["model", "term", "sample", "date"]).reset_index(drop=True)
    df_c.to_csv(OUT_DIR / "nightly_cohorts.csv", index=False)
    log(f"Saved nightly_cohorts.csv ({len(df_c):,} rows, {time.time() - t2:.1f}s)")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--repro-only", action="store_true")
    args = ap.parse_args()

    t0 = time.time()
    ensure_dirs()
    nodes = B0.load_nodes()
    run_reproduction_check(nodes)
    if not args.repro_only:
        run_nightly_series(nodes, workers=args.workers)
    print(f"Total elapsed: {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
