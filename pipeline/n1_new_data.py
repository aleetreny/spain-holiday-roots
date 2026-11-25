"""n1_new_data: Incorporate, validate, and test new open Spanish datasets for Question B.

Uses pipeline/config.py and pipeline/b00_common.py conventions (read-only on all other
pipeline/ and explore/ files). Run with:
    .venv/bin/python -I pipeline/n1_new_data.py [--stages=nodes,network,validate,roots_mun,roots_gen2,mediation]

Outputs (all in results/n1_*):
  n1_mun_nodes.parquet              Mainland municipality table (N=7,974) with pop, loss, area, Censo 2021 housing
  n1_mitma_zone_housing.csv         Censo 2021 housing aggregated to the 2,435 MITMA mainland zones
  n1_ine_od_inventory.csv           Summary of all 87 months (2019-07..2026-09) of INE municipal OD
  n1_ine_network_props.csv          Basic network properties of the true municipal INE OD network (density < 5%)
  n1_mun_network_node_tests.csv     Node-level emptied-village tests on true municipalities (INE OD)
  n1_mitma_validation.csv           Validation of INE OD vs MITMA across all 6 overlapping months of 2024
  n1_city_village_roots.csv         City x Village roots regressions (1st-gen & 2nd-gen strict) on INE municipal OD
  n1_city_village_qap.csv           Within-destination-province permutation null for City x Village roots
  n1_gen2_city_prov.csv             2nd-gen strict vs 1st-gen roots in R2's neighbouring-city design (MITMA & INE)
  n1_second_home_mediation.csv      Second-home stock mediation of emptied-village refill (municipal & MITMA zones)
  n1_ecepov_casa_del_pueblo.csv     Household-level 'casa del pueblo' match rates & permutation null (ECEPOV 2021)
"""
from __future__ import annotations

import argparse
import gzip
import random
import re
import sys
import time
import zlib
from collections import Counter
from pathlib import Path

import duckdb
import igraph as ig
import numpy as np
import pandas as pd
import pyfixest as pf
import scipy.sparse as sps
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C  # noqa: E402
import b00_common as B0  # noqa: E402

T0 = time.time()
EXT = C.ROOT / "_local" / "external"
OUT = C.RESULTS
NET = C.RESULTS / "b01_networks"

INE_OD_DIR = EXT / "ine_turismo_moviles_od"
INE_TM_RAW = EXT / "ine_turismo_moviles" / "raw"
ECEPOV_DIR = EXT / "ine_ecepov_2021"
VIV_RAW = EXT / "ine_censo_viviendas_2021" / "raw"
LAU_CSV = EXT / "eurostat_lau_2021" / "es_lau_2021.csv"
BIRTH_MUN_CSV = EXT / "ine_birth_residence" / "mun_residence_x_prov_birth_censo_2021_2025.csv"
ADRH_CSV = EXT / "ine_adrh" / "adrh_municipios.csv"

EXCL = set(C.EXCLUDED_PROVINCES)  # 07, 35, 38, 51, 52
WEBB = np.array([-np.sqrt(1.5), -1.0, -np.sqrt(0.5), np.sqrt(0.5), 1.0, np.sqrt(1.5)])

PROV2CCAA_NUM = {
    **{p: "01" for p in ["04", "11", "14", "18", "21", "23", "29", "41"]},
    **{p: "02" for p in ["22", "44", "50"]},
    "33": "03", "07": "04", "35": "05", "38": "05", "39": "06",
    **{p: "07" for p in ["05", "09", "24", "34", "37", "40", "42", "47", "49"]},
    **{p: "08" for p in ["02", "13", "16", "19", "45"]},
    **{p: "09" for p in ["08", "17", "25", "43"]},
    **{p: "10" for p in ["03", "12", "46"]},
    "06": "11", "10": "11",
    **{p: "12" for p in ["15", "27", "32", "36"]},
    "28": "13", "30": "14", "31": "15",
    **{p: "16" for p in ["01", "20", "48"]},
    "26": "17", "51": "18", "52": "19",
}
INLAND_CCAA_NUM = {"02", "07", "08", "11", "13", "15", "17"}

PROV_NAME_TO_CODE = {
    "Araba/Álava": "01", "Albacete": "02", "Alicante/Alacant": "03", "Almería": "04",
    "Ávila": "05", "Badajoz": "06", "Balears, Illes": "07", "Barcelona": "08",
    "Burgos": "09", "Cáceres": "10", "Cádiz": "11", "Castellón/Castelló": "12",
    "Ciudad Real": "13", "Córdoba": "14", "Coruña, A": "15", "Cuenca": "16",
    "Girona": "17", "Granada": "18", "Guadalajara": "19", "Gipuzkoa": "20",
    "Huelva": "21", "Huesca": "22", "Jaén": "23", "León": "24", "Lleida": "25",
    "Rioja, La": "26", "Lugo": "27", "Madrid": "28", "Málaga": "29", "Murcia": "30",
    "Navarra": "31", "Ourense": "32", "Asturias": "33", "Palencia": "34",
    "Palmas, Las": "35", "Pontevedra": "36", "Salamanca": "37",
    "Santa Cruz de Tenerife": "38", "Cantabria": "39", "Segovia": "40",
    "Sevilla": "41", "Soria": "42", "Tarragona": "43", "Teruel": "44",
    "Toledo": "45", "Valencia/València": "46", "Valladolid": "47", "Bizkaia": "48",
    "Zamora": "49", "Zaragoza": "50", "Ceuta": "51", "Melilla": "52",
}


def log(*a):
    print(f"[{time.time() - T0:7.1f}s]", *a, flush=True)


def seed_for(*parts) -> int:
    ints = [C.SEED] + [
        zlib.crc32(str(p).encode()) if not isinstance(p, (int, np.integer)) else int(p)
        for p in parts
    ]
    return int(np.random.SeedSequence(ints).generate_state(1)[0])


def parse_es_num(s: str) -> float | None:
    s = str(s).strip().strip('"')
    if not s or s in (".", "..", "-", "nan", "None"):
        return None
    try:
        return float(s.replace(".", "").replace(",", "."))
    except ValueError:
        return None


# =============================================================================
# Shared econometric helpers (FWL FE-OLS + Webb WCB; HDFE PPML + WCR score)
# =============================================================================
def fe_ols_wcb(
    y: np.ndarray,
    x_focus: np.ndarray,
    X_ctrl: np.ndarray,
    fe: np.ndarray,
    cl: np.ndarray,
    n_boot: int = 9999,
    seed: int = 0,
) -> dict:
    """One-way FE OLS with CRV1 (t with G-1 df) and Webb wild-cluster bootstrap."""
    y = np.asarray(y, float)
    xf = np.asarray(x_focus, float)
    Xc = np.asarray(X_ctrl, float) if X_ctrl.size else np.empty((len(y), 0))
    n = len(y)
    _, fe_i = np.unique(fe, return_inverse=True)
    _, cl_i = np.unique(cl, return_inverse=True)
    n_fe = int(fe_i.max()) + 1
    G = int(cl_i.max()) + 1

    def dm(A):
        A = np.asarray(A, float)
        one = A.ndim == 1
        M = A[:, None] if one else A.copy()
        cnt = np.bincount(fe_i, minlength=n_fe).astype(float)[:, None]
        s = np.zeros((n_fe, M.shape[1]))
        np.add.at(s, fe_i, M)
        out = M - (s / cnt)[fe_i]
        return out[:, 0] if one else out

    yd = dm(y)
    X = np.column_stack([xf, Xc])
    Xd = dm(X)
    k = Xd.shape[1]
    beta, *_ = np.linalg.lstsq(Xd, yd, rcond=None)
    resid = yd - Xd @ beta
    b = float(beta[0])

    # Partial out controls from x_focus (FWL) for fast scalar WCB
    if Xc.shape[1] > 0:
        Xc_d = Xd[:, 1:]
        g_x, *_ = np.linalg.lstsq(Xc_d, Xd[:, 0], rcond=None)
        xt = Xd[:, 0] - Xc_d @ g_x
        g_y, *_ = np.linalg.lstsq(Xc_d, yd, rcond=None)
        yt = yd - Xc_d @ g_y
    else:
        xt = Xd[:, 0]
        yt = yd
    xtx = float(xt @ xt)
    df_c = (G / (G - 1.0)) * ((n - 1.0) / (n - n_fe - k))

    def se_scalar(scores_mat):
        # scores_mat: (n, B) -> (B,)
        S = np.zeros((G, scores_mat.shape[1]))
        np.add.at(S, cl_i, scores_mat)
        return np.sqrt(df_c * np.sum(S ** 2, axis=0)) / xtx

    se = float(se_scalar((xt * resid)[:, None])[0])
    t_stat = b / se
    p_crv = float(2.0 * stats.t.sf(abs(t_stat), df=G - 1))
    tcrit = float(stats.t.ppf(0.975, df=G - 1))
    out = {
        "n": n, "G": G, "beta": b, "se_crv1": se, "t_crv1": t_stat, "p_crv1": p_crv,
        "ci_crv_lo": b - tcrit * se, "ci_crv_hi": b + tcrit * se,
    }
    if n_boot > 0:
        rng = np.random.default_rng(seed)
        # Restricted residual under H0: beta_0 = 0 -> e_r = yt
        er = yt
        # Unrestricted residual -> resid
        tr_all, tu_all = [], []
        bs = 1000
        for start in range(0, n_boot, bs):
            m = min(bs, n_boot - start)
            W = WEBB[rng.integers(0, len(WEBB), size=(G, m))][cl_i]  # (n, m)
            # Restricted bootstrap (WCR) for p-value
            yr_star = er[:, None] * W
            # Re-project out FE + controls? Since er is orthogonal to FE & controls, W*er needs FE/ctrl projection
            # For exactness when cluster == FE (province), W is constant within FE!
            # Check if cluster == FE (within-cluster constant weights):
            if np.array_equal(fe_i, cl_i):
                # W is constant within FE => W*er is already demeaned within FE!
                # Only project out Xc_d if needed:
                if Xc.shape[1] > 0:
                    br = np.linalg.lstsq(Xd, yr_star, rcond=None)[0]
                    ur = yr_star - Xd @ br
                    db_r = br[0]
                else:
                    db_r = (xt[:, None] * yr_star).sum(axis=0) / xtx
                    ur = yr_star - xt[:, None] * db_r[None, :]
                tr_all.append(db_r / se_scalar(xt[:, None] * ur))

                eu_star = resid[:, None] * W
                if Xc.shape[1] > 0:
                    bu = np.linalg.lstsq(Xd, eu_star, rcond=None)[0]
                    uu = eu_star - Xd @ bu
                    db_u = bu[0]
                else:
                    db_u = (xt[:, None] * eu_star).sum(axis=0) / xtx
                    uu = eu_star - xt[:, None] * db_u[None, :]
                tu_all.append(db_u / se_scalar(xt[:, None] * uu))
            else:
                db_r = (xt[:, None] * yr_star).sum(axis=0) / xtx
                ur = yr_star - xt[:, None] * db_r[None, :]
                tr_all.append(db_r / se_scalar(xt[:, None] * ur))
                eu_star = resid[:, None] * W
                db_u = (xt[:, None] * eu_star).sum(axis=0) / xtx
                uu = eu_star - xt[:, None] * db_u[None, :]
                tu_all.append(db_u / se_scalar(xt[:, None] * uu))
        tr = np.concatenate(tr_all)
        tu = np.concatenate(tu_all)
        out["p_wcr"] = float((np.sum(np.abs(tr) >= abs(t_stat)) + 1) / (n_boot + 1))
        out["ci_wcu_lo"] = float(b - np.quantile(tu, 0.975) * se)
        out["ci_wcu_hi"] = float(b - np.quantile(tu, 0.025) * se)
        out["n_boot"] = n_boot
    return out


class HDFE:
    """Weighted alternating projections over fixed-effect factors (from r2_roots.py)."""

    def __init__(self, codes_list):
        self.D = []
        for c in codes_list:
            c = np.asarray(c, dtype=np.int32)
            n, L = len(c), int(c.max()) + 1
            self.D.append(sps.csr_matrix((np.ones(n), (np.arange(n), c)), shape=(n, L)))

    def demean(self, X, w, tol=1e-9, maxit=2000):
        X = np.array(X, dtype=float, copy=True)
        one = X.ndim == 1
        if one:
            X = X[:, None]
        if X.shape[1] == 0:
            return X
        scale = max(1.0, float(np.abs(X).max()))
        ws = [D.T @ w for D in self.D]
        for _ in range(maxit):
            delta = 0.0
            for D, s in zip(self.D, ws):
                m = (D.T @ (w[:, None] * X)) / s[:, None]
                X -= D @ m
                delta = max(delta, float(np.abs(m).max()))
            if delta < tol * scale:
                return X[:, 0] if one else X
        raise RuntimeError("HDFE demean did not converge")


def poisson_dev(y, mu):
    t = np.where(y > 0, y * np.log(np.where(y > 0, y, 1.0) / mu), 0.0)
    return float(2.0 * np.sum(t - (y - mu)))


def drop_zero_groups(y, codes_list):
    keep = np.ones(len(y), bool)
    while True:
        bad = np.zeros(len(y), bool)
        for c in codes_list:
            s = np.bincount(c[keep], y[keep], minlength=int(c.max()) + 1)
            bad |= keep & (s[c] <= 0)
        if not bad.any():
            return keep
        keep &= ~bad


def ppml_hdfe(y, X, codes_list, offset=None, mu0=None, hd_pre=None, maxit=150):
    if hd_pre is None:
        keep = drop_zero_groups(y, codes_list)
        yk = y[keep]
        Xk = X[keep]
        offk = np.zeros(len(yk)) if offset is None else offset[keep]
        re_codes = [np.unique(c[keep], return_inverse=True)[1] for c in codes_list]
        hd = HDFE(re_codes)
    else:
        keep = np.ones(len(y), bool)
        yk, Xk = y, X
        offk = np.zeros(len(yk)) if offset is None else offset
        hd = hd_pre
    n, k = len(yk), Xk.shape[1]
    mu = (yk + yk.mean()) / 2.0 if mu0 is None else np.clip(mu0.copy(), 1e-12, 1e12)
    eta = np.log(mu)
    dev_old, beta = np.inf, np.zeros(k)
    for _ in range(maxit):
        z = eta - offk + (yk - mu) / mu
        Dm = hd.demean(np.column_stack([z, Xk]), mu)
        zt, Xt = Dm[:, 0], Dm[:, 1:]
        if k:
            XtW = Xt.T * mu
            b_new = np.linalg.solve(XtW @ Xt, XtW @ zt)
            e = zt - Xt @ b_new
        else:
            b_new, e = beta, zt
        eta = z - e + offk
        mu = np.exp(np.clip(eta, -30, 30))
        dev = poisson_dev(yk, mu)
        db = float(np.max(np.abs(b_new - beta))) if k else 0.0
        beta = b_new
        if abs(dev - dev_old) <= 1e-9 * (abs(dev) + 0.1) and db < 1e-7:
            return beta, mu, dev, keep, hd, Xt
        dev_old = dev
    raise RuntimeError("PPML did not converge")


def ppml_crv1_and_wcr(y, X, codes_list, cl, focus_idx=0, n_boot=9999, seed=0, do_ci=True):
    """PPML-HDFE with CRV1 (G/(G-1), t(G-1)) and Kline-Santos WCR score bootstrap for focus_idx."""
    beta, mu, dev, keep, hd, Xt = ppml_hdfe(y, X, codes_list)
    yk = y[keep]
    Xk = X[keep]
    clk = cl[keep]
    _, ci = np.unique(clk, return_inverse=True)
    G = int(ci.max()) + 1
    resid = yk - mu
    H = (Xt.T * mu) @ Xt
    Hinv = np.linalg.inv(H)
    scores = Xt * resid[:, None]
    S = np.zeros((G, Xk.shape[1]))
    np.add.at(S, ci, scores)
    meat = (G / (G - 1.0)) * (S.T @ S)
    V = Hinv @ meat @ Hinv
    se_all = np.sqrt(np.diag(V))
    b = float(beta[focus_idx])
    se = float(se_all[focus_idx])
    t_stat = b / se
    p_crv = float(2.0 * stats.t.sf(abs(t_stat), df=G - 1))
    tcrit = float(stats.t.ppf(0.975, df=G - 1))

    res = {
        "n": int(len(yk)), "G": G, "beta": b, "se_crv1": se, "t_crv1": t_stat, "p_crv1": p_crv,
        "ci_crv_lo": b - tcrit * se, "ci_crv_hi": b + tcrit * se,
        "beta_all": beta, "se_all": se_all,
    }
    if n_boot > 0:
        rng = np.random.default_rng(seed)
        W = WEBB[rng.integers(0, len(WEBB), size=(G, n_boot))]
        other_cols = [j for j in range(Xk.shape[1]) if j != focus_idx]
        X_o = Xk[:, other_cols]
        xf_raw = Xk[:, focus_idx]

        def wcr_pval(b0: float) -> float:
            off = b0 * xf_raw
            _, mu0, _, _, _, Xt_o = ppml_hdfe(yk, X_o, None, offset=off, mu0=mu, hd_pre=hd)
            xf_d = hd.demean(xf_raw, mu0)
            if Xt_o.shape[1] > 0:
                Hoo = (Xt_o.T * mu0) @ Xt_o
                g = np.linalg.solve(Hoo, (Xt_o.T * mu0) @ xf_d)
                x_eff = xf_d - Xt_o @ g
            else:
                x_eff = xf_d
            sc_i = x_eff * (yk - mu0)
            Sc = np.bincount(ci, weights=sc_i, minlength=G)
            obs = float(Sc.sum())
            se0 = float(np.sqrt((G / (G - 1.0)) * np.sum(Sc ** 2)))
            t0 = abs(obs) / se0 if se0 > 0 else 0.0
            ScW = Sc[:, None] * W
            num = np.abs(ScW.sum(axis=0))
            den = np.sqrt((G / (G - 1.0)) * np.sum((ScW - ScW.mean(axis=0, keepdims=True)) ** 2, axis=0))
            tb = np.where(den > 0, num / den, 0.0)
            return float((np.sum(tb >= t0) + 1) / (n_boot + 1))

        res["p_wcr"] = wcr_pval(0.0)
        if do_ci:
            def bisect_bound(direction: float) -> float:
                step = max(se, 0.02)
                lo_b = b
                hi_b = b + direction * 2.5 * step
                for _ in range(6):
                    if wcr_pval(hi_b) < 0.05:
                        break
                    hi_b += direction * 2.5 * step
                else:
                    return hi_b
                a, c = (lo_b, hi_b) if direction > 0 else (hi_b, lo_b)
                for _ in range(9):
                    mid = 0.5 * (a + c)
                    pm = wcr_pval(mid)
                    if direction > 0:
                        if pm >= 0.05:
                            a = mid
                        else:
                            c = mid
                    else:
                        if pm >= 0.05:
                            c = mid
                        else:
                            a = mid
                return 0.5 * (a + c)

            res["ci_wcr_lo"] = float(bisect_bound(-1.0))
            res["ci_wcr_hi"] = float(bisect_bound(+1.0))
        else:
            res["ci_wcr_lo"] = np.nan
            res["ci_wcr_hi"] = np.nan
        res["n_boot"] = n_boot
    return res


# =============================================================================
# STAGE 1: Municipality Node Table + Censo 2021 Housing + Eurostat LAU Area
# =============================================================================
def stage_nodes() -> tuple[pd.DataFrame, pd.DataFrame]:
    log("STAGE 1: Building true mainland municipality node table + Censo 2021 housing...")
    mitma_nodes = B0.load_nodes()
    xy = mitma_nodes[["x", "y"]].values / 1000.0
    Dm = np.hypot(xy[:, None, 0] - xy[None, :, 0], xy[:, None, 1] - xy[None, :, 1])
    big = mitma_nodes["population"].values >= 100_000
    coast = mitma_nodes["coastal_true"].values.astype(bool)
    mitma_nodes["dist_city_km"] = Dm[:, big].min(axis=1)
    mitma_nodes["dist_coast_km"] = Dm[:, coast].min(axis=1)
    zid = mitma_nodes.set_index("id")

    # Crosswalk INE municipality -> primary MITMA zone (by census section count, as R1)
    rel_full = pd.read_csv(B0.INE_MITMA_REL, sep="|", dtype=str)
    zone_of = (
        rel_full.groupby(["municipio_ine", "municipio_mitma"])
        .size()
        .reset_index(name="n")
        .sort_values(["municipio_ine", "n", "municipio_mitma"], ascending=[True, False, True])
        .drop_duplicates("municipio_ine")
        .set_index("municipio_ine")["municipio_mitma"]
    )

    # Historical population (1900-2025)
    H = pd.read_csv(B0.HIST_POP_CSV, dtype={"ine_code": str, "province_code": str})
    H = H[~H["province_code"].isin(EXCL)].copy()
    H["zone_id"] = H["ine_code"].map(zone_of)
    H = H[H["zone_id"].isin(zid.index)].copy()
    peak = H[B0._PEAK_COLS].max(axis=1)
    H["loss_peak_2021"] = 1.0 - H["pop_2021"] / peak
    H["m2_incomplete"] = H[B0._PEAK_COLS + ["pop_2021"]].isna().any(axis=1)

    # Eurostat GISCO LAU 2021 exact municipal area_km2
    lau = pd.read_csv(LAU_CSV, dtype=str)
    lau["ine_code"] = lau["LAU_ID"].str.zfill(5)
    lau["area_km2"] = pd.to_numeric(lau["AREA_KM2"], errors="coerce")
    H = H.merge(lau[["ine_code", "area_km2"]], on="ine_code", how="left")

    # Censo 2021 table 59525: Total, principal, no principal (all 8,131 municipalities)
    d5 = pd.read_csv(VIV_RAW / "censo2021_59525.csv", sep=";", dtype=str)
    d5 = d5[d5["Municipios"].notna() & d5["Municipios"].str.match(r"^\d{5}\s")].copy()
    d5["ine_code"] = d5["Municipios"].str[:5]
    d5["val"] = d5["Total"].map(parse_es_num)
    w5 = d5.pivot(index="ine_code", columns="Tipo de vivienda (principal o no)", values="val").reset_index()
    w5 = w5.rename(columns={
        "Total": "viv_total",
        "Vivienda principal": "viv_principal",
        "Vivienda no principal": "viv_no_principal",
    })
    assert np.allclose(w5["viv_total"], w5["viv_principal"] + w5["viv_no_principal"], equal_nan=True)
    w5["share_no_principal"] = w5["viv_no_principal"] / w5["viv_total"]

    # Censo 2021 table 59531: Electricity-based empty & sporadic use (municipalities >= ~500)
    d31 = pd.read_csv(VIV_RAW / "censo2021_59531.csv", sep=";", dtype=str)
    d31 = d31[d31["Municipios"].notna() & d31["Municipios"].str.match(r"^\d{5}\s")].copy()
    d31["ine_code"] = d31["Municipios"].str[:5]
    d31["val"] = d31["Total"].map(parse_es_num)
    keep31 = ["Viviendas totales", "Viviendas vacías", "Viviendas con bajo consumo", "Viviendas de uso esporádico"]
    w31 = d31[d31["Consumo eléctrico"].isin(keep31)].pivot(
        index="ine_code", columns="Consumo eléctrico", values="val"
    ).reset_index().rename(columns={
        "Viviendas totales": "viv_tot_31",
        "Viviendas vacías": "viv_vacias",
        "Viviendas con bajo consumo": "viv_bajo_consumo",
        "Viviendas de uso esporádico": "viv_esporadico",
    })

    H = H.merge(w5, on="ine_code", how="left").merge(w31, on="ine_code", how="left")
    H["viv_secundaria_est"] = (H["viv_no_principal"] - H["viv_vacias"]).clip(lower=0)
    H["share_secundaria_est"] = H["viv_secundaria_est"] / H["viv_total"]
    H["share_esporadico"] = H["viv_esporadico"] / H["viv_total"]
    H["share_vacias"] = H["viv_vacias"] / H["viv_total"]
    H["no_principal_per_res"] = H["viv_no_principal"] / H["padron_2024"].replace(0, np.nan)

    # Parent MITMA zone geography
    for col in ["coastal_true", "dist_city_km", "dist_coast_km", "x", "y", "is_am", "n_municipalities"]:
        H[f"zone_{col}"] = H["zone_id"].map(zid[col])
    H["coastal_true"] = H["zone_coastal_true"].astype(bool)
    H["zone_single"] = H["zone_n_municipalities"] == 1
    H["prov"] = H["ine_code"].str[:2]
    H["ccaa"] = H["prov"].map(B0.PROV2CCAA)

    # Census 137 mainland cities flag
    Bm = pd.read_csv(BIRTH_MUN_CSV, dtype={"ine_code": str})
    cities137 = {c for c in Bm["ine_code"].unique() if c[:2] not in EXCL}
    assert len(cities137) == 137
    H["is_census_city"] = H["ine_code"].isin(cities137)

    # Percentile rank of loss over valid municipalities
    ok = H["loss_peak_2021"].notna() & H["padron_2024"].notna() & (H["padron_2024"] > 0)
    H["loss_pct"] = np.nan
    H.loc[ok, "loss_pct"] = H.loc[ok, "loss_peak_2021"].rank(pct=True, method="average")

    cols_keep = [
        "ine_code", "prov", "ccaa", "name", "zone_id", "zone_single", "zone_is_am", "is_census_city",
        "pop_2021", "padron_2024", "loss_peak_2021", "loss_pct", "m2_incomplete", "area_km2",
        "coastal_true", "zone_dist_city_km", "zone_dist_coast_km", "zone_x", "zone_y",
        "viv_total", "viv_principal", "viv_no_principal", "share_no_principal", "no_principal_per_res",
        "viv_vacias", "viv_esporadico", "viv_bajo_consumo", "viv_secundaria_est",
        "share_vacias", "share_esporadico", "share_secundaria_est",
    ]
    mun = H[cols_keep].sort_values("ine_code").reset_index(drop=True)
    mun.to_parquet(OUT / "n1_mun_nodes.parquet", index=False)

    # Aggregate 59525 housing to the 2,435 MITMA zones
    zh = (
        mun.groupby("zone_id", as_index=False)[["viv_total", "viv_principal", "viv_no_principal"]]
        .sum(min_count=1)
    )
    zh["share_no_principal"] = zh["viv_no_principal"] / zh["viv_total"]
    zh = mitma_nodes[["id", "prov", "population", "coastal_true", "loss_peak_2021", "n_municipalities", "area_km2", "is_am"]].merge(
        zh, left_on="id", right_on="zone_id", how="left"
    )
    zh["no_principal_per_res"] = zh["viv_no_principal"] / zh["population"]
    zh.to_csv(OUT / "n1_mitma_zone_housing.csv", index=False)

    log(
        f"  mun nodes N={len(mun):,}, matched area={mun['area_km2'].notna().sum():,}, "
        f"matched 59525={mun['viv_total'].notna().sum():,}, matched 59531={mun['viv_vacias'].notna().sum():,}; "
        f"MITMA zones matched 59525={zh['viv_total'].notna().sum()}/{len(zh)}"
    )
    return mun, mitma_nodes


# =============================================================================
# STAGE 2: True Municipal INE OD Network Properties & Node-Level Emptied Tests
# =============================================================================
def rich_club_by_rank(g: ig.Graph, top_ns=(50, 100, 200), n_rewire=20, seed=0) -> dict:
    """Rank-based rich club (audit M1) on undirected graph g with rewired ceiling."""
    deg = np.array(g.degree(), dtype=int)
    order = np.lexsort((np.arange(g.vcount()), -deg))  # ties broken by node index ascending
    adj_sets = [set(g.neighbors(u)) for u in range(g.vcount())]

    def phi_of(graph_adj):
        out = {}
        for n_top in top_ns:
            club = set(order[:n_top].tolist())
            e_in = sum(len(graph_adj[u] & club) for u in club) // 2
            out[n_top] = (2.0 * e_in) / (n_top * (n_top - 1))
        return out

    phi_obs = phi_of(adj_sets)
    rew_phis = {n_top: [] for n_top in top_ns}
    rng = random.Random(seed)
    E = g.ecount()
    for _ in range(n_rewire):
        gr = g.copy()
        random.seed(rng.randrange(1, 2**31 - 1))
        gr.rewire(n=10 * E, mode="simple")
        r_adj = [set(gr.neighbors(u)) for u in range(gr.vcount())]
        pr = phi_of(r_adj)
        for n_top in top_ns:
            rew_phis[n_top].append(pr[n_top])
    res = {}
    for n_top in top_ns:
        m_rew = float(np.mean(rew_phis[n_top]))
        res[f"rc_phi_{n_top}"] = phi_obs[n_top]
        res[f"rc_phi_rew_{n_top}"] = m_rew
        res[f"rc_rho_{n_top}"] = phi_obs[n_top] / m_rew if m_rew > 0 else np.nan
        res[f"rc_ceil_{n_top}"] = 1.0 / m_rew if m_rew > 0 else np.nan
    return res


def stage_network(mun: pd.DataFrame):
    log("STAGE 2: Building true municipal INE OD network & computing topology + node tests...")
    # 1. Inventory across all 8 years (2019-2026)
    inv_rows = []
    for yr in range(2019, 2027):
        p = INE_OD_DIR / f"ine_od_mun_{yr}.parquet"
        if not p.exists():
            continue
        df = pd.read_parquet(p)
        for mes, g in df.groupby("mes"):
            ml = g[~g["orig"].str[:2].isin(EXCL) & ~g["dest"].str[:2].isin(EXCL)]
            inv_rows.append({
                "year": yr, "mes": mes,
                "edges_all": len(g), "turistas_all": int(g["turistas"].sum()),
                "edges_mainland": len(ml), "turistas_mainland": int(ml["turistas"].sum()),
                "n_orig_ml": ml["orig"].nunique(), "n_dest_ml": ml["dest"].nunique(),
                "n_union_ml": len(set(ml["orig"]) | set(ml["dest"])),
                "min_turistas": int(g["turistas"].min()),
            })
    inv = pd.DataFrame(inv_rows).sort_values("mes")
    inv.to_csv(OUT / "n1_ine_od_inventory.csv", index=False)
    log(f"  INE OD inventory: {len(inv)} months ({inv['mes'].min()}..{inv['mes'].max()}), total rows={inv['edges_all'].sum():,}")

    # 2. Fixed mainland municipality node set
    ml_ids = mun["ine_code"].values
    N_FIX = len(ml_ids)
    pos = {c: k for k, c in enumerate(ml_ids)}

    # Load 2022, 2023, 2024, 2025 target months
    target_months = [
        ("2024", "2024-02", "ordinary"),
        ("2024", "2024-03", "easter"),
        ("2024", "2024-08", "summer"),
        ("2024", "2024-10", "october"),
        ("2024", "2024-11", "allsaints"),
        ("2024", "2024-12", "christmas"),
        ("2022", "2022-02", "ordinary"),
        ("2022", "2022-08", "summer"),
        ("2023", "2023-02", "ordinary"),
        ("2023", "2023-08", "summer"),
        ("2025", "2025-02", "ordinary"),
        ("2025", "2025-08", "summer"),
    ]
    by_yr = {}
    for yr in ["2022", "2023", "2024", "2025"]:
        df = pd.read_parquet(INE_OD_DIR / f"ine_od_mun_{yr}.parquet")
        df = df[df["orig"].isin(pos) & df["dest"].isin(pos) & (df["orig"] != df["dest"])].copy()
        df["i"] = df["orig"].map(pos).astype(np.int32)
        df["j"] = df["dest"].map(pos).astype(np.int32)
        by_yr[yr] = df

    # Reference E for equal_edges = E_dir of 2024-02 published (w >= 30)
    feb24 = by_yr["2024"][by_yr["2024"]["mes"] == "2024-02"]
    E_EQ = len(feb24)
    log(f"  Fixed mainland N={N_FIX:,}; equal_edges E (2024-02 pub >=30) = {E_EQ:,}")

    prop_rows = []
    node_metrics_2024 = {}

    for yr, mes, grp in target_months:
        sub = by_yr[yr][by_yr[yr]["mes"] == mes][["i", "j", "turistas"]].rename(columns={"turistas": "w"})
        backbones = [
            ("pub_ge30", sub),
            ("w100", sub[sub["w"] >= 100]),
            ("equal_edges", sub.sort_values(["w", "i", "j"], ascending=[False, True, True], kind="mergesort").head(E_EQ)),
        ]
        for bname, ed in backbones:
            # Skip equal_edges/w100 rewiring on non-2024 replication years to save time (still compute exact graph metrics)
            do_nulls = (yr == "2024" and bname in ("pub_ge30", "equal_edges"))
            gd = B0.to_igraph(ed, N_FIX, directed=True)
            gu = B0.to_igraph(ed, N_FIX, directed=False)
            E_dir = gd.ecount()
            E_und = gu.ecount()
            dens_und = 2.0 * E_und / (N_FIX * (N_FIX - 1))
            assert dens_und < C.MAX_DENSITY, f"density {dens_und:.4f} >= 0.05"

            k_in = np.array(gd.indegree(), dtype=float)
            k_out = np.array(gd.outdegree(), dtype=float)
            k_und = np.array(gu.degree(), dtype=float)
            active = (k_und > 0)
            n_active = int(active.sum())
            dens_active = 2.0 * E_und / (n_active * (n_active - 1)) if n_active > 1 else 0.0

            # Spearman(k_in, k_out) over active nodes (and over all N_FIX)
            rho_io_act = float(stats.spearmanr(k_in[active], k_out[active]).statistic)
            rho_io_all = float(stats.spearmanr(k_in, k_out).statistic)
            r_und = float(gu.assortativity_degree(directed=False))
            c_avg = float(gu.transitivity_avglocal_undirected(mode="zero"))
            c_glob = float(gu.transitivity_undirected())
            recip = float(gd.reciprocity())
            core_u = np.array(gu.coreness(), dtype=int)
            core_in = np.array(gd.coreness(mode="in"), dtype=int)
            comps = gu.connected_components()
            gcc = comps.giant()
            gcc_n = gcc.vcount()
            # Approximate or exact L on gcc (sample 400 sources for speed on 6k nodes -> SE < 0.003)
            rng_l = np.random.default_rng(seed_for("ine_L", mes, bname))
            srcs = rng_l.choice(gcc_n, size=min(400, gcc_n), replace=False).tolist()
            dmat = np.array(gcc.distances(source=srcs), dtype=float)
            L_gcc = float(dmat[dmat > 0].mean())

            row = {
                "mes": mes, "group": grp, "backbone": bname,
                "N_fixed": N_FIX, "N_active": n_active,
                "E_dir": E_dir, "E_und": E_und,
                "density_und_fixed": dens_und, "density_und_active": dens_active,
                "total_tourists": float(ed["w"].sum()),
                "mean_k_und_active": 2.0 * E_und / n_active,
                "max_k_out": int(k_out.max()), "max_k_in": int(k_in.max()),
                "spearman_kin_kout_active": rho_io_act, "spearman_kin_kout_all": rho_io_all,
                "reciprocity": recip, "r_und": r_und,
                "C_avg_all": c_avg, "C_avg_active": float(np.mean(np.array(gu.transitivity_local_undirected(mode="zero"))[active])),
                "C_global": c_glob, "L_gcc_est": L_gcc, "gcc_n": gcc_n,
                "max_core_und": int(core_u.max()), "max_core_in": int(core_in.max()),
            }

            if do_nulls:
                # Degree-preserving undirected rewiring (N_REWIRE = 20 on 100k edges; fast & exact)
                n_rew = 20
                r_rews, c_rews = [], []
                rng_r = random.Random(seed_for("ine_rew", mes, bname))
                for _ in range(n_rew):
                    gr = gu.copy()
                    random.seed(rng_r.randrange(1, 2**31 - 1))
                    gr.rewire(n=10 * E_und, mode="simple")
                    r_rews.append(float(gr.assortativity_degree(directed=False)))
                    c_rews.append(float(gr.transitivity_avglocal_undirected(mode="zero")))
                row["r_rew_mean"] = float(np.mean(r_rews))
                row["r_rew_sd"] = float(np.std(r_rews, ddof=1))
                row["C_rew_mean"] = float(np.mean(c_rews))
                row["C_er_expected"] = dens_active
                rc = rich_club_by_rank(gu, top_ns=(50, 100), n_rewire=n_rew, seed=seed_for("ine_rc", mes, bname))
                row.update(rc)

            prop_rows.append(row)

            # Save 2024 pub_ge30 node metrics for emptied-village test
            if yr == "2024" and bname == "pub_ge30":
                s_in = np.bincount(ed["j"].values, weights=ed["w"].values, minlength=N_FIX)
                s_out = np.bincount(ed["i"].values, weights=ed["w"].values, minlength=N_FIX)
                node_metrics_2024[mes] = pd.DataFrame({
                    "ine_code": ml_ids,
                    "s_in": s_in, "s_out": s_out,
                    "k_in": k_in, "k_out": k_out,
                    "core_und": core_u.astype(float), "core_in": core_in.astype(float),
                    "net_bal": np.where((s_in + s_out) > 0, (s_in - s_out) / (s_in + s_out), np.nan),
                    "net_importer": np.where((s_in + s_out) > 0, (s_in > s_out).astype(float), np.nan),
                })

    props_df = pd.DataFrame(prop_rows)
    props_df.to_csv(OUT / "n1_ine_network_props.csv", index=False)
    log(f"  Saved n1_ine_network_props.csv ({len(props_df)} rows); max density_und={props_df['density_und_fixed'].max():.4%}")

    # 3. True municipal node-level emptied-village tests on INE OD (zero '_AM' aggregation!)
    base_m = node_metrics_2024["2024-02"].set_index("ine_code")
    nd = mun.copy()
    nd["lpop"] = np.log(nd["padron_2024"].clip(lower=1))
    nd["larea"] = np.log(nd["area_km2"].clip(lower=0.01))
    nd["coastal"] = nd["coastal_true"].astype(float)
    nd["ldcity"] = np.log1p(nd["zone_dist_city_km"])
    nd["ldcoast"] = np.log1p(nd["zone_dist_coast_km"])
    ok_base = nd["loss_pct"].notna() & nd["area_km2"].notna() & (nd["padron_2024"] > 0)
    ctrl_cols = ["lpop", "larea", "coastal", "ldcity", "ldcoast"]

    nt_rows = []
    for hol_mes, hol_lbl in [("2024-08", "Aug_vs_Feb"), ("2024-03", "Mar_vs_Feb"), ("2024-11", "Nov_vs_Feb"), ("2024-12", "Dec_vs_Feb")]:
        hm = node_metrics_2024[hol_mes].set_index("ine_code")
        outcomes = {
            "log_s_in_ratio": np.log((hm.loc[nd["ine_code"], "s_in"].values + 1.0) / (base_m.loc[nd["ine_code"], "s_in"].values + 1.0)),
            "log_k_in_ratio": np.log((hm.loc[nd["ine_code"], "k_in"].values + 1.0) / (base_m.loc[nd["ine_code"], "k_in"].values + 1.0)),
            "core_und_diff": hm.loc[nd["ine_code"], "core_und"].values - base_m.loc[nd["ine_code"], "core_und"].values,
            "core_in_diff": hm.loc[nd["ine_code"], "core_in"].values - base_m.loc[nd["ine_code"], "core_in"].values,
            "net_bal_diff": hm.loc[nd["ine_code"], "net_bal"].values - base_m.loc[nd["ine_code"], "net_bal"].values,
            "net_importer_level": hm.loc[nd["ine_code"], "net_importer"].values,
        }
        for oname, yv in outcomes.items():
            for samp_name, smask in [
                ("all_mainland_mun", ok_base),
                ("inland_mun", ok_base & ~nd["coastal_true"]),
                ("inside_AM_mun", ok_base & nd["zone_is_am"]),
                ("single_zone_mun", ok_base & nd["zone_single"]),
            ]:
                if samp_name != "all_mainland_mun" and hol_mes != "2024-08":
                    continue
                m = smask.values & np.isfinite(yv)
                d = nd.loc[m]
                c_use = [c for c in ctrl_cols if not (samp_name == "inland_mun" and c == "coastal")]
                res = fe_ols_wcb(
                    yv[m], d["loss_pct"].values, d[c_use].values,
                    d["prov"].values, d["prov"].values,
                    n_boot=1999, seed=seed_for("nt", hol_lbl, oname, samp_name),
                )
                res.update(contrast=hol_lbl, outcome=oname, sample=samp_name)
                nt_rows.append(res)
    pd.DataFrame(nt_rows).to_csv(OUT / "n1_mun_network_node_tests.csv", index=False)
    log(f"  Saved n1_mun_network_node_tests.csv ({len(nt_rows)} specs)")


# =============================================================================
# STAGE 3: Validation of INE OD vs MITMA across all 6 overlapping months of 2024
# =============================================================================
def stage_validate(mun: pd.DataFrame, mitma_nodes: pd.DataFrame):
    log("STAGE 3: Validating INE OD vs MITMA across all 6 overlapping months of 2024...")
    mun2zone = mun.set_index("ine_code")["zone_id"].to_dict()
    single_zones = set(mitma_nodes.loc[mitma_nodes["n_municipalities"] == 1, "id"])
    valid_zones = set(mitma_nodes["id"])

    # Load INE 52997 (prov x prov: Turistas and Pernoctaciones) for 2024
    d97 = pd.read_csv(INE_TM_RAW / "turismo_52997_prov_x_prov.csv.gz", sep=";", dtype=str)
    d97 = d97[d97["Periodo"].str.startswith("2024M")].copy()
    # Map exact INE province names in column .2 -> 2-digit province code
    d97["po"] = d97["CCAA y provincia de origen.2"].str.strip().map(PROV_NAME_TO_CODE)
    d97["pd"] = d97["CCAA y provincia de destino.2"].str.strip().map(PROV_NAME_TO_CODE)
    d97["val"] = d97["Total"].map(parse_es_num)
    p97 = d97[d97["po"].notna() & d97["pd"].notna() & (~d97["po"].isin(EXCL)) & (~d97["pd"].isin(EXCL)) & (d97["po"] != d97["pd"])].copy()
    assert p97["po"].nunique() == 47 and p97["pd"].nunique() == 47, f"expected 47 provs, got {p97['po'].nunique()}"

    # Load INE municipal OD 2024
    od24 = pd.read_parquet(INE_OD_DIR / "ine_od_mun_2024.parquet")
    od24 = od24[~od24["orig"].str[:2].isin(EXCL) & ~od24["dest"].str[:2].isin(EXCL)].copy()
    od24["zo"] = od24["orig"].map(mun2zone)
    od24["zd"] = od24["dest"].map(mun2zone)

    # Query MITMA monthly interprovincial totals via DuckDB
    con = duckdb.connect()
    mitma_months = ["02", "03", "08", "10", "11", "12"]
    val_rows = []

    for mm in mitma_months:
        mes = f"2024-{mm}"
        per = f"2024M{mm}"
        # MITMA monthly zone x zone interprovincial sum (days present in that month)
        m_df = con.execute(
            """
            SELECT r, s, SUM(p) AS mitma_nights, COUNT(DISTINCT date) AS n_days
            FROM read_parquet(?)
            WHERE strftime(date, '%m') = ? AND substr(r, 1, 2) <> substr(s, 1, 2)
            GROUP BY 1, 2
            """,
            [str(C.PERNOCT_2024), mm],
        ).df()
        n_days_month = int(con.execute(
            "SELECT COUNT(DISTINCT date) FROM read_parquet(?) WHERE strftime(date, '%m') = ?",
            [str(C.PERNOCT_2024), mm],
        ).fetchone()[0])
        m_df = m_df[m_df["r"].isin(valid_zones) & m_df["s"].isin(valid_zones)].copy()
        m_df["po"] = m_df["r"].str[:2]
        m_df["pd"] = m_df["s"].str[:2]

        # 1. Province x Province (47 x 46 = 2,162 off-diagonal mainland pairs)
        m_pp = m_df.groupby(["po", "pd"], as_index=False)["mitma_nights"].sum()
        ine97_m = p97[p97["Periodo"] == per]
        i97_tur = ine97_m[ine97_m["Concepto turístico"] == "Turistas"][["po", "pd", "val"]].rename(columns={"val": "ine52997_tur"})
        i97_per = ine97_m[ine97_m["Concepto turístico"] == "Pernoctaciones"][["po", "pd", "val"]].rename(columns={"val": "ine52997_pernoct"})
        od_m = od24[od24["mes"] == mes]
        od_pp = od_m.groupby([od_m["orig"].str[:2].rename("po"), od_m["dest"].str[:2].rename("pd")], as_index=False)["turistas"].sum().rename(columns={"turistas": "ine_od_tur"})

        pp = m_pp.merge(i97_tur, on=["po", "pd"], how="outer").merge(i97_per, on=["po", "pd"], how="outer").merge(od_pp, on=["po", "pd"], how="outer").fillna(0)
        cov_pct = pp["ine_od_tur"].sum() / pp["ine52997_tur"].sum()

        # 2. MITMA zone x zone interprovincial
        od_zz = od_m.dropna(subset=["zo", "zd"]).groupby(["zo", "zd"], as_index=False)["turistas"].sum()
        zz = m_df[["r", "s", "mitma_nights"]].merge(od_zz, left_on=["r", "s"], right_on=["zo", "zd"], how="inner")
        # Single-municipality x single-municipality (TRUE 1-to-1 municipal comparison!)
        zz_1to1 = zz[zz["r"].isin(single_zones) & zz["s"].isin(single_zones)]

        # 3. Destination zone total inflows (interprovincial)
        m_dest = m_df.groupby("s")["mitma_nights"].sum()
        i_dest = od_m.dropna(subset=["zd"]).groupby("zd")["turistas"].sum()
        dest_df = pd.DataFrame({"mitma": m_dest, "ine_od": i_dest}).fillna(0)

        val_rows.append({
            "mes": mes, "mitma_days": n_days_month,
            "mitma_interprov_nights": float(pp["mitma_nights"].sum()),
            "ine52997_pernoct": float(pp["ine52997_pernoct"].sum()),
            "ine52997_tur": float(pp["ine52997_tur"].sum()),
            "ine_od_tur": float(pp["ine_od_tur"].sum()),
            "ine_od_coverage_vs_52997": float(cov_pct),
            "pp_r_mitma_vs_52997pernoct": float(np.corrcoef(pp["mitma_nights"], pp["ine52997_pernoct"])[0, 1]),
            "pp_rlog_mitma_vs_52997pernoct": float(np.corrcoef(np.log1p(pp["mitma_nights"]), np.log1p(pp["ine52997_pernoct"]))[0, 1]),
            "pp_rho_mitma_vs_52997pernoct": float(stats.spearmanr(pp["mitma_nights"], pp["ine52997_pernoct"]).statistic),
            "pp_rlog_mitma_vs_52997tur": float(np.corrcoef(np.log1p(pp["mitma_nights"]), np.log1p(pp["ine52997_tur"]))[0, 1]),
            "pp_rlog_mitma_vs_od": float(np.corrcoef(np.log1p(pp["mitma_nights"]), np.log1p(pp["ine_od_tur"]))[0, 1]),
            "pp_rho_mitma_vs_od": float(stats.spearmanr(pp["mitma_nights"], pp["ine_od_tur"]).statistic),
            "zz_overlap_pairs": len(zz),
            "zz_rlog_mitma_vs_od": float(np.corrcoef(np.log1p(zz["mitma_nights"]), np.log1p(zz["turistas"]))[0, 1]),
            "zz_rho_mitma_vs_od": float(stats.spearmanr(zz["mitma_nights"], zz["turistas"]).statistic),
            "single1to1_pairs": len(zz_1to1),
            "single1to1_rlog": float(np.corrcoef(np.log1p(zz_1to1["mitma_nights"]), np.log1p(zz_1to1["turistas"]))[0, 1]),
            "single1to1_rho": float(stats.spearmanr(zz_1to1["mitma_nights"], zz_1to1["turistas"]).statistic),
            "dest_zones_n": len(dest_df),
            "dest_rlog": float(np.corrcoef(np.log1p(dest_df["mitma"]), np.log1p(dest_df["ine_od"]))[0, 1]),
            "dest_rho": float(stats.spearmanr(dest_df["mitma"], dest_df["ine_od"]).statistic),
        })
    vdf = pd.DataFrame(val_rows)
    vdf.to_csv(OUT / "n1_mitma_validation.csv", index=False)
    log("  Saved n1_mitma_validation.csv:\n" + vdf[["mes", "ine_od_coverage_vs_52997", "pp_rlog_mitma_vs_52997pernoct", "zz_rlog_mitma_vs_od", "single1to1_rlog", "dest_rho"]].to_string(index=False))


# =============================================================================
# STAGE 4: True City x Village Roots (1st & 2nd Gen Strict) on INE Municipal OD
# =============================================================================
def stage_roots_mun(mun: pd.DataFrame):
    log("STAGE 4: Building City x Village roots (ECEPOV 2021) and testing on INE municipal OD...")
    dh = pd.read_parquet(ECEPOV_DIR / "ecepov_hogar_2021.parquet")
    for c in ["IDQ_PV", "IDQ_MUN", "PROVN", "MUNIN", "PROVNP1", "MUNINP1", "PROVNP2", "MUNINP2"]:
        dh[c] = dh[c].astype(str).str.strip()

    ml_set = set(mun["ine_code"])
    cities = sorted(mun.loc[mun["is_census_city"], "ine_code"].unique())
    assert len(cities) == 137
    city_set = set(cities)

    ch = dh[dh["IDQ_MUN"].isin(city_set)].copy()
    # 1st-gen: born in a mainland municipality of ANOTHER province
    g1 = ch[
        (ch["PROVN"] != ch["IDQ_PV"]) & (~ch["PROVN"].isin(EXCL)) & ch["MUNIN"].isin(ml_set)
    ].groupby(["IDQ_MUN", "MUNIN"])["FACTOR"].sum().rename("born_g1")

    # 2nd-gen STRICT: person born in OWN province (NOT migrant themselves), parent born in OTHER mainland prov village
    own_born = ch[ch["PROVN"] == ch["IDQ_PV"]]
    p1 = own_born[
        (own_born["PROVNP1"] != own_born["IDQ_PV"]) & (~own_born["PROVNP1"].isin(EXCL)) & own_born["MUNINP1"].isin(ml_set)
    ].groupby(["IDQ_MUN", "MUNINP1"])["FACTOR"].sum() * 0.5
    p2 = own_born[
        (own_born["PROVNP2"] != own_born["IDQ_PV"]) & (~own_born["PROVNP2"].isin(EXCL)) & own_born["MUNINP2"].isin(ml_set)
    ].groupby(["IDQ_MUN", "MUNINP2"])["FACTOR"].sum() * 0.5
    p1.index.names = ["IDQ_MUN", "MUN"]
    p2.index.names = ["IDQ_MUN", "MUN"]
    g2 = p1.add(p2, fill_value=0.0).rename("born_g2")
    g1.index.names = ["IDQ_MUN", "MUN"]

    roots_cv = pd.concat([g1, g2], axis=1).fillna(0.0).reset_index()
    roots_cv.columns = ["orig", "dest", "born_g1", "born_g2"]
    roots_cv["born_any"] = roots_cv["born_g1"] + roots_cv["born_g2"]
    roots_cv.to_parquet(OUT / "n1_roots_city_village_ecepov2021.parquet", index=False)
    log(
        f"  ECEPOV City x Village roots: {len(roots_cv):,} non-zero (city, village) pairs; "
        f"g1>0={(roots_cv['born_g1']>0).sum():,}, g2>0={(roots_cv['born_g2']>0).sum():,}, "
        f"unique dest villages={roots_cv['dest'].nunique():,}"
    )

    # Build city x village panel for INE municipal OD across 2024 months (+ 2022, 2023, 2025)
    # Universe for FE estimation: all active destination villages that receive >= 1 published link from the 137 cities
    # in the analysed months, crossed with the 137 cities (off-own-province).
    mun_info = mun.set_index("ine_code")
    od24 = pd.read_parquet(INE_OD_DIR / "ine_od_mun_2024.parquet")
    od24 = od24[od24["orig"].isin(city_set) & od24["dest"].isin(ml_set)].copy()

    # Active destination villages across 2024 target months (02, 03, 08, 10, 11, 12)
    target_m = ["2024-02", "2024-03", "2024-08", "2024-10", "2024-11", "2024-12"]
    od_sub = od24[od24["mes"].isin(target_m)]
    active_dests = sorted(od_sub["dest"].unique())
    log(f"  Active destination municipalities from 137 cities in 2024 target months: {len(active_dests):,}")

    # Build balanced (city, active_dest) off-own-province skeleton with distances & roots
    c_arr = np.array(cities)
    d_arr = np.array(active_dests)
    c_prov = np.array([x[:2] for x in c_arr])
    d_prov = np.array([x[:2] for x in d_arr])
    cx = mun_info.loc[c_arr, "zone_x"].values
    cy = mun_info.loc[c_arr, "zone_y"].values
    dx = mun_info.loc[d_arr, "zone_x"].values
    dy = mun_info.loc[d_arr, "zone_y"].values

    ci_idx, di_idx = np.where(c_prov[:, None] != d_prov[None, :])
    dist_km = np.hypot(cx[ci_idx] - dx[di_idx], cy[ci_idx] - dy[di_idx]) / 1000.0
    base = pd.DataFrame({
        "orig": c_arr[ci_idx],
        "dest": d_arr[di_idx],
        "oprov": c_prov[ci_idx],
        "dprov": d_prov[di_idx],
        "lkm": np.log(np.maximum(dist_km, 1.0)),
    })
    base["cp_fe"] = base["orig"] + "_" + base["dprov"]
    base = base.merge(roots_cv, on=["orig", "dest"], how="left").fillna({"born_g1": 0.0, "born_g2": 0.0, "born_any": 0.0})
    base["lg1"] = np.log1p(base["born_g1"])
    base["lg2"] = np.log1p(base["born_g2"])
    base["lany"] = np.log1p(base["born_any"])
    base["has_g1"] = (base["born_g1"] > 0).astype(float)
    base["has_g2"] = (base["born_g2"] > 0).astype(float)
    base["loss_pct"] = base["dest"].map(mun_info["loss_pct"])
    base["share_np"] = base["dest"].map(mun_info["share_no_principal"])
    base["coastal"] = base["dest"].map(mun_info["coastal_true"]).astype(bool)
    base["dest_pop"] = base["dest"].map(mun_info["padron_2024"])
    log(f"  City x Village skeleton: {len(base):,} off-province cells ({len(cities)} cities x {len(active_dests)} dests)")

    # Helper to run pyfixest fepois / feols with (city x dest_prov) FE + village FE, clustered on dest_prov (G=47)
    def run_cv_spec(df_in, formula, method="ppml", vcov={"CRV1": "dprov"}):
        if method == "ppml":
            m = pf.fepois(formula, data=df_in, vcov=vcov, iwls_tol=1e-7, iwls_maxiter=40)
        else:
            m = pf.feols(formula, data=df_in, vcov=vcov)
        tb = m.tidy()
        return tb, int(m._N)

    # Pre-merge flows for each 2024 month + 2022/2023/2025 Aug & Feb
    month_flows = {}
    for mm in target_m:
        f = od24[od24["mes"] == mm].set_index(["orig", "dest"])["turistas"]
        month_flows[mm] = base.set_index(["orig", "dest"]).index.map(f).fillna(0.0).values
    for yr in ["2022", "2023", "2025"]:
        od_y = pd.read_parquet(INE_OD_DIR / f"ine_od_mun_{yr}.parquet")
        od_y = od_y[od_y["orig"].isin(city_set) & od_y["dest"].isin(ml_set)]
        for mm in [f"{yr}-02", f"{yr}-08"]:
            f = od_y[od_y["mes"] == mm].set_index(["orig", "dest"])["turistas"]
            month_flows[mm] = base.set_index(["orig", "dest"]).index.map(f).fillna(0.0).values

    cv_rows = []
    # Main monthly specs: (1) PPML all dests (lg1, lg2, joint lg1+lg2, lany),
    # (2) PPML inland dests only, (3) PPML small villages (pop < 5,000),
    # (4) Extensive margin LPM (1[flow >= 30]) to verify censoring robustness,
    # (5) Positive-only log-OLS (conditional on flow >= 30) to verify intensive margin!
    for mes, yvals in month_flows.items():
        base["y"] = yvals
        base["y_pos"] = (yvals >= 30).astype(float)
        is_main24 = mes.startswith("2024-")

        specs = [
            ("all_ppml_g1", base, "y ~ lg1 + lkm | cp_fe + dest", "ppml"),
            ("all_ppml_g2", base, "y ~ lg2 + lkm | cp_fe + dest", "ppml"),
            ("all_ppml_joint", base, "y ~ lg1 + lg2 + lkm | cp_fe + dest", "ppml"),
            ("all_ppml_any", base, "y ~ lany + lkm | cp_fe + dest", "ppml"),
            ("inland_ppml_any", base[~base["coastal"]], "y ~ lany + lkm | cp_fe + dest", "ppml"),
        ]
        if is_main24:
            specs += [
                ("inland_ppml_joint", base[~base["coastal"]], "y ~ lg1 + lg2 + lkm | cp_fe + dest", "ppml"),
                ("small5k_ppml_any", base[base["dest_pop"] < 5000], "y ~ lany + lkm | cp_fe + dest", "ppml"),
                ("ext_margin_lpm_any", base, "y_pos ~ lany + lkm | cp_fe + dest", "ols"),
                ("int_margin_pos_ols_any", base[base["y"] >= 30].assign(ly=lambda d: np.log(d["y"])), "ly ~ lany + lkm | cp_fe + dest", "ols"),
            ]
        if mes in ("2024-02", "2024-08", "2024-12"):
            # Interaction: do emptied municipalities & high second-home municipalities pull MORE from their diaspora cities?
            sub_i = base[base["loss_pct"].notna() & base["share_np"].notna() & ~base["coastal"]].copy()
            sub_i["lany_x_loss"] = sub_i["lany"] * sub_i["loss_pct"]
            sub_i["lany_x_shnp"] = sub_i["lany"] * sub_i["share_np"]
            specs += [
                ("inland_interact_loss", sub_i, "y ~ lany + lany_x_loss + lkm | cp_fe + dest", "ppml"),
                ("inland_interact_loss_shnp", sub_i, "y ~ lany + lany_x_loss + lany_x_shnp + lkm | cp_fe + dest", "ppml"),
            ]
        for sname, df_s, form, meth in specs:
            tb, n_obs = run_cv_spec(df_s, form, method=meth)
            for term in tb.index:
                if term == "lkm":
                    continue
                r = tb.loc[term]
                cv_rows.append({
                    "mes": mes, "spec": sname, "method": meth, "term": term, "n": n_obs,
                    "beta": float(r["Estimate"]), "se": float(r["Std. Error"]),
                    "t": float(r["t value"]), "p": float(r["Pr(>|t|)"]),
                    "ci_lo": float(r["2.5%"]), "ci_hi": float(r["97.5%"]),
                })

    # Stacked difference tests: holiday vs ordinary (Aug vs Feb, Aug vs Oct, Mar vs Feb, Nov vs Feb, Dec vs Feb)
    # Using exact stacked PPML with month-specific (cp_fe, dest, lkm) and roots x holiday interaction, clustered on dprov!
    for h_mes, o_mes, lbl in [
        ("2024-08", "2024-02", "Aug_vs_Feb"),
        ("2024-08", "2024-10", "Aug_vs_Oct"),
        ("2024-03", "2024-02", "Mar_vs_Feb"),
        ("2024-11", "2024-02", "Nov_vs_Feb"),
        ("2024-12", "2024-02", "Dec_vs_Feb"),
    ]:
        for subset_lbl, smask in [("all", np.ones(len(base), bool)), ("inland", ~base["coastal"].values)]:
            b_sub = base.loc[smask, ["orig", "dest", "dprov", "cp_fe", "lkm", "lg1", "lg2", "lany"]].copy()
            d0 = b_sub.copy()
            d0["y"] = month_flows[o_mes][smask]
            d0["hol"] = 0.0
            d0["cp_m"] = d0["cp_fe"] + "_0"
            d0["dest_m"] = d0["dest"] + "_0"
            d1 = b_sub.copy()
            d1["y"] = month_flows[h_mes][smask]
            d1["hol"] = 1.0
            d1["cp_m"] = d1["cp_fe"] + "_1"
            d1["dest_m"] = d1["dest"] + "_1"
            st = pd.concat([d0, d1], ignore_index=True)
            st["lkm_h"] = st["lkm"] * st["hol"]
            st["lany_h"] = st["lany"] * st["hol"]
            st["lg1_h"] = st["lg1"] * st["hol"]
            st["lg2_h"] = st["lg2"] * st["hol"]
            for sname, form in [
                (f"diff_{lbl}_{subset_lbl}_any", "y ~ lany + lany_h + lkm + lkm_h | cp_m + dest_m"),
                (f"diff_{lbl}_{subset_lbl}_g2", "y ~ lg2 + lg2_h + lkm + lkm_h | cp_m + dest_m"),
            ]:
                if lbl not in ("Aug_vs_Feb", "Aug_vs_Oct", "Dec_vs_Feb") and sname.endswith("_g2"):
                    continue
                tb, n_obs = run_cv_spec(st, form, method="ppml")
                for term in tb.index:
                    if term.endswith("_h"):
                        r = tb.loc[term]
                        cv_rows.append({
                            "mes": lbl, "spec": sname, "method": "ppml_diff", "term": term, "n": n_obs,
                            "beta": float(r["Estimate"]), "se": float(r["Std. Error"]),
                            "t": float(r["t value"]), "p": float(r["Pr(>|t|)"]),
                            "ci_lo": float(r["2.5%"]), "ci_hi": float(r["97.5%"]),
                        })
    cv_df = pd.DataFrame(cv_rows)
    cv_df.to_csv(OUT / "n1_city_village_roots.csv", index=False)
    log(f"  Saved n1_city_village_roots.csv ({len(cv_df)} rows)")

    # Within-destination-province permutation null (QAP-style) for City x Village roots!
    # Permute destination villages within each destination province -> preserves city x dest_prov margins & dest_prov structure!
    # Fast & exact on the demeaned log1p(flow) vs log1p(roots) residual after removing (cp_fe) and (dest) means & lkm
    rng_q = np.random.default_rng(seed_for("cv_qap"))
    n_perm = 1000
    # Index villages grouped by dprov (all active dests with >= 2 villages in their dprov)
    d_df = pd.DataFrame({"dest": d_arr, "dprov": d_prov, "di": np.arange(len(d_arr))})
    prov_groups = [g["di"].values for _, g in d_df.groupby("dprov") if len(g) >= 2]
    # Build dense 137 x N_dest matrices for fast vectorized permutation!
    M_any = np.zeros((len(c_arr), len(d_arr)))
    M_g2 = np.zeros((len(c_arr), len(d_arr)))
    M_lkm = np.log(np.maximum(np.hypot(cx[:, None] - dx[None, :], cy[:, None] - dy[None, :]) / 1000.0, 1.0))
    mask_off = c_prov[:, None] != d_prov[None, :]
    c_pos = {c: i for i, c in enumerate(c_arr)}
    d_pos = {d: i for i, d in enumerate(d_arr)}
    for r in roots_cv.itertuples(index=False):
        if r.orig in c_pos and r.dest in d_pos:
            M_any[c_pos[r.orig], d_pos[r.dest]] = np.log1p(r.born_any)
            M_g2[c_pos[r.orig], d_pos[r.dest]] = np.log1p(r.born_g2)

    def within_prov_demean(M):
        # Demean within each (city c, dest_prov p) AND subtract village mean (1 pass of alternating demeaning suffices to ~1e-3, do 5)
        Z = M.copy()
        Z[~mask_off] = 0.0
        for _ in range(8):
            # Subtract (city, dest_prov) mean
            for grp in prov_groups:
                sub = Z[:, grp]
                Z[:, grp] = sub - sub.mean(axis=1, keepdims=True)
            # Subtract column (village) mean over off-own-prov cities
            col_cnt = mask_off.sum(axis=0).clip(min=1)
            Z = np.where(mask_off, Z - (Z.sum(axis=0, keepdims=True) / col_cnt[None, :]), 0.0)
        return Z

    Z_lkm = within_prov_demean(M_lkm)
    Z_any = within_prov_demean(M_any)
    Z_g2 = within_prov_demean(M_g2)
    # Partial out Z_lkm from Z_any and Z_g2
    b_lkm_any = (Z_any[mask_off] @ Z_lkm[mask_off]) / (Z_lkm[mask_off] @ Z_lkm[mask_off])
    Z_any_r = np.where(mask_off, Z_any - b_lkm_any * Z_lkm, 0.0)
    b_lkm_g2 = (Z_g2[mask_off] @ Z_lkm[mask_off]) / (Z_lkm[mask_off] @ Z_lkm[mask_off])
    Z_g2_r = np.where(mask_off, Z_g2 - b_lkm_g2 * Z_lkm, 0.0)

    # Pre-build residualized flow matrices for 2024 months
    Z_Y = {}
    for mm in target_m:
        MY = np.zeros((len(c_arr), len(d_arr)))
        MY[ci_idx, di_idx] = np.log1p(month_flows[mm])
        zy = within_prov_demean(MY)
        by = (zy[mask_off] @ Z_lkm[mask_off]) / (Z_lkm[mask_off] @ Z_lkm[mask_off])
        Z_Y[mm] = np.where(mask_off, zy - by * Z_lkm, 0.0)

    qap_rows = []
    for r_lbl, ZR in [("any_g1g2", Z_any_r), ("strict_g2", Z_g2_r)]:
        rn = float(np.linalg.norm(ZR[mask_off]))
        # Test each month AND difference (holiday - Feb)
        tests = [(mm, Z_Y[mm]) for mm in target_m] + [
            (f"{hm}_minus_Feb", Z_Y[hm] - Z_Y["2024-02"]) for hm in ["2024-03", "2024-08", "2024-11", "2024-12"]
        ] + [("2024-08_minus_Oct", Z_Y["2024-08"] - Z_Y["2024-10"])]
        for tname, ZY in tests:
            yn = float(np.linalg.norm(ZY[mask_off]))
            r_obs = float(np.sum(ZY * ZR) / (yn * rn))
            null_r = np.empty(n_perm)
            for k in range(n_perm):
                perm = np.arange(len(d_arr))
                for grp in prov_groups:
                    perm[grp] = rng_q.permutation(grp)
                # Note: permuting columns within dprov keeps mask_off strictly invariant!
                ZP = ZR[:, perm]
                null_r[k] = float(np.sum(ZY * ZP) / (yn * rn))
            p_one = float((1 + np.sum(null_r >= r_obs)) / (n_perm + 1))
            qap_rows.append({
                "roots_type": r_lbl, "test": tname, "r_obs": r_obs,
                "null_mean": float(null_r.mean()), "null_sd": float(null_r.std(ddof=1)),
                "z_score": float((r_obs - null_r.mean()) / null_r.std(ddof=1)),
                "p_one_sided": p_one, "n_perm": n_perm,
            })
    pd.DataFrame(qap_rows).to_csv(OUT / "n1_city_village_qap.csv", index=False)
    log(f"  Saved n1_city_village_qap.csv ({len(qap_rows)} tests)")


# =============================================================================
# STAGE 5: 2nd-Generation Strict vs 1st-Gen in R2's Neighbouring-City Design
# =============================================================================
def stage_roots_gen2(mitma_nodes: pd.DataFrame):
    log("STAGE 5: Testing 2nd-generation strict roots in R2's exact neighbouring-city PPML design...")
    PROVS = np.array(sorted(mitma_nodes["prov"].unique()))
    P2I = {p: k for k, p in enumerate(PROVS)}
    NID = mitma_nodes["id"].values
    NPI = np.array([P2I[p] for p in mitma_nodes["prov"].values])
    NCOAST = mitma_nodes["coastal_true"].values.astype(bool)

    # Exact R2 city set (137 mainland cities, 112 identifying in 22 provinces)
    Bm = pd.read_csv(BIRTH_MUN_CSV, dtype={"ine_code": str, "prov_birth": str})
    Bm = Bm[Bm["year"] == 2024]
    CITIES = np.array(sorted(c for c in Bm["ine_code"].unique() if c[:2] not in EXCL))
    assert len(CITIES) == 137
    CPROV = np.array([c[:2] for c in CITIES])
    pcount = Counter(CPROV)
    IDROW = np.array([pcount[p] >= 2 for p in CPROV])
    assert IDROW.sum() == 112
    C2R = {c: i for i, c in enumerate(CITIES)}
    CITY_NODE = np.array([np.where(NID == c)[0][0] for c in CITIES])
    NODE2ROW = -np.ones(len(NID), dtype=int)
    NODE2ROW[CITY_NODE] = np.arange(len(CITIES))

    # 1st-gen roots (Censo 2024 68545)
    tot_res = Bm.groupby("ine_code")["persons"].sum()
    bm = Bm[Bm["prov_birth"].isin(PROVS)].set_index(["ine_code", "prov_birth"])["persons"]
    BORN1 = np.array([[bm.get((c, p), 0.0) for p in PROVS] for c in CITIES], dtype=float)
    TOT = tot_res.reindex(CITIES).values.astype(float)
    SHARE1 = BORN1 / TOT[:, None]

    # 2nd-gen STRICT roots (ECEPOV 2021: own-province born, parents born in province p != prov(c))
    dh = pd.read_parquet(ECEPOV_DIR / "ecepov_hogar_2021.parquet")
    for c in ["IDQ_PV", "IDQ_MUN", "PROVN", "PROVNP1", "PROVNP2"]:
        dh[c] = dh[c].astype(str).str.strip()
    ch = dh[dh["IDQ_MUN"].isin(set(CITIES))].copy()
    ec_tot = ch.groupby("IDQ_MUN")["FACTOR"].sum().reindex(CITIES).values
    own_b = ch[ch["PROVN"] == ch["IDQ_PV"]]
    s1 = own_b[own_b["PROVNP1"].isin(PROVS)].groupby(["IDQ_MUN", "PROVNP1"])["FACTOR"].sum() * 0.5
    s2 = own_b[own_b["PROVNP2"].isin(PROVS)].groupby(["IDQ_MUN", "PROVNP2"])["FACTOR"].sum() * 0.5
    s1.index.names = ["c", "p"]
    s2.index.names = ["c", "p"]
    g2_s = s1.add(s2, fill_value=0.0)
    BORN2 = np.array([[g2_s.get((c, p), 0.0) for p in PROVS] for c in CITIES], dtype=float)
    # Also ECEPOV 1st-gen (for apples-to-apples noise comparison!)
    g1_ec = ch[ch["PROVN"].isin(PROVS)].groupby(["IDQ_MUN", "PROVN"])["FACTOR"].sum()
    BORN1_EC = np.array([[g1_ec.get((c, p), 0.0) for p in PROVS] for c in CITIES], dtype=float)

    OWN = CPROV[:, None] == PROVS[None, :]
    MASK_C = (~OWN) & IDROW[:, None]
    ri_all, pj_all = np.where(MASK_C)

    # Distances & ADRH income (identical to R2)
    pw = mitma_nodes["population"].values.astype(float)
    PCX = np.array([np.average(mitma_nodes["x"].values[NPI == k], weights=pw[NPI == k]) for k in range(47)])
    PCY = np.array([np.average(mitma_nodes["y"].values[NPI == k], weights=pw[NPI == k]) for k in range(47)])
    CX, CY = mitma_nodes["x"].values[CITY_NODE], mitma_nodes["y"].values[CITY_NODE]
    KM_C = np.hypot(CX[:, None] - PCX[None, :], CY[:, None] - PCY[None, :]) / 1000.0
    lkm = np.log(KM_C[ri_all, pj_all])
    lkm_c = lkm - lkm.mean()

    A = pd.read_csv(ADRH_CSV, dtype={"ine_code": str})
    inc = A[(A["year"] == 2024) & (A["indicator"] == "Renta neta media por persona")].set_index("ine_code")["value"]
    LINC = np.log(inc.reindex(CITIES).values.astype(float))
    LINC = LINC - LINC.mean()
    # Income x destination-province interactions (46 cols, drop first prov as reference)
    inc_cols = []
    for k in range(1, 47):
        inc_cols.append(np.where(pj_all == k, LINC[ri_all], 0.0))
    X_ctrl = np.column_stack([lkm_c, lkm_c ** 2] + inc_cols)

    # Regressors (add tiny pseudocount to 2nd-gen share where ECEPOV sample has 0 cells: eps = min positive / 2)
    eps2 = float(BORN2[MASK_C & (BORN2 > 0)].min() / 2.0)
    eps1e = float(BORN1_EC[MASK_C & (BORN1_EC > 0)].min() / 2.0)
    x_g1 = np.log(SHARE1[ri_all, pj_all])
    x_g2 = np.log((BORN2[ri_all, pj_all] + eps2) / ec_tot[ri_all])
    x_g1_ec = np.log((BORN1_EC[ri_all, pj_all] + eps1e) / ec_tot[ri_all])

    # Fixed effect codes: city FE (ri_all) + (origin_prov x dest_prov) FE
    oprov_i = np.array([P2I[p] for p in CPROV[ri_all]])
    pair_fe = oprov_i * 47 + pj_all
    fe_codes = [ri_all, pair_fe]
    cl_oprov = oprov_i

    # Check collinearity between x_g1 and x_g2 (within FE!)
    hd_chk = HDFE([np.unique(c, return_inverse=True)[1] for c in fe_codes])
    w1 = np.ones(len(x_g1))
    r1_d = hd_chk.demean(x_g1, w1)
    r2_d = hd_chk.demean(x_g2, w1)
    log(f"  Within-FE corr(1st-gen Censo, 2nd-gen strict ECEPOV) = {np.corrcoef(r1_d, r2_d)[0,1]:.4f}; 2nd-gen zeros={(BORN2[MASK_C]==0).sum()}/{len(ri_all)}")

    # Load MITMA main/appendix nights (all + inland zones) AND INE 53000 monthly (2024M02, M03, M08, M10, M11, M12)
    datasets = {}
    for dlabel, ddate in [
        ("MITMA_Feb14", "2024-02-14"),
        ("MITMA_Aug15", "2024-08-15"),
        ("MITMA_Nov01", "2024-11-01"),
        ("MITMA_Dec25", "2024-12-25"),
        ("MITMA_Mar29", "2024-03-29"),
    ]:
        f = NET / f"edges_{ddate}__full.parquet"
        e = pd.read_parquet(f) if f.exists() else B0.build_edges(ddate, mitma_nodes)[0]
        i, j, w = e["i"].values, e["j"].values, e["w"].values
        r = NODE2ROW[i]
        pj = NPI[j]
        for kind, keep in [("all", np.ones(len(j), bool)), ("inland", ~NCOAST[j])]:
            M = np.zeros((137, 47))
            sel = (r >= 0) & keep
            np.add.at(M, (r[sel], pj[sel]), w[sel])
            datasets[f"{dlabel}_{kind}"] = M[ri_all, pj_all]

    # INE 53000 (137 cities x 47 provinces, monthly tourists)
    recs = []
    target_per = {"2024M02", "2024M03", "2024M08", "2024M10", "2024M11", "2024M12"}
    with gzip.open(INE_TM_RAW / "turismo_53000_munorig_x_provdest.csv.gz", "rt", encoding="utf-8") as f:
        f.readline()
        for line in f:
            p = line.rstrip("\n").split(";")
            if len(p) < 7 or p[5].strip() not in target_per:
                continue
            oc = p[0][:5]
            if oc not in C2R:
                continue
            dp = PROV_NAME_TO_CODE.get(p[3].strip())
            if dp not in P2I:
                continue
            v = parse_es_num(p[6])
            if v is not None:
                recs.append((p[5].strip(), C2R[oc], P2I[dp], v))
    d53 = pd.DataFrame(recs, columns=["per", "ri", "pj", "val"])
    assert len(d53) > 10000, f"53000 parse too small: {len(d53)}"
    for per in sorted(target_per):
        M = np.zeros((137, 47))
        sub = d53[d53["per"] == per]
        M[sub["ri"].values, sub["pj"].values] = sub["val"].values
        datasets[f"INE53000_{per}"] = M[ri_all, pj_all]

    rows = []
    for dname, yv in datasets.items():
        # Spec 1: 1st-gen only (replicates R2)
        # Spec 2: 2nd-gen strict only (ECEPOV)
        # Spec 3: Joint 1st-gen (Censo) + 2nd-gen strict (ECEPOV) -> focus on BOTH!
        # Spec 4: Joint 1st-gen (ECEPOV) + 2nd-gen strict (ECEPOV) -> equal measurement-error comparison!
        for sname, Xmat, fidx, term in [
            ("g1_censo_only", np.column_stack([x_g1, X_ctrl]), 0, "log_share_g1_censo"),
            ("g2_strict_only", np.column_stack([x_g2, X_ctrl]), 0, "log_share_g2_strict"),
            ("joint_censo_g1_plus_g2__g1", np.column_stack([x_g1, x_g2, X_ctrl]), 0, "log_share_g1_censo"),
            ("joint_censo_g1_plus_g2__g2", np.column_stack([x_g1, x_g2, X_ctrl]), 1, "log_share_g2_strict"),
            ("joint_ecepov_g1_plus_g2__g1", np.column_stack([x_g1_ec, x_g2, X_ctrl]), 0, "log_share_g1_ecepov"),
            ("joint_ecepov_g1_plus_g2__g2", np.column_stack([x_g1_ec, x_g2, X_ctrl]), 1, "log_share_g2_strict"),
        ]:
            # Full WCR + CI on headline specs; WCR p + CRV1 CI on secondary specs for speed
            is_main = dname in ("MITMA_Feb14_all", "MITMA_Aug15_all", "MITMA_Aug15_inland", "MITMA_Nov01_all", "MITMA_Dec25_all", "INE53000_2024M02", "INE53000_2024M08", "INE53000_2024M12")
            do_ci = is_main and sname in ("g2_strict_only", "joint_censo_g1_plus_g2__g1", "joint_censo_g1_plus_g2__g2", "joint_ecepov_g1_plus_g2__g2")
            nb = 9999 if do_ci else 1999
            res = ppml_crv1_and_wcr(yv, Xmat, fe_codes, cl_oprov, focus_idx=fidx, n_boot=nb, seed=seed_for("g2", dname, sname), do_ci=do_ci)
            rows.append({
                "dataset": dname, "spec": sname, "term": term,
                "n": res["n"], "G": res["G"],
                "beta": res["beta"], "se_crv1": res["se_crv1"], "p_crv1": res["p_crv1"],
                "ci_crv_lo": res["ci_crv_lo"], "ci_crv_hi": res["ci_crv_hi"],
                "p_wcr": res["p_wcr"], "ci_wcr_lo": res["ci_wcr_lo"], "ci_wcr_hi": res["ci_wcr_hi"],
                "n_boot": res["n_boot"],
            })
        log(f"    stage5 done {dname}")

    # Stacked difference tests for 2nd-gen strict (MITMA Aug15 - Feb14 all & inland; Dec25 - Feb14; INE M08 - M02, M08 - M10, M12 - M02)
    for h_key, o_key, lbl in [
        ("MITMA_Aug15_all", "MITMA_Feb14_all", "MITMA_Aug15_vs_Feb14_all"),
        ("MITMA_Aug15_inland", "MITMA_Feb14_inland", "MITMA_Aug15_vs_Feb14_inland"),
        ("MITMA_Dec25_all", "MITMA_Feb14_all", "MITMA_Dec25_vs_Feb14_all"),
        ("INE53000_2024M08", "INE53000_2024M02", "INE53000_Aug_vs_Feb"),
        ("INE53000_2024M08", "INE53000_2024M10", "INE53000_Aug_vs_Oct"),
        ("INE53000_2024M12", "INE53000_2024M02", "INE53000_Dec_vs_Feb"),
    ]:
        y_st = np.concatenate([datasets[o_key], datasets[h_key]])
        hol = np.concatenate([np.zeros(len(ri_all)), np.ones(len(ri_all))])
        fe_st = [
            np.concatenate([ri_all, ri_all + 200]),
            np.concatenate([pair_fe, pair_fe + 3000]),
        ]
        cl_st = np.concatenate([cl_oprov, cl_oprov])
        Xc_st = np.vstack([
            np.hstack([X_ctrl, np.zeros_like(X_ctrl)]),
            np.hstack([np.zeros_like(X_ctrl), X_ctrl]),
        ])
        for tname, xv in [("diff_g1_censo", x_g1), ("diff_g2_strict", x_g2)]:
            x_st = np.concatenate([xv, xv])
            X_st = np.column_stack([x_st * hol, x_st, Xc_st])
            res = ppml_crv1_and_wcr(y_st, X_st, fe_st, cl_st, focus_idx=0, n_boot=9999, seed=seed_for("g2diff", lbl, tname))
            rows.append({
                "dataset": lbl, "spec": tname, "term": f"{tname}_x_holiday",
                "n": res["n"], "G": res["G"],
                "beta": res["beta"], "se_crv1": res["se_crv1"], "p_crv1": res["p_crv1"],
                "p_wcr": res["p_wcr"], "ci_wcr_lo": res["ci_wcr_lo"], "ci_wcr_hi": res["ci_wcr_hi"],
                "n_boot": res["n_boot"],
            })

    out_df = pd.DataFrame(rows)
    out_df.to_csv(OUT / "n1_gen2_city_prov.csv", index=False)
    log(f"  Saved n1_gen2_city_prov.csv ({len(out_df)} rows)")


# =============================================================================
# STAGE 6: Second-Home Mediation (Censo 2021) + Household 'Casa del Pueblo' (ECEPOV)
# =============================================================================
def load_ine53001_totals(periods: set[str]) -> pd.DataFrame:
    recs = []
    with gzip.open(INE_TM_RAW / "turismo_53001_provorig_x_mundest.csv.gz", "rt", encoding="utf-8") as f:
        f.readline()
        for line in f:
            p = line.rstrip("\n").split(";")
            if len(p) < 7 or p[0].strip() != "Total Nacional" or p[1].strip() or p[2].strip():
                continue
            per = p[5].strip()
            if per not in periods:
                continue
            code = p[3][:5]
            if not re.match(r"^\d{5}$", code):
                continue
            v = parse_es_num(p[6])
            if v is not None:
                recs.append((code, per, v))
    return pd.DataFrame(recs, columns=["ine_code", "per", "tur"]).pivot(index="ine_code", columns="per", values="tur")


def load_ine53000_totals(periods: set[str]) -> pd.DataFrame:
    recs = []
    with gzip.open(INE_TM_RAW / "turismo_53000_munorig_x_provdest.csv.gz", "rt", encoding="utf-8") as f:
        f.readline()
        for line in f:
            p = line.rstrip("\n").split(";")
            if len(p) < 7 or p[1].strip() != "Total Nacional" or p[2].strip() or p[3].strip():
                continue
            per = p[5].strip()
            if per not in periods:
                continue
            code = p[0][:5]
            if not re.match(r"^\d{5}$", code):
                continue
            v = parse_es_num(p[6])
            if v is not None:
                recs.append((code, per, v))
    return pd.DataFrame(recs, columns=["ine_code", "per", "tur"]).pivot(index="ine_code", columns="per", values="tur")


def cluster_boot_mediation(y, x_loss, x_med, X_ctrl, fe, cl, n_boot=1999, seed=0):
    """Exact FE-OLS mediation (beta_base vs beta_with_med) + province-cluster pairs bootstrap for % mediated."""
    def fit_pair(y_, xl_, xm_, Xc_, fe_):
        _, fi = np.unique(fe_, return_inverse=True)
        n_fe = int(fi.max()) + 1
        cnt = np.bincount(fi, minlength=n_fe).astype(float)[:, None]

        def dm(A):
            s = np.zeros((n_fe, A.shape[1]))
            np.add.at(s, fi, A)
            return A - (s / cnt)[fi]

        M0 = dm(np.column_stack([y_, xl_, Xc_]))
        b0 = np.linalg.lstsq(M0[:, 1:], M0[:, 0], rcond=None)[0][0]
        M1 = dm(np.column_stack([y_, xl_, xm_, Xc_]))
        b1_vec = np.linalg.lstsq(M1[:, 1:], M1[:, 0], rcond=None)[0]
        b1, bm = b1_vec[0], b1_vec[1]
        return float(b0), float(b1), float(bm)

    # Point estimates + full WCB for base and mediated models
    r_base = fe_ols_wcb(y, x_loss, X_ctrl, fe, cl, n_boot=1999, seed=seed + 1)
    r_med = fe_ols_wcb(y, x_loss, np.column_stack([x_med, X_ctrl]), fe, cl, n_boot=1999, seed=seed + 2)
    r_sh = fe_ols_wcb(y, x_med, np.column_stack([x_loss, X_ctrl]), fe, cl, n_boot=1999, seed=seed + 3)

    b0, b1, bm = r_base["beta"], r_med["beta"], r_sh["beta"]
    assert abs(b0 - fit_pair(y, x_loss, x_med, X_ctrl, fe)[0]) < 1e-8
    med_diff = b0 - b1
    med_pct = med_diff / b0 if abs(b0) > 1e-9 else np.nan

    # Province-cluster pairs bootstrap for mediated share CI
    _, ci = np.unique(cl, return_inverse=True)
    G = int(ci.max()) + 1
    grp_idx = [np.where(ci == g)[0] for g in range(G)]
    rng = np.random.default_rng(seed)
    boot_diff = np.empty(n_boot)
    boot_pct = np.empty(n_boot)
    for k in range(n_boot):
        draw_g = rng.integers(0, G, size=G)
        idx = np.concatenate([grp_idx[g] for g in draw_g])
        # Assign distinct FE id per drawn cluster copy (since FE == province == cluster)
        fe_k = np.concatenate([np.full(len(grp_idx[g]), j, dtype=int) for j, g in enumerate(draw_g)])
        bb0, bb1, _ = fit_pair(y[idx], x_loss[idx], x_med[idx], X_ctrl[idx], fe_k)
        boot_diff[k] = bb0 - bb1
        boot_pct[k] = (bb0 - bb1) / bb0 if abs(bb0) > 1e-6 else np.nan

    return {
        "n": r_base["n"], "G": r_base["G"],
        "beta_base": b0, "se_base": r_base["se_crv1"], "ci_base_lo": r_base["ci_wcu_lo"], "ci_base_hi": r_base["ci_wcu_hi"], "p_base_wcr": r_base["p_wcr"],
        "beta_after_med": b1, "se_after": r_med["se_crv1"], "ci_after_lo": r_med["ci_wcu_lo"], "ci_after_hi": r_med["ci_wcu_hi"], "p_after_wcr": r_med["p_wcr"],
        "beta_mediator": bm, "se_med": r_sh["se_crv1"], "ci_med_lo": r_sh["ci_wcu_lo"], "ci_med_hi": r_sh["ci_wcu_hi"], "p_med_wcr": r_sh["p_wcr"],
        "mediated_diff": med_diff, "med_diff_ci_lo": float(np.nanquantile(boot_diff, 0.025)), "med_diff_ci_hi": float(np.nanquantile(boot_diff, 0.975)),
        "mediated_share": med_pct, "med_share_ci_lo": float(np.nanquantile(boot_pct, 0.025)), "med_share_ci_hi": float(np.nanquantile(boot_pct, 0.975)),
    }


def stage_mediation(mun: pd.DataFrame, mitma_nodes: pd.DataFrame):
    log("STAGE 6: Testing second-home mediation (Censo 2021) + household 'casa del pueblo' (ECEPOV)...")
    # --- Part A: Municipal stock mediation on INE 53001/53000 ---
    pers = {"2024M02", "2024M03", "2024M08", "2024M10", "2024M11", "2024M12"}
    T_in = load_ine53001_totals(pers)
    T_out = load_ine53000_totals(pers)

    nd = mun.copy()
    nd["lpop"] = np.log(nd["padron_2024"].clip(lower=1))
    nd["larea"] = np.log(nd["area_km2"].clip(lower=0.01))
    nd["coastal"] = nd["coastal_true"].astype(float)
    nd["ldcity"] = np.log1p(nd["zone_dist_city_km"])
    nd["ldcoast"] = np.log1p(nd["zone_dist_coast_km"])
    nd["l_np_per_res"] = np.log1p(nd["no_principal_per_res"])
    for p in pers:
        nd[f"in_{p}"] = nd["ine_code"].map(T_in[p]) if p in T_in else np.nan
        nd[f"out_{p}"] = nd["ine_code"].map(T_out[p]) if p in T_out else np.nan

    ok = nd["loss_pct"].notna() & nd["area_km2"].notna() & nd["share_no_principal"].notna() & (nd["padron_2024"] > 0)
    ctrl_m = ["lpop", "larea", "coastal", "ldcity", "ldcoast"]

    med_rows = []
    # Outcomes on INE 53001/53000
    ine_outcomes = [
        ("INE53001_log_Aug_Feb", np.log(nd["in_2024M08"] / nd["in_2024M02"]), ok),
        ("INE53001_log_Aug_Oct", np.log(nd["in_2024M08"] / nd["in_2024M10"]), ok),
        ("INE53001_log_Mar_Feb", np.log(nd["in_2024M03"] / nd["in_2024M02"]), ok),
        ("INE53001_log_Nov_Feb", np.log(nd["in_2024M11"] / nd["in_2024M02"]), ok),
        ("INE53001_log_Dec_Feb", np.log(nd["in_2024M12"] / nd["in_2024M02"]), ok),
        ("INE53001_Aug_minus_Feb_per_res_zero", (nd["in_2024M08"].fillna(0) - nd["in_2024M02"].fillna(0)) / nd["padron_2024"], ok),
        ("INE53001_Aug_minus_Feb_per_res_lb29", (nd["in_2024M08"].fillna(0) - nd["in_2024M02"].fillna(29)) / nd["padron_2024"], ok),
        ("INE53001_Aug_minus_Feb_per_res_ub29", (nd["in_2024M08"].fillna(29) - nd["in_2024M02"].fillna(0)) / nd["padron_2024"], ok),
        ("INE_net_Aug_minus_Feb_per_res", ((nd["in_2024M08"] - nd["out_2024M08"]) - (nd["in_2024M02"] - nd["out_2024M02"])) / nd["padron_2024"], ok),
    ]
    for oname, yv, smask in ine_outcomes:
        # Test share_no_principal (all municipalities, 59525)
        m = smask.values & np.isfinite(yv.values)
        d = nd.loc[m]
        res = cluster_boot_mediation(
            yv.values[m], d["loss_pct"].values, d["share_no_principal"].values,
            d[ctrl_m].values, d["prov"].values, d["prov"].values,
            n_boot=1999, seed=seed_for("med", oname, "shnp"),
        )
        res.update(geography="municipality_59525", outcome=oname, mediator="share_no_principal")
        med_rows.append(res)

        # Also test l_np_per_res on main Aug/Feb and Dec/Feb, and 59531 share_esporadico + share_secundaria_est on >=500 pop!
        if oname in ("INE53001_log_Aug_Feb", "INE53001_log_Dec_Feb", "INE_net_Aug_minus_Feb_per_res"):
            res2 = cluster_boot_mediation(
                yv.values[m], d["loss_pct"].values, d["l_np_per_res"].values,
                d[ctrl_m].values, d["prov"].values, d["prov"].values,
                n_boot=1999, seed=seed_for("med", oname, "lnpr"),
            )
            res2.update(geography="municipality_59525", outcome=oname, mediator="log1p_no_principal_per_res")
            med_rows.append(res2)

        if oname == "INE53001_log_Aug_Feb":
            for mcol in ["share_esporadico", "share_secundaria_est", "share_no_principal"]:
                m31 = ok.values & np.isfinite(yv.values) & nd[mcol].notna().values & nd["viv_vacias"].notna().values
                d31 = nd.loc[m31]
                res3 = cluster_boot_mediation(
                    yv.values[m31], d31["loss_pct"].values, d31[mcol].values,
                    d31[ctrl_m].values, d31["prov"].values, d31["prov"].values,
                    n_boot=1999, seed=seed_for("med31", oname, mcol),
                )
                res3.update(geography="municipality_59531_ge500", outcome=oname, mediator=mcol)
                med_rows.append(res3)

    # Also test second-home mediation on MITMA 2,435 zones (using R1's exact outcomes & controls!)
    zh = pd.read_csv(OUT / "n1_mitma_zone_housing.csv", dtype={"id": str, "prov": str})
    r1m = pd.read_parquet(OUT / "r1_node_metrics.parquet")
    r1m["id"] = r1m["id"].astype(str)
    r1m["date"] = r1m["date"].astype(str)
    xy = mitma_nodes[["x", "y"]].values / 1000.0
    Dm = np.hypot(xy[:, None, 0] - xy[None, :, 0], xy[:, None, 1] - xy[None, :, 1])
    big = mitma_nodes["population"].values >= 100_000
    coast = mitma_nodes["coastal_true"].values.astype(bool)
    zh["ldcity"] = np.log1p(Dm[:, big].min(axis=1))
    zh["ldcoast"] = np.log1p(Dm[:, coast].min(axis=1))
    zh["lpop"] = np.log(zh["population"])
    zh["lnmun"] = np.log(zh["n_municipalities"])
    zh["larea"] = np.log(zh["area_km2"])
    zh["coastal"] = zh["coastal_true"].astype(float)
    zh["artefact"] = zh["id"].isin(C.ARTEFACT_ZONES)
    ok_z = zh["loss_peak_2021"].notna() & (~zh["artefact"]) & zh["share_no_principal"].notna()
    zh["loss_pct"] = np.nan
    zh.loc[ok_z, "loss_pct"] = zh.loc[ok_z, "loss_peak_2021"].rank(pct=True, method="average")
    ctrl_z = ["lpop", "lnmun", "larea", "coastal", "ldcity", "ldcoast"]

    f_feb = r1m[r1m["date"] == "2024-02-14"].set_index("id")
    for ddate, dlbl in [("2024-08-15", "Aug15"), ("2024-11-01", "Nov01"), ("2024-12-25", "Dec25")]:
        f_h = r1m[r1m["date"] == ddate].set_index("id")
        if f_h.empty:
            continue
        for oname, yv in [
            (f"MITMA_{dlbl}_log_s_in_ratio", np.log((f_h.loc[zh["id"], "s_in"].values + 1) / (f_feb.loc[zh["id"], "s_in"].values + 1))),
            (f"MITMA_{dlbl}_net_bal_diff", f_h.loc[zh["id"], "net_bal"].values - f_feb.loc[zh["id"], "net_bal"].values),
            (f"MITMA_{dlbl}_net_per_res_diff", f_h.loc[zh["id"], "net_per_res"].values - f_feb.loc[zh["id"], "net_per_res"].values),
            (f"MITMA_{dlbl}_log_k_in_ratio", np.log((f_h.loc[zh["id"], "k_in_w10"].values + 1) / (f_feb.loc[zh["id"], "k_in_w10"].values + 1))),
        ]:
            m = ok_z.values & np.isfinite(yv)
            d = zh.loc[m]
            res = cluster_boot_mediation(
                yv[m], d["loss_pct"].values, d["share_no_principal"].values,
                d[ctrl_z].values, d["prov"].values, d["prov"].values,
                n_boot=1999, seed=seed_for("med_mitma", oname),
            )
            res.update(geography="mitma_zones_2435", outcome=oname, mediator="share_no_principal")
            med_rows.append(res)

    pd.DataFrame(med_rows).to_csv(OUT / "n1_second_home_mediation.csv", index=False)
    log(f"  Saved n1_second_home_mediation.csv ({len(med_rows)} specs)")

    # --- Part B: Household-level 'La Casa del Pueblo' direct test (ECEPOV 2021) ---
    dh = pd.read_parquet(ECEPOV_DIR / "ecepov_hogar_2021.parquet")
    dv = pd.read_parquet(ECEPOV_DIR / "ecepov_vivienda_2021.parquet")
    for c in ["IDQ_PV", "PROVN", "PROVNP1", "PROVNP2"]:
        dh[c] = dh[c].astype(str).str.strip()
    for c in ["IDQ_PV", "SEGUNRESI", "LUGSEGUNRESI", "OTRACCAA", "DIASUSA"]:
        dv[c] = dv[c].astype(str).str.strip()

    dh["ccaa_n"] = dh["PROVN"].map(PROV2CCAA_NUM)
    dh["ccaa_p1"] = dh["PROVNP1"].map(PROV2CCAA_NUM)
    dh["ccaa_p2"] = dh["PROVNP2"].map(PROV2CCAA_NUM)
    dv["ccaa_res"] = dv["IDQ_PV"].map(PROV2CCAA_NUM)

    # Households with second home in ANOTHER CCAA
    sh = dv[(dv["SEGUNRESI"] == "1") & (dv["LUGSEGUNRESI"] == "4") & (dv["OTRACCAA"].isin(set(PROV2CCAA_NUM.values())))].copy()
    # Collect per-household set of member birth CCAAs (g1) and parent birth CCAAs (g2)
    hm = dh[["IDEN", "ccaa_n", "ccaa_p1", "ccaa_p2"]].merge(sh[["IDEN", "OTRACCAA"]], on="IDEN", how="inner")
    hm["m_g1"] = hm["ccaa_n"] == hm["OTRACCAA"]
    hm["m_g2"] = (hm["ccaa_p1"] == hm["OTRACCAA"]) | (hm["ccaa_p2"] == hm["OTRACCAA"])
    hm["m_any"] = hm["m_g1"] | hm["m_g2"]
    hh_m = hm.groupby("IDEN")[["m_g1", "m_g2", "m_any"]].any().reset_index()
    sh = sh.merge(hh_m, on="IDEN", how="left")

    # Build lookup of each household's birth-CCAA bitmask (19 CCAAs -> uint32) for fast exact permutation null!
    ccaa_list = sorted(set(PROV2CCAA_NUM.values()))
    c2bit = {c: (1 << k) for k, c in enumerate(ccaa_list)}
    hm["b_g1"] = hm["ccaa_n"].map(c2bit).fillna(0).astype(np.int64)
    hm["b_g2"] = hm["ccaa_p1"].map(c2bit).fillna(0).astype(np.int64) | hm["ccaa_p2"].map(c2bit).fillna(0).astype(np.int64)
    hm["b_any"] = hm["b_g1"] | hm["b_g2"]
    hh_bits = hm.groupby("IDEN").agg(b_g1=("b_g1", lambda s: int(np.bitwise_or.reduce(s.values))),
                                     b_g2=("b_g2", lambda s: int(np.bitwise_or.reduce(s.values))),
                                     b_any=("b_any", lambda s: int(np.bitwise_or.reduce(s.values)))).reset_index()
    sh = sh.merge(hh_bits, on="IDEN", how="left")
    sh["dest_bit"] = sh["OTRACCAA"].map(c2bit).astype(np.int64)

    # Assert bitmask matches exact match columns
    assert ((sh["b_any"] & sh["dest_bit"]) > 0).equals(sh["m_any"].astype(bool))

    # Subgroups to summarise + province-cluster bootstrap 95% CIs + within-origin-province permutation null
    # (null: permute OTRACCAA across cross-CCAA second-home households of the SAME origin province -> preserves
    # origin province, exact destination-CCAA popularity from that origin province, and distance profile!)
    subgroups = [
        ("all_cross_ccaa", np.ones(len(sh), bool)),
        ("inland_dest_ccaa", sh["OTRACCAA"].isin(INLAND_CCAA_NUM).values),
        ("coastal_dest_ccaa", (~sh["OTRACCAA"].isin(INLAND_CCAA_NUM)).values),
        ("origin_Madrid", (sh["ccaa_res"] == "13").values),
        ("origin_Cataluna", (sh["ccaa_res"] == "09").values),
        ("origin_PaisVasco", (sh["ccaa_res"] == "16").values),
        ("origin_Madrid_to_inland", ((sh["ccaa_res"] == "13") & sh["OTRACCAA"].isin(INLAND_CCAA_NUM)).values),
        ("origin_Cataluna_to_inland", ((sh["ccaa_res"] == "09") & sh["OTRACCAA"].isin(INLAND_CCAA_NUM)).values),
        ("days_lt15", (sh["DIASUSA"] == "1").values),
        ("days_15_29", (sh["DIASUSA"] == "2").values),
        ("days_30_59", (sh["DIASUSA"] == "3").values),
        ("days_60plus", (sh["DIASUSA"] == "4").values),
    ]

    rng_ec = np.random.default_rng(seed_for("ecepov_sh"))
    n_boot = 1999
    n_perm = 1000
    ec_rows = []
    for sname, smask in subgroups:
        sub = sh.loc[smask].reset_index(drop=True)
        w = sub["FACTOR"].values
        wtot = float(w.sum())
        p_g1 = float((w * sub["m_g1"].values).sum() / wtot)
        p_g2 = float((w * sub["m_g2"].values).sum() / wtot)
        p_any = float((w * sub["m_any"].values).sum() / wtot)

        # Within-origin-province permutation null (preserves origin->destination CCAA flows & distances!)
        prov_grps = [g.index.values for _, g in sub.groupby("IDQ_PV") if len(g) >= 2]
        dbits = sub["dest_bit"].values
        bany = sub["b_any"].values
        bg1 = sub["b_g1"].values
        null_any = np.empty(n_perm)
        null_g1 = np.empty(n_perm)
        for k in range(n_perm):
            perm_d = dbits.copy()
            for grp in prov_grps:
                perm_d[grp] = rng_ec.permutation(dbits[grp])
            null_any[k] = float((w * ((bany & perm_d) > 0)).sum() / wtot)
            null_g1[k] = float((w * ((bg1 & perm_d) > 0)).sum() / wtot)

        # Bootstrap CI: province-cluster bootstrap when G >= 10, else household-level weighted bootstrap
        provs_u = sub["IDQ_PV"].unique()
        if len(provs_u) >= 10:
            ci_mode = f"cluster_prov_G{len(provs_u)}"
            pidx = [np.where(sub["IDQ_PV"].values == p)[0] for p in provs_u]
            G = len(provs_u)
            b_any_s = np.empty(n_boot)
            b_g1_s = np.empty(n_boot)
            for k in range(n_boot):
                draw = rng_ec.integers(0, G, size=G)
                ii = np.concatenate([pidx[g] for g in draw])
                wk = w[ii]
                b_any_s[k] = float((wk * sub["m_any"].values[ii]).sum() / wk.sum())
                b_g1_s[k] = float((wk * sub["m_g1"].values[ii]).sum() / wk.sum())
        else:
            ci_mode = f"hh_boot_n{len(sub)}"
            n_s = len(sub)
            b_any_s = np.empty(n_boot)
            b_g1_s = np.empty(n_boot)
            for k in range(n_boot):
                ii = rng_ec.integers(0, n_s, size=n_s)
                wk = w[ii]
                b_any_s[k] = float((wk * sub["m_any"].values[ii]).sum() / wk.sum())
                b_g1_s[k] = float((wk * sub["m_g1"].values[ii]).sum() / wk.sum())

        ec_rows.append({
            "subgroup": sname, "n_hh": len(sub), "weighted_hh": wtot, "ci_mode": ci_mode,
            "match_g1": p_g1, "g1_ci_lo": float(np.quantile(b_g1_s, 0.025)), "g1_ci_hi": float(np.quantile(b_g1_s, 0.975)),
            "match_parent": p_g2,
            "match_any": p_any, "any_ci_lo": float(np.quantile(b_any_s, 0.025)), "any_ci_hi": float(np.quantile(b_any_s, 0.975)),
            "null_any_mean": float(null_any.mean()), "null_any_sd": float(null_any.std(ddof=1)),
            "lift_any_over_null": float(p_any / null_any.mean()) if null_any.mean() > 0 else np.nan,
            "excess_any_pp": float((p_any - null_any.mean()) * 100.0),
            "p_perm_any": float((1 + np.sum(null_any >= p_any)) / (n_perm + 1)),
            "null_g1_mean": float(null_g1.mean()),
        })
    ec_df = pd.DataFrame(ec_rows)
    ec_df.to_csv(OUT / "n1_ecepov_casa_del_pueblo.csv", index=False)
    log("  Saved n1_ecepov_casa_del_pueblo.csv:\n" + ec_df[["subgroup", "n_hh", "match_g1", "match_parent", "match_any", "null_any_mean", "lift_any_over_null", "p_perm_any"]].to_string(index=False))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stages", default="nodes,network,validate,roots_mun,roots_gen2,mediation")
    args = ap.parse_args()
    stages = [s.strip() for s in args.stages.split(",") if s.strip()]
    log("Starting n1_new_data.py stages:", stages)

    mun, mitma_nodes = stage_nodes()
    if "network" in stages:
        stage_network(mun)
    if "validate" in stages:
        stage_validate(mun, mitma_nodes)
    if "roots_mun" in stages:
        stage_roots_mun(mun)
    if "roots_gen2" in stages:
        stage_roots_gen2(mitma_nodes)
    if "mediation" in stages:
        stage_mediation(mun, mitma_nodes)
    log("ALL REQUESTED STAGES COMPLETE.")


if __name__ == "__main__":
    main()
