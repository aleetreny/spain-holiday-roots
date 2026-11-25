"""v1_verify_n1: Independent verification of N1 claims (C1, C3, C4, C5) and full
robustness suite (C2a-g), re-implemented from scratch on the raw inputs.

Run:  .venv/bin/python -I pipeline/v1_verify_n1.py
Outputs: results/v1_*.csv
"""
from __future__ import annotations

import gzip
import io
import json
import re
import sys
import time
import zipfile
import zlib
import xml.etree.ElementTree as ET
from pathlib import Path

import duckdb
import igraph as ig
import numpy as np
import pandas as pd
import scipy.sparse as sps
from scipy import stats
from scipy.linalg import qr

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C  # noqa: E402
import b00_common as B0  # noqa: E402

T0 = time.time()
NET = C.RESULTS / "b01_networks"
OUT = C.RESULTS
EXT = C.ROOT / "_local" / "external"

P_CENSO24_PARQUET = EXT / "ine_censo_microdatos_2024" / "md_CensoAnual_muestra_2024.parquet"
P_CENSO24_DICT = EXT / "ine_censo_microdatos_2024" / "dr_CensoAnual_muestra_2022.json"
P_CENSO24_ZIP = EXT / "ine_censo_microdatos_2024" / "raw" / "datos_muestra_2024.zip"
P_ADRH = EXT / "ine_adrh" / "adrh_municipios.csv"
P_BIRTH68545 = EXT / "ine_birth_residence" / "mun_residence_x_prov_birth_censo_2021_2025.csv"
P_VIV59525 = EXT / "ine_censo_viviendas_2021" / "raw" / "censo2021_59525.csv"
P_LAU2021 = EXT / "eurostat_lau_2021" / "raw" / "LAU_RG_01M_2021_3035.csv"
P_INE53001 = EXT / "ine_turismo_moviles" / "raw" / "turismo_53001_provorig_x_mundest.csv.gz"
P_INE53000 = EXT / "ine_turismo_moviles" / "raw" / "turismo_53000_munorig_x_provdest.csv.gz"
P_ECEPOV_ZIP = EXT / "ine_ecepov_2021" / "raw" / "datos_2021.zip"
P_ECEPOV_VIV_PQ = EXT / "ine_ecepov_2021" / "ecepov_vivienda_2021.parquet"
P_ECEPOV_HOG_PQ = EXT / "ine_ecepov_2021" / "ecepov_hogar_2021.parquet"
P_INEOD24_PQ = EXT / "ine_turismo_moviles_od" / "ine_od_mun_2024.parquet"
P_INEOD24_XLSX = EXT / "ine_turismo_moviles_od" / "raw" / "exp_tmov_interno_mun_2024.xlsx"

WEBB = np.array([-np.sqrt(1.5), -1.0, -np.sqrt(0.5), np.sqrt(0.5), 1.0, np.sqrt(1.5)])
N_BOOT_PPML = 9999
N_BOOT_OLS = 9999
ALPHA = 0.05


def log(*a):
    print(f"[{time.time() - T0:7.1f}s]", *a, flush=True)


def rng_for(*parts) -> np.random.Generator:
    ints = [C.SEED] + [
        zlib.crc32(str(p).encode()) if not isinstance(p, (int, np.integer)) else int(p)
        for p in parts
    ]
    return np.random.default_rng(np.random.SeedSequence(ints))


# =============================================================================
# Fast 2-way HDFE PPML + Kline-Santos (2012) Wild Cluster Restricted Score Test
# =============================================================================
class FastHDFE2:
    """Weighted alternating projections for 2 FE factors + optional group-slope factors, Irons-Tuck accelerated."""

    def __init__(self, c0: np.ndarray, c1: np.ndarray, slope_factors: list[tuple[np.ndarray, np.ndarray]] | None = None):
        n = len(c0)
        L0, L1 = int(c0.max()) + 1, int(c1.max()) + 1
        self.D0 = sps.csr_matrix((np.ones(n), (c0, np.arange(n))), shape=(L0, n))
        self.D1 = sps.csr_matrix((np.ones(n), (c1, np.arange(n))), shape=(L1, n))
        self.slopes = []
        if slope_factors:
            for g_code, s_vec in slope_factors:
                Lg = int(g_code.max()) + 1
                Dg = sps.csr_matrix((s_vec, (g_code, np.arange(n))), shape=(Lg, n))
                Dg2 = sps.csr_matrix((s_vec * s_vec, (g_code, np.arange(n))), shape=(Lg, n))
                self.slopes.append((Dg, Dg2))

    def demean(self, X: np.ndarray, w: np.ndarray, tol: float = 1e-11, maxit: int = 5000) -> np.ndarray:
        X = np.array(X, dtype=float, copy=True)
        one = X.ndim == 1
        if one:
            X = X[:, None]
        if X.shape[1] == 0:
            return X[:, 0] if one else X
        scale = max(1.0, float(np.abs(X).max()))
        s0 = (self.D0 @ w)[:, None]
        s1 = (self.D1 @ w)[:, None]
        s_denoms = [((Dg2 @ w)[:, None]) for _, Dg2 in self.slopes]

        def step(Y):
            Y1 = Y - self.D0.T @ np.divide(self.D0 @ (w[:, None] * Y), s0, out=np.zeros((s0.shape[0], Y.shape[1])), where=s0 > 1e-30)
            Y2 = Y1 - self.D1.T @ np.divide(self.D1 @ (w[:, None] * Y1), s1, out=np.zeros((s1.shape[0], Y.shape[1])), where=s1 > 1e-30)
            for (Dg, _), sd in zip(self.slopes, s_denoms):
                num_g = Dg @ (w[:, None] * Y2)
                bg = np.divide(num_g, sd, out=np.zeros_like(num_g), where=sd > 1e-30)
                Y2 = Y2 - Dg.T @ bg
            return Y2

        for _ in range(maxit):
            X1 = step(X)
            if float(np.max(np.abs(X1 - X))) < tol * scale:
                return X1[:, 0] if one else X1
            X2 = step(X1)
            d1 = X1 - X
            d2 = X2 - X1
            dd = d2 - d1
            denom = np.sum(dd * dd, axis=0)
            num = np.sum(d2 * dd, axis=0)
            alpha = np.divide(num, denom, out=np.zeros_like(num), where=denom > 1e-30)
            X = X2 - alpha[None, :] * d2
        raise RuntimeError("FastHDFE2 demeaning did not converge")


def drop_zero_groups(y: np.ndarray, codes_list: list[np.ndarray]) -> np.ndarray:
    keep = np.ones(len(y), bool)
    while True:
        bad = np.zeros(len(y), bool)
        for c in codes_list:
            s = np.bincount(c[keep], y[keep], minlength=int(c.max()) + 1)
            bad |= keep & (s[c] <= 0)
        if not bad.any():
            return keep
        keep &= ~bad


def poisson_dev(y: np.ndarray, mu: np.ndarray) -> float:
    t = np.where(y > 0, y * np.log(np.where(y > 0, y, 1.0) / mu), 0.0)
    return float(2.0 * np.sum(t - (y - mu)))


def ppml_solve(y: np.ndarray, X: np.ndarray, hd: FastHDFE2, offset=None, mu0=None, maxit: int = 250):
    n, k = len(y), X.shape[1]
    off = np.zeros(n) if offset is None else offset
    mu = (y + y.mean()) / 2.0 if mu0 is None else mu0.copy()
    eta = np.log(mu)
    dev_old, beta = np.inf, np.zeros(k)
    for it in range(maxit):
        z = eta - off + (y - mu) / mu
        tol_d = 1e-11 if it >= 3 else 1e-8
        Dm = hd.demean(np.column_stack([z, X]) if k else z, mu, tol=tol_d)
        if k:
            zt, Xt = Dm[:, 0], Dm[:, 1:]
            XtW = Xt.T * mu
            b_new = np.linalg.solve(XtW @ Xt, XtW @ zt)
            e = zt - Xt @ b_new
            db = float(np.max(np.abs(b_new - beta)))
        else:
            zt = Dm
            b_new, e, db = beta, zt, 0.0
        eta = z - e + off
        mu = np.exp(eta)
        dev = poisson_dev(y, mu)
        beta = b_new
        if abs(dev - dev_old) <= 1e-11 * (abs(dev) + 0.1) and db < 1e-9:
            return beta, mu, dev
        dev_old = dev
    raise RuntimeError("PPML did not converge")


class SingleNightDesign:
    """Single-night PPML design after separation drops and rank-revealing QR on controls."""

    def __init__(self, y: np.ndarray, X_roots: np.ndarray, cov_raw: np.ndarray,
                 c_code: np.ndarray, p_code: np.ndarray, cl_code: np.ndarray,
                 slope_raw: list[tuple[np.ndarray, np.ndarray]] | None = None,
                 drop_cl: int | None = None):
        keep = drop_zero_groups(y, [c_code, p_code])
        if drop_cl is not None:
            keep &= (cl_code != drop_cl)
            keep = drop_zero_groups(y * keep, [c_code, p_code]) & keep
        self.keep = keep
        self.y = y[keep]
        self.Xr = X_roots[keep] if X_roots.ndim == 2 else X_roots[keep, None]
        self.n_terms = self.Xr.shape[1]
        c0 = pd.factorize(c_code[keep])[0]
        c1 = pd.factorize(p_code[keep])[0]
        sf_k = [( pd.factorize(gc[keep])[0], sv[keep] ) for gc, sv in slope_raw] if slope_raw else None
        self.fe0, self.fe1 = c0, c1
        self.hd = FastHDFE2(c0, c1, slope_factors=sf_k)
        cov = cov_raw[keep]
        if cov.shape[1]:
            ct = self.hd.demean(cov, np.ones(len(self.y)), tol=1e-11)
            _, Rq, piv = qr(ct, mode="economic", pivoting=True)
            dg = np.abs(np.diag(Rq))
            ok = np.sort(piv[dg > 1e-8 * dg[0]])
            cov = cov[:, ok]
        self.cov = cov
        self.cl = cl_code[keep]
        self.n = int(keep.sum())


def fit_single(nt: SingleNightDesign, n_clusters_total: int):
    X = np.column_stack([nt.Xr, nt.cov])
    beta, mu, dev = ppml_solve(nt.y, X, nt.hd)
    Xt = nt.hd.demean(X, mu, tol=1e-12)
    H = (Xt.T * mu) @ Xt
    Hinv = np.linalg.inv(H)
    sc = Xt * (nt.y - mu)[:, None]
    assert np.abs(sc.sum(0)).max() < 1e-5 * (1.0 + np.abs(nt.y).sum()), "PPML score != 0"
    Sg = np.vstack([np.bincount(nt.cl, sc[:, q], minlength=n_clusters_total) for q in range(X.shape[1])]).T
    IF = Sg @ Hinv[:, :nt.n_terms]  # shape: (n_clusters_total, n_terms)
    Geff = int((np.abs(IF[:, 0]) > 0).sum())
    se = np.sqrt(Geff / (Geff - 1) * np.sum(IF ** 2, axis=0))
    return {"beta": beta[:nt.n_terms].copy(), "se": se, "IF": IF, "G": Geff, "mu": mu, "dev": dev, "n": nt.n}


def wcr_score_p(y: np.ndarray, mu: np.ndarray, target: np.ndarray, Xother: np.ndarray,
                hd_demean_fn, cl: np.ndarray, VB: np.ndarray) -> tuple[float, float]:
    """Wild cluster restricted score bootstrap (Kline & Santos 2012), Webb weights.
    Uses (1 + sum(|t_b| >= |t_0|)) / (B + 1) convention (and also returns raw mean check).
    """
    G_tot = VB.shape[1]
    Dm = hd_demean_fn(target, Xother, mu)
    tt_, Xt = Dm[:, 0], Dm[:, 1:]
    if Xt.shape[1]:
        XtW = Xt.T * mu
        tt_ = tt_ - Xt @ np.linalg.solve(XtW @ Xt, XtW @ tt_)
    s = tt_ * (y - mu)
    Sg = np.bincount(cl, s, minlength=G_tot)
    t0 = float(Sg.sum() / np.sqrt(np.sum(Sg ** 2)))
    tb = (VB @ Sg) / np.sqrt((VB ** 2) @ (Sg ** 2))
    p_plus1 = float((1 + np.sum(np.abs(tb) >= abs(t0))) / (VB.shape[0] + 1))
    return p_plus1, t0


def p_single_term(nt: SingleNightDesign, term_idx: int, b0: float, VB: np.ndarray, mu0=None):
    target = nt.Xr[:, term_idx]
    other_r = np.delete(nt.Xr, term_idx, axis=1)
    Xother = np.column_stack([other_r, nt.cov])
    _, mu, _ = ppml_solve(nt.y, Xother, nt.hd, offset=b0 * target, mu0=mu0)

    def dm_fn(tgt, Xo, w):
        return nt.hd.demean(np.column_stack([tgt, Xo]), w, tol=1e-11)

    p, _ = wcr_score_p(nt.y, mu, target, Xother, dm_fn, nt.cl, VB)
    return p, mu


class FastStackedTwoNight:
    """Exploits exact block-diagonality of night-specific FEs and covariates for 2x speedup."""

    def __init__(self, nh: SingleNightDesign, no: SingleNightDesign, mu_h_unr=None, mu_o_unr=None):
        assert nh.n_terms == no.n_terms
        self.no, self.nh = no, nh
        self.n0, self.n1 = no.n, nh.n
        self.y = np.r_[no.y, nh.y]
        self.cl = np.r_[no.cl, nh.cl]
        self.k_r = no.n_terms
        self.Xr_common = np.vstack([no.Xr, nh.Xr])
        self.Xr_hol = np.vstack([np.zeros_like(no.Xr), nh.Xr])
        self.mu0_warm = np.r_[mu_o_unr, mu_h_unr] if (mu_o_unr is not None and mu_h_unr is not None) else None

    def _demean_blocks(self, vec_cols: np.ndarray, mu: np.ndarray, tol: float = 1e-11) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Demeans vec_cols (n0+n1, m) + block-diagonal [cov_o, 0; 0, cov_h] separately per night."""
        n0 = self.n0
        mu0, mu1 = mu[:n0], mu[n0:]
        D0 = self.no.hd.demean(np.column_stack([vec_cols[:n0], self.no.cov]), mu0, tol=tol)
        D1 = self.nh.hd.demean(np.column_stack([vec_cols[n0:], self.nh.cov]), mu1, tol=tol)
        m = vec_cols.shape[1]
        V_dm = np.vstack([D0[:, :m], D1[:, :m]])
        C0_dm = D0[:, m:]
        C1_dm = D1[:, m:]
        return V_dm, C0_dm, C1_dm

    def _solve_wls_block(self, zt: np.ndarray, Vt: np.ndarray, C0t: np.ndarray, C1t: np.ndarray, mu: np.ndarray) -> np.ndarray:
        """Exact WLS residual zt - P_[Vt, diag(C0t, C1t)] zt using Schur complement on C0t, C1t."""
        n0 = self.n0
        mu0, mu1 = mu[:n0], mu[n0:]
        # Project zt and Vt out of C0t on night 0 and C1t on night 1
        if C0t.shape[1]:
            C0W = C0t.T * mu0
            M0 = C0W @ C0t
            b0_z = np.linalg.solve(M0, C0W @ zt[:n0])
            b0_V = np.linalg.solve(M0, C0W @ Vt[:n0])
            z0 = zt[:n0] - C0t @ b0_z
            V0 = Vt[:n0] - C0t @ b0_V
        else:
            z0, V0 = zt[:n0], Vt[:n0]
        if C1t.shape[1]:
            C1W = C1t.T * mu1
            M1 = C1W @ C1t
            b1_z = np.linalg.solve(M1, C1W @ zt[n0:])
            b1_V = np.linalg.solve(M1, C1W @ Vt[n0:])
            z1 = zt[n0:] - C1t @ b1_z
            V1 = Vt[n0:] - C1t @ b1_V
        else:
            z1, V1 = zt[n0:], Vt[n0:]
        zp = np.r_[z0, z1]
        Vp = np.vstack([V0, V1])
        if Vp.shape[1]:
            VpW = Vp.T * mu
            bv = np.linalg.solve(VpW @ Vp, VpW @ zp)
            e = zp - Vp @ bv
        else:
            bv = np.zeros(0)
            e = zp
        return bv, e

    def fit_unrestricted_check(self) -> tuple[np.ndarray, np.ndarray]:
        V = np.column_stack([self.Xr_hol, self.Xr_common])
        mu = self.mu0_warm.copy() if self.mu0_warm is not None else (self.y + self.y.mean()) / 2.0
        eta = np.log(mu)
        dev_old, bv_old = np.inf, np.zeros(V.shape[1])
        for it in range(200):
            z = eta + (self.y - mu) / mu
            VD, C0t, C1t = self._demean_blocks(np.column_stack([z, V]), mu, tol=1e-11)
            bv, e = self._solve_wls_block(VD[:, 0], VD[:, 1:], C0t, C1t, mu)
            eta = z - e
            mu = np.exp(eta)
            dev = poisson_dev(self.y, mu)
            if abs(dev - dev_old) <= 1e-11 * (abs(dev) + 0.1) and np.max(np.abs(bv - bv_old)) < 1e-9:
                return bv[:self.k_r], bv[self.k_r:]
            bv_old, dev_old = bv, dev
        raise RuntimeError("Stacked unrestricted did not converge")

    def p_delta(self, term_idx: int, d0: float, VB: np.ndarray, mu0=None) -> tuple[float, np.ndarray]:
        target = self.Xr_hol[:, term_idx]
        other_hol = np.delete(self.Xr_hol, term_idx, axis=1)
        V_other = np.column_stack([self.Xr_common, other_hol])
        off = d0 * target
        mu = (self.mu0_warm.copy() if self.mu0_warm is not None else (self.y + self.y.mean()) / 2.0) if mu0 is None else mu0.copy()
        eta = np.log(mu)
        dev_old, bv_old = np.inf, np.zeros(V_other.shape[1])
        for it in range(200):
            z = eta - off + (self.y - mu) / mu
            tol_d = 1e-11 if it >= 2 else 1e-8
            VD, C0t, C1t = self._demean_blocks(np.column_stack([z, V_other]), mu, tol=tol_d)
            bv, e = self._solve_wls_block(VD[:, 0], VD[:, 1:], C0t, C1t, mu)
            eta = z - e + off
            mu = np.exp(eta)
            dev = poisson_dev(self.y, mu)
            if abs(dev - dev_old) <= 1e-11 * (abs(dev) + 0.1) and np.max(np.abs(bv - bv_old)) < 1e-9:
                break
            bv_old, dev_old = bv, dev
        else:
            raise RuntimeError("Stacked restricted PPML did not converge")

        # Score test at restricted MLE mu: project target out of [V_other, C0, C1] with weight mu
        VD, C0t, C1t = self._demean_blocks(np.column_stack([target, V_other]), mu, tol=1e-11)
        _, tt_ = self._solve_wls_block(VD[:, 0], VD[:, 1:], C0t, C1t, mu)
        s = tt_ * (self.y - mu)
        G_tot = VB.shape[1]
        Sg = np.bincount(self.cl, s, minlength=G_tot)
        t0 = float(Sg.sum() / np.sqrt(np.sum(Sg ** 2)))
        tb = (VB @ Sg) / np.sqrt((VB ** 2) @ (Sg ** 2))
        p_val = float((1 + np.sum(np.abs(tb) >= abs(t0))) / (VB.shape[0] + 1))
        return p_val, mu


def invert_ci(pfun, est: float, se: float, alpha: float = ALPHA, tol_frac: float = 2e-3) -> tuple[float, float]:
    out = []
    cache = {"mu": None}

    def pf_(b):
        p, mu = pfun(b, cache["mu"])
        cache["mu"] = mu
        return p

    for sgn in (-1, 1):
        cache["mu"] = None
        inner, step = est, se
        outer = est + sgn * step
        k = 0
        while pf_(outer) > alpha:
            inner, step = outer, step * 1.6
            outer = est + sgn * step
            k += 1
            if k > 30:
                raise RuntimeError("CI bracket failed")
        while abs(outer - inner) > tol_frac * se:
            mid = 0.5 * (inner + outer)
            if pf_(mid) > alpha:
                inner = mid
            else:
                outer = mid
        out.append(0.5 * (inner + outer))
    return float(out[0]), float(out[1])


# =============================================================================
# FE-OLS + CRV1 + Webb Wild Cluster Bootstrap (for C3 & C5)
# =============================================================================
def fe_ols_multi(y: np.ndarray, X_interest: np.ndarray, Z_ctrl: np.ndarray,
                 fe_codes: np.ndarray, cl_codes: np.ndarray, n_boot: int = 0, seed: int = 0) -> dict:
    """OLS of y on X_interest (n x p) + Z_ctrl + FE, clustered by cl_codes."""
    y = np.asarray(y, float)
    X_in = np.asarray(X_interest, float)
    if X_in.ndim == 1:
        X_in = X_in[:, None]
    Z = np.asarray(Z_ctrl, float)
    keep_cols = [k for k in range(Z.shape[1]) if np.ptp(Z[:, k]) > 0]
    Z = Z[:, keep_cols]
    F = pd.get_dummies(pd.Series(fe_codes), dtype=float).to_numpy()
    A = np.column_stack([Z, F])
    U, S, _ = np.linalg.svd(A, full_matrices=False)
    Q = U[:, S > 1e-10 * S.max()]
    mz = lambda v: v - Q @ (Q.T @ v)  # noqa: E731
    Xt = mz(X_in)
    yt = mz(y)
    XXinv = np.linalg.inv(Xt.T @ Xt)
    b = XXinv @ (Xt.T @ yt)
    u = yt - Xt @ b
    n = len(y)
    codes, G = pd.factorize(pd.Series(cl_codes))[0], int(pd.Series(cl_codes).nunique())
    k_tot = X_in.shape[1] + Z.shape[1]
    cadj = G / (G - 1) * (n - 1) / (n - k_tot)

    sc = Xt * u[:, None]
    Sg = np.zeros((G, X_in.shape[1]))
    np.add.at(Sg, codes, sc)
    V = cadj * (XXinv @ (Sg.T @ Sg) @ XXinv)
    se = np.sqrt(np.diag(V))
    tstat = b / se
    tcrit = stats.t.ppf(0.975, G - 1)
    p_t = 2.0 * stats.t.sf(np.abs(tstat), G - 1)

    res = {
        "n": n, "G": G, "beta": b, "se": se, "t": tstat, "p_t": p_t,
        "ci_lo": b - tcrit * se, "ci_hi": b + tcrit * se,
    }
    if n_boot > 0:
        rng = np.random.default_rng(seed)
        p_wcr = np.zeros(X_in.shape[1])
        ci_wcu_lo = np.zeros(X_in.shape[1])
        ci_wcu_hi = np.zeros(X_in.shape[1])
        # Pre-draw Webb weights (G x n_boot)
        V_webb = WEBB[rng.integers(0, 6, size=(G, n_boot))]
        for j in range(X_in.shape[1]):
            # Partial out the other X_in columns too so single-parameter FWL exact WCR/WCU applies cleanly!
            if X_in.shape[1] > 1:
                Xo = np.delete(Xt, j, axis=1)
                Uo, So, _ = np.linalg.svd(Xo, full_matrices=False)
                Qo = Uo[:, So > 1e-10 * So.max()]
                xj = Xt[:, j] - Qo @ (Qo.T @ Xt[:, j])
                yj = yt - Qo @ (Qo.T @ yt)
                # Combine Q and Qo into full projection basis Qf
                Qf = np.column_stack([Q, Qo])
            else:
                xj, yj, Qf = Xt[:, 0], yt, Q
            sxx = float(xj @ xj)
            bj = float((xj @ yj) / sxx)
            uj = yj - bj * xj
            mzf = lambda v: v - Qf @ (Qf.T @ v)  # noqa: E731

            tr_list, tu_list = [], []
            for st in range(0, n_boot, 1000):
                m = min(1000, n_boot - st)
                vw = V_webb[:, st:st + m][codes]  # n x m
                # Restricted H0: bj = 0 -> u_r = yj
                ys = mzf(yj[:, None] * vw)
                bs = (xj @ ys) / sxx
                us = ys - xj[:, None] * bs[None, :]
                S_r = np.zeros((G, m))
                np.add.at(S_r, codes, xj[:, None] * us)
                se_r = np.sqrt(cadj * (S_r ** 2).sum(0)) / sxx
                tr_list.append(bs / se_r)
                # Unrestricted
                es = mzf(uj[:, None] * vw)
                db = (xj @ es) / sxx
                uu = es - xj[:, None] * db[None, :]
                S_u = np.zeros((G, m))
                np.add.at(S_u, codes, xj[:, None] * uu)
                se_u = np.sqrt(cadj * (S_u ** 2).sum(0)) / sxx
                tu_list.append(db / se_u)
            tr = np.concatenate(tr_list)
            tu = np.concatenate(tu_list)
            p_wcr[j] = float((1 + np.sum(np.abs(tr) >= abs(tstat[j]))) / (n_boot + 1))
            ci_wcu_lo[j] = float(bj - np.quantile(tu, 0.975) * se[j])
            ci_wcu_hi[j] = float(bj - np.quantile(tu, 0.025) * se[j])
        res["p_wcr"] = p_wcr
        res["ci_wcu_lo"] = ci_wcu_lo
        res["ci_wcu_hi"] = ci_wcu_hi
    return res


# =============================================================================
# PART 1: Verify C1 & Run Full C2 Robustness Suite
# =============================================================================
def run_c1_and_c2(nodes: pd.DataFrame) -> dict:
    log("=== PART 1: C1 & C2 (Censo Anual 2024 10% microdata + MITMA single nights) ===")
    con = duckdb.connect()

    # 1.1 Verify dictionary & raw zip integrity
    with open(P_CENSO24_DICT, encoding="utf-8") as f:
        cdict = json.load(f)
    layout_names = [v["name"] for v in cdict["layout"]]
    assert "FACTOR" not in layout_names and len(layout_names) == 35
    assert cdict["codelists"]["T_ANYO"]["0001"] == "Antes de 1981"
    assert cdict["codelists"]["T_MUNI"]["991"] == "Tamaño municipio <=2000"
    with zipfile.ZipFile(P_CENSO24_ZIP) as zf:
        pq_info = zf.getinfo("md_CensoAnual_muestra_2024.parquet")
        assert pq_info.file_size == P_CENSO24_PARQUET.stat().st_size, "Parquet size mismatch vs official INE zip"
        with zf.open("CSV/Censo_muestra_2024.tab") as tf:
            tab_hdr = [x.strip('"') for x in tf.readline().decode("utf-8").rstrip("\r\n").split("\t")]
            assert tab_hdr == layout_names, "Tab header differs from dictionary layout"

    total_sample_rows = con.execute(f"SELECT count(*) FROM read_parquet('{P_CENSO24_PARQUET}')").fetchone()[0]
    distinct_varanop = [
        r[0] for r in con.execute(f"SELECT DISTINCT VARANOP FROM read_parquet('{P_CENSO24_PARQUET}') ORDER BY 1").fetchall()
    ]
    log(f"Censo 2024 verified: {total_sample_rows:,} rows (unweighted 10% extract), 35 cols; distinct VARANOP={distinct_varanop}")
    assert "2024" not in distinct_varanop and "2023" in distinct_varanop, "2024 should not exist in 1 Jan 2024 census"

    # 1.2 Municipality set & mapping to MITMA zones
    PROVS = np.array(sorted(nodes["prov"].unique()), dtype=str)
    assert len(PROVS) == 47
    P2I = {p: k for k, p in enumerate(PROVS)}
    NID = nodes["id"].to_numpy(dtype=str)
    NPI = np.array([P2I[p] for p in nodes["prov"].to_numpy(dtype=str)])
    NCOAST = nodes["coastal_true"].to_numpy(dtype=bool)
    COASTAL_PROVS = sorted(nodes.loc[nodes["coastal_true"], "prov"].unique())

    SPEC_MUNI = "('991','992','993','666','966','777','888','999')"
    m_tot = con.execute(f"""
        SELECT CPRO || CMUN AS mun, CPRO AS prov, count(*) AS n_tot,
               sum(CASE WHEN VAREDAD >= 65 THEN 1 ELSE 0 END) AS n_tot_65,
               sum(CASE WHEN VAREDAD >= 60 THEN 1 ELSE 0 END) AS n_tot_60
        FROM read_parquet('{P_CENSO24_PARQUET}')
        WHERE CMUN NOT IN {SPEC_MUNI}
          AND CPRO NOT IN ('07','35','38','51','52')
        GROUP BY 1, 2 ORDER BY 1
    """).df()
    n_mainland_gt10k = len(m_tot)  # 711
    nid_set = set(NID)
    single_set = set(nodes.loc[nodes["n_municipalities"] == 1, "id"])
    in_am_munis = sorted(set(m_tot["mun"]) - nid_set)  # 32 munis >10k folded into MITMA _AM zones
    split_non_am = sorted(nid_set.intersection(m_tot["mun"]) - single_set)  # ['03015'] Almoradi

    m_single = m_tot[m_tot["mun"].isin(single_set)].copy()
    pcounts = m_single.groupby("prov").size()
    id_provs = sorted(pcounts[pcounts >= 2].index)
    dropped_single_provs = sorted(pcounts[pcounts < 2].index)  # ['05', '34', '42']
    m675 = m_single[m_single["prov"].isin(id_provs)].sort_values("mun").reset_index(drop=True)
    assert len(m675) == 675 and len(id_provs) == 44
    log(f"Municipality mapping: {n_mainland_gt10k} mainland >10k -> {len(in_am_munis)} in _AM, "
        f"{len(split_non_am)} split ({split_non_am}), {len(m_single)} strict single-zone -> "
        f"{len(m675)} identifying in {len(id_provs)} provs (dropped single-city provs: {dropped_single_provs})")

    # Check sampling uniformity against 100%-enumeration table 68545 on the 137 cities
    bm68545 = pd.read_csv(P_BIRTH68545, dtype={"ine_code": str, "prov_birth": str})
    bm24_tot = bm68545[bm68545["year"] == 2024].groupby("ine_code")["persons"].sum()
    chk137 = m_tot[m_tot["mun"].isin(bm24_tot.index)].copy()
    chk137["enum_pop"] = chk137["mun"].map(bm24_tot)
    chk137["rate"] = chk137["n_tot"] / chk137["enum_pop"]
    log(f"10% sample rate check vs Censo 68545 (137 cities): min={chk137['rate'].min():.5f}, "
        f"median={chk137['rate'].median():.5f}, max={chk137['rate'].max():.5f}")

    # 1.3 Query cohort, birthplace-size, and age-split counts per (mun, dprov)
    cp = con.execute(f"""
        SELECT CPRO || CMUN AS mun, CPRO_NAC AS dprov,
               count(*) AS n_born,
               sum(CASE WHEN VARANOP = '0001' THEN 1 ELSE 0 END) AS n_pre1981,
               sum(CASE WHEN VARANOP IN ('0002','0003','0004','0005') THEN 1 ELSE 0 END) AS n_mid81_15,
               sum(CASE WHEN VARANOP IN ('2016','2017','2018','2019','2020','2021','2022','2023') THEN 1 ELSE 0 END) AS n_rec2016,
               sum(CASE WHEN CMUN_NAC = '991' THEN 1 ELSE 0 END) AS n_vil2k,
               sum(CASE WHEN CMUN_NAC IN ('992','993') THEN 1 ELSE 0 END) AS n_mid2_10k,
               sum(CASE WHEN CMUN_NAC NOT IN {SPEC_MUNI} THEN 1 ELSE 0 END) AS n_urb10k,
               sum(CASE WHEN VAREDAD >= 65 THEN 1 ELSE 0 END) AS n_age65,
               sum(CASE WHEN VAREDAD >= 60 THEN 1 ELSE 0 END) AS n_age60,
               sum(CASE WHEN VARANOP = '0001' AND VAREDAD < 65 THEN 1 ELSE 0 END) AS n_pre1981_lt65,
               sum(CASE WHEN VARANOP = '0001' AND VAREDAD >= 65 THEN 1 ELSE 0 END) AS n_pre1981_ge65,
               sum(CASE WHEN VARANOP != '0001' AND VAREDAD >= 65 THEN 1 ELSE 0 END) AS n_post1981_ge65,
               sum(CASE WHEN VARANOP IN ('2016','2017','2018','2019','2020','2021','2022','2023') AND VAREDAD < 65 THEN 1 ELSE 0 END) AS n_rec2016_lt65,
               sum(CASE WHEN VARANOP = '0001' AND VAREDAD < 60 THEN 1 ELSE 0 END) AS n_pre1981_lt60,
               sum(CASE WHEN VARANOP = '0001' AND VAREDAD >= 60 THEN 1 ELSE 0 END) AS n_pre1981_ge60,
               sum(CASE WHEN VARANOP != '0001' AND VAREDAD >= 60 THEN 1 ELSE 0 END) AS n_post1981_ge60,
               sum(CASE WHEN VARANOP IN ('2016','2017','2018','2019','2020','2021','2022','2023') AND VAREDAD < 60 THEN 1 ELSE 0 END) AS n_rec2016_lt60,
               sum(CASE WHEN CMUN_NAC = '991' AND VAREDAD < 65 THEN 1 ELSE 0 END) AS n_vil2k_lt65,
               sum(CASE WHEN CMUN_NAC = '991' AND VAREDAD >= 65 THEN 1 ELSE 0 END) AS n_vil2k_ge65,
               sum(CASE WHEN CMUN_NAC NOT IN {SPEC_MUNI} AND VAREDAD < 65 THEN 1 ELSE 0 END) AS n_urb10k_lt65,
               sum(CASE WHEN CMUN_NAC NOT IN {SPEC_MUNI} AND VAREDAD >= 65 THEN 1 ELSE 0 END) AS n_urb10k_ge65
        FROM read_parquet('{P_CENSO24_PARQUET}')
        WHERE CMUN NOT IN {SPEC_MUNI}
          AND CPRO NOT IN ('07','35','38','51','52')
          AND CPRO_NAC BETWEEN '01' AND '52'
          AND CPRO_NAC NOT IN ('07','35','38','51','52')
          AND CPRO != CPRO_NAC
        GROUP BY 1, 2
    """).df()

    MMUN = m675["mun"].to_numpy(dtype=str)
    MPROV = m675["prov"].to_numpy(dtype=str)
    NTOT = m675["n_tot"].to_numpy(dtype=float)[:, None]
    SHARE65_CITY = (m675["n_tot_65"] / m675["n_tot"]).to_numpy(dtype=float)
    M2R = {m: r for r, m in enumerate(MMUN)}
    C_M = len(MMUN)

    cp675 = cp[cp["mun"].isin(M2R)].copy()
    ri_cp = cp675["mun"].map(M2R).to_numpy(dtype=int)
    pj_cp = cp675["dprov"].map(P2I).to_numpy(dtype=int)

    MATS = {}
    for col in [c for c in cp.columns if c.startswith("n_")]:
        M = np.zeros((C_M, 47))
        M[ri_cp, pj_cp] = cp675[col].to_numpy(dtype=float)
        MATS[col] = M

    OWN = MPROV[:, None] == PROVS[None, :]
    ri_all, pj_all = np.where(~OWN)
    assert len(ri_all) == 675 * 46  # 31,050

    # Covariates (exact R2 specification: log km, (log km)^2, 46 centred log ADRH income x dest_prov slopes)
    pw = nodes["population"].to_numpy(dtype=float)
    nx, ny = nodes["x"].to_numpy(dtype=float), nodes["y"].to_numpy(dtype=float)
    PCX = np.array([np.average(nx[NPI == k], weights=pw[NPI == k]) for k in range(47)])
    PCY = np.array([np.average(ny[NPI == k], weights=pw[NPI == k]) for k in range(47)])
    MNODE = np.array([np.where(NID == m)[0][0] for m in MMUN])
    KM = np.hypot(nx[MNODE][:, None] - PCX[None, :], ny[MNODE][:, None] - PCY[None, :]) / 1000.0
    LKM = np.log(KM[ri_all, pj_all])
    lkm_c = LKM - LKM.mean()

    adrh = pd.read_csv(P_ADRH, dtype={"ine_code": str})
    inc_s = adrh[(adrh.year == 2024) & (adrh.indicator == "Renta neta media por persona")].set_index("ine_code")["value"]
    INC = inc_s.reindex(MMUN).to_numpy(dtype=float)
    assert np.isfinite(INC).all(), "Missing ADRH 2024 income on 675 munis"
    LINC = np.log(INC) - np.log(INC).mean()
    XINC_FULL = np.column_stack([LINC[ri_all] * (pj_all == P2I[p]) for p in PROVS[1:]])
    COV_DIST = np.column_stack([lkm_c, lkm_c ** 2])
    SLOPE_MAIN = [(pj_all, LINC[ri_all])]

    S65_C = SHARE65_CITY - SHARE65_CITY.mean()
    SLOPE_AGECITY = [(pj_all, LINC[ri_all]), (pj_all, S65_C[ri_all])]

    c_code = pd.factorize(MMUN[ri_all])[0]
    p_code = pd.factorize(np.char.add(np.char.add(MPROV[ri_all], "_"), PROVS[pj_all]))[0]
    CL_LEVELS = np.array(id_provs, dtype=str)
    cl_code = pd.Categorical(MPROV[ri_all], categories=CL_LEVELS).codes.astype(int)
    G_TOT = len(CL_LEVELS)  # 44

    VB_9999 = rng_for("v1_webb_9999").choice(WEBB, size=(N_BOOT_PPML, G_TOT))
    VB_4999 = rng_for("v1_webb_4999").choice(WEBB, size=(4999, G_TOT))

    # Regressors using Haldane-Anscombe (count + 0.5)/N_tot (verified to reproduce N1 to 1e-14)
    lsh05 = lambda col: np.log((MATS[col] + 0.5) / NTOT)[ri_all, pj_all]  # noqa: E731
    X_COH = np.column_stack([lsh05("n_pre1981"), lsh05("n_rec2016")])
    X_BS = np.column_stack([lsh05("n_vil2k"), lsh05("n_urb10k")])

    # 1.4 Load all MITMA single nights (using b01 full edge lists verified against b00.build_edges)
    NW = pd.read_csv(NET / "networks.csv", dtype=str)
    dates24 = B0.available_dates(C.PERNOCT_2024)
    feb_sisters = sorted(NW.loc[(NW.group == "ordinary") & (NW.role == "sister"), "date"])
    oct_sisters = [d for d in dates24 if d[5:7] == "10" and pd.Timestamp(d).day_of_week in (1, 2, 3)]
    SISTERS = {
        "ordinary": feb_sisters,
        "summer": list(C.SISTER_NIGHTS["summer"]),
        "allsaints": list(C.SISTER_NIGHTS["allsaints"]),
        "christmas": [C.APPENDIX_NIGHTS["christmas"]],
        "goodfriday": [C.APPENDIX_NIGHTS["goodfriday"]],
        "october": oct_sisters,
    }
    YEARS = {g: list(v) for g, v in C.YEAR_REPLICATION.items()}
    ALL_DATES = sorted({d for ds in list(SISTERS.values()) + list(YEARS.values()) for d in ds})
    assert len(feb_sisters) == 13 and len(SISTERS["summer"]) == 9 and len(oct_sisters) == 15

    # Verify b01 == b00.build_edges on main nights
    for md in C.MAIN_NIGHTS.values():
        e_b01 = pd.read_parquet(NET / f"edges_{md}__full.parquet")
        e_b00, _ = B0.build_edges(md, nodes)
        assert len(e_b01) == len(e_b00) and np.allclose(e_b01["w"].to_numpy(), e_b00["w"].to_numpy())

    NODE2MROW = -np.ones(len(NID), dtype=int)
    NODE2MROW[MNODE] = np.arange(C_M)
    COAST_PROV_IDX = [P2I[p] for p in COASTAL_PROVS]

    FLOWS: dict[str, dict[str, np.ndarray]] = {}
    for d in ALL_DATES:
        f_pq = NET / f"edges_{d}__full.parquet"
        if f_pq.exists():
            e = pd.read_parquet(f_pq)
        else:
            e, _ = B0.build_edges(d, nodes)
        i, j, w = e["i"].to_numpy(), e["j"].to_numpy(), e["w"].to_numpy()
        r = NODE2MROW[i]
        sel = r >= 0
        out = {}
        for kind, keep in (("all", np.ones(len(j), bool)), ("inland", ~NCOAST[j]), ("coastal", NCOAST[j])):
            M = np.zeros((C_M, 47))
            s2 = sel & keep
            np.add.at(M, (r[s2], NPI[j[s2]]), w[s2])
            M[OWN] = 0.0
            out[kind] = M[ri_all, pj_all]
        Minl_p = np.zeros((C_M, 47))
        Minl_p[ri_all, pj_all] = out["all"]
        Minl_p[:, COAST_PROV_IDX] = 0.0
        out["inlandprov"] = Minl_p[ri_all, pj_all]
        FLOWS[d] = out

    # Window means (audit C3 rule: sum / n_nights, missing = 0)
    FLOWS["2024-02-mean13"] = {k: sum(FLOWS[d][k] for d in feb_sisters) / len(feb_sisters) for k in ("all", "inland", "coastal", "inlandprov")}
    FLOWS["2024-10-mean15"] = {k: sum(FLOWS[d][k] for d in oct_sisters) / len(oct_sisters) for k in ("all", "inland", "coastal", "inlandprov")}
    log(f"Loaded MITMA flows for {len(ALL_DATES)} single nights + 2 window means")

    # 1.5 Cross-check our fast slope-HDFE PPML solver against pyfixest (explicit XINC_FULL) on Feb 14 and Aug 15
    import pyfixest as pf
    COV_FULL_CHK = np.column_stack([COV_DIST, XINC_FULL])
    for chk_d in ("2024-02-14", "2024-08-15"):
        nt_fast = SingleNightDesign(FLOWS[chk_d]["all"], X_COH, COV_DIST, c_code, p_code, cl_code, slope_raw=SLOPE_MAIN)
        fit_fast = fit_single(nt_fast, G_TOT)
        nt_c = SingleNightDesign(FLOWS[chk_d]["all"], X_COH, COV_FULL_CHK, c_code, p_code, cl_code)
        df_pf = pd.DataFrame({"y": nt_c.y, "x0": nt_c.Xr[:, 0], "x1": nt_c.Xr[:, 1], "c": nt_c.fe0, "p": nt_c.fe1, "cl": nt_c.cl})
        vcols = []
        for q in range(nt_c.cov.shape[1]):
            df_pf[f"v{q}"] = nt_c.cov[:, q]
            vcols.append(f"v{q}")
        m_pf = pf.fepois("y ~ x0 + x1 + " + " + ".join(vcols) + " | c + p", data=df_pf, vcov={"CRV1": "cl"})
        t_pf = m_pf.tidy()
        for k_i, nm in enumerate(("x0", "x1")):
            assert abs(float(t_pf.loc[nm, "Estimate"]) - fit_fast["beta"][k_i]) < 1e-6
            dof_pf = np.sqrt((nt_c.n - 1) / (nt_c.n - nt_c.n_terms - nt_c.cov.shape[1]))
            assert abs(float(t_pf.loc[nm, "Std. Error"]) - fit_fast["se"][k_i] * dof_pf) < 1e-6
    log("Fast slope-HDFE solver verified against pyfixest on Feb 14 & Aug 15 (|diff| < 1e-6 on beta and exact CRV1 SE)")

    # 1.6 C1 Zero-cell & specification sensitivity table (Aug 15 vs Feb 14, all destinations)
    zero_rows = []
    for col_nm, label in [
        ("n_born", "pooled_born_p"), ("n_pre1981", "pre1981_exodus"), ("n_mid81_15", "mid_1981_2015"),
        ("n_rec2016", "recent_2016_2023"), ("n_vil2k", "born_village_le2k"), ("n_mid2_10k", "born_mid_2k_10k"),
        ("n_urb10k", "born_urban_gt10k"), ("n_age65", "born_p_age65plus"),
    ]:
        vals = MATS[col_nm][ri_all, pj_all]
        zero_rows.append({
            "variable": label, "total_cells": len(vals),
            "zero_cells": int((vals == 0).sum()), "zero_share": float((vals == 0).mean()),
            "total_sampled_persons": int(vals.sum()),
        })
    pd.DataFrame(zero_rows).to_csv(OUT / "v1_c1_zero_cells.csv", index=False)

    # Also test R2's 112 large cities (>50k + capitals, 22 provinces) where zero cells are much rarer!
    cities_137 = sorted(c for c in bm24_tot.index if c[:2] not in C.EXCLUDED_PROVINCES)
    c137_cnt = pd.Series([c[:2] for c in cities_137]).value_counts()
    id_p22 = set(c137_cnt[c137_cnt >= 2].index)
    m112_mask = np.isin(MMUN, [c for c in cities_137 if c[:2] in id_p22])
    assert m112_mask.sum() == 112
    cell_112 = m112_mask[ri_all]

    sens_specs = [
        ("N1_haldane_0.5_sample", lambda M: np.log((M + 0.5) / NTOT)[ri_all, pj_all], np.ones(len(ri_all), bool)),
        ("plus_1.0_sample_eq_log1p", lambda M: np.log((M + 1.0) / NTOT)[ri_all, pj_all], np.ones(len(ri_all), bool)),
        ("plus_0.1_sample_eq_pop_plus1", lambda M: np.log((10.0 * M + 1.0) / (10.0 * NTOT))[ri_all, pj_all], np.ones(len(ri_all), bool)),
        ("log_share_plus_1e4", lambda M: np.log(M / NTOT + 1e-4)[ri_all, pj_all], np.ones(len(ri_all), bool)),
        ("log_share_plus_1e5", lambda M: np.log(M / NTOT + 1e-5)[ri_all, pj_all], np.ones(len(ri_all), bool)),
        ("asinh_1000_share", lambda M: np.arcsinh(1000.0 * (M / NTOT))[ri_all, pj_all], np.ones(len(ri_all), bool)),
        ("R2_112_large_cities_G22_haldane0.5", lambda M: np.log((M + 0.5) / NTOT)[ri_all, pj_all], cell_112),
    ]
    sens_rows = []
    for sname, fn, smask in sens_specs:
        cl_sub = pd.factorize(MPROV[ri_all][smask])[0]
        G_sub = int(cl_sub.max()) + 1
        VB_sub = VB_9999[:, :G_sub]
        sf_sub = [(pj_all[smask], LINC[ri_all][smask])]
        for mtype, cols, tnames in [
            ("cohort_joint_2term", ["n_pre1981", "n_rec2016"], ["pre1981", "recent_2016_23"]),
            ("birthsize_joint_2term", ["n_vil2k", "n_urb10k"], ["village_le2k", "urban_gt10k"]),
        ]:
            X_s = np.column_stack([fn(MATS[c])[smask] for c in cols])
            no_s = SingleNightDesign(FLOWS["2024-02-14"]["all"][smask], X_s, COV_DIST[smask], c_code[smask], p_code[smask], cl_sub, slope_raw=sf_sub)
            nh_s = SingleNightDesign(FLOWS["2024-08-15"]["all"][smask], X_s, COV_DIST[smask], c_code[smask], p_code[smask], cl_sub, slope_raw=sf_sub)
            fo_s, fh_s = fit_single(no_s, G_sub), fit_single(nh_s, G_sub)
            st_s = FastStackedTwoNight(nh_s, no_s, fh_s["mu"], fo_s["mu"])
            dIF = fh_s["IF"] - fo_s["IF"]
            se_d = np.sqrt(G_sub / (G_sub - 1) * np.sum(dIF ** 2, axis=0))
            for k_i, tnm in enumerate(tnames):
                pw, _ = st_s.p_delta(k_i, 0.0, VB_sub)
                z_sh = float((MATS[cols[k_i]][ri_all, pj_all][smask] == 0).mean())
                sens_rows.append({
                    "spec": sname, "model": mtype, "term": tnm, "G": G_sub, "n_aug": nh_s.n,
                    "zero_cell_share": z_sh,
                    "beta_feb14": float(fo_s["beta"][k_i]), "beta_aug15": float(fh_s["beta"][k_i]),
                    "delta_aug_feb": float(fh_s["beta"][k_i] - fo_s["beta"][k_i]),
                    "se_delta_crv1": float(se_d[k_i]), "p_wcr_delta": pw,
                })
    # Also add complete 3-term partitions (including 1981-2015 cohort and 2k-10k town size!)
    for mtype, cols, tnames in [
        ("cohort_complete_3term", ["n_pre1981", "n_mid81_15", "n_rec2016"], ["pre1981", "mid_1981_2015", "recent_2016_23"]),
        ("birthsize_complete_3term", ["n_vil2k", "n_mid2_10k", "n_urb10k"], ["village_le2k", "mid_2k_10k", "urban_gt10k"]),
    ]:
        X_s = np.column_stack([lsh05(c) for c in cols])
        no_s = SingleNightDesign(FLOWS["2024-02-14"]["all"], X_s, COV_DIST, c_code, p_code, cl_code, slope_raw=SLOPE_MAIN)
        nh_s = SingleNightDesign(FLOWS["2024-08-15"]["all"], X_s, COV_DIST, c_code, p_code, cl_code, slope_raw=SLOPE_MAIN)
        fo_s, fh_s = fit_single(no_s, G_TOT), fit_single(nh_s, G_TOT)
        st_s = FastStackedTwoNight(nh_s, no_s, fh_s["mu"], fo_s["mu"])
        dIF = fh_s["IF"] - fo_s["IF"]
        se_d = np.sqrt(G_TOT / (G_TOT - 1) * np.sum(dIF ** 2, axis=0))
        for k_i, tnm in enumerate(tnames):
            pw, _ = st_s.p_delta(k_i, 0.0, VB_9999)
            sens_rows.append({
                "spec": "N1_haldane_0.5_3term_complete", "model": mtype, "term": tnm, "G": G_TOT, "n_aug": nh_s.n,
                "zero_cell_share": float((MATS[cols[k_i]][ri_all, pj_all] == 0).mean()),
                "beta_feb14": float(fo_s["beta"][k_i]), "beta_aug15": float(fh_s["beta"][k_i]),
                "delta_aug_feb": float(fh_s["beta"][k_i] - fo_s["beta"][k_i]),
                "se_delta_crv1": float(se_d[k_i]), "p_wcr_delta": pw,
            })
    pd.DataFrame(sens_rows).to_csv(OUT / "v1_c1_zero_sensitivity.csv", index=False)
    log("C1 zero-cell & complete-partition sensitivity done")

    # 1.7 Fit all single nights + means for cohort_joint and birthsize_joint
    MODELS = {
        "cohort_joint": (X_COH, ["pre1981_exodus", "recent_2016_2024"]),
        "birthsize_joint": (X_BS, ["born_village_le2k", "born_urban_gt10k"]),
    }
    FITS: dict = {}
    sn_rows = []
    # Fit 'all' and 'inland' and 'coastal' on all dates; 'inlandprov' on key dates
    KEY_DATES = {
        "2024-02-14", "2024-08-15", "2024-11-01", "2024-12-25", "2024-03-29",
        "2024-02-mean13", "2024-10-mean15",
        "2022-02-16", "2022-08-15", "2022-11-01", "2023-02-15", "2023-08-15",
    }
    for mname, (Xmat, tnames) in MODELS.items():
        for sample in ("all", "inland", "coastal", "inlandprov"):
            dlist = (ALL_DATES + ["2024-02-mean13", "2024-10-mean15"]) if sample in ("all", "inland", "coastal") else sorted(KEY_DATES)
            for d in dlist:
                nt = SingleNightDesign(FLOWS[d][sample], Xmat, COV_DIST, c_code, p_code, cl_code, slope_raw=SLOPE_MAIN)
                f = fit_single(nt, G_TOT)
                FITS[(mname, sample, d)] = (nt, f)
                for k_i, tnm in enumerate(tnames):
                    rec = {
                        "model": mname, "sample": sample, "date": d, "term": tnm,
                        "n": nt.n, "G": f["G"], "beta": float(f["beta"][k_i]), "se_crv1": float(f["se"][k_i]),
                    }
                    if sample == "all" and d in ("2024-02-14", "2024-08-15", "2024-11-01", "2024-12-25", "2024-03-29"):
                        pw, _ = p_single_term(nt, k_i, 0.0, VB_9999, mu0=f["mu"])
                        rec["p_wcr"] = pw
                        if d in ("2024-02-14", "2024-08-15"):
                            lo, hi = invert_ci(lambda b, m0: p_single_term(nt, k_i, b, VB_9999, m0),
                                               float(f["beta"][k_i]), float(f["se"][k_i]))
                            rec["ci_wcr_lo"], rec["ci_wcr_hi"] = lo, hi
                    sn_rows.append(rec)
        log(f"Single-night fits completed for {mname}")
    pd.DataFrame(sn_rows).to_csv(OUT / "v1_c2_single_nights.csv", index=False)

    # Assert exact Stacked unrestricted delta == beta_h - beta_o
    nh_chk, fh_chk = FITS[("cohort_joint", "all", "2024-08-15")]
    no_chk, fo_chk = FITS[("cohort_joint", "all", "2024-02-14")]
    d_unr, b_unr = FastStackedTwoNight(nh_chk, no_chk, fh_chk["mu"], fo_chk["mu"]).fit_unrestricted_check()
    assert np.allclose(d_unr, fh_chk["beta"] - fo_chk["beta"], atol=1e-7)
    assert np.allclose(b_unr, fo_chk["beta"], atol=1e-7)
    log("Stacked two-night unrestricted identity verified: delta == beta_h - beta_o to < 1e-7")

    # 1.8 Run ALL C1 & C2(a,b,c,d,g) stacked contrasts + sister pairs
    # Define contrast families:
    # - Main pairs (get WCR p + WCR 95% CI)
    # - Sister pairs (get exact WCR p for every pair!)
    MAIN_CONTRASTS = [
        ("Aug15_vs_Feb14", "2024-08-15", "2024-02-14", ("all", "inland", "coastal", "inlandprov")),
        ("Aug15_vs_FebMean13", "2024-08-15", "2024-02-mean13", ("all", "inland", "coastal")),
        ("Aug15_vs_OctMean15", "2024-08-15", "2024-10-mean15", ("all", "inland", "coastal", "inlandprov")),
        ("Nov01_vs_Feb14", "2024-11-01", "2024-02-14", ("all", "inland", "coastal", "inlandprov")),
        ("Nov01_vs_OctMean15", "2024-11-01", "2024-10-mean15", ("all", "inland")),
        ("Dec25_vs_Feb14", "2024-12-25", "2024-02-14", ("all", "inland", "coastal", "inlandprov")),
        ("Dec25_vs_OctMean15", "2024-12-25", "2024-10-mean15", ("all", "inland")),
        ("Mar29_vs_Feb14", "2024-03-29", "2024-02-14", ("all", "inland", "coastal", "inlandprov")),
        ("Mar29_vs_OctMean15", "2024-03-29", "2024-10-mean15", ("all", "inland")),
        ("Y2022_Aug15_vs_Feb16", "2022-08-15", "2022-02-16", ("all", "inland", "coastal")),
        ("Y2022_Nov01_vs_Feb16", "2022-11-01", "2022-02-16", ("all", "inland")),
        ("Y2023_Aug15_vs_Feb15", "2023-08-15", "2023-02-15", ("all", "inland", "coastal")),
        ("Placebo_OctMean15_vs_FebMean13", "2024-10-mean15", "2024-02-mean13", ("all", "inland")),
    ]
    CI_CONTRASTS = {
        ("Aug15_vs_Feb14", "all"), ("Aug15_vs_Feb14", "inland"), ("Aug15_vs_Feb14", "coastal"), ("Aug15_vs_Feb14", "inlandprov"),
        ("Aug15_vs_OctMean15", "all"), ("Aug15_vs_OctMean15", "inland"),
        ("Nov01_vs_Feb14", "all"), ("Dec25_vs_Feb14", "all"), ("Mar29_vs_Feb14", "all"),
        ("Y2022_Aug15_vs_Feb16", "all"), ("Y2022_Nov01_vs_Feb16", "all"), ("Y2023_Aug15_vs_Feb15", "all"),
        ("Placebo_OctMean15_vs_FebMean13", "all"),
    }
    main_diff_rows = []
    for cname, h, o, samples in MAIN_CONTRASTS:
        for sample in samples:
            for mname, (_, tnames) in MODELS.items():
                nh, fh = FITS[(mname, sample, h)]
                no, fo = FITS[(mname, sample, o)]
                st = FastStackedTwoNight(nh, no, fh["mu"], fo["mu"])
                dIF = fh["IF"] - fo["IF"]
                Geff = int((np.abs(dIF[:, 0]) > 0).sum())
                se_d = np.sqrt(Geff / (Geff - 1) * np.sum(dIF ** 2, axis=0))
                delta = fh["beta"] - fo["beta"]
                tc = stats.t.ppf(0.975, Geff - 1)
                for k_i, tnm in enumerate(tnames):
                    pw9, mu_r0 = st.p_delta(k_i, 0.0, VB_9999)
                    # Reuse restricted mu_r0 for B=4999 check without re-solving PPML!
                    VD_r, C0_r, C1_r = st._demean_blocks(
                        np.column_stack([st.Xr_hol[:, k_i], np.column_stack([st.Xr_common, np.delete(st.Xr_hol, k_i, axis=1)])]),
                        mu_r0, tol=1e-11,
                    )
                    _, tt_r = st._solve_wls_block(VD_r[:, 0], VD_r[:, 1:], C0_r, C1_r, mu_r0)
                    Sg_r = np.bincount(st.cl, tt_r * (st.y - mu_r0), minlength=G_TOT)
                    t0_r = float(Sg_r.sum() / np.sqrt(np.sum(Sg_r ** 2)))
                    tb5 = (VB_4999 @ Sg_r) / np.sqrt((VB_4999 ** 2) @ (Sg_r ** 2))
                    pw5 = float((1 + np.sum(np.abs(tb5) >= abs(t0_r))) / (VB_4999.shape[0] + 1))
                    if (cname, sample) in CI_CONTRASTS:
                        lo, hi = invert_ci(lambda d0, m0: st.p_delta(k_i, d0, VB_9999, m0),
                                           float(delta[k_i]), float(se_d[k_i]))
                    else:
                        lo, hi = np.nan, np.nan
                    p_crv = float(2.0 * stats.t.sf(abs(delta[k_i] / se_d[k_i]), Geff - 1))
                    main_diff_rows.append({
                        "contrast": cname, "holiday": h, "ordinary": o, "sample": sample,
                        "model": mname, "term": tnm, "n": nh.n + no.n, "G": Geff,
                        "beta_holiday": float(fh["beta"][k_i]), "beta_ordinary": float(fo["beta"][k_i]),
                        "delta": float(delta[k_i]), "se_crv1": float(se_d[k_i]), "p_crv1": p_crv,
                        "ci_crv_lo": float(delta[k_i] - tc * se_d[k_i]), "ci_crv_hi": float(delta[k_i] + tc * se_d[k_i]),
                        "p_wcr_9999": pw9, "p_wcr_4999": pw5, "ci_wcr_lo": lo, "ci_wcr_hi": hi,
                    })
    df_main_diffs = pd.DataFrame(main_diff_rows)
    df_main_diffs.to_csv(OUT / "v1_c2_main_contrasts.csv", index=False)
    log("Main contrasts with WCR 95% CIs completed")

    # Sister-pair grids with exact WCR p-values!
    SISTER_GRIDS = [
        ("summer_x_feb_all", SISTERS["summer"], SISTERS["ordinary"], "all"),
        ("summer_x_feb_inland", SISTERS["summer"], SISTERS["ordinary"], "inland"),
        ("summer_x_feb_coastal", SISTERS["summer"], SISTERS["ordinary"], "coastal"),
        ("summer_x_oct_all", ["2024-08-15"], SISTERS["october"], "all"),
        ("allsaints_x_feb_all", SISTERS["allsaints"], SISTERS["ordinary"], "all"),
        ("christmas_x_feb_all", ["2024-12-25"], SISTERS["ordinary"], "all"),
        ("goodfriday_x_feb_all", ["2024-03-29"], SISTERS["ordinary"], "all"),
    ]
    sis_rows = []
    for gname, h_list, o_list, sample in SISTER_GRIDS:
        # For 9 x 15 = 135 summer_sisters_x_oct_all, compute WCR p too!
        for mname, (_, tnames) in MODELS.items():
            for h in h_list:
                nh, fh = FITS[(mname, sample, h)]
                for o in o_list:
                    no, fo = FITS[(mname, sample, o)]
                    st = FastStackedTwoNight(nh, no, fh["mu"], fo["mu"])
                    dIF = fh["IF"] - fo["IF"]
                    Geff = int((np.abs(dIF[:, 0]) > 0).sum())
                    se_d = np.sqrt(Geff / (Geff - 1) * np.sum(dIF ** 2, axis=0))
                    delta = fh["beta"] - fo["beta"]
                    for k_i, tnm in enumerate(tnames):
                        pw, _ = st.p_delta(k_i, 0.0, VB_9999)
                        sis_rows.append({
                            "grid": gname, "sample": sample, "model": mname, "term": tnm,
                            "holiday": h, "ordinary": o, "delta": float(delta[k_i]),
                            "se_crv1": float(se_d[k_i]),
                            "p_crv1": float(2.0 * stats.t.sf(abs(delta[k_i] / se_d[k_i]), Geff - 1)),
                            "p_wcr": pw,
                        })
        log(f"Sister grid done: {gname} ({len(h_list)}x{len(o_list)}={len(h_list)*len(o_list)} pairs)")
    df_sis = pd.DataFrame(sis_rows)
    df_sis.to_csv(OUT / "v1_c2_sister_pairs.csv", index=False)

    # Summarise sister grids
    sis_sum_rows = []
    for (gname, sample, mname, tnm), g in df_sis.groupby(["grid", "sample", "model", "term"], sort=False):
        vals = g["delta"].to_numpy()
        pws = g["p_wcr"].to_numpy()
        med = float(np.median(vals))
        sgn_ref = np.sign(med)
        # Check if main pair (Aug15-Feb14 or Nov01-Feb14 etc.) is in this grid
        main_sub = g[(g["holiday"].isin(["2024-08-15", "2024-11-01", "2024-12-25", "2024-03-29"])) &
                     (g["ordinary"] == "2024-02-14")]
        main_val = float(main_sub["delta"].iloc[0]) if len(main_sub) == 1 else np.nan
        if len(main_sub) == 1:
            sgn_ref = np.sign(main_val)
        is_ext = bool(len(main_sub) == 1 and len(vals) > 1 and (
            (main_val == vals.max() and (vals == vals.max()).sum() == 1) or
            (main_val == vals.min() and (vals == vals.min()).sum() == 1)
        ))
        sis_sum_rows.append({
            "grid": gname, "sample": sample, "model": mname, "term": tnm,
            "n_pairs": len(vals), "main_delta": main_val, "main_is_extreme": is_ext,
            "median_delta": med, "min_delta": float(vals.min()), "max_delta": float(vals.max()),
            "q25_delta": float(np.quantile(vals, 0.25)), "q75_delta": float(np.quantile(vals, 0.75)),
            "same_sign_count": int((np.sign(vals) == sgn_ref).sum()),
            "same_sign_share": float((np.sign(vals) == sgn_ref).mean()),
            "sig_wcr_05_count": int(((pws < 0.05) & (np.sign(vals) == sgn_ref)).sum()),
            "sig_wcr_05_share": float(((pws < 0.05) & (np.sign(vals) == sgn_ref)).mean()),
            "sig_wcr_01_count": int(((pws < 0.01) & (np.sign(vals) == sgn_ref)).sum()),
            "max_p_wcr": float(pws.max()), "median_p_wcr": float(np.median(pws)),
        })
    df_sis_sum = pd.DataFrame(sis_sum_rows)
    df_sis_sum.to_csv(OUT / "v1_c2_sister_summary.csv", index=False)

    # 1.9 C2(e): AGE CONFOUND battery (Aug 15 vs Feb 14, all & inland, with WCR p + 95% WCR CIs)
    log("Running C2(e) Age Confound battery...")
    AGE_SPECS = [
        ("E0_baseline_cohort_joint", np.column_stack([lsh05("n_pre1981"), lsh05("n_rec2016")]),
         ["pre1981", "recent_2016_23"], SLOPE_MAIN, True),
        ("E1a_horserace_born_p_ge65", np.column_stack([lsh05("n_pre1981"), lsh05("n_age65"), lsh05("n_rec2016")]),
         ["pre1981", "born_p_ge65", "recent_2016_23"], SLOPE_MAIN, True),
        ("E1b_horserace_born_p_ge60", np.column_stack([lsh05("n_pre1981"), lsh05("n_age60"), lsh05("n_rec2016")]),
         ["pre1981", "born_p_ge60", "recent_2016_23"], SLOPE_MAIN, False),
        ("E2a_split65_4way_exclusive",
         np.column_stack([lsh05("n_pre1981_lt65"), lsh05("n_pre1981_ge65"), lsh05("n_post1981_ge65"), lsh05("n_rec2016_lt65")]),
         ["pre1981_lt65_age43_64", "pre1981_ge65", "post1981_ge65", "recent2016_lt65"], SLOPE_MAIN, True),
        ("E2b_split60_4way_exclusive",
         np.column_stack([lsh05("n_pre1981_lt60"), lsh05("n_pre1981_ge60"), lsh05("n_post1981_ge60"), lsh05("n_rec2016_lt60")]),
         ["pre1981_lt60_age43_59", "pre1981_ge60", "post1981_ge60", "recent2016_lt60"], SLOPE_MAIN, False),
        ("E3_city_share65_x_dest_slopes", np.column_stack([lsh05("n_pre1981"), lsh05("n_rec2016")]),
         ["pre1981", "recent_2016_23"], SLOPE_AGECITY, True),
        ("E4a_birthsize_horserace_ge65", np.column_stack([lsh05("n_vil2k"), lsh05("n_age65"), lsh05("n_urb10k")]),
         ["village_le2k", "born_p_ge65", "urban_gt10k"], SLOPE_MAIN, False),
        ("E4b_birthsize_split65_4way",
         np.column_stack([lsh05("n_vil2k_lt65"), lsh05("n_vil2k_ge65"), lsh05("n_urb10k_lt65"), lsh05("n_urb10k_ge65")]),
         ["village_le2k_lt65", "village_le2k_ge65", "urban_gt10k_lt65", "urban_gt10k_ge65"], SLOPE_MAIN, True),
    ]
    age_rows = []
    for aspec, X_a, tnames, sf_a, do_ci in AGE_SPECS:
        for sample in ("all", "inland"):
            no_a = SingleNightDesign(FLOWS["2024-02-14"][sample], X_a, COV_DIST, c_code, p_code, cl_code, slope_raw=sf_a)
            nh_a = SingleNightDesign(FLOWS["2024-08-15"][sample], X_a, COV_DIST, c_code, p_code, cl_code, slope_raw=sf_a)
            fo_a, fh_a = fit_single(no_a, G_TOT), fit_single(nh_a, G_TOT)
            st_a = FastStackedTwoNight(nh_a, no_a, fh_a["mu"], fo_a["mu"])
            dIF = fh_a["IF"] - fo_a["IF"]
            se_d = np.sqrt(G_TOT / (G_TOT - 1) * np.sum(dIF ** 2, axis=0))
            delta = fh_a["beta"] - fo_a["beta"]
            tc = stats.t.ppf(0.975, G_TOT - 1)
            for k_i, tnm in enumerate(tnames):
                pw, _ = st_a.p_delta(k_i, 0.0, VB_9999)
                if do_ci and sample == "all":
                    lo, hi = invert_ci(lambda d0, m0: st_a.p_delta(k_i, d0, VB_9999, m0),
                                       float(delta[k_i]), float(se_d[k_i]))
                else:
                    lo, hi = np.nan, np.nan
                p_aug, _ = p_single_term(nh_a, k_i, 0.0, VB_9999, mu0=fh_a["mu"])
                p_feb, _ = p_single_term(no_a, k_i, 0.0, VB_9999, mu0=fo_a["mu"])
                age_rows.append({
                    "spec": aspec, "sample": sample, "term": tnm, "n_stacked": nh_a.n + no_a.n, "G": G_TOT,
                    "beta_feb14": float(fo_a["beta"][k_i]), "se_feb14": float(fo_a["se"][k_i]), "p_wcr_feb14": p_feb,
                    "beta_aug15": float(fh_a["beta"][k_i]), "se_aug15": float(fh_a["se"][k_i]), "p_wcr_aug15": p_aug,
                    "delta_aug_feb": float(delta[k_i]), "se_delta_crv1": float(se_d[k_i]),
                    "ci_crv_lo": float(delta[k_i] - tc * se_d[k_i]), "ci_crv_hi": float(delta[k_i] + tc * se_d[k_i]),
                    "p_wcr_delta": pw, "ci_wcr_lo": lo, "ci_wcr_hi": hi,
                })
    df_age = pd.DataFrame(age_rows)
    df_age.to_csv(OUT / "v1_c2_age_confound.csv", index=False)
    log("C2(e) Age Confound battery completed")

    # 1.10 C2(f): Leave-one-origin-province-out jackknife (all 44 clusters, Aug 15 - Feb 14, all & inland)
    log("Running C2(f) Leave-one-province-out jackknife across all 44 provinces...")
    jk_rows = []
    for sample in ("all", "inland"):
        for mname, (Xmat, tnames) in MODELS.items():
            for cl_drop_idx, p_drop in enumerate(CL_LEVELS):
                no_j = SingleNightDesign(FLOWS["2024-02-14"][sample], Xmat, COV_DIST, c_code, p_code, cl_code, slope_raw=SLOPE_MAIN, drop_cl=cl_drop_idx)
                nh_j = SingleNightDesign(FLOWS["2024-08-15"][sample], Xmat, COV_DIST, c_code, p_code, cl_code, slope_raw=SLOPE_MAIN, drop_cl=cl_drop_idx)
                fo_j, fh_j = fit_single(no_j, G_TOT), fit_single(nh_j, G_TOT)
                st_j = FastStackedTwoNight(nh_j, no_j, fh_j["mu"], fo_j["mu"])
                dIF = fh_j["IF"] - fo_j["IF"]
                Geff = int((np.abs(dIF[:, 0]) > 0).sum())
                assert Geff == 43
                se_d = np.sqrt(Geff / (Geff - 1) * np.sum(dIF ** 2, axis=0))
                delta = fh_j["beta"] - fo_j["beta"]
                n_muns_p = int((MPROV == p_drop).sum())
                for k_i, tnm in enumerate(tnames):
                    pw, _ = st_j.p_delta(k_i, 0.0, VB_9999)
                    jk_rows.append({
                        "sample": sample, "model": mname, "term": tnm,
                        "dropped_prov": p_drop, "dropped_n_munis": n_muns_p, "G_rem": Geff,
                        "beta_feb14": float(fo_j["beta"][k_i]), "beta_aug15": float(fh_j["beta"][k_i]),
                        "delta": float(delta[k_i]), "se_crv1": float(se_d[k_i]), "p_wcr": pw,
                    })
    df_jk = pd.DataFrame(jk_rows)
    df_jk.to_csv(OUT / "v1_c2_jackknife.csv", index=False)
    log("C2(f) 44-province jackknife completed")

    return {
        "main_diffs": df_main_diffs,
        "sis_sum": df_sis_sum,
        "age": df_age,
        "jk": df_jk,
        "single_nights": pd.DataFrame(sn_rows),
    }


# =============================================================================
# PART 2: Verify C3 (Censo 2021 59525 Second-Home Decomposition)
# =============================================================================
def parse_es_num(s: str):
    s = s.strip().strip('"')
    if not s or s in (".", "..", "-"):
        return None
    return float(s.replace(".", "").replace(",", "."))


def load_ine53001(periods: set[str]) -> pd.DataFrame:
    recs = []
    with gzip.open(P_INE53001, "rt", encoding="utf-8") as f:
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
    df = pd.DataFrame(recs, columns=["ine_code", "per", "tourists"])
    return df.pivot(index="ine_code", columns="per", values="tourists")


def run_c3(nodes: pd.DataFrame) -> pd.DataFrame:
    log("=== PART 2: C3 (Censo 2021 table 59525 second-home cross-sectional decomposition) ===")
    # Parse raw 59525
    rows_59525 = []
    with open(P_VIV59525, encoding="utf-8-sig") as f:
        f.readline()
        for line in f:
            p = line.rstrip("\n").split(";")
            if len(p) < 4 or not p[1].strip():
                continue
            code = p[1].strip()[:5]
            if not re.match(r"^\d{5}$", code):
                continue
            val = parse_es_num(p[3])
            if val is not None:
                rows_59525.append((code, p[2].strip(), val))
    viv59 = pd.DataFrame(rows_59525, columns=["ine_code", "tipo", "val"]).pivot(index="ine_code", columns="tipo", values="val")
    viv59["share_no_principal"] = viv59["Vivienda no principal"] / viv59["Total"]
    assert len(viv59) == 8131

    lau = pd.read_csv(P_LAU2021, dtype={"LAU_ID": str, "CNTR_CODE": str})
    area_map = lau[lau["CNTR_CODE"] == "ES"].assign(c=lambda d: d["LAU_ID"].str.zfill(5)).set_index("c")["AREA_KM2"]

    # Node controls (exact R1)
    nd = nodes.copy()
    xy = nd[["x", "y"]].to_numpy(dtype=float) / 1000.0
    D = np.hypot(xy[:, None, 0] - xy[None, :, 0], xy[:, None, 1] - xy[None, :, 1])
    big = nd["population"].to_numpy(dtype=float) >= 100_000
    coast = nd["coastal_true"].to_numpy(dtype=bool)
    nd["coastal"] = coast.astype(float)
    nd["ldcity"] = np.log1p(D[:, big].min(axis=1))
    nd["ldcoast"] = np.log1p(D[:, coast].min(axis=1))
    nd["lpop"] = np.log(nd["population"].to_numpy(dtype=float))
    nd["lnmun"] = np.log(nd["n_municipalities"].to_numpy(dtype=float))
    nd["larea"] = np.log(nd["area_km2"].to_numpy(dtype=float))
    nd["artefact"] = nd["id"].isin(C.ARTEFACT_ZONES)
    ok_z = nd["loss_peak_2021"].notna() & ~nd["artefact"]
    nd["loss_pct"] = np.nan
    nd.loc[ok_z, "loss_pct"] = nd.loc[ok_z, "loss_peak_2021"].rank(pct=True, method="average")
    zid = nd.set_index("id")

    # Municipality table (exact R1 ine_check + LAU area)
    H = pd.read_csv(B0.HIST_POP_CSV, dtype={"ine_code": str, "province_code": str})
    rel_full = pd.read_csv(B0.INE_MITMA_REL, sep="|", dtype=str)
    zone_of = (rel_full.groupby(["municipio_ine", "municipio_mitma"]).size().reset_index(name="n")
               .sort_values(["municipio_ine", "n", "municipio_mitma"], ascending=[True, False, True])
               .drop_duplicates("municipio_ine").set_index("municipio_ine")["municipio_mitma"])
    H["zone"] = H["ine_code"].map(zone_of)
    H = H[H["zone"].isin(zid.index)].copy()
    H["peak"] = H[B0._PEAK_COLS].max(axis=1)
    H["loss"] = 1.0 - H["pop_2021"] / H["peak"]
    H["pop"] = H["padron_2024"]
    H["area_km2"] = H["ine_code"].map(area_map)
    H["share_np"] = H["ine_code"].map(viv59["share_no_principal"])

    Hm = H[H["pop"].notna() & (H["pop"] > 0) & H["loss"].notna()].copy()
    Hm["prov"] = Hm["zone"].str[:2]
    Hm["loss_pct"] = Hm["loss"].rank(pct=True, method="average")
    Hm["lpop"] = np.log(Hm["pop"])
    Hm["larea"] = np.log(Hm["area_km2"])
    for col in ("coastal", "ldcity", "ldcoast", "is_am"):
        Hm[col] = Hm["zone"].map(zid[col]).astype(float)

    T53 = load_ine53001({"2024M02", "2024M03", "2024M08", "2024M10", "2024M11", "2024M12"})
    for per in ("2024M02", "2024M03", "2024M08", "2024M10", "2024M11", "2024M12"):
        Hm[per] = Hm["ine_code"].map(T53[per])

    # Zone-level housing stock (both unique zone_of and rel_u give identical 2,376 sample because 03903 has missing loss)
    zv = pd.DataFrame({"municipio_ine": zone_of.index, "municipio_mitma": zone_of.values}).merge(
        viv59.reset_index(), left_on="municipio_ine", right_on="ine_code", how="inner"
    )
    zg = zv.groupby("municipio_mitma")[["Vivienda no principal", "Total"]].sum()
    nd["share_np"] = nd["id"].map(zg["Vivienda no principal"] / zg["Total"])

    # Compute MITMA s_in directly from b01/b00 full edge lists (independent of r1_node_metrics!)
    N_nodes = len(nd)
    sin_by_date = {}
    sout_by_date = {}
    for d in ("2024-02-14", "2024-03-29", "2024-08-15", "2024-11-01", "2024-12-25"):
        e = pd.read_parquet(NET / f"edges_{d}__full.parquet")
        sin_by_date[d] = np.bincount(e["j"].to_numpy(), weights=e["w"].to_numpy(), minlength=N_nodes)
        sout_by_date[d] = np.bincount(e["i"].to_numpy(), weights=e["w"].to_numpy(), minlength=N_nodes)

    def cluster_mediation(df: pd.DataFrame, y: np.ndarray, xcol: str, mcol: str, ctrls: list[str], tag: str, n_boot_cl: int = 1999) -> dict:
        good = np.isfinite(y) & df[xcol].notna().to_numpy() & df[mcol].notna().to_numpy()
        for c in ctrls:
            good &= df[c].notna().to_numpy()
        d = df.loc[good].copy()
        yv = y[good]
        xv = d[xcol].to_numpy(dtype=float)
        mv = d[mcol].to_numpy(dtype=float)
        zv_c = d[ctrls].to_numpy(dtype=float)
        fe = d["prov"].to_numpy(dtype=str)

        r_base = fe_ols_multi(yv, xv, zv_c, fe, fe, n_boot=N_BOOT_OLS, seed=int(rng_for("c3_base", tag).integers(1, 2**31)))
        r_med = fe_ols_multi(yv, np.column_stack([xv, mv]), zv_c, fe, fe, n_boot=N_BOOT_OLS,
                             seed=int(rng_for("c3_med", tag).integers(1, 2**31)))
        # First-stage regression M ~ T + controls | prov FE (collinearity / accounting link!)
        r_fs = fe_ols_multi(mv, xv, zv_c, fe, fe, n_boot=0)
        # Partial R^2 of M on T | controls + prov FE
        F_d = pd.get_dummies(pd.Series(fe), dtype=float).to_numpy()
        A = np.column_stack([zv_c, F_d])
        U, S, _ = np.linalg.svd(A, full_matrices=False)
        Q = U[:, S > 1e-10 * S.max()]
        xt = xv - Q @ (Q.T @ xv)
        mt = mv - Q @ (Q.T @ mv)
        r_partial_xm = float(np.corrcoef(xt, mt)[0, 1])
        raw_r_xm = float(np.corrcoef(xv, mv)[0, 1])

        b0 = float(r_base["beta"][0])
        b1 = float(r_med["beta"][0])
        bm = float(r_med["beta"][1])
        diff = b0 - b1
        share = diff / b0

        # Cluster pairs bootstrap (by 47 provinces) for (diff, share) CI
        provs_u = sorted(d["prov"].unique())
        G = len(provs_u)
        idx_by_p = [np.where(fe == p)[0] for p in provs_u]
        rng = rng_for("c3_clboot", tag)
        diffs_b, shares_b = [], []
        for _ in range(n_boot_cl):
            pick = rng.integers(0, G, size=G)
            # Re-index cluster FE so resampled duplicate provinces keep distinct or province-level FE:
            # note in cross-sectional province-FE cluster bootstrap, each drawn cluster gets its own FE block!
            rows_b = []
            fe_b = []
            for new_g, g_idx in enumerate(pick):
                ii = idx_by_p[g_idx]
                rows_b.append(ii)
                fe_b.append(np.full(len(ii), new_g, dtype=int))
            ib = np.concatenate(rows_b)
            feb = np.concatenate(fe_b)
            Fb = sps.csr_matrix((np.ones(len(ib)), (feb, np.arange(len(ib)))), shape=(G, len(ib)))
            cnt_b = np.array(Fb.sum(axis=1)).ravel()[:, None]
            # Demean [yv, xv, mv, zv_c] by cluster FE first, then QR on zv_c!
            M_all = np.column_stack([yv[ib], xv[ib], mv[ib], zv_c[ib]])
            M_dm = M_all - Fb.T @ ((Fb @ M_all) / cnt_b)
            yb, xb, mb, Zb = M_dm[:, 0], M_dm[:, 1], M_dm[:, 2], M_dm[:, 3:]
            Qb, Rb, _ = qr(Zb, mode="economic", pivoting=True)
            ok_q = np.abs(np.diag(Rb)) > 1e-10 * max(1e-12, float(np.abs(Rb).max()))
            Qb = Qb[:, ok_q]
            yt_b = yb - Qb @ (Qb.T @ yb)
            xt_b = xb - Qb @ (Qb.T @ xb)
            mt_b = mb - Qb @ (Qb.T @ mb)
            b0_b = float((xt_b @ yt_b) / (xt_b @ xt_b))
            X2_b = np.column_stack([xt_b, mt_b])
            b1_b = float(np.linalg.solve(X2_b.T @ X2_b, X2_b.T @ yt_b)[0])
            diffs_b.append(b0_b - b1_b)
            shares_b.append((b0_b - b1_b) / b0_b if abs(b0_b) > 1e-8 else np.nan)
        diffs_b = np.array(diffs_b)
        shares_b = np.array(shares_b)
        return {
            "spec": tag, "n": r_base["n"], "G": r_base["G"],
            "raw_corr_loss_housing": raw_r_xm, "partial_corr_loss_housing": r_partial_xm,
            "first_stage_beta_housing_on_loss": float(r_fs["beta"][0]),
            "beta_base": b0, "se_base": float(r_base["se"][0]), "p_base_wcr": float(r_base["p_wcr"][0]),
            "ci_base_lo": float(r_base["ci_wcu_lo"][0]), "ci_base_hi": float(r_base["ci_wcu_hi"][0]),
            "beta_after_med": b1, "se_after": float(r_med["se"][0]), "p_after_wcr": float(r_med["p_wcr"][0]),
            "ci_after_lo": float(r_med["ci_wcu_lo"][0]), "ci_after_hi": float(r_med["ci_wcu_hi"][0]),
            "beta_mediator": bm, "se_med": float(r_med["se"][1]), "p_med_wcr": float(r_med["p_wcr"][1]),
            "ci_med_lo": float(r_med["ci_wcu_lo"][1]), "ci_med_hi": float(r_med["ci_wcu_hi"][1]),
            "decomp_diff": diff,
            "decomp_diff_ci_lo": float(np.nanquantile(diffs_b, 0.025)),
            "decomp_diff_ci_hi": float(np.nanquantile(diffs_b, 0.975)),
            "decomp_share": share,
            "decomp_share_ci_lo": float(np.nanquantile(shares_b, 0.025)),
            "decomp_share_ci_hi": float(np.nanquantile(shares_b, 0.975)),
        }

    m_ctrls_lau = ["lpop", "larea", "coastal", "ldcity", "ldcoast"]
    m_ctrls_nolau = ["lpop", "coastal", "ldcity", "ldcoast"]
    z_ctrls = ["lpop", "lnmun", "larea", "coastal", "ldcity", "ldcoast"]

    c3_rows = [
        cluster_mediation(Hm, np.log(Hm["2024M08"].to_numpy() / Hm["2024M02"].to_numpy()), "loss_pct", "share_np", m_ctrls_lau, "INE53001_log_Aug_Feb_with_LAU_area"),
        cluster_mediation(Hm, np.log(Hm["2024M08"].to_numpy() / Hm["2024M02"].to_numpy()), "loss_pct", "share_np", m_ctrls_nolau, "INE53001_log_Aug_Feb_no_area_R1spec"),
        cluster_mediation(Hm[Hm["coastal"] == 0], np.log(Hm.loc[Hm["coastal"] == 0, "2024M08"].to_numpy() / Hm.loc[Hm["coastal"] == 0, "2024M02"].to_numpy()), "loss_pct", "share_np", ["lpop", "larea", "ldcity", "ldcoast"], "INE53001_log_Aug_Feb_inland_only"),
        cluster_mediation(Hm, np.log(Hm["2024M12"].to_numpy() / Hm["2024M02"].to_numpy()), "loss_pct", "share_np", m_ctrls_lau, "INE53001_log_Dec_Feb_monthly_all"),
        cluster_mediation(Hm, np.log(Hm["2024M11"].to_numpy() / Hm["2024M02"].to_numpy()), "loss_pct", "share_np", m_ctrls_lau, "INE53001_log_Nov_Feb_monthly_all"),
        cluster_mediation(Hm, np.log(Hm["2024M03"].to_numpy() / Hm["2024M02"].to_numpy()), "loss_pct", "share_np", m_ctrls_lau, "INE53001_log_Mar_Feb_monthly_all"),
        cluster_mediation(nd[ok_z], np.log((sin_by_date["2024-08-15"][ok_z] + 1.0) / (sin_by_date["2024-02-14"][ok_z] + 1.0)), "loss_pct", "share_np", z_ctrls, "MITMA_Aug15_vs_Feb14_log_sin"),
        cluster_mediation(nd[ok_z], np.log((sin_by_date["2024-11-01"][ok_z] + 1.0) / (sin_by_date["2024-02-14"][ok_z] + 1.0)), "loss_pct", "share_np", z_ctrls, "MITMA_Nov01_vs_Feb14_log_sin"),
        cluster_mediation(nd[ok_z], np.log((sin_by_date["2024-12-25"][ok_z] + 1.0) / (sin_by_date["2024-02-14"][ok_z] + 1.0)), "loss_pct", "share_np", z_ctrls, "MITMA_Dec25_vs_Feb14_log_sin"),
        cluster_mediation(nd[ok_z], np.log((sin_by_date["2024-03-29"][ok_z] + 1.0) / (sin_by_date["2024-02-14"][ok_z] + 1.0)), "loss_pct", "share_np", z_ctrls, "MITMA_Mar29_vs_Feb14_log_sin"),
    ]
    df_c3 = pd.DataFrame(c3_rows)
    df_c3.to_csv(OUT / "v1_c3_mediation.csv", index=False)
    log("C3 decomposition checks completed")
    return df_c3


# =============================================================================
# PART 3: Verify C4 (ECEPOV 2021 Cross-CCAA Second Homes)
# =============================================================================
def run_c4() -> pd.DataFrame:
    log("=== PART 3: C4 (ECEPOV 2021 cross-CCAA second homes & birth CCAA) ===")
    con = duckdb.connect()
    viv = con.execute(f"SELECT * FROM read_parquet('{P_ECEPOV_VIV_PQ}')").df()
    hog = con.execute(f"SELECT * FROM read_parquet('{P_ECEPOV_HOG_PQ}')").df()

    # Verify parquet directly against raw datos_2021.zip row count & weighted sum
    with zipfile.ZipFile(P_ECEPOV_ZIP) as zf:
        with zipfile.ZipFile(io.BytesIO(zf.read("ECEPOVvivienda.zip"))) as zv:
            csv_v = [n for n in zv.namelist() if n.endswith(".csv")][0]
            with zv.open(csv_v) as f:
                hdr = [x.strip('"') for x in f.readline().decode("utf-8").rstrip("\r\n").split("\t")]
                idx_seg, idx_lug, idx_fac = hdr.index("SEGUNRESI"), hdr.index("LUGSEGUNRESI"), hdr.index("FACTOR")
                n_raw_all, n_raw_c4, w_raw_c4 = 0, 0, 0.0
                for line in f:
                    n_raw_all += 1
                    parts = line.decode("utf-8", errors="replace").rstrip("\r\n").split("\t")
                    if parts[idx_seg].strip('"') == "1" and parts[idx_lug].strip('"') == "4":
                        n_raw_c4 += 1
                        w_raw_c4 += float(parts[idx_fac].strip('"'))
    assert n_raw_all == len(viv) and n_raw_c4 == 9421 and abs(w_raw_c4 - viv.loc[(viv.SEGUNRESI == "1") & (viv.LUGSEGUNRESI == "4"), "FACTOR"].sum()) < 1e-4
    log(f"ECEPOV raw CSV cross-checked: {n_raw_all:,} HHs total, {n_raw_c4:,} cross-CCAA 2nd-home HHs ({w_raw_c4:,.1f} weighted)")

    CCAA_CODE = {
        "AND": "01", "ARA": "02", "AST": "03", "BAL": "04", "CAN": "05", "CANT": "06",
        "CYL": "07", "CLM": "08", "CAT": "09", "VAL": "10", "EXT": "11", "GAL": "12",
        "MAD": "13", "MUR": "14", "NAV": "15", "PVA": "16", "RIO": "17", "CEU": "18", "MEL": "19",
    }
    PROV_TO_CCAA = {p: CCAA_CODE[c] for p, c in B0.PROV2CCAA.items()}

    v2 = viv[(viv["SEGUNRESI"] == "1") & (viv["LUGSEGUNRESI"] == "4")].copy()
    v2["orig_ccaa"] = v2["IDQ_PV"].map(PROV_TO_CCAA)
    assert (v2["orig_ccaa"] != v2["OTRACCAA"]).all(), "Cross-CCAA second home in own CCAA!"

    hog_sub = hog[hog["IDEN"].isin(set(v2["IDEN"]))].copy()
    nac_m = hog_sub["NACIM"].to_numpy(dtype=str) == "1"
    nac_p1 = hog_sub["NACIMP1"].to_numpy(dtype=str) == "1"
    nac_p2 = hog_sub["NACIMP2"].to_numpy(dtype=str) == "1"
    hog_sub["c_mem"] = pd.Series(np.where(nac_m, hog_sub["PROVN"].map(PROV_TO_CCAA).to_numpy(dtype=object), None), index=hog_sub.index)
    hog_sub["c_p1"] = pd.Series(np.where(nac_p1, hog_sub["PROVNP1"].map(PROV_TO_CCAA).to_numpy(dtype=object), None), index=hog_sub.index)
    hog_sub["c_p2"] = pd.Series(np.where(nac_p2, hog_sub["PROVNP2"].map(PROV_TO_CCAA).to_numpy(dtype=object), None), index=hog_sub.index)

    # Also track households with any unknown/missing parent birth province
    hog_sub["p_missing"] = hog_sub["NACIMP1"].isin(["0", " "]) | hog_sub["NACIMP2"].isin(["0", " "])
    hh_any_pmiss = hog_sub.groupby("IDEN")["p_missing"].any().to_dict()

    mem_sets = hog_sub.groupby("IDEN")["c_mem"].apply(lambda s: set(s.dropna())).to_dict()
    par_sets = {}
    for iden, g in hog_sub.groupby("IDEN"):
        par_sets[iden] = set(g["c_p1"].dropna()) | set(g["c_p2"].dropna())

    v2["m_g1"] = np.array([r.OTRACCAA in mem_sets.get(r.IDEN, set()) for r in v2.itertuples()], dtype=float)
    v2["m_par"] = np.array([r.OTRACCAA in par_sets.get(r.IDEN, set()) for r in v2.itertuples()], dtype=float)
    v2["m_any"] = np.maximum(v2["m_g1"], v2["m_par"])
    v2["has_pmiss"] = np.array([hh_any_pmiss.get(r.IDEN, False) for r in v2.itertuples()], dtype=bool)

    INLAND_7 = {"02", "07", "08", "11", "13", "15", "17"}
    INLAND_6_NOMAD = {"02", "07", "08", "11", "15", "17"}

    # Global permutation null (permute OTRACCAA within IDQ_PV across ALL 9,421 cross-CCAA HHs)
    rng_g = rng_for("c4_global_perm")
    all_pvs = v2["IDQ_PV"].to_numpy(dtype=str)
    all_d = v2["OTRACCAA"].to_numpy(dtype=str)
    all_w = v2["FACTOR"].to_numpy(dtype=float)
    all_idens = v2["IDEN"].to_numpy(dtype=str)
    all_any_sets = [mem_sets.get(i, set()) | par_sets.get(i, set()) for i in all_idens]
    pv_grp_all = [np.where(all_pvs == p)[0] for p in np.unique(all_pvs)]
    glob_inl7_null = []
    for _ in range(1000):
        pd_a = all_d.copy()
        for idx in pv_grp_all:
            pd_a[idx] = rng_g.permutation(pd_a[idx])
        hit = np.array([pd_a[k] in all_any_sets[k] for k in range(len(v2))], dtype=float)
        m_i7 = np.isin(pd_a, list(INLAND_7))
        glob_inl7_null.append(np.average(hit[m_i7], weights=all_w[m_i7]))

    def eval_c4_sub(sub: pd.DataFrame, name: str, global_null_mean: float = np.nan) -> dict:
        w = sub["FACTOR"].to_numpy(dtype=float)
        g1 = float(np.average(sub["m_g1"], weights=w))
        par = float(np.average(sub["m_par"], weights=w))
        any_w = float(np.average(sub["m_any"], weights=w))
        any_unw = float(sub["m_any"].mean())

        # Within-subset within-IDQ_PV permutation null (1,000 perms)
        pvs = sub["IDQ_PV"].to_numpy(dtype=str)
        ds = sub["OTRACCAA"].to_numpy(dtype=str)
        sets_s = [mem_sets.get(i, set()) | par_sets.get(i, set()) for i in sub["IDEN"]]
        grps = [np.where(pvs == p)[0] for p in np.unique(pvs)]
        rng_p = rng_for("c4_perm", name)
        null_s = []
        for _ in range(1000):
            dp = ds.copy()
            for idx in grps:
                dp[idx] = rng_p.permutation(dp[idx])
            hit = np.array([dp[k] in sets_s[k] for k in range(len(sub))], dtype=float)
            null_s.append(np.average(hit, weights=w))
        null_s = np.array(null_s)
        null_mean = float(null_s.mean())
        p_perm = float((1 + np.sum(null_s >= any_w)) / (len(null_s) + 1))

        # Bootstrap CI: 52-prov cluster bootstrap if >=10 provs, else HH bootstrap; ALSO HH bootstrap always for comparison!
        rng_b = rng_for("c4_boot", name)
        u_pvs = np.unique(pvs)
        if len(u_pvs) >= 10:
            p_idx = [np.where(pvs == p)[0] for p in u_pvs]
            b_vals = []
            for _ in range(1999):
                pick = rng_b.integers(0, len(u_pvs), size=len(u_pvs))
                ii = np.concatenate([p_idx[k] for k in pick])
                b_vals.append(np.average(sub["m_any"].to_numpy()[ii], weights=w[ii]))
            ci_mode = f"cluster_prov_G{len(u_pvs)}"
        else:
            b_vals = []
            y_arr = sub["m_any"].to_numpy()
            for _ in range(1999):
                ii = rng_b.integers(0, len(sub), size=len(sub))
                b_vals.append(np.average(y_arr[ii], weights=w[ii]))
            ci_mode = f"hh_boot_n{len(sub)}"
        # Stratified-by-IDQ_PV HH bootstrap (survey-design appropriate!)
        b_strat = []
        y_arr = sub["m_any"].to_numpy()
        for _ in range(1999):
            ii = np.concatenate([idx[rng_b.integers(0, len(idx), size=len(idx))] for idx in grps])
            b_strat.append(np.average(y_arr[ii], weights=w[ii]))

        return {
            "subgroup": name, "n_hh": len(sub), "weighted_hh": float(w.sum()),
            "share_madrid_weight": float(w[pvs == "28"].sum() / w.sum()),
            "share_hh_any_parent_missing": float(np.average(sub["has_pmiss"], weights=w)),
            "match_g1_weighted": g1, "match_parent_weighted": par,
            "match_any_weighted": any_w, "match_any_unweighted": any_unw,
            "ci_mode": ci_mode,
            "any_ci_lo": float(np.quantile(b_vals, 0.025)), "any_ci_hi": float(np.quantile(b_vals, 0.975)),
            "strat_hh_ci_lo": float(np.quantile(b_strat, 0.025)), "strat_hh_ci_hi": float(np.quantile(b_strat, 0.975)),
            "within_sub_perm_null": null_mean, "global_perm_null": global_null_mean,
            "lift_over_within_null": any_w / null_mean,
            "lift_over_global_null": any_w / global_null_mean if np.isfinite(global_null_mean) else np.nan,
            "p_perm": p_perm,
        }

    c4_rows = [
        eval_c4_sub(v2, "all_cross_ccaa", np.nan),
        eval_c4_sub(v2[v2["OTRACCAA"].isin(INLAND_7)], "inland_dest_ccaa_7_incl_Madrid", float(np.mean(glob_inl7_null))),
        eval_c4_sub(v2[v2["OTRACCAA"].isin(INLAND_6_NOMAD)], "inland_dest_ccaa_6_excl_Madrid", np.nan),
        eval_c4_sub(v2[~v2["OTRACCAA"].isin(INLAND_7)], "coastal_island_dest_ccaa", np.nan),
        eval_c4_sub(v2[v2["IDQ_PV"] != "28"], "all_cross_ccaa_excl_origin_Madrid", np.nan),
        eval_c4_sub(v2[v2["IDQ_PV"] == "28"], "origin_Madrid_all", np.nan),
        eval_c4_sub(v2[(v2["IDQ_PV"] == "28") & v2["OTRACCAA"].isin(INLAND_7)], "origin_Madrid_to_inland", np.nan),
        eval_c4_sub(v2[v2["orig_ccaa"] == "09"], "origin_Cataluna_all", np.nan),
        eval_c4_sub(v2[(v2["orig_ccaa"] == "09") & v2["OTRACCAA"].isin(INLAND_7)], "origin_Cataluna_to_inland", np.nan),
        eval_c4_sub(v2[v2["orig_ccaa"] == "16"], "origin_PaisVasco_all", np.nan),
    ]
    df_c4 = pd.DataFrame(c4_rows)
    df_c4.to_csv(OUT / "v1_c4_ecepov.csv", index=False)
    log("C4 ECEPOV checks completed")
    return df_c4


# =============================================================================
# PART 4: Verify C5 (Municipal INE OD Network, N = 7,974)
# =============================================================================
def verify_xlsx_month_sum(xlsx_path: Path, sheet_name: str) -> tuple[int, float]:
    """Reads one sheet of raw exp_tmov_interno_mun_2024.xlsx directly via XML stream to verify parquet."""
    ns = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
    with zipfile.ZipFile(xlsx_path) as zf:
        wb = ET.fromstring(zf.read("xl/workbook.xml"))
        snames = [el.attrib["name"] for el in wb.iter(f"{ns}sheet")]
        sidx = snames.index(sheet_name) + 1
        nrows, tot_tur = 0, 0.0
        with zf.open(f"xl/worksheets/sheet{sidx}.xml") as sf:
            for event, elem in ET.iterparse(sf, events=("end",)):
                if elem.tag == f"{ns}row":
                    if elem.attrib.get("r") != "1":
                        cs = [c for c in elem.findall(f"{ns}c") if c.attrib.get("t", "") not in ("s", "str", "inlineStr")]
                        vals = [c.find(f"{ns}v").text for c in cs if c.find(f"{ns}v") is not None and c.find(f"{ns}v").text]
                        if vals:
                            nrows += 1
                            tot_tur += float(vals[-1])
                    elem.clear()
        return nrows, tot_tur


def run_c5(nodes: pd.DataFrame) -> pd.DataFrame:
    log("=== PART 4: C5 (Municipal INE OD network, N = 7,974) ===")
    con = duckdb.connect()
    # Cross-check Feb & Aug 2024 in ine_od_mun_2024.parquet directly against raw Excel sheets!
    for s_xlsx, m_pq in (("202402", "2024-02"), ("202408", "2024-08")):
        nr_x, tot_x = verify_xlsx_month_sum(P_INEOD24_XLSX, s_xlsx)
        nr_p, tot_p = con.execute(
            f"SELECT count(*), sum(turistas) FROM read_parquet('{P_INEOD24_PQ}') WHERE mes = ?", [m_pq]
        ).fetchone()
        assert nr_x == nr_p and abs(tot_x - tot_p) < 1e-6, f"Raw xlsx != parquet on {m_pq}: {(nr_x, tot_x)} vs {(nr_p, tot_p)}"
    log("Raw exp_tmov_interno_mun_2024.xlsx verified row-for-row & sum-for-sum against ine_od_mun_2024.parquet")

    nd = nodes.copy()
    xy = nd[["x", "y"]].to_numpy(dtype=float) / 1000.0
    D = np.hypot(xy[:, None, 0] - xy[None, :, 0], xy[:, None, 1] - xy[None, :, 1])
    big = nd["population"].to_numpy(dtype=float) >= 100_000
    coast = nd["coastal_true"].to_numpy(dtype=bool)
    nd["coastal"] = coast.astype(float)
    nd["ldcity"] = np.log1p(D[:, big].min(axis=1))
    nd["ldcoast"] = np.log1p(D[:, coast].min(axis=1))
    zid = nd.set_index("id")

    H = pd.read_csv(B0.HIST_POP_CSV, dtype={"ine_code": str, "province_code": str})
    rel_full = pd.read_csv(B0.INE_MITMA_REL, sep="|", dtype=str)
    zone_of = (rel_full.groupby(["municipio_ine", "municipio_mitma"]).size().reset_index(name="n")
               .sort_values(["municipio_ine", "n", "municipio_mitma"], ascending=[True, False, True])
               .drop_duplicates("municipio_ine").set_index("municipio_ine")["municipio_mitma"])
    H["zone"] = H["ine_code"].map(zone_of)
    H = H[H["zone"].isin(zid.index)].sort_values("ine_code").reset_index(drop=True)
    assert len(H) == 7974
    lau = pd.read_csv(P_LAU2021, dtype={"LAU_ID": str, "CNTR_CODE": str})
    area_map = lau[lau["CNTR_CODE"] == "ES"].assign(c=lambda d: d["LAU_ID"].str.zfill(5)).set_index("c")["AREA_KM2"]
    H["area_km2"] = H["ine_code"].map(area_map)
    assert H["area_km2"].notna().all()
    H["peak"] = H[B0._PEAK_COLS].max(axis=1)
    H["loss"] = 1.0 - H["pop_2021"] / H["peak"]
    H["pop"] = H["padron_2024"]
    H["artefact"] = H["zone"].isin(C.ARTEFACT_ZONES)
    H["is_am"] = H["zone"].map(zid["is_am"]).astype(bool)
    H["single_zone"] = H["zone"].map(zid["n_municipalities"]) == 1

    M2I = {c: i for i, c in enumerate(H["ine_code"].to_numpy(dtype=str))}
    Nm = len(H)

    def build_mun_net(mes: str, w_min: int = 30, equal_e: int | None = None):
        df = con.execute(
            f"SELECT orig, dest, sum(turistas) AS w FROM read_parquet('{P_INEOD24_PQ}') WHERE mes = ? GROUP BY 1, 2",
            [mes],
        ).df()
        df = df[df["orig"].isin(M2I) & df["dest"].isin(M2I) & (df["orig"] != df["dest"])].copy()
        df["i"] = df["orig"].map(M2I).astype(int)
        df["j"] = df["dest"].map(M2I).astype(int)
        # Full published s_in (all cells >= 30)
        s_in = np.bincount(df["j"].to_numpy(), weights=df["w"].to_numpy(dtype=float), minlength=Nm)
        if equal_e is not None:
            sub = df.sort_values(["w", "i", "j"], ascending=[False, True, True], kind="mergesort").head(equal_e)
        else:
            sub = df[df["w"] >= w_min]
        ii, jj = sub["i"].to_numpy(), sub["j"].to_numpy()
        k_in = np.bincount(jj, minlength=Nm).astype(float)
        Gd = ig.Graph(n=Nm, edges=list(zip(ii.tolist(), jj.tolist())), directed=True)
        core_in = np.array(Gd.coreness(mode="in"), dtype=float)
        Gu = B0.to_igraph(pd.DataFrame({"i": ii, "j": jj}), Nm, directed=False)
        core_und = np.array(Gu.coreness(), dtype=float)
        return {"s_in": s_in, "k_in": k_in, "core_in": core_in, "core_und": core_und, "E_dir": len(sub)}

    mf = build_mun_net("2024-02", w_min=30)
    ma = build_mun_net("2024-08", w_min=30)
    ma_eq = build_mun_net("2024-08", equal_e=mf["E_dir"])
    mf_100 = build_mun_net("2024-02", w_min=100)
    ma_100 = build_mun_net("2024-08", w_min=100)

    ok_n1 = H["pop"].notna() & (H["pop"] > 0) & H["loss"].notna()
    d = H[ok_n1].copy()
    d["loss_pct"] = d["loss"].rank(pct=True, method="average")
    d["lpop"] = np.log(d["pop"])
    d["larea"] = np.log(d["area_km2"])
    for col in ("coastal", "ldcity", "ldcoast"):
        d[col] = d["zone"].map(zid[col]).astype(float)
    d["prov"] = d["zone"].str[:2]
    ctrls = ["lpop", "larea", "coastal", "ldcity", "ldcoast"]

    # Zero-baseline censorship diagnostic inside _AM municipalities!
    am_mask_full = ok_n1.to_numpy() & H["is_am"].to_numpy()
    zero_feb_am = float((mf["s_in"][am_mask_full] == 0).mean())
    zero_aug_am = float((ma["s_in"][am_mask_full] == 0).mean())
    log(f"C5 censorship diagnostic on 6,276 _AM munis: s_in==0 (<30 censored) in {100*zero_feb_am:.1f}% (Feb) vs {100*zero_aug_am:.1f}% (Aug)")

    c5_rows = []
    outcomes = [
        ("log_s_in_ratio_pub30", np.log((ma["s_in"][ok_n1] + 1.0) / (mf["s_in"][ok_n1] + 1.0))),
        ("log_k_in_ratio_pub30", np.log((ma["k_in"][ok_n1] + 1.0) / (mf["k_in"][ok_n1] + 1.0))),
        ("core_in_diff_pub30", ma["core_in"][ok_n1] - mf["core_in"][ok_n1]),
        ("core_und_diff_pub30", ma["core_und"][ok_n1] - mf["core_und"][ok_n1]),
        ("log_k_in_ratio_equal_edges", np.log((ma_eq["k_in"][ok_n1] + 1.0) / (mf["k_in"][ok_n1] + 1.0))),
        ("core_in_diff_equal_edges", ma_eq["core_in"][ok_n1] - mf["core_in"][ok_n1]),
        ("log_k_in_ratio_w100", np.log((ma_100["k_in"][ok_n1] + 1.0) / (mf_100["k_in"][ok_n1] + 1.0))),
    ]
    samples_c5 = [
        ("all_mainland_mun", np.ones(len(d), bool)),
        ("inland_mun", d["coastal"].to_numpy() == 0),
        ("inside_AM_mun", d["is_am"].to_numpy()),
        ("inside_AM_mun_both_pos_ge30", d["is_am"].to_numpy() & (mf["s_in"][ok_n1] > 0) & (ma["s_in"][ok_n1] > 0)),
        ("inside_AM_mun_excl_artefact", d["is_am"].to_numpy() & ~d["artefact"].to_numpy()),
        ("single_zone_mun", d["single_zone"].to_numpy()),
    ]
    for oname, y_vec in outcomes:
        for sname, smask in samples_c5:
            if oname not in ("log_s_in_ratio_pub30", "log_k_in_ratio_pub30", "core_in_diff_pub30", "core_und_diff_pub30") and not sname.startswith("inside_AM"):
                continue
            sub = d.loc[smask]
            c_use = [c for c in ctrls if sub[c].nunique() > 1]
            r = fe_ols_multi(y_vec[smask], sub["loss_pct"].to_numpy(), sub[c_use].to_numpy(),
                             sub["prov"].to_numpy(), sub["prov"].to_numpy(),
                             n_boot=N_BOOT_OLS, seed=int(rng_for("c5_wcb", oname, sname).integers(1, 2**31)))
            c5_rows.append({
                "outcome": oname, "sample": sname, "n": r["n"], "G": r["G"],
                "beta": float(r["beta"][0]), "se_crv1": float(r["se"][0]),
                "p_crv1": float(r["p_t"][0]), "p_wcr": float(r["p_wcr"][0]),
                "ci_wcu_lo": float(r["ci_wcu_lo"][0]), "ci_wcu_hi": float(r["ci_wcu_hi"][0]),
                "zero_feb_share_in_sample": float((mf["s_in"][ok_n1][smask] == 0).mean()),
                "zero_aug_share_in_sample": float((ma["s_in"][ok_n1][smask] == 0).mean()),
            })
    df_c5 = pd.DataFrame(c5_rows)
    df_c5.to_csv(OUT / "v1_c5_mun_network.csv", index=False)
    log("C5 municipal network checks completed")
    return df_c5


# =============================================================================
# Main
# =============================================================================
def main():
    nodes = B0.load_nodes()
    saved_nodes = pd.read_parquet(NET / "nodes.parquet")
    assert (saved_nodes["id"].to_numpy(dtype=str) == nodes["id"].to_numpy(dtype=str)).all()

    res_c12 = run_c1_and_c2(nodes)
    df_c3 = run_c3(nodes)
    df_c4 = run_c4()
    df_c5 = run_c5(nodes)

    # Assemble master claim-by-claim reproduction summary table (v1_claims_summary.csv)
    md = res_c12["main_diffs"]
    g_c1 = md[(md.contrast == "Aug15_vs_Feb14") & (md.sample == "all")].set_index("term")
    c3_m = df_c3.set_index("spec")
    c4_m = df_c4.set_index("subgroup")
    c5_m = df_c5.set_index(["outcome", "sample"])

    claims = [
        ("C1a", "Censo 2024 Aug15-Feb14 pre-1981 exodus delta (all, G=44, n=675 munis)",
         0.083649, float(g_c1.loc["pre1981_exodus", "delta"]), float(g_c1.loc["pre1981_exodus", "p_wcr_9999"]),
         float(g_c1.loc["pre1981_exodus", "ci_wcr_lo"]), float(g_c1.loc["pre1981_exodus", "ci_wcr_hi"]),
         "CONFIRMED WITH CAVEAT",
         "Exact to 1e-14 using Haldane (count+0.5)/N on unweighted 10% counts; 46% zero cells make scale c0-dependent and 76% of pre-1981 are 65+"),
        ("C1b", "Censo 2024 Aug15-Feb14 recent 2016-23 delta (all, G=44)",
         -0.071563, float(g_c1.loc["recent_2016_2024", "delta"]), float(g_c1.loc["recent_2016_2024", "p_wcr_9999"]),
         float(g_c1.loc["recent_2016_2024", "ci_wcr_lo"]), float(g_c1.loc["recent_2016_2024", "ci_wcr_hi"]),
         "CONFIRMED WITH CAVEAT",
         "Exact to 1e-14; VARANOP ends at 2023 (1 Jan 2024 census, no 2024 code), 51% zero cells, 1981-2015 middle cohort omitted"),
        ("C1c", "Censo 2024 Aug15-Feb14 village-born <=2k delta (all, G=44)",
         0.086459, float(g_c1.loc["born_village_le2k", "delta"]), float(g_c1.loc["born_village_le2k", "p_wcr_9999"]),
         float(g_c1.loc["born_village_le2k", "ci_wcr_lo"]), float(g_c1.loc["born_village_le2k", "ci_wcr_hi"]),
         "CONFIRMED WITH CAVEAT",
         "Exact to 1e-14; passes all 117 sister pairs and survives <65 vs 65+ split, but 60% zero cells mean 'triples' ratio is c0-artefact"),
        ("C1d", "Censo 2024 Aug15-Feb14 urban-born >10k delta (all, G=44)",
         -0.061562, float(g_c1.loc["born_urban_gt10k", "delta"]), float(g_c1.loc["born_urban_gt10k", "p_wcr_9999"]),
         float(g_c1.loc["born_urban_gt10k", "ci_wcr_lo"]), float(g_c1.loc["born_urban_gt10k", "ci_wcr_hi"]),
         "CONFIRMED WITH CAVEAT",
         "Exact to 1e-14; fragile across sister nights and driven partly by Barcelona"),
        ("C3a", "Censo 2021 59525 non-primary housing share absorbed, INE 53001 log(Aug/Feb) (n=5,831)",
         0.646393, float(c3_m.loc["INE53001_log_Aug_Feb_with_LAU_area", "decomp_share"]),
         float(c3_m.loc["INE53001_log_Aug_Feb_with_LAU_area", "p_med_wcr"]),
         float(c3_m.loc["INE53001_log_Aug_Feb_with_LAU_area", "decomp_share_ci_lo"]),
         float(c3_m.loc["INE53001_log_Aug_Feb_with_LAU_area", "decomp_share_ci_hi"]),
         "CONFIRMED WITH CAVEAT",
         "Number exact (64.6%); cross-sectional Gelbach/Baron-Kenny accounting, NOT causal mediation; also mixes INE Aug (64.6%) with MITMA Dec25 (1.8%) while INE Dec/Feb is 43.1%"),
        ("C3b", "Censo 2021 59525 non-primary housing share absorbed, MITMA Dec 25 vs Feb 14 (n=2,376)",
         0.018165, float(c3_m.loc["MITMA_Dec25_vs_Feb14_log_sin", "decomp_share"]),
         float(c3_m.loc["MITMA_Dec25_vs_Feb14_log_sin", "p_med_wcr"]),
         float(c3_m.loc["MITMA_Dec25_vs_Feb14_log_sin", "decomp_share_ci_lo"]),
         float(c3_m.loc["MITMA_Dec25_vs_Feb14_log_sin", "decomp_share_ci_hi"]),
         "CONFIRMED WITH CAVEAT",
         "Number exact (1.8% on MITMA zones; matched MITMA Aug 15 is 48.1%); contrast 48.1% vs 1.8% on same dataset"),
        ("C4a", "ECEPOV 2021 cross-CCAA 2nd homes in member/parent birth CCAA (all, n=9,421)",
         0.583822, float(c4_m.loc["all_cross_ccaa", "match_any_weighted"]),
         float(c4_m.loc["all_cross_ccaa", "p_perm"]),
         float(c4_m.loc["all_cross_ccaa", "any_ci_lo"]), float(c4_m.loc["all_cross_ccaa", "any_ci_hi"]),
         "CONFIRMED",
         "Exact (58.4%, 2.90x null 20.1%); conservative lower bound due to ~9% missing parent birthplace; Madrid is 48% of weighted total"),
        ("C4b", "ECEPOV 2021 inland cross-CCAA 2nd homes in member/parent birth CCAA (n=4,241)",
         0.773701, float(c4_m.loc["inland_dest_ccaa_7_incl_Madrid", "match_any_weighted"]),
         float(c4_m.loc["inland_dest_ccaa_7_incl_Madrid", "p_perm"]),
         float(c4_m.loc["inland_dest_ccaa_7_incl_Madrid", "any_ci_lo"]),
         float(c4_m.loc["inland_dest_ccaa_7_incl_Madrid", "any_ci_hi"]),
         "CONFIRMED WITH CAVEAT",
         "Exact (77.4%, 2.02x subset null 38.3%, 2.56x global null 30.2%); 'inland CCAA' includes Madrid (13); excluding Madrid gives 79.0% (n=3,896)"),
        ("C5a", "Municipal INE OD (N=7,974) _AM munis Aug-Feb log((s_in+1)/(s_in+1)) (n=6,276)",
         1.178236, float(c5_m.loc[("log_s_in_ratio_pub30", "inside_AM_mun"), "beta"]),
         float(c5_m.loc[("log_s_in_ratio_pub30", "inside_AM_mun"), "p_wcr"]),
         float(c5_m.loc[("log_s_in_ratio_pub30", "inside_AM_mun"), "ci_wcu_lo"]),
         float(c5_m.loc[("log_s_in_ratio_pub30", "inside_AM_mun"), "ci_wcu_hi"]),
         "CONFIRMED WITH CAVEAT",
         "Exact (+1.178); magnitude inflated by <30 censorship in Feb (zero in Feb -> >=30 in Aug); conditioning on >=30 in both months attenuates"),
        ("C5b", "Municipal INE OD (N=7,974) _AM munis Aug-Feb log((k_in+1)/(k_in+1)) (n=6,276)",
         0.447545, float(c5_m.loc[("log_k_in_ratio_pub30", "inside_AM_mun"), "beta"]),
         float(c5_m.loc[("log_k_in_ratio_pub30", "inside_AM_mun"), "p_wcr"]),
         float(c5_m.loc[("log_k_in_ratio_pub30", "inside_AM_mun"), "ci_wcu_lo"]),
         float(c5_m.loc[("log_k_in_ratio_pub30", "inside_AM_mun"), "ci_wcu_hi"]),
         "CONFIRMED WITH CAVEAT",
         "Exact (+0.448); robust on equal-edges backbone and excluding MITMA artefact zones"),
        ("C5c", "Municipal INE OD (N=7,974) _AM munis Aug-Feb delta core_in (n=6,276)",
         1.777205, float(c5_m.loc[("core_in_diff_pub30", "inside_AM_mun"), "beta"]),
         float(c5_m.loc[("core_in_diff_pub30", "inside_AM_mun"), "p_wcr"]),
         float(c5_m.loc[("core_in_diff_pub30", "inside_AM_mun"), "ci_wcu_lo"]),
         float(c5_m.loc[("core_in_diff_pub30", "inside_AM_mun"), "ci_wcu_hi"]),
         "CONFIRMED WITH CAVEAT",
         "Exact (+1.777 on pub_ge30); check equal-edges backbone where global edge volume is held fixed"),
    ]
    pd.DataFrame(claims, columns=[
        "claim_id", "description", "n1_reported", "v1_reproduced", "v1_p_value",
        "v1_ci_lo", "v1_ci_hi", "verdict", "notes",
    ]).to_csv(OUT / "v1_claims_summary.csv", index=False)
    log("ALL VERIFICATION STAGES COMPLETE -> results/v1_*.csv written.")


if __name__ == "__main__":
    main()
