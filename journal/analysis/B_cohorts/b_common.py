"""Common data loaders and multi-block HDFE PPML + WCR score bootstrap for Analysis B.

Imports frozen definitions and solvers from the coursework pipeline (read-only) and
reads heavy data by absolute path from the main checkout without touching the held-out
confirmation files in _local/external_journal/mitma_pernoct/.
"""
from __future__ import annotations

import sys
import time
import zlib
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
from scipy import stats

MAIN_ROOT = Path("/Users/alejandrotreny/Proyectos/spain-holiday-roots")
JOURNAL_DATA_BD = Path("/Users/alejandrotreny/Proyectos/spain-holiday-roots-journal/_local/journal_data/BD")

# Reuse coursework modules by import only (never edit pipeline/)
sys.path.insert(0, str(MAIN_ROOT / "pipeline"))
import config as C  # noqa: E402
import v1_verify_n1 as V1  # noqa: E402

PERNOCT_2024 = MAIN_ROOT / "data" / "derived" / "pernoct_mun_2024_6m.parquet"
PERNOCT_2022_23 = MAIN_ROOT / "data" / "derived" / "pernoct_mun_2022_2023.parquet"
NODES_PQ = MAIN_ROOT / "results" / "b01_networks" / "nodes.parquet"
NET_DIR = MAIN_ROOT / "results" / "b01_networks"
EXT_DIR = MAIN_ROOT / "_local" / "external"
CENSO24_PQ = EXT_DIR / "ine_censo_microdatos_2024" / "md_CensoAnual_muestra_2024.parquet"
ADRH_CSV = EXT_DIR / "ine_adrh" / "adrh_municipios.csv"
BIRTH68545_CSV = EXT_DIR / "ine_birth_residence" / "mun_residence_x_prov_birth_censo_2021_2025.csv"

WEBB = V1.WEBB
N_BOOT_WCR = 9999
ALPHA = 0.05

# Binding pre-registration §1 arrival-code sets (audit correction: drop 7777 / special / missing)
DATED_ARRIVAL_CODES = (
    "0001", "0002", "0003", "0004", "0005",
    "2016", "2017", "2018", "2019", "2020", "2021", "2022", "2023",
)
POST1981_DATED_CODES = (
    "0002", "0003", "0004", "0005",
    "2016", "2017", "2018", "2019", "2020", "2021", "2022", "2023",
)
MID_ARRIVAL_CODES = ("0002", "0003", "0004", "0005")
RECENT_ARRIVAL_CODES = ("2016", "2017", "2018", "2019", "2020", "2021", "2022", "2023")
SPEC_MUNI_CODES = ("991", "992", "993", "666", "966", "777", "888", "999")
DANA_DATES = ("2024-10-29", "2024-10-30", "2024-10-31")
VAL_HOL_OCT_DATES = ("2024-10-09", "2024-10-10")  # evenings 8-9 Oct (Valencia Day 9 Oct)


def rng_for(*parts) -> np.random.Generator:
    ints = [C.SEED] + [
        zlib.crc32(str(p).encode()) if not isinstance(p, (int, np.integer)) else int(p)
        for p in parts
    ]
    return np.random.default_rng(np.random.SeedSequence(ints))


def _sql_in(vals: tuple[str, ...]) -> str:
    return "(" + ",".join(f"'{v}'" for v in vals) + ")"


# =============================================================================
# Generalized Multi-Block Stacked HDFE PPML + Kline-Santos WCR Score Bootstrap
# =============================================================================
class FastStackedMultiBlock:
    """Exact block-diagonal stacked HDFE PPML over B blocks (nights and/or strata).

    Each block b has its own SingleNightDesign (block-specific city FE, orig_prov x dest_prov FE,
    block-specific income x dest_prov slopes, and block-specific distance covariates).
    V_full (sum_b n_b, m) holds the shared/interacted roots regressors across all blocks.
    """

    def __init__(self, designs: list[V1.SingleNightDesign], V_full: np.ndarray, mu0_warm: np.ndarray | None = None):
        self.designs = designs
        self.B = len(designs)
        self.ns = [d.n for d in designs]
        self.slices = []
        pos = 0
        for nb in self.ns:
            self.slices.append((pos, pos + nb))
            pos += nb
        self.N = pos
        assert V_full.shape[0] == self.N
        self.V = np.asarray(V_full, dtype=float)
        self.y = np.concatenate([d.y for d in designs])
        self.cl = np.concatenate([d.cl for d in designs])
        self.mu0_warm = None if mu0_warm is None else np.asarray(mu0_warm, dtype=float).copy()

    def _demean_blocks(self, mat: np.ndarray, mu: np.ndarray, tol: float = 1e-11):
        m = mat.shape[1]
        V_parts = []
        C_parts = []
        for d, (s, e) in zip(self.designs, self.slices):
            if d.cov.shape[1]:
                Dm = d.hd.demean(np.column_stack([mat[s:e], d.cov]), mu[s:e], tol=tol)
                V_parts.append(Dm[:, :m])
                C_parts.append(Dm[:, m:])
            else:
                Dm = d.hd.demean(mat[s:e], mu[s:e], tol=tol)
                V_parts.append(Dm if Dm.ndim == 2 else Dm[:, None])
                C_parts.append(np.empty((e - s, 0)))
        return np.vstack(V_parts), C_parts

    def _solve_wls_blocks(self, zt: np.ndarray, Vt: np.ndarray, C_parts: list[np.ndarray], mu: np.ndarray):
        z_parts, V_parts = [], []
        for C_b, (s, e) in zip(C_parts, self.slices):
            zb, Vb, mub = zt[s:e], Vt[s:e], mu[s:e]
            if C_b.shape[1]:
                CW = C_b.T * mub
                M = CW @ C_b
                bz = np.linalg.solve(M, CW @ zb)
                bV = np.linalg.solve(M, CW @ Vb)
                z_parts.append(zb - C_b @ bz)
                V_parts.append(Vb - C_b @ bV)
            else:
                z_parts.append(zb)
                V_parts.append(Vb)
        zp = np.concatenate(z_parts)
        Vp = np.vstack(V_parts)
        if Vp.shape[1]:
            VW = Vp.T * mu
            bv = np.linalg.solve(VW @ Vp, VW @ zp)
            e_res = zp - Vp @ bv
        else:
            bv = np.zeros(0)
            e_res = zp
        return bv, e_res

    def fit_unrestricted(self, n_clusters_total: int, maxit: int = 200) -> dict:
        mu = self.mu0_warm.copy() if self.mu0_warm is not None else (self.y + self.y.mean()) / 2.0
        eta = np.log(mu)
        dev_old, bv_old = np.inf, np.zeros(self.V.shape[1])
        for it in range(maxit):
            z = eta + (self.y - mu) / mu
            tol_d = 1e-11 if it >= 2 else 1e-8
            VD, Cs = self._demean_blocks(np.column_stack([z, self.V]), mu, tol=tol_d)
            bv, e_res = self._solve_wls_blocks(VD[:, 0], VD[:, 1:], Cs, mu)
            eta = z - e_res
            mu = np.exp(eta)
            dev = V1.poisson_dev(self.y, mu)
            if abs(dev - dev_old) <= 1e-11 * (abs(dev) + 0.1) and np.max(np.abs(bv - bv_old)) < 1e-9:
                break
            bv_old, dev_old = bv, dev
        else:
            raise RuntimeError("FastStackedMultiBlock unrestricted PPML did not converge")

        # Residualise V against block FEs + block covariates at converged mu for exact CRV1 IF
        VD, Cs = self._demean_blocks(self.V, mu, tol=1e-12)
        Vp_parts = []
        for C_b, (s, e) in zip(Cs, self.slices):
            Vb, mub = VD[s:e], mu[s:e]
            if C_b.shape[1]:
                CW = C_b.T * mub
                Vp_parts.append(Vb - C_b @ np.linalg.solve(CW @ C_b, CW @ Vb))
            else:
                Vp_parts.append(Vb)
        Vp = np.vstack(Vp_parts)
        H = (Vp.T * mu) @ Vp
        Hinv = np.linalg.inv(H)
        sc = Vp * (self.y - mu)[:, None]
        assert np.abs(sc.sum(0)).max() < 1e-4 * (1.0 + np.abs(self.y).sum()), f"MultiBlock score != 0: {np.abs(sc.sum(0)).max()}"
        Sg = np.vstack([np.bincount(self.cl, sc[:, q], minlength=n_clusters_total) for q in range(self.V.shape[1])]).T
        IF = Sg @ Hinv
        Geff = int((np.abs(IF[:, 0]) > 0).sum())
        se = np.sqrt(Geff / (Geff - 1) * np.sum(IF ** 2, axis=0))
        return {"beta": bv.copy(), "se": se, "IF": IF, "G": Geff, "mu": mu, "dev": dev, "n": self.N}

    def p_param(
        self, target_idx: int, val0: float, VB: np.ndarray, mu0=None, maxit: int = 200,
        return_stats: bool = False,
    ):
        target = self.V[:, target_idx]
        V_other = np.delete(self.V, target_idx, axis=1)
        off = val0 * target
        if mu0 is not None:
            mu = mu0.copy()
        elif self.mu0_warm is not None:
            mu = self.mu0_warm.copy()
        else:
            mu = (self.y + self.y.mean()) / 2.0
        eta = np.log(mu)
        dev_old, bv_old = np.inf, np.zeros(V_other.shape[1])
        db = 1.0
        for it in range(maxit):
            z = eta - off + (self.y - mu) / mu
            tol_d = 1e-6 if it == 0 else max(1e-9, min(1e-6, 1e-3 * db))
            VD, Cs = self._demean_blocks(np.column_stack([z, V_other]), mu, tol=tol_d)
            bv, e_res = self._solve_wls_blocks(VD[:, 0], VD[:, 1:], Cs, mu)
            eta = z - e_res + off
            mu = np.exp(eta)
            dev = V1.poisson_dev(self.y, mu)
            db = float(np.max(np.abs(bv - bv_old))) if len(bv) else 0.0
            if abs(dev - dev_old) <= 1e-9 * (abs(dev) + 0.1) and db < 1e-6:
                break
            bv_old, dev_old = bv, dev
        else:
            raise RuntimeError("FastStackedMultiBlock restricted PPML did not converge")

        VD, Cs = self._demean_blocks(np.column_stack([target, V_other]), mu, tol=1e-9)
        _, tt_ = self._solve_wls_blocks(VD[:, 0], VD[:, 1:], Cs, mu)
        s = tt_ * (self.y - mu)
        G_tot = VB.shape[1]
        Sg = np.bincount(self.cl, s, minlength=G_tot)
        denom0 = float(np.sqrt(np.sum(Sg ** 2)))
        if denom0 <= 1e-30:
            return (1.0, mu, 0.0, 1.96) if return_stats else (1.0, mu)
        t0 = float(Sg.sum() / denom0)
        tb = (VB @ Sg) / np.sqrt((VB ** 2) @ (Sg ** 2))
        abs_tb = np.abs(tb)
        p_val = float((1 + np.sum(abs_tb >= abs(t0))) / (VB.shape[0] + 1))
        if return_stats:
            crit = float(np.quantile(abs_tb, 1.0 - ALPHA))
            return p_val, mu, t0, crit
        return p_val, mu


def build_two_night_stack(nh: V1.SingleNightDesign, no: V1.SingleNightDesign,
                          mu_h=None, mu_o=None) -> FastStackedMultiBlock:
    """Two-night stack: columns are [Xr_hol (k), Xr_common (k)]."""
    Xr_common = np.vstack([no.Xr, nh.Xr])
    Xr_hol = np.vstack([np.zeros_like(no.Xr), nh.Xr])
    V = np.column_stack([Xr_hol, Xr_common])
    mu0 = np.r_[mu_o, mu_h] if (mu_o is not None and mu_h is not None) else None
    return FastStackedMultiBlock([no, nh], V, mu0_warm=mu0)


def build_two_night_diff_stack(nh: V1.SingleNightDesign, no: V1.SingleNightDesign,
                               idx_pos: int, idx_neg: int,
                               mu_h=None, mu_o=None) -> FastStackedMultiBlock:
    """Two-night stack reparameterized so column 0 is delta[idx_pos] - delta[idx_neg].

    Uses identity: d_pos * x_pos + d_neg * x_neg = (d_pos - d_neg) * x_pos + d_neg * (x_pos + x_neg).
    """
    Xr_common = np.vstack([no.Xr, nh.Xr])
    Xr_hol = np.vstack([np.zeros_like(no.Xr), nh.Xr])
    k = no.n_terms
    v_diff = Xr_hol[:, idx_pos]
    v_sum = Xr_hol[:, idx_pos] + Xr_hol[:, idx_neg]
    other_cols = [Xr_hol[:, j] for j in range(k) if j not in (idx_pos, idx_neg)]
    V = np.column_stack([v_diff, v_sum] + other_cols + [Xr_common])
    mu0 = np.r_[mu_o, mu_h] if (mu_o is not None and mu_h is not None) else None
    return FastStackedMultiBlock([no, nh], V, mu0_warm=mu0)


def build_four_block_inland_coast_stack(
    no_c: V1.SingleNightDesign, nh_c: V1.SingleNightDesign,
    no_i: V1.SingleNightDesign, nh_i: V1.SingleNightDesign,
    mu_oc=None, mu_hc=None, mu_oi=None, mu_hi=None,
) -> FastStackedMultiBlock:
    """4-block stack [ord_coast, hol_coast, ord_inland, hol_inland].

    Columns are [gamma (k): roots x hol x inland, delta_coast (k): roots x hol,
                 alpha_inland (k): roots x inland, beta_base (k): roots].
    Unrestricted MLE satisfies gamma == (beta_hi - beta_oi) - (beta_hc - beta_oc).
    """
    z_oc = np.zeros_like(no_c.Xr)
    z_hc = np.zeros_like(nh_c.Xr)
    z_oi = np.zeros_like(no_i.Xr)
    V_base = np.vstack([no_c.Xr, nh_c.Xr, no_i.Xr, nh_i.Xr])
    V_hol = np.vstack([z_oc, nh_c.Xr, z_oi, nh_i.Xr])
    V_inl = np.vstack([z_oc, z_hc, no_i.Xr, nh_i.Xr])
    V_gam = np.vstack([z_oc, z_hc, z_oi, nh_i.Xr])
    V = np.column_stack([V_gam, V_hol, V_inl, V_base])
    if all(m is not None for m in (mu_oc, mu_hc, mu_oi, mu_hi)):
        mu0 = np.r_[mu_oc, mu_hc, mu_oi, mu_hi]
    else:
        mu0 = None
    return FastStackedMultiBlock([no_c, nh_c, no_i, nh_i], V, mu0_warm=mu0)


def build_four_block_diff_stack(
    no_c: V1.SingleNightDesign, nh_c: V1.SingleNightDesign,
    no_i: V1.SingleNightDesign, nh_i: V1.SingleNightDesign,
    idx_pos: int, idx_neg: int,
    mu_oc=None, mu_hc=None, mu_oi=None, mu_hi=None,
) -> FastStackedMultiBlock:
    """4-block stack where column 0 is gamma[idx_pos] - gamma[idx_neg]."""
    st = build_four_block_inland_coast_stack(no_c, nh_c, no_i, nh_i, mu_oc, mu_hc, mu_oi, mu_hi)
    k = no_c.n_terms
    V_gam = st.V[:, :k]
    V_rest = st.V[:, k:]
    v_diff = V_gam[:, idx_pos]
    v_sum = V_gam[:, idx_pos] + V_gam[:, idx_neg]
    other = [V_gam[:, j] for j in range(k) if j not in (idx_pos, idx_neg)]
    V_new = np.column_stack([v_diff, v_sum] + other + [V_rest])
    return FastStackedMultiBlock(st.designs, V_new, mu0_warm=st.mu0_warm)


def build_multinight_stack(
    designs_ord: list[V1.SingleNightDesign],
    designs_hol: list[V1.SingleNightDesign],
    mu0_list: list[np.ndarray] | None = None,
) -> FastStackedMultiBlock:
    """PREREGISTRATION §3 multi-night stack: all O nights + all H nights with night-interacted FEs/covariates.

    Columns are [Xr_hol (k), Xr_common (k)].
    """
    designs = list(designs_ord) + list(designs_hol)
    Xr_common = np.vstack([d.Xr for d in designs])
    Xr_hol = np.vstack([np.zeros_like(d.Xr) for d in designs_ord] + [d.Xr for d in designs_hol])
    V = np.column_stack([Xr_hol, Xr_common])
    mu0 = np.concatenate(mu0_list) if mu0_list is not None else None
    return FastStackedMultiBlock(designs, V, mu0_warm=mu0)


def invert_ci_safe(
    pfun, est: float, se: float, alpha: float = ALPHA, tol_frac: float = 5e-3
) -> tuple[float, float]:
    """Invert two-sided WCR score test for 95% CI using bracketed secant on g(b) = |t0(b)| - crit(b).

    pfun(b, mu0) may return (p, mu, t0, crit) for secant acceleration, or (p, mu) for bisection fallback.
    Guarantees an exact sign-bracketing interval [inner, outer] with p(inner) > alpha >= p(outer)
    and |outer - inner| <= tol_frac * se, converging in ~4-5 evals per endpoint instead of ~11.
    """
    est = float(est)
    se = max(float(se), 1e-4)
    tol = tol_frac * se
    out = []

    for sgn in (-1.0, 1.0):
        cache = {"mu": None}

        def ev(b: float) -> tuple[float, float | None, float | None]:
            res = pfun(b, cache["mu"])
            if len(res) == 4:
                p, mu, t0, crit = res
                cache["mu"] = mu
                return float(p), abs(float(t0)), float(crit)
            p, mu = res
            cache["mu"] = mu
            return float(p), None, None

        # Evaluate probe 1 at 1.96 * se
        b1 = est + sgn * 1.96 * se
        p1, at1, c1 = ev(b1)
        g1 = (at1 - c1) if (at1 is not None and c1 is not None) else None

        # Use secant jump from (est, -c1) through (b1, at1 - c1) if score stats available
        pts = [(est, 1.0, -c1 if c1 is not None else None), (b1, p1, g1)]
        if at1 is not None and c1 is not None and at1 > 1e-8:
            b2 = est + (b1 - est) * (c1 / at1)
            b2 = est + sgn * float(np.clip(abs(b2 - est), 0.35 * se, 6.0 * se))
            if abs(b2 - b1) > 0.5 * tol:
                p2, at2, c2 = ev(b2)
                g2 = (at2 - c2) if (at2 is not None and c2 is not None) else None
                pts.append((b2, p2, g2))

        # Sort points by distance from est to extract tightest valid bracket
        pts.sort(key=lambda x: abs(x[0] - est))
        inner, g_in = est, pts[0][2]
        outer, g_out = None, None
        for bp, pp, gp in pts[1:]:
            if pp > alpha:
                inner, g_in = bp, gp
            else:
                outer, g_out = bp, gp
                break

        # Expand outward if no rejected point yet
        k = 0
        while outer is None:
            k += 1
            if k > 30:
                raise RuntimeError(f"CI bracket failed at est={est}, se={se}")
            step = max(abs(inner - est), 1.96 * se)
            cand = est + sgn * step * 1.5
            pc, atc, cc = ev(cand)
            gc = (atc - cc) if (atc is not None and cc is not None) else None
            if pc > alpha:
                inner, g_in = cand, gc
            else:
                outer, g_out = cand, gc

        # Bracketed secant / regula falsi with midpoint safeguard
        # If secant jump lands very close to root, test a tight +/- 0.45*tol bracket to terminate fast!
        it_b = 0
        while abs(outer - inner) > tol and it_b < 25:
            it_b += 1
            if (g_in is not None and g_out is not None and g_in < 0.0 < g_out and (g_out - g_in) > 1e-9):
                # Secant estimate inside (inner, outer)
                frac = -g_in / (g_out - g_in)
                # Keep inside [0.05, 0.95] of bracket unless bracket is already narrow; nudge slightly toward
                # inner/outer on alternating steps to ensure rapid 2-sided convergence
                frac = float(np.clip(frac, 0.05, 0.95))
                cand = inner + frac * (outer - inner)
                # If secant predicts very close to boundary, probe just past it by 0.35*tol to close bracket immediately
                if abs(cand - inner) < 0.35 * tol:
                    cand = inner + np.sign(outer - inner) * 0.45 * tol
                elif abs(outer - cand) < 0.35 * tol:
                    cand = outer - np.sign(outer - inner) * 0.45 * tol
            else:
                cand = 0.5 * (inner + outer)
            pc, atc, cc = ev(cand)
            gc = (atc - cc) if (atc is not None and cc is not None) else None
            if pc > alpha:
                inner, g_in = cand, gc
            else:
                outer, g_out = cand, gc
            # If secant converged tightly on one side, test the opposite side at distance 0.9*tol to verify bracket!
            if abs(outer - inner) > tol and it_b >= 2 and g_in is not None and g_out is not None:
                # Re-estimate root
                if g_in < 0.0 < g_out and (g_out - g_in) > 1e-9:
                    b_hat = inner + (-g_in / (g_out - g_in)) * (outer - inner)
                    if abs(b_hat - inner) < 0.4 * tol:
                        probe_hi = inner + np.sign(outer - inner) * 0.95 * tol
                        if abs(probe_hi - inner) < abs(outer - inner):
                            ph, ath, ch = ev(probe_hi)
                            gh = (ath - ch) if (ath is not None and ch is not None) else None
                            if ph <= alpha:
                                outer, g_out = probe_hi, gh
                            else:
                                inner, g_in = probe_hi, gh
                    elif abs(outer - b_hat) < 0.4 * tol:
                        probe_lo = outer - np.sign(outer - inner) * 0.95 * tol
                        if abs(outer - probe_lo) < abs(outer - inner):
                            pl, atl, cl = ev(probe_lo)
                            gl = (atl - cl) if (atl is not None and cl is not None) else None
                            if pl > alpha:
                                inner, g_in = probe_lo, gl
                            else:
                                outer, g_out = probe_lo, gl

        out.append(0.5 * (inner + outer))
    return float(out[0]), float(out[1])


# =============================================================================
# Data loading: 675 cities, census matrices, covariates, and cached MITMA flows
# =============================================================================
def load_census_and_design() -> dict:
    """Build the frozen 675-city x 46-destination-province frame (31,050 dyads, G=44)."""
    nodes = pd.read_parquet(NODES_PQ).sort_values("idx").reset_index(drop=True)
    assert len(nodes) == C.EXPECTED_N_NODES
    PROVS = np.array(sorted(nodes["prov"].unique()), dtype=str)
    assert len(PROVS) == 47
    P2I = {p: k for k, p in enumerate(PROVS)}
    NID = nodes["id"].to_numpy(dtype=str)
    NPI = np.array([P2I[p] for p in nodes["prov"].to_numpy(dtype=str)])
    NCOAST = nodes["coastal_true"].to_numpy(dtype=bool)
    NARTEFACT = nodes["id"].isin(C.ARTEFACT_ZONES).to_numpy(dtype=bool)
    COASTAL_PROVS = sorted(nodes.loc[nodes["coastal_true"], "prov"].unique())
    assert len(COASTAL_PROVS) == 19

    con = duckdb.connect()
    spec_sql = _sql_in(SPEC_MUNI_CODES)
    dated_sql = _sql_in(DATED_ARRIVAL_CODES)
    post81_sql = _sql_in(POST1981_DATED_CODES)
    mid_sql = _sql_in(MID_ARRIVAL_CODES)
    rec_sql = _sql_in(RECENT_ARRIVAL_CODES)

    m_tot = con.execute(f"""
        SELECT CPRO || CMUN AS mun, CPRO AS prov,
               count(*) AS n_tot,
               sum(CASE WHEN VAREDAD >= 65 THEN 1 ELSE 0 END) AS n_tot_65
        FROM read_parquet('{CENSO24_PQ}')
        WHERE CMUN NOT IN {spec_sql}
          AND CPRO NOT IN ('07','35','38','51','52')
        GROUP BY 1, 2 ORDER BY 1
    """).df()
    single_set = set(nodes.loc[nodes["n_municipalities"] == 1, "id"])
    m_single = m_tot[m_tot["mun"].isin(single_set)].copy()
    pcounts = m_single.groupby("prov").size()
    id_provs = sorted(pcounts[pcounts >= 2].index)
    m675 = m_single[m_single["prov"].isin(id_provs)].sort_values("mun").reset_index(drop=True)
    assert len(m675) == 675 and len(id_provs) == 44

    # Query ALL needed cell counts (both dated-only per pre-reg §1 AND raw V1 for audit diff)
    cp = con.execute(f"""
        SELECT CPRO || CMUN AS mun, CPRO_NAC AS dprov,
               count(*) AS n_born_raw,
               sum(CASE WHEN VARANOP IN {dated_sql} THEN 1 ELSE 0 END) AS n_born_dated,
               sum(CASE WHEN VARANOP = '7777' THEN 1 ELSE 0 END) AS n_special_7777,

               -- Arrival cohorts (explicit codes only; 0001 / mid / 2016-2023)
               sum(CASE WHEN VARANOP = '0001' THEN 1 ELSE 0 END) AS n_pre1981,
               sum(CASE WHEN VARANOP IN {mid_sql} THEN 1 ELSE 0 END) AS n_mid81_15,
               sum(CASE WHEN VARANOP IN {rec_sql} THEN 1 ELSE 0 END) AS n_rec2016,
               sum(CASE WHEN VARANOP IN {post81_sql} THEN 1 ELSE 0 END) AS n_post1981_dated,

               -- Birthplace size: dated arrivals only (pre-reg §1) AND raw V1 (including 7777)
               sum(CASE WHEN CMUN_NAC = '991' AND VARANOP IN {dated_sql} THEN 1 ELSE 0 END) AS n_vil2k_dated,
               sum(CASE WHEN CMUN_NAC IN ('992','993') AND VARANOP IN {dated_sql} THEN 1 ELSE 0 END) AS n_mid2_10k_dated,
               sum(CASE WHEN CMUN_NAC NOT IN {spec_sql} AND VARANOP IN {dated_sql} THEN 1 ELSE 0 END) AS n_urb10k_dated,
               sum(CASE WHEN CMUN_NAC = '991' THEN 1 ELSE 0 END) AS n_vil2k_raw,
               sum(CASE WHEN CMUN_NAC IN ('992','993') THEN 1 ELSE 0 END) AS n_mid2_10k_raw,
               sum(CASE WHEN CMUN_NAC NOT IN {spec_sql} THEN 1 ELSE 0 END) AS n_urb10k_raw,

               -- Age-matched cohorts (explicit dated arrival codes only!)
               sum(CASE WHEN VARANOP IN {dated_sql} AND VAREDAD >= 65 THEN 1 ELSE 0 END) AS n_age65_dated,
               sum(CASE WHEN VARANOP = '0001' AND VAREDAD BETWEEN 43 AND 64 THEN 1 ELSE 0 END) AS n_pre1981_43_64,
               sum(CASE WHEN VARANOP IN {post81_sql} AND VAREDAD BETWEEN 43 AND 64 THEN 1 ELSE 0 END) AS n_post1981_43_64_dated,
               sum(CASE WHEN VARANOP IN {mid_sql} AND VAREDAD BETWEEN 43 AND 64 THEN 1 ELSE 0 END) AS n_mid81_15_43_64,
               sum(CASE WHEN VARANOP IN {rec_sql} AND VAREDAD BETWEEN 43 AND 64 THEN 1 ELSE 0 END) AS n_rec2016_43_64,

               sum(CASE WHEN VARANOP = '0001' AND VAREDAD >= 65 THEN 1 ELSE 0 END) AS n_pre1981_ge65,
               sum(CASE WHEN VARANOP IN {post81_sql} AND VAREDAD >= 65 THEN 1 ELSE 0 END) AS n_post1981_ge65_dated,
               sum(CASE WHEN VARANOP != '0001' AND VAREDAD >= 65 THEN 1 ELSE 0 END) AS n_post1981_ge65_raw_v1,
               sum(CASE WHEN VARANOP IN {mid_sql} AND VAREDAD >= 65 THEN 1 ELSE 0 END) AS n_mid81_15_ge65,
               sum(CASE WHEN VARANOP IN {rec_sql} AND VAREDAD >= 65 THEN 1 ELSE 0 END) AS n_rec2016_ge65,

               sum(CASE WHEN VARANOP IN {post81_sql} AND VAREDAD < 43 THEN 1 ELSE 0 END) AS n_post1981_lt43_dated,
               sum(CASE WHEN VARANOP IN {rec_sql} AND VAREDAD < 43 THEN 1 ELSE 0 END) AS n_rec2016_lt43,
               sum(CASE WHEN VARANOP IN {rec_sql} AND VAREDAD < 65 THEN 1 ELSE 0 END) AS n_rec2016_lt65,

               -- Age-matched 60+ sensitivity (43-59 vs 60+)
               sum(CASE WHEN VARANOP = '0001' AND VAREDAD BETWEEN 43 AND 59 THEN 1 ELSE 0 END) AS n_pre1981_43_59,
               sum(CASE WHEN VARANOP IN {post81_sql} AND VAREDAD BETWEEN 43 AND 59 THEN 1 ELSE 0 END) AS n_post1981_43_59_dated,
               sum(CASE WHEN VARANOP = '0001' AND VAREDAD >= 60 THEN 1 ELSE 0 END) AS n_pre1981_ge60,
               sum(CASE WHEN VARANOP IN {post81_sql} AND VAREDAD >= 60 THEN 1 ELSE 0 END) AS n_post1981_ge60_dated,
               sum(CASE WHEN VARANOP IN {rec_sql} AND VAREDAD < 60 THEN 1 ELSE 0 END) AS n_rec2016_lt60,

               -- Age-matched birthplace size (dated arrivals only)
               sum(CASE WHEN CMUN_NAC = '991' AND VARANOP IN {dated_sql} AND VAREDAD < 43 THEN 1 ELSE 0 END) AS n_vil2k_lt43_dated,
               sum(CASE WHEN CMUN_NAC = '991' AND VARANOP IN {dated_sql} AND VAREDAD BETWEEN 43 AND 64 THEN 1 ELSE 0 END) AS n_vil2k_43_64_dated,
               sum(CASE WHEN CMUN_NAC = '991' AND VARANOP IN {dated_sql} AND VAREDAD < 65 THEN 1 ELSE 0 END) AS n_vil2k_lt65_dated,
               sum(CASE WHEN CMUN_NAC = '991' AND VARANOP IN {dated_sql} AND VAREDAD >= 65 THEN 1 ELSE 0 END) AS n_vil2k_ge65_dated,
               sum(CASE WHEN CMUN_NAC NOT IN {spec_sql} AND VARANOP IN {dated_sql} AND VAREDAD < 43 THEN 1 ELSE 0 END) AS n_urb10k_lt43_dated,
               sum(CASE WHEN CMUN_NAC NOT IN {spec_sql} AND VARANOP IN {dated_sql} AND VAREDAD BETWEEN 43 AND 64 THEN 1 ELSE 0 END) AS n_urb10k_43_64_dated,
               sum(CASE WHEN CMUN_NAC NOT IN {spec_sql} AND VARANOP IN {dated_sql} AND VAREDAD < 65 THEN 1 ELSE 0 END) AS n_urb10k_lt65_dated,
               sum(CASE WHEN CMUN_NAC NOT IN {spec_sql} AND VARANOP IN {dated_sql} AND VAREDAD >= 65 THEN 1 ELSE 0 END) AS n_urb10k_ge65_dated
        FROM read_parquet('{CENSO24_PQ}')
        WHERE CMUN NOT IN {spec_sql}
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
        M = np.zeros((C_M, 47), dtype=float)
        M[ri_cp, pj_cp] = cp675[col].to_numpy(dtype=float)
        MATS[col] = M

    # Assert exact identity of special-code drop
    assert int(MATS["n_special_7777"].sum()) == 11350
    assert int(MATS["n_post1981_ge65_raw_v1"].sum() - MATS["n_post1981_ge65_dated"].sum()) == 3245
    assert np.allclose(MATS["n_pre1981_43_64"] + MATS["n_pre1981_ge65"], MATS["n_pre1981"])

    OWN = MPROV[:, None] == PROVS[None, :]
    ri_all, pj_all = np.where(~OWN)
    assert len(ri_all) == 675 * 46

    # Distance to population-weighted province centroids (exact R2/V1 spec)
    pw = nodes["population"].to_numpy(dtype=float)
    nx, ny = nodes["x"].to_numpy(dtype=float), nodes["y"].to_numpy(dtype=float)
    PCX = np.array([np.average(nx[NPI == k], weights=pw[NPI == k]) for k in range(47)])
    PCY = np.array([np.average(ny[NPI == k], weights=pw[NPI == k]) for k in range(47)])
    MNODE = np.array([np.where(NID == m)[0][0] for m in MMUN])
    KM = np.hypot(nx[MNODE][:, None] - PCX[None, :], ny[MNODE][:, None] - PCY[None, :]) / 1000.0
    LKM = np.log(KM[ri_all, pj_all])
    lkm_c = LKM - LKM.mean()
    COV_DIST = np.column_stack([lkm_c, lkm_c ** 2])

    # Stratum-specific centroids (inland centroid vs coastal centroid) for B3 sensitivity
    PCX_inl = np.array([np.average(nx[(NPI == k) & (~NCOAST)], weights=pw[(NPI == k) & (~NCOAST)]) for k in range(47)])
    PCY_inl = np.array([np.average(ny[(NPI == k) & (~NCOAST)], weights=pw[(NPI == k) & (~NCOAST)]) for k in range(47)])
    PCX_cst = np.array([
        np.average(nx[(NPI == k) & NCOAST], weights=pw[(NPI == k) & NCOAST]) if np.any((NPI == k) & NCOAST) else np.nan
        for k in range(47)
    ])
    PCY_cst = np.array([
        np.average(ny[(NPI == k) & NCOAST], weights=pw[(NPI == k) & NCOAST]) if np.any((NPI == k) & NCOAST) else np.nan
        for k in range(47)
    ])
    KM_inl = np.hypot(nx[MNODE][:, None] - PCX_inl[None, :], ny[MNODE][:, None] - PCY_inl[None, :]) / 1000.0
    KM_cst = np.hypot(nx[MNODE][:, None] - PCX_cst[None, :], ny[MNODE][:, None] - PCY_cst[None, :]) / 1000.0
    lkm_inl = np.log(KM_inl[ri_all, pj_all])
    lkm_inl_c = lkm_inl - np.nanmean(lkm_inl)
    COV_DIST_INL = np.column_stack([lkm_inl_c, lkm_inl_c ** 2])
    lkm_cst = np.log(KM_cst[ri_all, pj_all])
    lkm_cst_c = lkm_cst - np.nanmean(lkm_cst)
    COV_DIST_CST = np.column_stack([lkm_cst_c, lkm_cst_c ** 2])

    # ADRH 2024 income
    adrh = pd.read_csv(ADRH_CSV, dtype={"ine_code": str})
    inc_s = adrh[(adrh.year == 2024) & (adrh.indicator == "Renta neta media por persona")].set_index("ine_code")["value"]
    INC = inc_s.reindex(MMUN).to_numpy(dtype=float)
    assert np.isfinite(INC).all()
    LINC = np.log(INC) - np.log(INC).mean()
    SLOPE_MAIN = [(pj_all, LINC[ri_all])]
    S65_C = SHARE65_CITY - SHARE65_CITY.mean()
    SLOPE_AGECITY = [(pj_all, LINC[ri_all]), (pj_all, S65_C[ri_all])]

    # R2 112 large cities mask (G = 22)
    bm68545 = pd.read_csv(BIRTH68545_CSV, dtype={"ine_code": str, "prov_birth": str})
    bm24_tot = bm68545[bm68545["year"] == 2024].groupby("ine_code")["persons"].sum()
    cities_137 = sorted(c for c in bm24_tot.index if c[:2] not in C.EXCLUDED_PROVINCES)
    c137_cnt = pd.Series([c[:2] for c in cities_137]).value_counts()
    id_p22 = set(c137_cnt[c137_cnt >= 2].index)
    m112_mask = np.isin(MMUN, [c for c in cities_137 if c[:2] in id_p22])
    assert m112_mask.sum() == 112
    cell_112 = m112_mask[ri_all]

    c_code = pd.factorize(MMUN[ri_all])[0]
    p_code = pd.factorize(np.char.add(np.char.add(MPROV[ri_all], "_"), PROVS[pj_all]))[0]
    CL_LEVELS = np.array(id_provs, dtype=str)
    cl_code = pd.Categorical(MPROV[ri_all], categories=CL_LEVELS).codes.astype(int)
    G_TOT = len(CL_LEVELS)  # 44

    COAST_PROV_IDX = np.array([P2I[p] for p in COASTAL_PROVS], dtype=int)
    cell_coastal_prov = np.isin(pj_all, COAST_PROV_IDX)

    return {
        "nodes": nodes, "PROVS": PROVS, "P2I": P2I, "NID": NID, "NPI": NPI,
        "NCOAST": NCOAST, "NARTEFACT": NARTEFACT,
        "COASTAL_PROVS": COASTAL_PROVS, "COAST_PROV_IDX": COAST_PROV_IDX,
        "m675": m675, "MMUN": MMUN, "MPROV": MPROV, "NTOT": NTOT,
        "M2R": M2R, "MNODE": MNODE, "OWN": OWN, "ri_all": ri_all, "pj_all": pj_all,
        "MATS": MATS, "COV_DIST": COV_DIST,
        "COV_DIST_INL": COV_DIST_INL, "COV_DIST_CST": COV_DIST_CST,
        "LINC": LINC, "S65_C": S65_C,
        "SLOPE_MAIN": SLOPE_MAIN, "SLOPE_AGECITY": SLOPE_AGECITY,
        "c_code": c_code, "p_code": p_code, "cl_code": cl_code,
        "CL_LEVELS": CL_LEVELS, "G_TOT": G_TOT,
        "cell_112": cell_112, "cell_coastal_prov": cell_coastal_prov,
    }


def load_mitma_flows(ctx: dict, dates: list[str], exclude_artefact_dest: bool = True) -> dict[str, dict[str, np.ndarray]]:
    """Aggregate cached MITMA parquets to (city_675, dest_prov) flows per night & stratum.

    Only reads coursework cache parquets (PERNOCT_2024, PERNOCT_2022_23).
    Never opens _local/external_journal/mitma_pernoct/.
    """
    nodes = ctx["nodes"]
    PROVS = ctx["PROVS"]
    P2I = ctx["P2I"]
    M2R = ctx["M2R"]
    OWN = ctx["OWN"]
    ri_all, pj_all = ctx["ri_all"], ctx["pj_all"]
    COAST_PROV_IDX = ctx["COAST_PROV_IDX"]
    C_M = len(ctx["MMUN"])

    valid_dest_df = nodes[["id", "prov", "coastal_true"]].copy()
    if exclude_artefact_dest:
        valid_dest_df = valid_dest_df[~valid_dest_df["id"].isin(C.ARTEFACT_ZONES)].copy()
    valid_orig_df = pd.DataFrame({"r": ctx["MMUN"]})

    dates_24 = sorted(d for d in dates if d.startswith("2024"))
    dates_22_23 = sorted(d for d in dates if d[:4] in ("2022", "2023"))
    assert len(dates_24) + len(dates_22_23) == len(dates), "Unexpected year outside 2022-2024 cache"

    con = duckdb.connect()
    con.register("valid_dest", valid_dest_df)
    con.register("valid_orig", valid_orig_df)

    dfs = []
    for path, dlist in ((PERNOCT_2024, dates_24), (PERNOCT_2022_23, dates_22_23)):
        if not dlist:
            continue
        con.register("req_dates", pd.DataFrame({"d": dlist}))
        q = con.execute(f"""
            SELECT CAST(t.date AS VARCHAR) AS date,
                   t.r AS city,
                   vd.prov AS dprov,
                   vd.coastal_true AS is_coast,
                   sum(t.p) AS w
            FROM read_parquet('{path}') t
            JOIN req_dates rd ON CAST(t.date AS VARCHAR) = rd.d
            JOIN valid_orig vo ON t.r = vo.r
            JOIN valid_dest vd ON t.s = vd.id
            WHERE t.r != t.s
            GROUP BY 1, 2, 3, 4
        """).df()
        dfs.append(q)
    agg = pd.concat(dfs, ignore_index=True)

    FLOWS: dict[str, dict[str, np.ndarray]] = {}
    for d in dates:
        sub = agg[agg["date"] == d]
        assert not sub.empty, f"Date {d} missing from coursework cache"
        r_idx = sub["city"].map(M2R).to_numpy(dtype=int)
        p_idx = sub["dprov"].map(P2I).to_numpy(dtype=int)
        w_val = sub["w"].to_numpy(dtype=float)
        is_c = sub["is_coast"].to_numpy(dtype=bool)

        M_all = np.zeros((C_M, 47), dtype=float)
        M_inl = np.zeros((C_M, 47), dtype=float)
        M_cst = np.zeros((C_M, 47), dtype=float)
        np.add.at(M_all, (r_idx, p_idx), w_val)
        np.add.at(M_inl, (r_idx[~is_c], p_idx[~is_c]), w_val[~is_c])
        np.add.at(M_cst, (r_idx[is_c], p_idx[is_c]), w_val[is_c])
        M_all[OWN] = 0.0
        M_inl[OWN] = 0.0
        M_cst[OWN] = 0.0
        M_inlp = M_all.copy()
        M_inlp[:, COAST_PROV_IDX] = 0.0
        FLOWS[d] = {
            "all": M_all[ri_all, pj_all],
            "inland": M_inl[ri_all, pj_all],
            "coastal": M_cst[ri_all, pj_all],
            "inlandprov": M_inlp[ri_all, pj_all],
        }
    return FLOWS


def add_window_mean(FLOWS: dict, key: str, dlist: list[str]) -> None:
    """Audit C3 window mean: sum(flow) / len(dlist), missing nights = 0."""
    n_d = float(len(dlist))
    FLOWS[key] = {
        k: sum(FLOWS[d][k] for d in dlist) / n_d
        for k in ("all", "inland", "coastal", "inlandprov")
    }


def cached_night_sets() -> dict[str, list[str]]:
    """Return canonical night lists strictly from coursework cache parquets."""
    con = duckdb.connect()
    d24 = [str(r[0]) for r in con.execute(f"SELECT DISTINCT date FROM read_parquet('{PERNOCT_2024}') ORDER BY 1").fetchall()]
    d22_23 = [str(r[0]) for r in con.execute(f"SELECT DISTINCT date FROM read_parquet('{PERNOCT_2022_23}') ORDER BY 1").fetchall()]

    NW = pd.read_csv(NET_DIR / "networks.csv", dtype=str)
    feb24_sis = sorted(NW.loc[(NW.group == "ordinary") & (NW.role == "sister"), "date"])
    aug24_sis = list(C.SISTER_NIGHTS["summer"])
    nov24_sis = list(C.SISTER_NIGHTS["allsaints"])
    oct24_all15 = [d for d in d24 if d[5:7] == "10" and pd.Timestamp(d).day_of_week in (1, 2, 3)]
    oct24_clean12 = [d for d in oct24_all15 if d not in DANA_DATES]
    oct24_no_val10 = [d for d in oct24_clean12 if d not in VAL_HOL_OCT_DATES]

    feb22_sis = [d for d in d22_23 if d.startswith("2022-02") and pd.Timestamp(d).day_of_week in (1, 2, 3)]
    aug22_sis = [d for d in d22_23 if d.startswith("2022-08") and 10 <= int(d[8:10]) <= 18]
    feb23_sis = [d for d in d22_23 if d.startswith("2023-02") and pd.Timestamp(d).day_of_week in (1, 2, 3)]
    aug23_sis = [d for d in d22_23 if d.startswith("2023-08") and 10 <= int(d[8:10]) <= 18]

    assert len(feb24_sis) == 13 and len(aug24_sis) == 9
    assert len(oct24_all15) == 15 and len(oct24_clean12) == 12 and len(oct24_no_val10) == 10
    assert len(feb22_sis) == 12 and len(aug22_sis) == 9
    assert len(feb23_sis) == 12 and len(aug23_sis) == 9

    all_needed = sorted(
        set(feb24_sis + aug24_sis + nov24_sis + oct24_all15 + feb22_sis + aug22_sis + feb23_sis + aug23_sis)
        | {C.APPENDIX_NIGHTS["christmas"], C.APPENDIX_NIGHTS["goodfriday"], "2022-11-01"}
    )
    return {
        "feb24_sis": feb24_sis,
        "aug24_sis": aug24_sis,
        "nov24_sis": nov24_sis,
        "oct24_all15": oct24_all15,
        "oct24_clean12": oct24_clean12,
        "oct24_no_val10": oct24_no_val10,
        "feb22_sis": feb22_sis,
        "aug22_sis": aug22_sis,
        "feb23_sis": feb23_sis,
        "aug23_sis": aug23_sis,
        "all_needed": all_needed,
    }


def make_design(
    ctx: dict,
    y: np.ndarray,
    X_roots: np.ndarray,
    mask: np.ndarray | None = None,
    cov_dist: np.ndarray | None = None,
    slope_raw=None,
    drop_cl: int | None = None,
) -> tuple[V1.SingleNightDesign, np.ndarray, int]:
    """Build a SingleNightDesign with optional row subsetting (re-indexing cluster codes)."""
    cov = ctx["COV_DIST"] if cov_dist is None else cov_dist
    slp = ctx["SLOPE_MAIN"] if slope_raw is None else slope_raw
    if mask is None:
        nt = V1.SingleNightDesign(y, X_roots, cov, ctx["c_code"], ctx["p_code"], ctx["cl_code"],
                                  slope_raw=slp, drop_cl=drop_cl)
        return nt, ctx["cl_code"], ctx["G_TOT"]
    m = np.asarray(mask, dtype=bool)
    cl_sub, cl_uniques = pd.factorize(ctx["MPROV"][ctx["ri_all"]][m], sort=True)
    G_sub = len(cl_uniques)
    slp_sub = [(gc[m], sv[m]) for gc, sv in slp] if slp else None
    nt = V1.SingleNightDesign(y[m], X_roots[m], cov[m], ctx["c_code"][m], ctx["p_code"][m], cl_sub,
                              slope_raw=slp_sub, drop_cl=drop_cl)
    return nt, cl_sub, G_sub
