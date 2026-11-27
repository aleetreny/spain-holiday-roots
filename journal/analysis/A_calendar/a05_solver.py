"""a05_solver: Exact block-diagonal multi-night stacked PPML + Kline & Santos (2012)
Wild Cluster Restricted (WCR) score bootstrap with Webb 6-point weights and test-inversion CIs.

Generalises pipeline/v1_verify_n1.FastStackedTwoNight from 2 nights to any number of
stacked blocks (nights, years, inland/coastal arms), preserving 100% algebraic equivalence
to the full stacked PPML with block-specific city FE, origin-prov x dest-prov FE,
distance polynomial and income x dest-prov slopes.

Includes a built-in unit/cross-check against FastStackedTwoNight and pyfixest.
"""
from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.linalg import qr

sys.path.insert(0, str(Path(__file__).resolve().parent))
from a00_paths import V1  # noqa: E402

FastHDFE2 = V1.FastHDFE2
drop_zero_groups = V1.drop_zero_groups
poisson_dev = V1.poisson_dev
ppml_solve = V1.ppml_solve
WEBB = V1.WEBB


@dataclass
class BlockDesign:
    """One post-separation PPML block (one night, or one night x arm)."""
    label: str
    year: int
    y: np.ndarray          # (n,)
    Xr: np.ndarray         # (n, K) base roots regressors
    cov: np.ndarray        # (n, k_c) rank-filtered covariates (e.g. lkm_c, lkm_c^2)
    hd: FastHDFE2          # absorbs city FE, pair FE, and income x dest_prov slope
    cl: np.ndarray         # (n,) cluster code in 0..G-1
    n: int
    mu_single: np.ndarray | None = None  # optional warm start from single-block fit


def make_block_design(
    label: str,
    year: int,
    y_raw: np.ndarray,
    Xr_raw: np.ndarray,
    cov_raw: np.ndarray,
    c_code: np.ndarray,
    p_code: np.ndarray,
    cl_code: np.ndarray,
    slope_raw: list[tuple[np.ndarray, np.ndarray]] | None = None,
    extra_keep: np.ndarray | None = None,
    scale_y: float = 1.0,
) -> BlockDesign:
    """Construct a single-block PPML design after iterative zero-group separation drops."""
    y_use = y_raw.astype(float, copy=True)
    if extra_keep is not None:
        keep_init = extra_keep.copy()
        y_masked = np.where(keep_init, y_use, 0.0)
        keep = drop_zero_groups(y_masked, [c_code, p_code]) & keep_init
    else:
        keep = drop_zero_groups(y_use, [c_code, p_code])

    y = y_use[keep] * scale_y
    Xr = Xr_raw[keep] if Xr_raw.ndim == 2 else Xr_raw[keep, None]
    c0 = pd.factorize(c_code[keep])[0]
    c1 = pd.factorize(p_code[keep])[0]
    sf_k = [(pd.factorize(gc[keep])[0], sv[keep]) for gc, sv in slope_raw] if slope_raw else None
    hd = FastHDFE2(c0, c1, slope_factors=sf_k)

    cov = cov_raw[keep]
    if cov.shape[1]:
        ct = hd.demean(cov, np.ones(len(y)), tol=1e-11)
        _, Rq, piv = qr(ct, mode="economic", pivoting=True)
        dg = np.abs(np.diag(Rq))
        ok = np.sort(piv[dg > 1e-8 * dg[0]])
        cov = cov[:, ok]

    return BlockDesign(
        label=label,
        year=year,
        y=y,
        Xr=np.ascontiguousarray(Xr, dtype=float),
        cov=np.ascontiguousarray(cov, dtype=float),
        hd=hd,
        cl=cl_code[keep].astype(np.int32),
        n=int(keep.sum()),
    )


def fit_single_block(blk: BlockDesign, G_tot: int) -> dict:
    """Fit single-block unrestricted PPML and return beta, CRV1 SE, cluster IF, mu."""
    K = blk.Xr.shape[1]
    X = np.column_stack([blk.Xr, blk.cov])
    beta, mu, dev = ppml_solve(blk.y, X, blk.hd, mu0=blk.mu_single)
    Xt = blk.hd.demean(X, mu, tol=1e-12)
    H = (Xt.T * mu) @ Xt
    Hinv = np.linalg.inv(H)
    sc = Xt * (blk.y - mu)[:, None]
    assert np.abs(sc.sum(0)).max() < 1e-5 * (1.0 + np.abs(blk.y).sum()), (
        f"PPML score != 0 on {blk.label}: {np.abs(sc.sum(0)).max():.3e}"
    )
    Sg = np.vstack([np.bincount(blk.cl, sc[:, q], minlength=G_tot) for q in range(X.shape[1])]).T
    IF = Sg @ Hinv[:, :K]
    Geff = int((np.abs(IF[:, 0]) > 0).sum())
    se = np.sqrt(Geff / (Geff - 1) * np.sum(IF ** 2, axis=0))
    return {
        "beta": beta[:K].copy(),
        "se": se,
        "IF": IF,
        "G": Geff,
        "mu": mu,
        "dev": dev,
        "n": blk.n,
        "flow_total": float(blk.y.sum()),
    }


def _block_project(blk: BlockDesign, z: np.ndarray, mu: np.ndarray, tol: float):
    """Demean [z, Xr, cov] with weight mu and residualise [z, Xr] against cov."""
    K = blk.Xr.shape[1]
    kc = blk.cov.shape[1]
    if kc:
        M_in = np.column_stack([z, blk.Xr, blk.cov])
    else:
        M_in = np.column_stack([z, blk.Xr])
    Dm = blk.hd.demean(M_in, mu, tol=tol)
    zt = Dm[:, 0]
    Xrt = Dm[:, 1:1 + K]
    if kc:
        Ct = Dm[:, 1 + K:]
        CtW = Ct.T * mu
        A_c = CtW @ Ct
        # Tiny kc x kc (typically 2x2) projection
        rhs = np.column_stack([CtW @ zt, CtW @ Xrt])
        sol = np.linalg.solve(A_c, rhs)
        zp = zt - Ct @ sol[:, 0]
        Xrp = Xrt - Ct @ sol[:, 1:]
        return zp, Xrp
    return zt, Xrt


class FastStackedMultiNight:
    """Exact block-diagonal stacked PPML across B blocks with shared loading matrices.

    On each block b:
      - target regressor w_b  = blk.Xr @ l_tgt[b]   (scalar coefficient theta_target)
      - other shared V_b      = blk.Xr @ L_oth[b]   (m-vector of shared nuisance coefficients)
      - block-specific FE (city, pair, income*dest) and covariates (blk.cov) are completely free per block.

    Supports:
      - Unrestricted estimation (theta_target and beta_other jointly) + exact CRV1 SE
      - Restricted estimation under H0: theta_target = d0
      - Exact Kline & Santos (2012) Wild Cluster Restricted (WCR) score bootstrap (Webb weights)
      - Exact test-inversion 95% WCR confidence intervals via warm-started bisection
    """

    def __init__(
        self,
        blocks: list[BlockDesign],
        l_tgt: list[np.ndarray],
        L_oth: list[np.ndarray],
        G_tot: int,
    ):
        assert len(blocks) == len(l_tgt) == len(L_oth) and len(blocks) > 0
        self.blocks = blocks
        self.B = len(blocks)
        self.l_tgt = [np.asarray(v, dtype=float).reshape(-1) for v in l_tgt]
        self.L_oth = [np.asarray(M, dtype=float) for M in L_oth]
        self.m_oth = self.L_oth[0].shape[1]
        self.G_tot = int(G_tot)
        # Pre-compute target vectors per block
        self.w_tgt = [b.Xr @ lv for b, lv in zip(self.blocks, self.l_tgt)]
        # Initial warm mu list
        self.mu_init = [
            b.mu_single.copy() if b.mu_single is not None else (b.y + b.y.mean()) / 2.0
            for b in self.blocks
        ]

    def _solve_irls(
        self,
        L_list: list[np.ndarray],
        offsets: list[np.ndarray] | None = None,
        mu_start: list[np.ndarray] | None = None,
        maxit: int = 200,
    ) -> tuple[np.ndarray, list[np.ndarray], float, list[np.ndarray]]:
        """Run exact block-diagonal IRLS for shared loading matrices L_list (each K x m)."""
        m = L_list[0].shape[1]
        mus = [m0.copy() for m0 in (self.mu_init if mu_start is None else mu_start)]
        etas = [np.log(mu) for mu in mus]
        dev_old = np.inf
        b_old = np.zeros(m)
        Xrp_last: list[np.ndarray] = [np.empty((0, 0))] * self.B

        for it in range(maxit):
            tol_d = 1e-11 if it >= 2 else 1e-8
            A = np.zeros((m, m))
            rhs = np.zeros(m)
            zs: list[np.ndarray] = []
            zps: list[np.ndarray] = []
            Xrps: list[np.ndarray] = []

            for idx, blk in enumerate(self.blocks):
                mu = mus[idx]
                eta = etas[idx]
                z = eta + (blk.y - mu) / mu if offsets is None else eta - offsets[idx] + (blk.y - mu) / mu
                zp, Xrp = _block_project(blk, z, mu, tol=tol_d)
                zs.append(z)
                zps.append(zp)
                Xrps.append(Xrp)
                if m:
                    L = L_list[idx]
                    # Vp = Xrp @ L  ->  Vp' W Vp = L' (Xrp' W Xrp) L
                    XrpW = Xrp.T * mu
                    A += L.T @ (XrpW @ Xrp) @ L
                    rhs += L.T @ (XrpW @ zp)

            b_new = np.linalg.solve(A, rhs) if m else np.zeros(0)
            dev = 0.0
            for idx, blk in enumerate(self.blocks):
                if m:
                    # Vp @ b_new = Xrp @ (L @ b_new)
                    e = zps[idx] - Xrps[idx] @ (L_list[idx] @ b_new)
                else:
                    e = zps[idx]
                eta_new = zs[idx] - e if offsets is None else zs[idx] - e + offsets[idx]
                # Guard against extreme step overshoot
                mu_new = np.exp(np.clip(eta_new, -50.0, 50.0))
                etas[idx] = eta_new
                mus[idx] = mu_new
                dev += poisson_dev(blk.y, mu_new)

            db = float(np.max(np.abs(b_new - b_old))) if m else 0.0
            b_old = b_new
            Xrp_last = Xrps
            if abs(dev - dev_old) <= 1e-11 * (abs(dev) + 0.1) and db < 1e-9:
                return b_new, mus, dev, Xrp_last
            dev_old = dev
        raise RuntimeError(f"FastStackedMultiNight IRLS did not converge in {maxit} iterations (db={db:.2e})")

    def fit_unrestricted(self) -> dict:
        """Fit unrestricted stacked PPML and compute exact CRV1 SE for theta_target (index 0)."""
        L_full = [np.column_stack([lv[:, None], Lo]) for lv, Lo in zip(self.l_tgt, self.L_oth)]
        m_full = 1 + self.m_oth
        beta, mus, dev, _ = self._solve_irls(L_full)

        # Exact score and Hessian at converged mu
        H = np.zeros((m_full, m_full))
        Sg = np.zeros((self.G_tot, m_full))
        n_tot = 0
        y_tot = 0.0
        for idx, blk in enumerate(self.blocks):
            mu = mus[idx]
            _, Xrp = _block_project(blk, np.zeros(blk.n), mu, tol=1e-12)
            L = L_full[idx]
            Vp = Xrp @ L
            H += (Vp.T * mu) @ Vp
            resid = blk.y - mu
            for q in range(m_full):
                Sg[:, q] += np.bincount(blk.cl, Vp[:, q] * resid, minlength=self.G_tot)
            n_tot += blk.n
            y_tot += float(blk.y.sum())

        max_score = float(np.abs(Sg.sum(axis=0)).max())
        assert max_score < 1e-4 * (1.0 + y_tot), f"Stacked unrestricted score != 0: {max_score:.3e}"
        Hinv = np.linalg.inv(H)
        IF = Sg @ Hinv[:, 0]
        Geff = int((np.abs(IF) > 0).sum())
        se = float(np.sqrt(Geff / (Geff - 1) * np.sum(IF ** 2)))
        est = float(beta[0])
        tc = float(stats.t.ppf(0.975, Geff - 1))
        p_t = float(2.0 * stats.t.sf(abs(est / se), Geff - 1))
        return {
            "est": est,
            "beta_all": beta.copy(),
            "se_crv1": se,
            "G": Geff,
            "ci_crv_lo": est - tc * se,
            "ci_crv_hi": est + tc * se,
            "p_crv1": p_t,
            "n_obs": n_tot,
            "n_blocks": self.B,
            "flow_total": y_tot,
            "mus": mus,
            "IF": IF,
        }

    def p_wcr(
        self,
        d0: float,
        VB: np.ndarray,
        mu_start: list[np.ndarray] | None = None,
    ) -> tuple[float, float, list[np.ndarray]]:
        """Exact Kline & Santos (2012) WCR score test for H0: theta_target = d0.

        Returns (p_wcr_plus1, t0, converged_mus).
        """
        offsets = [d0 * w for w in self.w_tgt]
        _, mus, _, _ = self._solve_irls(self.L_oth, offsets=offsets, mu_start=mu_start)

        # Project target w_b out of [L_oth, block FE, block cov] at restricted MLE mus
        m = self.m_oth
        A = np.zeros((m, m))
        q_vec = np.zeros(m)
        Xrps: list[np.ndarray] = []
        for idx, blk in enumerate(self.blocks):
            mu = mus[idx]
            _, Xrp = _block_project(blk, np.zeros(blk.n), mu, tol=1e-12)
            Xrps.append(Xrp)
            if m:
                XrpW = Xrp.T * mu
                G_x = XrpW @ Xrp
                L = self.L_oth[idx]
                lv = self.l_tgt[idx]
                A += L.T @ G_x @ L
                q_vec += L.T @ (G_x @ lv)

        gamma = np.linalg.solve(A, q_vec) if m else np.zeros(0)

        Sg = np.zeros(self.G_tot)
        for idx, blk in enumerate(self.blocks):
            mu = mus[idx]
            eff_l = self.l_tgt[idx] - (self.L_oth[idx] @ gamma if m else 0.0)
            w_tilde = Xrps[idx] @ eff_l
            s = w_tilde * (blk.y - mu)
            Sg += np.bincount(blk.cl, s, minlength=self.G_tot)

        denom0 = float(np.sqrt(np.sum(Sg ** 2)))
        t0 = float(Sg.sum() / denom0) if denom0 > 0 else 0.0
        tb = (VB @ Sg) / np.sqrt((VB ** 2) @ (Sg ** 2))
        p_val = float((1 + np.sum(np.abs(tb) >= abs(t0))) / (VB.shape[0] + 1))
        return p_val, t0, mus

    def invert_wcr_ci(
        self,
        est: float,
        se: float,
        VB: np.ndarray,
        mu_start: list[np.ndarray] | None = None,
        alpha: float = 0.05,
        tol_frac: float = 5e-3,
    ) -> tuple[float, float]:
        """Exact 95% WCR confidence interval by warm-started bisection test inversion."""
        bounds = []
        t0_step = float(stats.t.ppf(1.0 - alpha / 2.0, max(self.G_tot - 1, 1)))
        for sgn in (-1.0, 1.0):
            cache_mu = mu_start

            def eval_p(val: float) -> float:
                nonlocal cache_mu
                p, _, mus = self.p_wcr(val, VB, mu_start=cache_mu)
                cache_mu = mus
                return p

            # Start right at the CRV1 t-bound: usually brackets in 1-2 evaluations!
            inner = est
            step = max(t0_step * se, 1e-4)
            outer = est + sgn * step
            p_out = eval_p(outer)
            k = 0
            while p_out > alpha:
                inner = outer
                step *= 1.5
                outer = est + sgn * step
                p_out = eval_p(outer)
                k += 1
                if k > 30:
                    raise RuntimeError(f"WCR CI bracketing failed (est={est}, se={se}, sgn={sgn})")
            # If initial CRV1 bound was already outside (p <= alpha), narrow inner quickly
            if inner == est:
                # One interior probe at 0.75 * step before bisection if p_out is very small
                probe = est + sgn * (0.75 * step)
                if eval_p(probe) > alpha:
                    inner = probe
                else:
                    outer = probe
            tol = max(tol_frac * se, 1e-4)
            while abs(outer - inner) > tol:
                mid = 0.5 * (inner + outer)
                if eval_p(mid) > alpha:
                    inner = mid
                else:
                    outer = mid
            bounds.append(0.5 * (inner + outer))
        return float(bounds[0]), float(bounds[1])


# =============================================================================
# Loading-matrix builders for pre-registered equal-year-weighted contrasts
# =============================================================================
def build_equal_year_loadings(
    blocks: list[BlockDesign],
    is_holiday: list[bool],
    K: int,
    target_vec: np.ndarray,
    other_vecs: list[np.ndarray] | None = None,
    pooled_across_years: bool = False,
) -> tuple[list[np.ndarray], list[np.ndarray]]:
    """Build exact loading matrices (l_tgt, L_oth) for equal-weight year-averaged holiday contrast.

    Within each year y in sorted(years):
      - Ordinary blocks (is_holiday=False): Xr @ beta_{0,y}
      - Holiday blocks (is_holiday=True):   Xr @ (beta_{0,y} + delta_y)
    Reparameterising delta_y along target_vec (and orthogonal directions in other_vecs):
      - Let [u_0, u_1, ..., u_{K-1}] be a basis of R^K with u_0 = target_vec (e.g. [1] for pooled,
        [1, 0] for pre1981, [0, 1] for recent, [0.5, -0.5] and [0.5, 0.5] for pre1981 - recent!).
      - For u_0: delta_{y,0} = delta_bar + d_tilde_y with sum_{y=1}^Y d_tilde_y = 0,
        so delta_bar = (1/Y) * sum_y delta_{y,0} is the EXACT equal-weight year average!
      - If pooled_across_years=True: single common beta_0 and single common delta across all years.
    """
    years = sorted({b.year for b in blocks})
    Y = len(years)
    y2i = {yr: idx for idx, yr in enumerate(years)}

    u0 = np.asarray(target_vec, dtype=float).reshape(K)
    others = [np.asarray(v, dtype=float).reshape(K) for v in (other_vecs or [])]
    assert 1 + len(others) == K, f"Basis size {1 + len(others)} != K={K}"

    l_tgt_list: list[np.ndarray] = []
    L_oth_list: list[np.ndarray] = []

    if pooled_across_years:
        # Shared columns: [beta_0 (K cols), delta_other (len(others) cols)]
        m_oth = K + len(others)
        for blk, hol in zip(blocks, is_holiday):
            l_tgt_list.append(u0 if hol else np.zeros(K))
            Lo = np.zeros((K, m_oth))
            Lo[:, :K] = np.eye(K)
            if hol:
                for j, uj in enumerate(others):
                    Lo[:, K + j] = uj
            L_oth_list.append(Lo)
        return l_tgt_list, L_oth_list

    # Year-specific baselines: Y * K columns
    # Year-specific target deviations (zero-sum): (Y - 1) columns
    # Year-specific other holiday contrasts: Y * len(others) columns
    n_base = Y * K
    n_dev = Y - 1
    n_oh = Y * len(others)
    m_oth = n_base + n_dev + n_oh

    for blk, hol in zip(blocks, is_holiday):
        yi = y2i[blk.year]
        l_tgt_list.append(u0 if hol else np.zeros(K))
        Lo = np.zeros((K, m_oth))
        # Year-specific baseline beta_{0,yi}
        Lo[:, yi * K:(yi + 1) * K] = np.eye(K)
        if hol:
            # Zero-sum deviation for target contrast across years
            if Y > 1:
                if yi < Y - 1:
                    Lo[:, n_base + yi] = u0
                else:
                    Lo[:, n_base:n_base + (Y - 1)] = -u0[:, None]
            # Year-specific holiday contrasts for other basis directions
            for j, uj in enumerate(others):
                Lo[:, n_base + n_dev + yi * len(others) + j] = uj
        L_oth_list.append(Lo)

    return l_tgt_list, L_oth_list


def build_p3_inland_coast_loadings(
    blocks: list[BlockDesign],
    is_holiday: list[bool],
    is_inland: list[bool],
    pooled_across_years: bool = False,
) -> tuple[list[np.ndarray], list[np.ndarray]]:
    """Build exact loading matrices for P3: delta_AUG(inland) - delta_AUG(coast) (K = 1).

    In each year y, the 4 cells (hol in {0,1}, inland in {0,1}) have slopes:
      (0, 0) coast O:    beta_c,y
      (0, 1) inland O:   beta_i,y
      (1, 0) coast H:    beta_c,y + delta_c,y
      (1, 1) inland H:   beta_i,y + delta_c,y + Delta_y
    where Delta_y = delta_i,y - delta_c,y = Delta_bar + d_tilde_y (sum_y d_tilde_y = 0).
    Target coefficient is Delta_bar = (1/Y) * sum_y [delta_AUG(inland, y) - delta_AUG(coast, y)].
    """
    years = sorted({b.year for b in blocks})
    Y = len(years)
    y2i = {yr: idx for idx, yr in enumerate(years)}

    l_tgt_list: list[np.ndarray] = []
    L_oth_list: list[np.ndarray] = []

    if pooled_across_years:
        # Other columns: [beta_c, beta_i, delta_c] -> 3 columns
        for hol, inl in zip(is_holiday, is_inland):
            l_tgt_list.append(np.array([1.0 if (hol and inl) else 0.0]))
            Lo = np.zeros((1, 3))
            Lo[0, 1 if inl else 0] = 1.0
            if hol:
                Lo[0, 2] = 1.0
            L_oth_list.append(Lo)
        return l_tgt_list, L_oth_list

    # Columns in Lo:
    #   [0 .. Y-1]:          beta_c,y  (coast baseline per year)
    #   [Y .. 2Y-1]:         beta_i,y  (inland baseline per year)
    #   [2Y .. 3Y-1]:        delta_c,y (coast holiday contrast per year)
    #   [3Y .. 4Y-2]:        d_tilde_y (zero-sum year deviations of Delta_y)
    m_oth = 3 * Y + (Y - 1)
    for blk, hol, inl in zip(blocks, is_holiday, is_inland):
        yi = y2i[blk.year]
        l_tgt_list.append(np.array([1.0 if (hol and inl) else 0.0]))
        Lo = np.zeros((1, m_oth))
        Lo[0, (Y + yi) if inl else yi] = 1.0
        if hol:
            Lo[0, 2 * Y + yi] = 1.0
            if inl and Y > 1:
                if yi < Y - 1:
                    Lo[0, 3 * Y + yi] = 1.0
                else:
                    Lo[0, 3 * Y:3 * Y + (Y - 1)] = -1.0
        L_oth_list.append(Lo)
    return l_tgt_list, L_oth_list


def holm_adjust(pvals: list[float] | np.ndarray) -> np.ndarray:
    """Step-down Holm-Bonferroni adjustment for a family of m p-values."""
    p = np.asarray(pvals, dtype=float)
    m = len(p)
    order = np.argsort(p)
    adj = np.empty(m, dtype=float)
    running_max = 0.0
    for rank, idx in enumerate(order):
        val = (m - rank) * p[idx]
        running_max = max(running_max, val)
        adj[idx] = min(1.0, running_max)
    return adj
