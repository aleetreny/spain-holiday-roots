"""r1_villages: Result 1 at node level. Do emptied villages gain a place in the holiday network
beyond what size, zone aggregation, distance to cities and coast predict?

Inputs (all produced by the clean core; nothing is rebuilt from raw MITMA files):
  results/b01_networks/nodes.parquet                node table written by b00_common.load_nodes()
  results/b01_networks/networks.csv                 (date, group, role) of every network in config
  results/b01_networks/edges_<date>__full.parquet   full directed weighted matrix, w = persons
  results/b01_networks/edges_<date>__w10.parquet    main binary backbone (w >= 10)
  results/b01_networks/edges_<date>__equal_edges.parquet   robustness backbone
  results/b02_coreness.parquet                      coreness of the undirected projection of each backbone
  _local/external/ine_turismo_moviles/raw/turismo_53001_provorig_x_mundest.csv.gz   INE table 53001
  _local/external/ine_historical_pop/municipal_population_1900_2025.csv   (b00_common.HIST_POP_CSV)

Definitions (frozen here; see .planning/2026-11-24-consolidation/decisions.md, lines tagged [R1]):
  loss        = loss_peak_2021 of the node table = 1 - pop_2021 / max(pop_1900, ..., pop_1960), census
                columns summed over the INE municipalities of the zone (b00_common.historical_loss).
                Zones with missing loss are excluded and counted (results/r1_sample.csv).
  loss_pct    = percentile rank of loss among the zones of the "all" sample (0 = grew most, 1 = emptied most).
                The regression coefficient is the change in the outcome from the least to the most emptied zone.
  s_in        = in-strength on the full matrix = non-resident persons sleeping in the zone that night.
  k_in        = in-degree on the backbone (w10 main, equal_edges check).
  core        = coreness of the undirected projection of the backbone (b02_coreness.parquet).
  n_prov      = number of origin provinces OTHER than the zone's own whose summed flow into the zone is >= 10 persons.
  eff_prov    = inverse Simpson index 1 / sum_p q_p^2 over ALL origin provinces (own included), q_p = share of s_in.
  share_far   = share of s_in coming from origin zones more than 100 km away.
  k_out       = out-degree on w10; core_in = in-coreness of the directed w10 backbone (igraph mode "in").
  s_out       = residents of the zone sleeping in another mainland zone; away_share = s_out / population.
  net_per_res = (s_in - s_out) / population; net_bal = (s_in - s_out) / (s_in + s_out);
  net_importer = 1[s_in > s_out]. The gravity null preserves s_in AND s_out, so it reproduces net balance exactly.
  Nights: config MAIN_NIGHTS, SISTER_NIGHTS (ordinary resolved in networks.csv), YEAR_REPLICATION,
  APPENDIX_NIGHTS (vs the ordinary main night). File date d = night ending on d. Nights b01 has not built are
  built with b00_common.build_edges + backbones and cached in results/r1_cache/.
  config.ARTEFACT_ZONES are excluded from every sample (and from the loss ranking).
  Volume curve (pattern #3): the loss gradient of each 2024 night (306 nights from results/p_zone_daily.parquet,
  cross-checked against this script's s_in, s_out) regressed on log(total away persons) over non-holiday 2024
  nights; residuals of holiday nights in sd units.
  Distances: Euclidean between zone representative points (EPSG:3042), km.
  dist_city   = distance to the nearest zone with population >= 100,000 (0 for those zones); 58 such zones.
  dist_coast  = distance to the nearest coastal_true zone (0 for coastal zones).
  Outcome of a (holiday h, ordinary o) pair: log ratio log((m_h + 1) / (m_o + 1)) for s_in, k_in, core, n_prov;
  log(eff_h / eff_o) for eff_prov; difference m_h - m_o for core, n_prov and share_far;
  (s_h - s_o) / population for s_in.
  Model: y ~ loss_pct + log(pop) + log(n_mun) + log(area) + coastal + log1p(dist_city) + log1p(dist_coast)
         | province FE; CRV1 clustered by province, t with G - 1 df; wild-cluster bootstrap (Webb weights,
         B = N_WCB) restricted p-value and unrestricted percentile-t 95% CI.

Gravity null: doubly constrained (origin and destination totals) with a non-parametric distance deterrence
f(bin) fitted jointly by 3-way IPF on the full weighted matrix of EACH night (diagonal = 0). The fitted mean
T_ij preserves every node's s_in exactly, so the null can only be tested on what s_in does not fix:
k_in, core, n_prov, eff_prov, share_far. Binary metrics under the null are computed from negative-binomial
draws around T_ij (variance T + alpha T^2, alpha by moments per night): analytic E[k_in] = sum_i P(W_ij >= 10)
for every night; coreness = mean over N_SPATIAL_DRAWS draws for the main nights.

Rewired null for coreness: N_REWIRE undirected degree-preserving rewirings of the w10 projection of each main night.

Run: .venv/bin/python -I pipeline/r1_villages.py [--recompute | --recompute=metrics,gravity,rewire]
Cached stage outputs (r1_node_metrics, r1_gravity_*, r1_rewired_*) are reused unless recomputed.
"""
from __future__ import annotations

import gzip
import random
import re
import sys
import time
import zlib
from pathlib import Path

import igraph as ig
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C  # noqa: E402
import b00_common as B  # noqa: E402

NET = C.RESULTS / "b01_networks"
OUT = C.RESULTS
DIAG = C.RESULTS / "r1_diag"
INE53001 = C.ROOT / "_local" / "external" / "ine_turismo_moviles" / "raw" / "turismo_53001_provorig_x_mundest.csv.gz"

BIG_CITY_POP = 100_000
FAR_KM = 100.0
N_WCB = 9999
DIST_BIN_EDGES_KM = [0, 5, 10, 15, 20, 25, 30, 40, 50, 60, 80, 100, 125, 150, 200, 250, 300, 400, 500, 600, 800, 1200]
IPF_TOL = 1e-8          # max relative error on every non-zero margin
IPF_MAXIT = 2000
W_MIN = C.BACKBONE_MAIN["w_min"]
INE_MIN_PUBLISHED = 30  # smallest published 53001 total (checked in data); suppressed cells are < 30
STREAMS = ("nb_draw", "rewire", "wcb")

CONTROLS = ["lpop", "lnmun", "larea", "coastal", "ldcity", "ldcoast"]


def seed_for(*parts) -> int:
    ints = [C.SEED] + [zlib.crc32(str(p).encode()) if not isinstance(p, (int, np.integer)) else int(p) for p in parts]
    return int(np.random.SeedSequence(ints).generate_state(1)[0])


# ---------------------------------------------------------------------------
# Nodes and controls
# ---------------------------------------------------------------------------
def load_nodes() -> pd.DataFrame:
    nd = pd.read_parquet(NET / "nodes.parquet")
    assert len(nd) == C.EXPECTED_N_NODES and (nd.idx.values == np.arange(len(nd))).all()
    xy = nd[["x", "y"]].values / 1000.0
    D = np.hypot(xy[:, None, 0] - xy[None, :, 0], xy[:, None, 1] - xy[None, :, 1])
    big = nd.population.values >= BIG_CITY_POP
    coast = nd.coastal_true.values.astype(bool)
    nd["dist_city_km"] = D[:, big].min(axis=1)
    nd["dist_coast_km"] = D[:, coast].min(axis=1)
    nd["lpop"] = np.log(nd.population)
    nd["lnmun"] = np.log(nd.n_municipalities)
    nd["larea"] = np.log(nd.area_km2)
    nd["coastal"] = coast.astype(float)
    nd["ldcity"] = np.log1p(nd.dist_city_km)
    nd["ldcoast"] = np.log1p(nd.dist_coast_km)
    nd["artefact"] = nd.id.isin(C.ARTEFACT_ZONES)
    assert nd.artefact.sum() == len(C.ARTEFACT_ZONES)
    ok = nd.loss_peak_2021.notna() & ~nd.artefact
    nd["loss_pct"] = np.nan
    nd.loc[ok, "loss_pct"] = nd.loc[ok, "loss_peak_2021"].rank(pct=True, method="average")
    nd["loss_q"] = np.nan
    nd.loc[ok, "loss_q"] = pd.qcut(nd.loc[ok, "loss_peak_2021"], 5, labels=False) + 1
    nd["m2_incomplete"] = nd.id.map(m2_incomplete_flags()).fillna(True)
    return nd, D


def m2_incomplete_flags() -> pd.Series:
    """Audit m2: True if any INE municipality of the zone lacks a census value in 1900-1960 or 2021.

    Uses the same two files and the same inner join as b00_common.historical_loss.
    """
    H = pd.read_csv(B.HIST_POP_CSV, dtype={"ine_code": str, "province_code": str})
    rel = pd.read_csv(B.INE_MITMA_REL, sep="|", dtype=str)[["municipio_ine", "municipio_mitma"]].drop_duplicates()
    H = H.merge(rel, left_on="ine_code", right_on="municipio_ine", how="inner")
    cols = B._PEAK_COLS + ["pop_2021"]
    H["miss"] = H[cols].isna().any(axis=1)
    return H.groupby("municipio_mitma")["miss"].any()


# ---------------------------------------------------------------------------
# Node metrics of one network
# ---------------------------------------------------------------------------
def node_metrics_from_matrix(W: np.ndarray, prov_idx: np.ndarray, n_prov: int, D: np.ndarray) -> dict:
    """Weighted node metrics from a dense N x N matrix (rows = residence, cols = sleep zone)."""
    s_in = W.sum(axis=0)
    P = np.zeros((n_prov, W.shape[1]))
    np.add.at(P, prov_idx, W)                       # origin province x destination
    own = np.zeros_like(P, dtype=bool)
    own[prov_idx, np.arange(W.shape[1])] = True
    n_other = ((P >= W_MIN) & ~own).sum(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        q = P / s_in[None, :]
        eff = 1.0 / (q ** 2).sum(axis=0)
        far = (W * (D > FAR_KM)).sum(axis=0) / s_in
    return {"s_in": s_in, "n_prov": n_other.astype(float), "eff_prov": eff, "share_far": far}


CACHE = C.RESULTS / "r1_cache"


def get_edges(date: str, kind: str) -> pd.DataFrame:
    """Edges of one night: the b01 file if it exists, otherwise built with the b00 loaders and cached.

    kind: full | w10 | equal_edges. equal_edges uses E = E_dir of the w10 backbone of the ordinary main
    night, read from results/b01_networks/manifest.csv (as b01 does).
    """
    f = NET / f"edges_{date}__{kind}.parquet"
    if f.exists():
        return pd.read_parquet(f)
    CACHE.mkdir(parents=True, exist_ok=True)
    g = CACHE / f"edges_{date}__{kind}.parquet"
    if g.exists():
        return pd.read_parquet(g)
    full_f = CACHE / f"edges_{date}__full.parquet"
    if full_f.exists():
        full = pd.read_parquet(full_f)
    else:
        nodes = pd.read_parquet(NET / "nodes.parquet")
        full, _ = B.build_edges(date, nodes)
        full.to_parquet(full_f, index=False)
    if kind == "full":
        return full
    if kind == "w10":
        e = B.backbone_threshold(full, W_MIN)
    elif kind == "equal_edges":
        man = pd.read_csv(NET / "manifest.csv", dtype={"date": str})
        E = int(man[(man.date == C.MAIN_NIGHTS["ordinary"]) & (man.backbone == "w10")].E.iloc[0])
        e = B.backbone_equal_edges(full, E)
    else:
        raise ValueError(kind)
    e.to_parquet(g, index=False)
    return e


def dense(e: pd.DataFrame, N: int) -> np.ndarray:
    W = np.zeros((N, N))
    W[e.i.values, e.j.values] = e.w.values
    return W


def undirected_coreness(i: np.ndarray, j: np.ndarray, N: int) -> np.ndarray:
    G = B.to_igraph(pd.DataFrame({"i": i, "j": j}), N, directed=False)
    return np.array(G.coreness(), dtype=float)


def in_coreness(i: np.ndarray, j: np.ndarray, N: int) -> np.ndarray:
    """In-coreness of the directed backbone (k-core on in-degree only): depends on who VISITS the zone."""
    G = ig.Graph(n=N, edges=list(zip(i.tolist(), j.tolist())), directed=True)
    return np.array(G.coreness(mode="in"), dtype=float)


def build_node_metrics(nd: pd.DataFrame, D: np.ndarray, dates: list[str]) -> pd.DataFrame:
    N = len(nd)
    provs = sorted(nd.prov.unique())
    pidx = nd.prov.map({p: k for k, p in enumerate(provs)}).values
    core = pd.read_parquet(C.RESULTS / "b02_coreness.parquet").set_index("id").loc[nd.id]
    rows = []
    for d in dates:
        W = dense(get_edges(d, "full"), N)
        m = node_metrics_from_matrix(W, pidx, len(provs), D)
        for bb in ("w10", "equal_edges"):
            e = get_edges(d, bb)
            m[f"k_in_{bb}"] = np.bincount(e.j.values, minlength=N).astype(float)
            m[f"core_{bb}"] = undirected_coreness(e.i.values, e.j.values, N)
            col = f"{d}__{bb}"
            if col in core.columns:   # consistency with b02 where b02 has the night
                assert np.array_equal(core[col].values.astype(float), m[f"core_{bb}"]), col
        e = get_edges(d, "w10")
        m["k_out_w10"] = np.bincount(e.i.values, minlength=N).astype(float)
        m["core_in_w10"] = in_coreness(e.i.values, e.j.values, N)
        s_out = W.sum(axis=1)
        m["s_out"] = s_out
        m["away_share"] = s_out / nd.population.values   # residents sleeping in another zone, per resident
        m["net_per_res"] = (m["s_in"] - s_out) / nd.population.values
        m["net_bal"] = (m["s_in"] - s_out) / (m["s_in"] + s_out)
        m["net_importer"] = (m["s_in"] > s_out).astype(float)
        m["share_of_away"] = m["s_in"] / W.sum()          # share of all away persons of the night (pattern #12)
        df = pd.DataFrame(m)
        df.insert(0, "id", nd.id.values)
        df.insert(0, "date", d)
        rows.append(df)
        print("metrics", d, flush=True)
    return pd.concat(rows, ignore_index=True)


# ---------------------------------------------------------------------------
# Gravity null
# ---------------------------------------------------------------------------
def dist_bins(D: np.ndarray) -> np.ndarray:
    b = np.digitize(D, DIST_BIN_EDGES_KM) - 1
    assert b.max() < len(DIST_BIN_EDGES_KM) - 1, "distance beyond last bin edge"
    np.fill_diagonal(b, -1)
    return b


def fit_gravity(W: np.ndarray, bins: np.ndarray) -> tuple[np.ndarray, dict]:
    """3-way IPF: T_ij = a_i b_j f(bin_ij), T_ii = 0, matching row sums, column sums and bin totals of W."""
    nb = bins.max() + 1
    off = bins >= 0
    O, Dd = W.sum(1), W.sum(0)
    Bt = np.bincount(bins[off], weights=W[off], minlength=nb)
    T = off.astype(float)
    T[:, Dd == 0] = 0.0
    T[O == 0, :] = 0.0
    errs = (np.inf, np.inf, np.inf)
    for it in range(1, IPF_MAXIT + 1):
        rs = T.sum(1)
        T *= np.divide(O, rs, out=np.zeros_like(O), where=rs > 0)[:, None]
        cs = T.sum(0)
        T *= np.divide(Dd, cs, out=np.zeros_like(Dd), where=cs > 0)[None, :]
        bt = np.bincount(bins[off], weights=T[off], minlength=nb)
        f = np.divide(Bt, bt, out=np.zeros_like(Bt), where=bt > 0)
        T[off] *= f[bins[off]]
        rs, cs = T.sum(1), T.sum(0)
        bt = np.bincount(bins[off], weights=T[off], minlength=nb)
        errs = (np.max(np.abs(rs - O)[O > 0] / O[O > 0]),
                np.max(np.abs(cs - Dd)[Dd > 0] / Dd[Dd > 0]),
                np.max(np.abs(bt - Bt)[Bt > 0] / Bt[Bt > 0]))
        if max(errs) < IPF_TOL:
            break
    assert max(errs) < IPF_TOL, f"IPF did not converge: {errs}"
    m = off & (O > 0)[:, None] & (Dd > 0)[None, :]
    alpha = float((((W[m] - T[m]) ** 2) - T[m]).sum() / (T[m] ** 2).sum())
    disp = float((((W[m] - T[m]) ** 2) / np.maximum(T[m], 1e-12)).sum() / m.sum())
    info = dict(ipf_iterations=it, max_rel_err_rows=errs[0], max_rel_err_cols=errs[1], max_rel_err_bins=errs[2],
                nb_alpha=alpha, pearson_dispersion=disp)
    return T, info


def nb_sf(T: np.ndarray, alpha: float, k: int) -> np.ndarray:
    """P(X >= k), X ~ NB(mean T, var T + alpha T^2)."""
    n = 1.0 / alpha
    return stats.nbinom.sf(k - 1, n, n / (n + T))


def nb_draw(T: np.ndarray, alpha: float, rng: np.random.Generator) -> np.ndarray:
    n = 1.0 / alpha
    lam = rng.gamma(shape=n, scale=T / n)
    return rng.poisson(lam)


def build_gravity(nd, D, dates, draw_dates) -> tuple[pd.DataFrame, pd.DataFrame]:
    N = len(nd)
    provs = sorted(nd.prov.unique())
    pidx = nd.prov.map({p: k for k, p in enumerate(provs)}).values
    bins = dist_bins(D)
    rows, fits = [], []
    for d in dates:
        t0 = time.time()
        W = dense(get_edges(d, "full"), N)
        T, info = fit_gravity(W, bins)
        m = node_metrics_from_matrix(T, pidx, len(provs), D)
        assert np.allclose(m["s_in"], W.sum(0), rtol=1e-6), "gravity must preserve s_in"
        sf = nb_sf(T, info["nb_alpha"], W_MIN)
        np.fill_diagonal(sf, 0.0)
        m["k_in_w10_nb"] = sf.sum(0)
        m["k_in_w10_poisson"] = stats.poisson.sf(W_MIN - 1, T).sum(0)
        ii, jj = np.nonzero(T >= W_MIN)
        m["k_in_w10_det"] = np.bincount(jj, minlength=N).astype(float)
        m["core_w10_det"] = undirected_coreness(ii, jj, N)
        info.update(date=d, E_obs_w10=int((W >= W_MIN).sum()), E_exp_nb=float(sf.sum()), E_det=int(len(ii)))
        if d in draw_dates:
            cs, ks, ci, es = np.zeros(N), np.zeros(N), np.zeros(N), []
            for rep in range(C.N_SPATIAL_DRAWS):
                rng = np.random.default_rng(seed_for("nb_draw", d, rep))
                X = nb_draw(T, info["nb_alpha"], rng)
                np.fill_diagonal(X, 0)
                ii, jj = np.nonzero(X >= W_MIN)
                cs += undirected_coreness(ii, jj, N)
                ci += in_coreness(ii, jj, N)
                ks += np.bincount(jj, minlength=N)
                es.append(len(ii))
            m["core_w10_nb"] = cs / C.N_SPATIAL_DRAWS
            m["core_in_w10_nb"] = ci / C.N_SPATIAL_DRAWS
            m["k_in_w10_nbdraw"] = ks / C.N_SPATIAL_DRAWS
            info.update(E_draw_mean=float(np.mean(es)), E_draw_sd=float(np.std(es)))
        df = pd.DataFrame(m)
        df.insert(0, "id", nd.id.values)
        df.insert(0, "date", d)
        rows.append(df)
        info["seconds"] = round(time.time() - t0, 1)
        fits.append(info)
        print("gravity", info, flush=True)
    return pd.concat(rows, ignore_index=True), pd.DataFrame(fits)


# ---------------------------------------------------------------------------
# Degree-preserving rewiring null for coreness
# ---------------------------------------------------------------------------
def build_rewired_coreness(nd, dates) -> tuple[pd.DataFrame, pd.DataFrame]:
    N = len(nd)
    rows, summ = [], []
    for d in dates:
        e = get_edges(d, "w10")
        G = B.to_igraph(e, N, directed=False)
        deg = np.array(G.degree(), float)
        core = np.array(G.coreness(), float)
        acc = np.zeros(N)
        acc2 = np.zeros(N)
        rho_rew = []
        for rep in range(C.N_REWIRE):
            ig.set_random_number_generator(random.Random(seed_for("rewire", d, rep)))
            H = G.copy()
            H.rewire(n=C.REWIRE_SWAPS_PER_EDGE * G.ecount(), allowed_edge_types="simple")
            assert np.array_equal(np.array(H.degree()), deg.astype(int))
            c = np.array(H.coreness(), float)
            acc += c
            acc2 += c ** 2
            rho_rew.append(stats.spearmanr(deg, c).statistic)
        mean = acc / C.N_REWIRE
        sd = np.sqrt(np.maximum(acc2 / C.N_REWIRE - mean ** 2, 0))
        rows.append(pd.DataFrame({"date": d, "id": nd.id.values, "degree": deg, "core": core,
                                  "core_rew_mean": mean, "core_rew_sd": sd}))
        # how much of coreness is a function of degree
        r2_real = r2_on_degree(core, deg)
        r2_rew = r2_on_degree(mean, deg)
        summ.append(dict(date=d, spearman_core_degree_real=stats.spearmanr(deg, core).statistic,
                         spearman_core_degree_rew_median=float(np.median(rho_rew)),
                         r2_core_on_degree_real=r2_real, r2_rewmean_on_degree=r2_rew,
                         corr_core_real_vs_rewmean=float(np.corrcoef(core, mean)[0, 1]),
                         mean_abs_core_minus_rewmean=float(np.mean(np.abs(core - mean))),
                         median_rew_sd=float(np.median(sd))))
        print("rewire", summ[-1], flush=True)
    return pd.concat(rows, ignore_index=True), pd.DataFrame(summ)


def r2_on_degree(c: np.ndarray, deg: np.ndarray) -> float:
    """R^2 of c on a flexible function of degree (mean of c within each distinct degree value)."""
    s = pd.Series(c).groupby(deg).transform("mean").values
    return float(1 - ((c - s) ** 2).sum() / ((c - c.mean()) ** 2).sum())


# ---------------------------------------------------------------------------
# Regression: FE-OLS, CRV1 by cluster, wild-cluster bootstrap
# ---------------------------------------------------------------------------
WEBB = np.array([-np.sqrt(1.5), -1.0, -np.sqrt(0.5), np.sqrt(0.5), 1.0, np.sqrt(1.5)])


def fe_ols(y, x, Z, fe, cl, n_boot=0, seed=0) -> dict:
    """OLS of y on x + Z + FE(fe), clustered by cl. FWL: everything residualised on [Z, FE dummies]."""
    y, x = np.asarray(y, float), np.asarray(x, float)
    Z = np.asarray(Z, float)
    keep_cols = [k for k in range(Z.shape[1]) if np.ptp(Z[:, k]) > 0]
    Z = Z[:, keep_cols]
    F = pd.get_dummies(pd.Series(fe)).values.astype(float)
    A = np.column_stack([Z, F])
    Q, R = np.linalg.qr(A)
    rank = np.sum(np.abs(np.diag(R)) > 1e-10 * np.abs(R).max())
    if rank < A.shape[1]:
        U, S, _ = np.linalg.svd(A, full_matrices=False)
        Q = U[:, S > 1e-10 * S.max()]
    mz = lambda v: v - Q @ (Q.T @ v)  # noqa: E731
    xt, yt = mz(x), mz(y)
    sxx = xt @ xt
    b = (xt @ yt) / sxx
    u = yt - b * xt
    n = len(y)
    codes, G = pd.factorize(pd.Series(cl))[0], pd.Series(cl).nunique()
    k = 1 + Z.shape[1]                      # FE nested in clusters: not counted (pyfixest/reghdfe convention)
    cadj = G / (G - 1) * (n - 1) / (n - k)

    def se_of(xu):  # xu: n x m array of x_i u_i
        S = np.zeros((G, xu.shape[1]))
        np.add.at(S, codes, xu)
        return np.sqrt(cadj * (S ** 2).sum(0)) / sxx

    se = float(se_of((xt * u)[:, None])[0])
    tcrit = stats.t.ppf(0.975, G - 1)
    t = b / se
    out = dict(coef=b, se=se, t=t, df=G - 1, p_t=float(2 * stats.t.sf(abs(t), G - 1)),
               ci_lo=b - tcrit * se, ci_hi=b + tcrit * se, n=n, G=G, k=k)
    if n_boot:
        rng = np.random.default_rng(seed)
        tr, tu = [], []
        u_r = yt                                   # restricted residuals (b = 0)
        for start in range(0, n_boot, 1000):
            m = min(1000, n_boot - start)
            v = WEBB[rng.integers(0, 6, size=(G, m))][codes]   # n x m
            # restricted (H0: b = 0)
            ys = mz(u_r[:, None] * v)
            bs = (xt @ ys) / sxx
            us = ys - xt[:, None] * bs[None, :]
            tr.append(bs / se_of(xt[:, None] * us))
            # unrestricted (percentile-t CI)
            es = mz(u[:, None] * v)
            db = (xt @ es) / sxx
            us = es - xt[:, None] * db[None, :]
            tu.append(db / se_of(xt[:, None] * us))
        tr, tu = np.concatenate(tr), np.concatenate(tu)
        out.update(p_wcr=float((np.sum(np.abs(tr) >= abs(t)) + 1) / (n_boot + 1)),
                   ci_wcu_lo=b - np.quantile(tu, 0.975) * se, ci_wcu_hi=b - np.quantile(tu, 0.025) * se,
                   n_boot=n_boot)
    return out


# ---------------------------------------------------------------------------
# Outcomes of a (holiday, ordinary) pair
# ---------------------------------------------------------------------------
OUTCOMES = {
    # name: (metric, transform)
    "s_in_logratio": ("s_in", "lr1"),
    "s_in_diff_per_res": ("s_in", "dpop"),
    "k_in_w10_logratio": ("k_in_w10", "lr1"),
    "k_in_w10_diff": ("k_in_w10", "diff"),
    "core_w10_diff": ("core_w10", "diff"),
    "core_w10_logratio": ("core_w10", "lr1"),
    "k_in_eq_logratio": ("k_in_equal_edges", "lr1"),
    "core_eq_diff": ("core_equal_edges", "diff"),
    "n_prov_diff": ("n_prov", "diff"),
    "n_prov_logratio": ("n_prov", "lr1"),
    "eff_prov_logratio": ("eff_prov", "lr"),
    "eff_prov_diff": ("eff_prov", "diff"),
    "share_far_diff": ("share_far", "diff"),
    "core_in_w10_diff": ("core_in_w10", "diff"),
    "core_in_w10_logratio": ("core_in_w10", "lr1"),
    "k_out_w10_logratio": ("k_out_w10", "lr1"),
    "away_share_logratio": ("away_share", "lr"),
    "net_per_res_diff": ("net_per_res", "diff"),
    "net_bal_diff": ("net_bal", "diff"),
    "net_importer_diff": ("net_importer", "diff"),
}
PRIMARY = ["s_in_logratio", "k_in_w10_logratio", "core_w10_diff", "n_prov_diff", "eff_prov_logratio", "share_far_diff",
           "net_per_res_diff", "net_bal_diff"]


def transform(mh, mo, kind, pop):
    if kind == "lr1":
        return np.log((mh + 1) / (mo + 1))
    if kind == "lr":
        return np.log(mh / mo)
    if kind == "diff":
        return mh - mo
    if kind == "dpop":
        return (mh - mo) / pop
    raise ValueError(kind)


def samples(nd) -> dict[str, np.ndarray]:
    ok = nd.loss_pct.notna().values
    return {
        "all": ok,
        "single_mun": ok & (nd.n_municipalities.values == 1),
        "inland": ok & ~nd.coastal_true.values.astype(bool),
        "m2_complete": ok & ~nd.m2_incomplete.values.astype(bool),
    }


def run_model(nd, y, mask, n_boot=0, seed=0):
    d = nd.loc[mask]
    yy = y[mask]
    good = np.isfinite(yy)
    d, yy = d[good], yy[good]
    return fe_ols(yy, d.loss_pct.values, d[CONTROLS].values, d.prov.values, d.prov.values, n_boot, seed)


# ---------------------------------------------------------------------------
# INE 53001, municipality level
# ---------------------------------------------------------------------------
def parse_es_num(s: str):
    s = s.strip().strip('"')
    if not s or s in (".", "..", "-"):
        return None
    return float(s.replace(".", "").replace(",", "."))


def load_ine53001(periods: set[str]) -> pd.DataFrame:
    """Total Nacional origin rows of INE 53001 (tourists arriving in each destination municipality per month).

    Same reading as explore/p89 part B: columns 0-2 = origin (Total Nacional, CCAA, province), column 3 =
    '<5-digit code> <name>', column 5 = period, column 6 = value in Spanish number format.
    Rows with an empty CCAA and province are the national-origin totals.
    """
    recs = []
    with gzip.open(INE53001, "rt", encoding="utf-8") as f:
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
    assert not df.duplicated(["ine_code", "per"]).any()
    return df.pivot(index="ine_code", columns="per", values="tourists")


INE53000 = C.ROOT / "_local" / "external" / "ine_turismo_moviles" / "raw" / "turismo_53000_munorig_x_provdest.csv.gz"


def load_ine53000_departures(periods: set[str]) -> pd.DataFrame:
    """INE 53000 rows '<code> <name>;Total Nacional;;;Turistas;<period>;<value>': tourists RESIDENT in the
    origin municipality travelling anywhere in Spain that month (the outbound counterpart of 53001)."""
    recs = []
    with gzip.open(INE53000, "rt", encoding="utf-8") as f:
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
    df = pd.DataFrame(recs, columns=["ine_code", "per", "tourists"])
    assert not df.duplicated(["ine_code", "per"]).any()
    return df.pivot(index="ine_code", columns="per", values="tourists")


def ine_check(nd: pd.DataFrame, D: np.ndarray) -> tuple[pd.DataFrame, pd.DataFrame]:
    pairs = [("2024M08", "2024M02"), ("2024M11", "2024M02"), ("2022M08", "2022M02"), ("2023M08", "2023M02"),
             ("2025M08", "2025M02"), ("2022M11", "2022M02"), ("2023M11", "2023M02")]
    periods = {p for pr in pairs for p in pr}
    T = load_ine53001(periods)
    assert T.min().min() >= INE_MIN_PUBLISHED, T.min().min()
    H = pd.read_csv(B.HIST_POP_CSV, dtype={"ine_code": str, "province_code": str})
    rel = pd.read_csv(B.INE_MITMA_REL, sep="|", dtype=str)[["municipio_ine", "municipio_mitma"]].drop_duplicates()
    # a municipality split across two zones (audit m2: 03903) is assigned to the zone holding most sections
    rel_full = pd.read_csv(B.INE_MITMA_REL, sep="|", dtype=str)
    zone_of = (rel_full.groupby(["municipio_ine", "municipio_mitma"]).size().reset_index(name="n")
               .sort_values(["municipio_ine", "n", "municipio_mitma"], ascending=[True, False, True])
               .drop_duplicates("municipio_ine").set_index("municipio_ine").municipio_mitma)
    n_split = int(rel.municipio_ine.duplicated().sum())
    H["zone"] = H.ine_code.map(zone_of)
    zid = nd.set_index("id")
    H = H[H.zone.isin(zid.index)].copy()                        # mainland zones of the node set
    H["peak"] = H[B._PEAK_COLS].max(axis=1)
    H["loss"] = 1 - H.pop_2021 / H.peak
    H["pop"] = H.padron_2024
    counts = {"mainland_municipalities_in_node_zones": len(H), "municipalities_split_across_zones": n_split}
    H = H[H["pop"].notna() & (H["pop"] > 0)]
    counts["with_padron_2024"] = len(H)
    counts["missing_loss_excluded"] = int(H.loss.isna().sum())
    H = H[H.loss.notna()].copy()
    counts["analysed"] = len(H)
    H["prov"] = H.zone.str[:2]
    H["loss_pct"] = H.loss.rank(pct=True)
    H["loss_q"] = pd.qcut(H.loss, 5, labels=[1, 2, 3, 4, 5]).astype(int)
    H["lpop"] = np.log(H["pop"])
    for c in ("coastal", "ldcity", "ldcoast", "is_am"):
        H[c] = H.zone.map(zid[c]).astype(float)
    H["zone_single"] = H.zone.map(zid.n_municipalities) == 1
    for p in periods:
        H[p] = H.ine_code.map(T[p]) if p in T else np.nan
    ctrl = ["lpop", "coastal", "ldcity", "ldcoast"]
    rows, qrows = [], []
    for h, o in pairs:
        both = H[h].notna() & H[o].notna()
        out = {
            "logratio_published": np.where(both, np.log(H[h] / H[o]), np.nan),
            "diff_per_res_lo": (H[h].fillna(0) - H[o].fillna(INE_MIN_PUBLISHED - 1)) / H["pop"],
            "diff_per_res_hi": (H[h].fillna(INE_MIN_PUBLISHED - 1) - H[o].fillna(0)) / H["pop"],
            "diff_per_res_zero": (H[h].fillna(0) - H[o].fillna(0)) / H["pop"],
        }
        subs = {"all": np.ones(len(H), bool), "in_AM_zone": H.is_am.values == 1,
                "single_mun_zone": H.zone_single.values, "inland": H.coastal.values == 0,
                "pop_lt_1000": H["pop"].values < 1000}
        for oname, yv in out.items():
            yv = np.asarray(yv, float)
            for sname, sm in subs.items():
                m = sm & np.isfinite(yv)
                d = H[m]
                nb = N_WCB if (h, o) == ("2024M08", "2024M02") and oname in ("logratio_published", "diff_per_res_zero") else 0
                r = fe_ols(yv[m], d.loss_pct.values, d[ctrl].values, d.prov.values, d.prov.values, nb,
                           seed_for("wcb", "ine", h, o, oname, sname))
                r.update(holiday=h, ordinary=o, outcome=oname, sample=sname,
                         n_published_both=int((both & sm).sum()))
                rows.append(r)
            if oname in ("logratio_published", "diff_per_res_zero"):
                g = pd.DataFrame({"q": H.loss_q.values, "y": yv, "pop": H["pop"].values, "am": H.is_am.values})
                for q, gg in g.groupby("q"):
                    qrows.append(dict(holiday=h, ordinary=o, outcome=oname, loss_q=q, n=int(gg.y.notna().sum()),
                                      median_y=float(gg.y.median()), median_pop=float(gg["pop"].median()),
                                      share_in_AM=float(gg.am.mean())))
    res = pd.DataFrame(rows)
    # levels (published cells) and ANCOVA, main year
    extra = []
    for per in ("2024M02", "2024M08", "2024M11"):
        yv = np.log(H[per].values / H["pop"].values)
        for sname, sm in {"all": np.ones(len(H), bool), "in_AM_zone": H.is_am.values == 1,
                          "single_mun_zone": H.zone_single.values}.items():
            m = sm & np.isfinite(yv)
            d = H[m]
            r = fe_ols(yv[m], d.loss_pct.values, d[ctrl].values, d.prov.values, d.prov.values, N_WCB,
                       seed_for("wcb", "ine_level", per, sname))
            r.update(holiday=per, ordinary="", outcome="level_log_per_res_published", sample=sname,
                     n_published_both=int(m.sum()))
            extra.append(r)
    for h in ("2024M08", "2024M11"):
        yh, yo = np.log(H[h].values / H["pop"].values), np.log(H["2024M02"].values / H["pop"].values)
        for sname, sm in {"all": np.ones(len(H), bool), "in_AM_zone": H.is_am.values == 1,
                          "single_mun_zone": H.zone_single.values}.items():
            m = sm & np.isfinite(yh) & np.isfinite(yo)
            d = H[m]
            r = fe_ols(yh[m], d.loss_pct.values, np.column_stack([d[ctrl].values, yo[m]]), d.prov.values,
                       d.prov.values, N_WCB, seed_for("wcb", "ine_ancova", h, sname))
            r.update(holiday=h, ordinary="2024M02", outcome="ancova_log_per_res_published", sample=sname,
                     n_published_both=int(m.sum()))
            extra.append(r)
    # net balance (arrivals 53001 - departures 53000) per resident, published both sides
    Dp = load_ine53000_departures({"2024M02", "2024M08", "2024M11", "2022M02", "2022M08", "2023M02", "2023M08"})
    counts["min_published_53000"] = float(Dp.min().min())
    for per in Dp.columns:
        H[f"dep_{per}"] = H.ine_code.map(Dp[per])
    net = {}
    for per in Dp.columns:
        if per in H:
            net[per] = (H[per] - H[f"dep_{per}"]) / H["pop"]
    subs3 = {"all": np.ones(len(H), bool), "in_AM_zone": H.is_am.values == 1, "single_mun_zone": H.zone_single.values,
             "inland": H.coastal.values == 0}
    for per, yv in net.items():
        yv = yv.values
        for sname, sm in subs3.items():
            m = sm & np.isfinite(yv)
            d = H[m]
            r = fe_ols(yv[m], d.loss_pct.values, d[ctrl].values, d.prov.values, d.prov.values,
                       N_WCB if per in ("2024M08", "2024M11", "2024M02") else 0, seed_for("wcb", "ine_net", per, sname))
            r.update(holiday=per, ordinary="", outcome="net_per_res_level_published", sample=sname,
                     n_published_both=int(m.sum()), share_net_importer=float((yv[m] > 0).mean()))
            extra.append(r)
            if sname == "all":
                g = pd.DataFrame({"q": H.loss_q.values[m], "y": yv[m]})
                for q, gg in g.groupby("q"):
                    qrows.append(dict(holiday=per, ordinary="", outcome="net_per_res_level_published", loss_q=q,
                                      n=len(gg), median_y=float(gg.y.median()),
                                      share_net_importer=float((gg.y > 0).mean())))
    for h, o in (("2024M08", "2024M02"), ("2024M11", "2024M02"), ("2022M08", "2022M02"), ("2023M08", "2023M02")):
        yv = (net[h] - net[o]).values
        for sname, sm in subs3.items():
            m = sm & np.isfinite(yv)
            d = H[m]
            r = fe_ols(yv[m], d.loss_pct.values, d[ctrl].values, d.prov.values, d.prov.values,
                       N_WCB if h.startswith("2024") else 0, seed_for("wcb", "ine_netdiff", h, sname))
            r.update(holiday=h, ordinary=o, outcome="net_per_res_diff_published", sample=sname,
                     n_published_both=int(m.sum()))
            extra.append(r)
    res = pd.concat([res, pd.DataFrame(extra)], ignore_index=True)
    res.attrs["counts"] = counts
    return res, pd.DataFrame(qrows), counts


# ---------------------------------------------------------------------------
# Volume curve (pattern #3)
# ---------------------------------------------------------------------------
P_ZONE_DAILY = C.RESULTS / "p_zone_daily.parquet"   # pattern agent: per-zone in/out away persons, 306 nights


def curve_fit_nights(dates: list[str]) -> list[str]:
    """2024 nights with no national holiday, long weekend or holiday season (fit set of the curve).

    Feb; Mar 1-22 (before Holy Week); Oct except 11-13 (National Day weekend); Nov 4-30; Dec 1-20 except 5-9
    (Constitution/Immaculate). Weekends are INCLUDED (they supply the high-volume ordinary nights).
    """
    keep = []
    for d in dates:
        t = pd.Timestamp(d)
        if t.year != 2024:
            continue
        m, dd = t.month, t.day
        if (m == 2 or (m == 3 and dd <= 22) or (m == 10 and dd not in (11, 12, 13)) or (m == 11 and dd >= 4)
                or (m == 12 and dd <= 20 and dd not in (5, 6, 7, 8, 9))):
            keep.append(d)
    return keep


def vol_curve(nd, smp, NM) -> None:
    Z = pd.read_parquet(P_ZONE_DAILY)
    Z["date"] = Z.date.astype(str)
    # cross-check against this script's node metrics on every overlapping night
    chk = Z.merge(NM[["date", "id", "s_in", "s_out"]], on=["date", "id"])
    assert len(chk) > 0
    assert np.allclose(chk.in_away_ml, chk.s_in, rtol=1e-9, atol=1e-6), "p_zone_daily in_away_ml != r1 s_in"
    assert np.allclose(chk.out_away_ml, chk.s_out, rtol=1e-9, atol=1e-6), "p_zone_daily out_away_ml != r1 s_out"
    pos = nd.set_index("id").idx
    rows = []
    for d, g in Z.groupby("date"):
        g = g.set_index("id").reindex(nd.id)
        s_in, s_out = g.in_away_ml.values, g.out_away_ml.values
        V = np.nansum(s_in)
        outs = {"log_share_of_away": np.log(s_in / V), "net_bal": (s_in - s_out) / (s_in + s_out),
                "net_per_res": (s_in - s_out) / nd.population.values}
        for oname, y in outs.items():
            for sname in ("all", "single_mun"):
                r = run_model(nd, y, smp[sname])
                rows.append(dict(date=d, away_persons=V, outcome=oname, sample=sname, coef=r["coef"],
                                 ci_lo=r["ci_lo"], ci_hi=r["ci_hi"], n=r["n"]))
    T = pd.DataFrame(rows)
    fitset = set(curve_fit_nights(sorted(T.date.unique())))
    T["in_fit_set"] = T.date.isin(fitset)
    out = []
    for (oname, sname), d in T.groupby(["outcome", "sample"]):
        f = d[d.in_fit_set]
        X = np.column_stack([np.ones(len(f)), np.log(f.away_persons)])
        b = np.linalg.lstsq(X, f.coef.values, rcond=None)[0]
        res_f = f.coef.values - X @ b
        sd = res_f.std(ddof=2)
        r2 = 1 - (res_f ** 2).sum() / ((f.coef - f.coef.mean()) ** 2).sum()
        vmax = f.away_persons.max()
        dd = d.assign(pred=b[0] + b[1] * np.log(d.away_persons))
        dd["resid_sd"] = (dd.coef - dd.pred) / sd
        dd["extrapolated"] = dd.away_persons > vmax
        dd["fit_r2"], dd["fit_slope"], dd["fit_resid_sd"], dd["fit_n"] = r2, b[1], sd, len(f)
        out.append(dd)
    T = pd.concat(out)
    T.to_csv(OUT / "r1_volume_curve_nights.csv", index=False)
    labels = {"summer_sisters": C.SISTER_NIGHTS["summer"], "allsaints_sisters": C.SISTER_NIGHTS["allsaints"],
              "summer_years": C.YEAR_REPLICATION["summer"], "allsaints_years": C.YEAR_REPLICATION["allsaints"],
              **{k: [v] for k, v in C.APPENDIX_NIGHTS.items()},
              "ordinary_main": [C.MAIN_NIGHTS["ordinary"]], "summer_main": [C.MAIN_NIGHTS["summer"]],
              "allsaints_main": [C.MAIN_NIGHTS["allsaints"]]}
    srows = []
    for (oname, sname), d in T.groupby(["outcome", "sample"]):
        for lab, ds in labels.items():
            x = d[d.date.isin(ds)]
            if x.empty:
                continue
            srows.append(dict(outcome=oname, sample=sname, nights=lab, n_nights=len(x),
                              fit_r2=x.fit_r2.iloc[0], fit_n=int(x.fit_n.iloc[0]), fit_slope=x.fit_slope.iloc[0],
                              coef_median=x.coef.median(), pred_median=x.pred.median(),
                              resid_sd_median=x.resid_sd.median(), resid_sd_min=x.resid_sd.min(),
                              resid_sd_max=x.resid_sd.max(), share_extrapolated=x.extrapolated.mean(),
                              away_persons_median=x.away_persons.median()))
    pd.DataFrame(srows).to_csv(OUT / "r1_volume_curve_summary.csv", index=False)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main(recompute: set[str] = frozenset()) -> None:
    DIAG.mkdir(parents=True, exist_ok=True)
    nd, D = load_nodes()
    nets = pd.read_csv(NET / "networks.csv", dtype=str)
    ord_sis = sorted(nets[(nets.group == "ordinary") & (nets.role == "sister")].date)
    main_nights = C.MAIN_NIGHTS
    sis = {g: (ord_sis if g == "ordinary" else sorted(C.SISTER_NIGHTS[g])) for g in main_nights}
    yrs = {g: sorted(C.YEAR_REPLICATION[g]) for g in main_nights}
    app = dict(C.APPENDIX_NIGHTS)
    dates = sorted(set(main_nights.values()) | {d for v in sis.values() for d in v} |
                   {d for v in yrs.values() for d in v} | set(app.values()))

    # sample bookkeeping
    smp = samples(nd)
    srows = [dict(sample=k, n=int(v.sum()), n_provinces=int(nd.loc[v, "prov"].nunique())) for k, v in smp.items()]
    miss = nd.loss_peak_2021.isna()
    srows.append(dict(sample="excluded_missing_loss", n=int(miss.sum()), n_provinces=int(nd.loc[miss, "prov"].nunique())))
    srows.append(dict(sample="excluded_missing_loss_single_mun", n=int((miss & (nd.n_municipalities == 1)).sum()),
                      n_provinces=np.nan))
    srows.append(dict(sample="excluded_artefact_zones", n=int(nd.artefact.sum()), n_provinces=np.nan))
    srows.append(dict(sample="big_city_zones_pop_ge_100k", n=int((nd.population >= BIG_CITY_POP).sum()), n_provinces=np.nan))
    pd.DataFrame(srows).to_csv(OUT / "r1_sample.csv", index=False)

    f_nm = OUT / "r1_node_metrics.parquet"
    if "metrics" in recompute or not f_nm.exists():
        build_node_metrics(nd, D, dates).to_parquet(f_nm, index=False)
    NM = pd.read_parquet(f_nm)

    f_gr, f_gf = OUT / "r1_gravity_node_metrics.parquet", OUT / "r1_gravity_fit.csv"
    if "gravity" in recompute or not f_gr.exists():
        g, gf = build_gravity(nd, D, dates, draw_dates=set(main_nights.values()))
        g.to_parquet(f_gr, index=False)
        gf.to_csv(f_gf, index=False)
    GM = pd.read_parquet(f_gr)

    f_rw, f_rs = OUT / "r1_rewired_coreness_nodes.parquet", OUT / "r1_rewired_coreness_summary.csv"
    if "rewire" in recompute or not f_rw.exists():
        rw, rs = build_rewired_coreness(nd, list(main_nights.values()))
        rw.to_parquet(f_rw, index=False)
        rs.to_csv(f_rs, index=False)
    RW = pd.read_parquet(f_rw)

    def metric(df, d, col):
        s = df[df.date == d].set_index("id")[col]
        return s.reindex(nd.id).values

    pop = nd.population.values

    # ---- pair list
    pairs = []
    o_main = main_nights["ordinary"]
    for g in ("summer", "allsaints"):
        pairs.append((g, "main", main_nights[g], o_main))
        for h in sis[g]:
            for o in sis["ordinary"]:
                pairs.append((g, "sister", h, o))
    for g in ("summer", "allsaints"):
        for h in yrs[g]:
            o = [x for x in yrs["ordinary"] if x[:4] == h[:4]][0]
            pairs.append((g, "year", h, o))
    for o in sis["ordinary"]:
        if o != o_main:
            pairs.append(("placebo", "ordinary_vs_ordinary", o, o_main))
    for g, h in app.items():
        pairs.append((g, "main", h, o_main))   # appendix nights: one night each, vs the ordinary main night

    # ---- observed models
    rows = []
    for g, role, h, o in pairs:
        for oname, (mcol, kind) in OUTCOMES.items():
            y = transform(metric(NM, h, mcol), metric(NM, o, mcol), kind, pop)
            for sname, mask in smp.items():
                if role != "main" and sname not in ("all", "single_mun"):
                    continue
                nb = N_WCB if role == "main" else 0
                r = run_model(nd, y, mask, nb, seed_for("wcb", h, o, oname, sname))
                r.update(group=g, role=role, holiday=h, ordinary=o, outcome=oname, sample=sname)
                rows.append(r)
    M = pd.DataFrame(rows)
    M.to_csv(OUT / "r1_models_all_pairs.csv", index=False)

    # ---- summary across sister nights and years
    summ = []
    for (g, oname, sname), d in M[M.group != "placebo"].groupby(["group", "outcome", "sample"]):
        mn = d[d.role == "main"].iloc[0]
        s = d[d.role == "sister"]
        yv = d[d.role == "year"]
        sg = np.sign(mn.coef)
        row = dict(group=g, outcome=oname, sample=sname, n=mn.n, G=mn.G, coef_main=mn.coef, se_main=mn.se,
                   ci_t_lo=mn.ci_lo, ci_t_hi=mn.ci_hi, p_t=mn.p_t, p_wcr=mn.p_wcr,
                   ci_wcu_lo=mn.ci_wcu_lo, ci_wcu_hi=mn.ci_wcu_hi)
        if len(s):
            row.update(sister_n_pairs=len(s), sister_median=s.coef.median(), sister_min=s.coef.min(),
                       sister_max=s.coef.max(), sister_same_sign=int((np.sign(s.coef) == sg).sum()),
                       sister_ci_excl0_same_sign=int(((np.sign(s.coef) == sg) & (np.sign(s.ci_lo) == np.sign(s.ci_hi))).sum()))
        for _, yr in yv.iterrows():
            row[f"year_{yr.holiday[:4]}"] = yr.coef
            row[f"year_{yr.holiday[:4]}_ci_lo"] = yr.ci_lo
            row[f"year_{yr.holiday[:4]}_ci_hi"] = yr.ci_hi
        pl = M[(M.group == "placebo") & (M.outcome == oname) & (M["sample"] == sname)]
        if len(pl):
            row.update(placebo_n=len(pl), placebo_min=pl.coef.min(), placebo_max=pl.coef.max())
        summ.append(row)
    S = pd.DataFrame(summ)
    S.to_csv(OUT / "r1_models_summary.csv", index=False)

    # ---- gravity null: observed vs expected gains
    GCOL = {  # observed metric -> expected metric (same definition under the null)
        "k_in_w10": "k_in_w10_nb", "core_w10": "core_w10_det", "n_prov": "n_prov", "eff_prov": "eff_prov",
        "share_far": "share_far", "core_in_w10": "core_in_w10_nb",
    }
    GOUT = {"k_in_w10_logratio": ("k_in_w10", "lr1"), "core_w10_diff": ("core_w10", "diff"),
            "n_prov_diff": ("n_prov", "diff"), "eff_prov_logratio": ("eff_prov", "lr"),
            "share_far_diff": ("share_far", "diff"), "core_in_w10_diff": ("core_in_w10", "diff")}
    grows, gq = [], []
    gpairs = [p for p in pairs if p[1] in ("main", "sister", "year")]
    for g, role, h, o in gpairs:
        for oname, (mcol, kind) in GOUT.items():
            variants = [GCOL[mcol]]
            drawn = role == "main" and g in ("summer", "allsaints")
            if drawn and mcol == "core_w10":
                variants.append("core_w10_nb")
            if mcol == "core_in_w10" and not drawn:
                continue
            if drawn and mcol == "k_in_w10":
                variants += ["k_in_w10_nbdraw", "k_in_w10_det"]
            yo = transform(metric(NM, h, mcol), metric(NM, o, mcol), kind, pop)
            for ecol in variants:
                ye = transform(metric(GM, h, ecol), metric(GM, o, ecol), kind, pop)
                for sname in ("all", "single_mun"):
                    mask = smp[sname]
                    for which, y in (("expected", ye), ("excess", yo - ye)):
                        nb = N_WCB if role == "main" and which == "excess" else 0
                        r = run_model(nd, y, mask, nb, seed_for("wcb", "grav", h, o, oname, ecol, sname, which))
                        r.update(group=g, role=role, holiday=h, ordinary=o, outcome=oname, null_metric=ecol,
                                 sample=sname, which=which)
                        grows.append(r)
                    if drawn:
                        q = pd.DataFrame({"q": nd.loss_q.values, "obs": yo, "exp": ye, "exc": yo - ye})[mask]
                        for qq, d in q.groupby("q"):
                            gq.append(dict(group=g, outcome=oname, null_metric=ecol, sample=sname, loss_q=int(qq),
                                           n=len(d), median_obs_gain=d.obs.median(), median_exp_gain=d.exp.median(),
                                           median_excess=d.exc.median(), mean_obs_gain=d.obs.mean(),
                                           mean_exp_gain=d.exp.mean()))
    GR = pd.DataFrame(grows)
    GR.to_csv(OUT / "r1_gravity_models_all_pairs.csv", index=False)
    pd.DataFrame(gq).to_csv(OUT / "r1_gravity_quintiles.csv", index=False)
    gs = []
    for (g, oname, ecol, sname, which), d in GR.groupby(["group", "outcome", "null_metric", "sample", "which"]):
        mn = d[d.role == "main"]
        if mn.empty:
            continue
        mn = mn.iloc[0]
        s = d[d.role == "sister"]
        row = dict(group=g, outcome=oname, null_metric=ecol, sample=sname, which=which, n=mn.n,
                   coef_main=mn.coef, ci_t_lo=mn.ci_lo, ci_t_hi=mn.ci_hi, p_wcr=mn.get("p_wcr", np.nan),
                   ci_wcu_lo=mn.get("ci_wcu_lo", np.nan), ci_wcu_hi=mn.get("ci_wcu_hi", np.nan))
        if len(s):
            row.update(sister_median=s.coef.median(), sister_min=s.coef.min(), sister_max=s.coef.max(),
                       sister_same_sign=int((np.sign(s.coef) == np.sign(mn.coef)).sum()), sister_n_pairs=len(s))
        for _, yr in d[d.role == "year"].iterrows():
            row[f"year_{yr.holiday[:4]}"] = yr.coef
        gs.append(row)
    pd.DataFrame(gs).to_csv(OUT / "r1_gravity_summary.csv", index=False)

    # ---- rewired null for coreness: excess coreness gain
    rrows = []
    for g in ("summer", "allsaints"):
        h = main_nights[g]
        exc = lambda d: (metric(RW, d, "core") - metric(RW, d, "core_rew_mean"))  # noqa: E731
        y_exc = exc(h) - exc(o_main)
        y_rew = metric(RW, h, "core_rew_mean") - metric(RW, o_main, "core_rew_mean")
        for sname in ("all", "single_mun"):
            for which, y in (("rewired_expected_gain", y_rew), ("excess_over_rewired", y_exc)):
                r = run_model(nd, y, smp[sname], N_WCB, seed_for("wcb", "rew", h, which, sname))
                r.update(group=g, holiday=h, ordinary=o_main, outcome="core_w10_diff", sample=sname, which=which)
                rrows.append(r)
    pd.DataFrame(rrows).to_csv(OUT / "r1_rewired_coreness_models.csv", index=False)

    # ---- quintile description of the observed gains (main pairs)
    qrows = []
    for g in ("summer", "allsaints"):
        h = main_nights[g]
        for oname in PRIMARY:
            mcol, kind = OUTCOMES[oname]
            y = transform(metric(NM, h, mcol), metric(NM, o_main, mcol), kind, pop)
            for sname in ("all", "single_mun"):
                mask = smp[sname]
                d = nd[mask].assign(y=y[mask])
                for q, gg in d.groupby("loss_q"):
                    qrows.append(dict(group=g, outcome=oname, sample=sname, loss_q=int(q), n=len(gg),
                                      median_loss=gg.loss_peak_2021.median(), median_pop=gg.population.median(),
                                      median_n_mun=gg.n_municipalities.median(), share_AM=gg.is_am.mean(),
                                      share_coastal=gg.coastal_true.mean(), median_y=gg.y.median(),
                                      q25_y=gg.y.quantile(.25), q75_y=gg.y.quantile(.75)))
                # adjusted Q5 - Q1 and every quintile vs Q1
                dd = d.copy()
                for q in (2, 3, 4, 5):
                    dd[f"q{q}"] = (dd.loss_q == q).astype(float)
                for q in (2, 3, 4, 5):
                    others = [f"q{k}" for k in (2, 3, 4, 5) if k != q]
                    r = fe_ols(dd.y.values, dd[f"q{q}"].values, dd[CONTROLS + others].values, dd.prov.values,
                               dd.prov.values)
                    qrows.append(dict(group=g, outcome=oname, sample=sname, loss_q=q, n=r["n"],
                                      adj_vs_Q1=r["coef"], adj_ci_lo=r["ci_lo"], adj_ci_hi=r["ci_hi"]))
    pd.DataFrame(qrows).to_csv(OUT / "r1_quintiles.csv", index=False)

    # ---- levels per night and ANCOVA (holiday level controlling the ordinary level of the same metric)
    LEVELS = {"s_in": np.log1p, "k_in_w10": np.log1p, "core_w10": lambda v: v, "core_in_w10": lambda v: v,
              "eff_prov": np.log, "n_prov": lambda v: v, "share_far": lambda v: v, "k_out_w10": np.log1p,
              "away_share": np.log, "net_per_res": lambda v: v, "net_bal": lambda v: v,
              "net_importer": lambda v: v}
    lrows = []
    level_dates = [(g, "main", d) for g, d in main_nights.items()] + \
                  [(g, "sister", d) for g in main_nights for d in sis[g]] + \
                  [(g, "year", d) for g in main_nights for d in yrs[g]] + \
                  [(g, "main", d) for g, d in app.items()]
    for g, role, d in level_dates:
        for mcol, f in LEVELS.items():
            y = f(metric(NM, d, mcol))
            for sname in ("all", "single_mun"):
                nb = N_WCB if role == "main" else 0
                r = run_model(nd, y, smp[sname], nb, seed_for("wcb", "level", d, mcol, sname))
                r.update(spec="level", group=g, role=role, holiday=d, ordinary="", metric=mcol, sample=sname,
                         median_metric=float(np.nanmedian(metric(NM, d, mcol)[smp[sname]])))
                lrows.append(r)
    for g, role, h, o in pairs:
        if g == "placebo":
            continue
        for mcol, f in LEVELS.items():
            yh, yo = f(metric(NM, h, mcol)), f(metric(NM, o, mcol))
            for sname in ("all", "single_mun"):
                m = smp[sname] & np.isfinite(yh) & np.isfinite(yo)
                dd = nd[m]
                nb = N_WCB if role == "main" else 0
                r = fe_ols(yh[m], dd.loss_pct.values, np.column_stack([dd[CONTROLS].values, yo[m]]),
                           dd.prov.values, dd.prov.values, nb, seed_for("wcb", "ancova", h, o, mcol, sname))
                r.update(spec="ancova", group=g, role=role, holiday=h, ordinary=o, metric=mcol, sample=sname)
                lrows.append(r)
    L = pd.DataFrame(lrows)
    L.to_csv(OUT / "r1_models_levels_ancova_all.csv", index=False)
    ls = []
    for (spec, g, mcol, sname), d in L.groupby(["spec", "group", "metric", "sample"]):
        mn = d[d.role == "main"].iloc[0]
        s = d[d.role == "sister"]
        row = dict(spec=spec, group=g, metric=mcol, sample=sname, n=mn.n, coef_main=mn.coef,
                   ci_wcu_lo=mn.ci_wcu_lo, ci_wcu_hi=mn.ci_wcu_hi, p_wcr=mn.p_wcr, ci_t_lo=mn.ci_lo, ci_t_hi=mn.ci_hi,
                   sister_median=s.coef.median(), sister_min=s.coef.min(), sister_max=s.coef.max(),
                   sister_same_sign=int((np.sign(s.coef) == np.sign(mn.coef)).sum()), sister_n=len(s))
        for _, yr in d[d.role == "year"].iterrows():
            row[f"year_{yr.holiday[:4]}"] = yr.coef
        ls.append(row)
    pd.DataFrame(ls).to_csv(OUT / "r1_models_levels_ancova_summary.csv", index=False)

    # ---- net balance description: share of net importers by loss quintile and for "emptied inland"
    nrows = []
    emptied_inland = (nd.loss_peak_2021 >= 0.5) & ~nd.coastal_true & ~nd.artefact
    for g, role, d in level_dates:
        v = metric(NM, d, "net_per_res")
        imp = metric(NM, d, "net_importer")
        nb = metric(NM, d, "net_bal")
        for sname in ("all", "single_mun"):
            mask = smp[sname]
            for q in (1, 2, 3, 4, 5):
                mm = mask & (nd.loss_q.values == q)
                nrows.append(dict(group=g, role=role, date=d, sample=sname, cls=f"Q{q}", n=int(mm.sum()),
                                  share_net_importer=float(imp[mm].mean()), median_net_per_res=float(np.median(v[mm])),
                                  median_net_bal=float(np.median(nb[mm]))))
            mm = mask & emptied_inland.values
            nrows.append(dict(group=g, role=role, date=d, sample=sname, cls="emptied_inland_loss_ge_0.5", n=int(mm.sum()),
                              share_net_importer=float(imp[mm].mean()), median_net_per_res=float(np.median(v[mm])),
                              median_net_bal=float(np.median(nb[mm]))))
    pd.DataFrame(nrows).to_csv(OUT / "r1_net_balance_by_quintile.csv", index=False)

    # ---- volume curve (pattern #3): is the loss gradient a holiday effect or a "more people away" effect?
    vol_curve(nd, smp, NM)

    # ---- INE 53001
    ine, ineq, counts = ine_check(nd, D)
    ine.to_csv(OUT / "r1_ine53001_models.csv", index=False)
    ineq.to_csv(OUT / "r1_ine53001_quintiles.csv", index=False)
    pd.Series(counts).to_csv(OUT / "r1_ine53001_sample.csv", header=["n"])

    # ---- cross-check of fe_ols against pyfixest on one main model
    try:
        import pyfixest as pf
        y = transform(metric(NM, main_nights["summer"], "k_in_w10"), metric(NM, o_main, "k_in_w10"), "lr1", pop)
        d = nd[smp["all"]].assign(y=y[smp["all"]])
        fit = pf.feols("y ~ loss_pct + lpop + lnmun + larea + coastal + ldcity + ldcoast | prov", data=d,
                       vcov={"CRV1": "prov"})
        tid = fit.tidy().loc["loss_pct"]
        mine = run_model(nd, y, smp["all"])
        chk = dict(pyfixest_coef=tid["Estimate"], pyfixest_se=tid["Std. Error"], mine_coef=mine["coef"], mine_se=mine["se"])
        pd.Series(chk).to_csv(DIAG / "check_fe_ols_vs_pyfixest.csv", header=["value"])
        assert abs(chk["pyfixest_coef"] - chk["mine_coef"]) < 1e-8, chk
        print("pyfixest check", chk)
    except ImportError:
        pass

    diagnostics(nd, NM, GM, main_nights, smp, pop)
    print("done")


def diagnostics(nd, NM, GM, main_nights, smp, pop):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    def metric(df, d, col):
        return df[df.date == d].set_index("id")[col].reindex(nd.id).values

    o = main_nights["ordinary"]
    fig, axes = plt.subplots(2, 5, figsize=(20, 7))
    specs = [("s_in", "lr1", None), ("k_in_w10", "lr1", "k_in_w10_nb"), ("core_w10", "diff", "core_w10_nb"),
             ("n_prov", "diff", "n_prov"), ("share_far", "diff", "share_far")]
    for r, g in enumerate(("summer", "allsaints")):
        h = main_nights[g]
        for c, (m, kind, e) in enumerate(specs):
            ax = axes[r, c]
            for sname, col in (("all", "k"), ("single_mun", "tab:blue")):
                msk = smp[sname]
                y = transform(metric(NM, h, m), metric(NM, o, m), kind, pop)
                b = pd.Series(y[msk]).groupby(pd.qcut(nd.loss_pct[msk].values, 10, labels=False)).median()
                ax.plot(b.index, b.values, "o-", color=col, label=f"obs {sname}")
                if e is not None:
                    ye = transform(metric(GM, h, e), metric(GM, o, e), kind, pop)
                    be = pd.Series(ye[msk]).groupby(pd.qcut(nd.loss_pct[msk].values, 10, labels=False)).median()
                    ax.plot(be.index, be.values, "x--", color=col, label=f"gravity {sname}")
            ax.set_title(f"{g}: {m} ({kind})", fontsize=9)
            ax.set_xlabel("loss decile (0 grew, 9 emptied)")
            if r == 0 and c == 0:
                ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(DIAG / "binned_gain_by_loss_decile.png", dpi=110)
    plt.close(fig)


if __name__ == "__main__":
    arg = [a for a in sys.argv[1:] if a.startswith("--recompute")]
    stages = {"metrics", "gravity", "rewire"} if arg == ["--recompute"] else (
        set(arg[0].split("=", 1)[1].split(",")) if arg else set())
    main(recompute=stages)
