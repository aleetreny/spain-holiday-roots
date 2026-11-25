"""b02_basic_properties: basic network properties, null models, Zhou methods, and representativeness.

Input: results/b01_networks/ (run b01_build_networks.py first).

Graphs per (network, backbone)
  Gd : directed binary graph of the backbone (edge i->j kept by the backbone).
  Gu : undirected simple projection of Gd (edge {i, j} if i->j or j->i).
Directed metrics on Gd: E, density, giant strong component, in/out degree stats,
  reciprocity (share of directed edges whose reverse exists), r_out_in = Pearson
  correlation over directed edges of k_out(source) and k_in(target).
Undirected metrics on Gu: E_und, density_und, giant weak (= connected) component,
  degree stats, average shortest path length L and diameter on the giant component
  (undirected, unweighted), average local clustering C_avg (nodes with k < 2 count as 0)
  and global transitivity C_global, degree assortativity r_und, rich-club
  coefficient phi(k) = 2 E_k / (N_k (N_k - 1)) over nodes with degree > k, rich club by rank
  (top-n nodes by degree), rich-club clique size n_clique, Zhou (2004) top-5% node-node
  link distribution (l_top5_any, l_top5_both), second-order assortativity R_avg and R_max
  (exact Zhou, Cox & Hansen excess-degree excluding-other-endpoint definition, plus the
  pre-change inclusive version R_avg_incl for traceability), coreness (k-core index) per
  node and its maximum.

Null models (seeds derived from config.SEED)
  rew_und : N_REWIRE undirected degree-preserving rewirings of Gu (igraph edge swaps,
            REWIRE_SWAPS_PER_EDGE * E_und swap trials, no self-loops or multi-edges).
            Used for every undirected metric and for rho(k) = phi(k) / <phi_rew(k)>.
  rew_dir : N_REWIRE directed rewirings of Gd preserving every k_in and k_out
            (REWIRE_SWAPS_PER_EDGE * E swap trials). Used for the directed metrics.
  er_und  : N_REWIRE Erdos-Renyi G(N, E_und) graphs (undirected metrics).
  er_dir  : N_REWIRE directed G(N, E) graphs (directed metrics).
  perm100 : 100 random permutations of one of the two directed-link sequences (length 2L)
            for r_und, R_avg, and R_max (Zhou, Cox & Hansen Sec. 2.2.2).
  Why two rewiring families: directed in/out-preserving swaps do not preserve the
  undirected degree of Gu (they destroy reciprocated pairs, so E_und grows), which would
  confound clustering, rich club and k-core with a change of density. The undirected
  metrics are therefore normalised by the null that fixes the degree sequence of the
  graph they are computed on. See decisions.md.

Outputs
  results/b02_basic_properties.csv  long: date, backbone, metric, null, real, null_mean, null_sd, z, n_null
  results/b02_rich_club_curves.csv  date, backbone, k, N_k, E_k, phi_real, phi_rew_mean, phi_rew_sd, rho, z, phi_er_mean
  results/b02_rich_club_rank.csv    rich club by rank: top-n nodes (config.RICH_CLUB_TOP_N), phi real, rewired
                                    mean/sd/min/max, rho = phi / <phi_rew>, ceiling 1 / <phi_rew> (max possible rho)
  results/b02_rich_club_zhou.csv    main nights: Zhou (2004, 2007) raw phi(r) vs normalised rank r/N (log-spaced),
                                    normalised rho(r) secondary, top-5% node-node link distribution, clique size
  results/b02_zhou_table.csv        main nights: Zhou & Mondragon (2007) Table 1 properties on giant component
  results/b02_rewiring_convergence.csv  convergence of r and phi(top-50) vs swaps/edge (up to 1000 x L on ordinary)
  results/b02_coreness.parquet      node x network coreness (columns id, <date>__<backbone>)
  results/b02_degree_ccdf.parquet   date, backbone, kind (in/out/und), k, ccdf = P(K >= k)
  results/b02_structural_constraints.csv  main nights: [r_min, r_max] from connected max-(dis)assortative rewiring,
                                    connected maximal-random Case I, and connected P(k)+phi(k)-preserving Case II
  results/b02_bosam_order_<group>_<date>__<backbone>.csv / b02_bosam_edges_*.parquet  exact iterative BOSAM
  results/b02_powerlaw.csv          main nights: CSN discrete power-law fit of k_in, k_out vs lognormal etc.
  results/b02_representativeness.csv
  results/b02_reconciliation.md

Run: .venv/bin/python -I pipeline/b02_basic_properties.py [--workers 12]
"""
from __future__ import annotations

import argparse
import random
import warnings
import sys
import zlib
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import igraph as ig
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C  # noqa: E402
import b00_common as B  # noqa: E402
import b02_zhou as Z  # noqa: E402

NET_DIR = C.RESULTS / "b01_networks"
BACKBONES = ("w10", "equal_edges")
NULLS = ("rew_und", "rew_dir", "er_und", "er_dir")
NULL_DESCRIPTION = {
    "none": "observed value, no null",
    "rew_und": "UNDIRECTED degree-preserving rewiring of the undirected projection (every node degree fixed)",
    "rew_dir": "DIRECTED degree-preserving rewiring of the directed backbone (every k_in and k_out fixed)",
    "er_und": "undirected Erdos-Renyi G(N, E_und)",
    "er_dir": "directed Erdos-Renyi G(N, E)",
    "perm100": "100 random permutations of one directed-link sequence of length 2L (Zhou, Cox & Hansen Sec. 2.2.2)",
}
# Append new stream names at the end so all existing null seeds stay bit-for-bit identical.
SEED_STREAMS = NULLS + ("max_assortative", "max_disassortative", "rc_preserving", "perm100", "maxrand_conn", "convergence")
# Rich-club k grid: about 50 log-spaced integers plus the k used in the legacy table.
K_GRID = np.array(sorted(set(np.unique(np.round(np.geomspace(1, 2500, 50)).astype(int)).tolist())
                         | {10, 20, 30, 50, 75, 100, 150, 200, 300, 400}))
RHO_REPORT_K = (50, 100, 200)

# Metrics preserved exactly by a rewiring null (no null row is written for them).
_PRESERVED = {
    "rew_und": {"E_und", "density_und", "k_mean", "k_median", "k_max"},
    "rew_dir": {"E", "density", "kin_mean", "kin_median", "kin_max", "kout_mean", "kout_median", "kout_max"},
    "er_und": {"E_und", "density_und", "k_mean"},
    "er_dir": {"E", "density", "kin_mean", "kout_mean"},
}


# ---------------------------------------------------------------------------
# Metric functions (pure; unit-tested in tests/test_b_metrics.py)
# ---------------------------------------------------------------------------
def density_directed(n: int, m: int) -> float:
    return m / (n * (n - 1))


def density_undirected(n: int, m: int) -> float:
    return 2 * m / (n * (n - 1))


def path_stats(Gu: ig.Graph) -> tuple[float, int]:
    """(average shortest path length, diameter) on the giant connected component, unweighted."""
    giant = Gu.connected_components().giant()
    if giant.vcount() < 2:
        return float("nan"), 0
    h = giant.path_length_hist(directed=False)
    tot = s = 0
    diam = 0
    for start, _end, cnt in h.bins():
        if cnt:
            tot += cnt
            s += start * cnt
            diam = max(diam, int(start))
    return s / tot, diam


def pearson(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    a0, b0 = a - a.mean(), b - b.mean()
    d = np.sqrt((a0 * a0).sum() * (b0 * b0).sum())
    return float((a0 * b0).sum() / d) if d > 0 else float("nan")


def r_out_in(Gd: ig.Graph) -> float:
    """Directed degree assortativity: Pearson(k_out(u), k_in(v)) over edges u->v."""
    el = np.array(Gd.get_edgelist(), dtype=np.int64).reshape(-1, 2)
    kout = np.array(Gd.outdegree())
    kin = np.array(Gd.indegree())
    return pearson(kout[el[:, 0]], kin[el[:, 1]])


def r_undirected(Gu: ig.Graph) -> float:
    """Undirected degree assortativity (Newman 2002): Pearson over both ends of every edge."""
    el = np.array(Gu.get_edgelist(), dtype=np.int64).reshape(-1, 2)
    k = np.array(Gu.degree())
    a = np.concatenate([k[el[:, 0]], k[el[:, 1]]])
    b = np.concatenate([k[el[:, 1]], k[el[:, 0]]])
    return pearson(a, b)


def reciprocity(Gd: ig.Graph) -> float:
    """Share of directed edges u->v whose reverse v->u also exists."""
    if Gd.ecount() == 0:
        return float("nan")
    return float(Gd.reciprocity(ignore_loops=True, mode="default"))


def rich_club_phi(Gu: ig.Graph, ks) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """phi(k) = 2 E_k / (N_k (N_k - 1)), nodes with degree > k. Returns (N_k, E_k, phi); phi = nan if N_k < 2."""
    ks = np.asarray(ks)
    deg = np.array(Gu.degree())
    el = np.array(Gu.get_edgelist(), dtype=np.int64).reshape(-1, 2)
    mind = np.sort(np.minimum(deg[el[:, 0]], deg[el[:, 1]])) if len(el) else np.array([], int)
    sd = np.sort(deg)
    Nk = len(sd) - np.searchsorted(sd, ks, side="right")
    Ek = len(mind) - np.searchsorted(mind, ks, side="right")
    with np.errstate(invalid="ignore", divide="ignore"):
        phi = np.where(Nk >= 2, 2.0 * Ek / (Nk * (Nk - 1.0)), np.nan)
    return Nk, Ek, phi


def second_order_r_incl(Gu: ig.Graph) -> float:
    """Pre-change inclusive second-order assortativity R_avg_incl (kept for traceability):
    Pearson over both ends of every edge of the plain average neighbour degree k_nn
    INCLUDING the other endpoint of the edge."""
    el = np.array(Gu.get_edgelist(), dtype=np.int64).reshape(-1, 2)
    deg = np.array(Gu.degree(), float)
    s = np.zeros(len(deg))
    np.add.at(s, el[:, 0], deg[el[:, 1]])
    np.add.at(s, el[:, 1], deg[el[:, 0]])
    knn = np.divide(s, deg, out=np.zeros_like(deg), where=deg > 0)
    a = np.concatenate([knn[el[:, 0]], knn[el[:, 1]]])
    b = np.concatenate([knn[el[:, 1]], knn[el[:, 0]]])
    return pearson(a, b)


def second_order_r(Gu: ig.Graph) -> float:
    """Exact Zhou, Cox & Hansen (0903.0687) second-order assortativity R_avg:
    excess degrees (k - 1), neighbours of a excluding b, Pearson over both directions."""
    return Z.second_order_r_avg(Gu)


def second_order_r_max(Gu: ig.Graph) -> float:
    """Exact Zhou, Cox & Hansen (0903.0687) second-order assortativity R_max:
    excess degrees (k - 1), max over neighbours of a excluding b, Pearson over both directions."""
    return Z.second_order_r_max(Gu)


def rich_club_top(Gu: ig.Graph, ns) -> tuple[np.ndarray, np.ndarray]:
    """Rich club by rank: phi among the top-n nodes by degree (ties broken by node index ascending).
    Returns (E_n, phi_n). Degree-preserving rewiring keeps every node's degree, so the club is the
    same node set in the real and rewired graphs."""
    deg = np.array(Gu.degree())
    order = np.lexsort((np.arange(len(deg)), -deg))
    rank = np.empty(len(deg), dtype=np.int64)
    rank[order] = np.arange(len(deg))
    el = np.array(Gu.get_edgelist(), dtype=np.int64).reshape(-1, 2)
    mr = np.sort(np.maximum(rank[el[:, 0]], rank[el[:, 1]])) if len(el) else np.array([], int)
    ns = np.asarray(ns)
    En = np.searchsorted(mr, ns, side="left")  # edges whose two ends both have rank < n
    return En, 2.0 * En / (ns * (ns - 1.0))


def und_metrics(Gu: ig.Graph) -> dict:
    n, m = Gu.vcount(), Gu.ecount()
    deg = np.array(Gu.degree())
    L, diam = path_stats(Gu)
    R_avg, R_max = Z.second_order_mixing(Gu)
    nld = Z.node_node_link_distribution(Gu, frac=0.05)
    return dict(
        N=n, E_und=m, density_und=density_undirected(n, m),
        giant_weak=int(max(Gu.connected_components().sizes())),
        k_mean=float(deg.mean()), k_median=float(np.median(deg)), k_max=int(deg.max()),
        L=L, diameter=diam,
        C_avg=float(Gu.transitivity_avglocal_undirected(mode="zero")),
        C_global=float(Gu.transitivity_undirected()),
        r_und=r_undirected(Gu),
        R_avg=R_avg,
        R_max=R_max,
        R_avg_incl=second_order_r_incl(Gu),
        kcore_max=int(max(Gu.coreness())),
        n_clique=int(Z.rich_club_clique_size(Gu)),
        l_top5_any=int(nld["l_top_any"]),
        l_top5_both=int(nld["l_top_both"]),
    )


def dir_metrics(Gd: ig.Graph) -> dict:
    n, m = Gd.vcount(), Gd.ecount()
    kin = np.array(Gd.indegree())
    kout = np.array(Gd.outdegree())
    return dict(
        E=m, density=density_directed(n, m),
        giant_strong=int(max(Gd.connected_components(mode="strong").sizes())),
        kin_mean=float(kin.mean()), kin_median=float(np.median(kin)), kin_max=int(kin.max()),
        kout_mean=float(kout.mean()), kout_median=float(np.median(kout)), kout_max=int(kout.max()),
        reciprocity=reciprocity(Gd), r_out_in=r_out_in(Gd),
    )


def ccdf(values) -> tuple[np.ndarray, np.ndarray]:
    """Unique k and P(K >= k) over all nodes (zeros included)."""
    v = np.sort(np.asarray(values))
    k = np.unique(v)
    return k, 1.0 - np.searchsorted(v, k, side="left") / len(v)


def zscore(real: float, null: np.ndarray) -> tuple[float, float, float]:
    null = np.asarray(null, float)
    null = null[~np.isnan(null)]
    if len(null) < 2:
        return (float(null.mean()) if len(null) else float("nan")), float("nan"), float("nan")
    mu, sd = float(null.mean()), float(null.std(ddof=1))
    # sd below 1e-12 relative is floating-point noise around a value the null preserves exactly.
    return mu, sd, ((real - mu) / sd if sd > 1e-12 * max(1.0, abs(mu)) else float("nan"))


# ---------------------------------------------------------------------------
# Null models
# ---------------------------------------------------------------------------
def null_seed(date: str, backbone: str, null: str, rep: int) -> int:
    ss = np.random.SeedSequence([C.SEED, zlib.crc32(date.encode()), BACKBONES.index(backbone), SEED_STREAMS.index(null), rep])
    return int(ss.generate_state(1)[0])


def rewire(G: ig.Graph, seed: int, swaps_per_edge: int = C.REWIRE_SWAPS_PER_EDGE, keep_connected: bool = False) -> ig.Graph:
    """Degree-preserving rewiring (simple graphs only).

    Undirected G: preserves every degree; when ``keep_connected=True``, only accepts swaps
    that keep all non-isolated nodes in one connected component (Zhou & Mondragon 2007).
    Directed G: preserves every k_in and k_out (igraph edge swaps).
    """
    if keep_connected:
        if G.is_directed():
            raise ValueError("keep_connected=True is defined for undirected graphs (Zhou 2007)")
        H, _ = Z.random_rewire(G, seed=seed, swaps_per_edge=swaps_per_edge, keep_connected=True)
        return H
    ig.set_random_number_generator(random.Random(seed))
    H = G.copy()
    H.rewire(n=swaps_per_edge * G.ecount(), allowed_edge_types="simple")
    return H


def erdos_renyi(n: int, m: int, directed: bool, seed: int) -> ig.Graph:
    ig.set_random_number_generator(random.Random(seed))
    return ig.Graph.Erdos_Renyi(n=n, m=m, directed=directed, loops=False)


# ---------------------------------------------------------------------------
# One (network, backbone) task
# ---------------------------------------------------------------------------
def analyse(date: str, backbone: str, n_null: int = C.N_REWIRE, swaps_per_edge: int = C.REWIRE_SWAPS_PER_EDGE) -> dict:
    warnings.filterwarnings("ignore", category=RuntimeWarning)  # nanmean over k with N_k < 2
    e = pd.read_parquet(NET_DIR / f"edges_{date}__{backbone}.parquet")
    N = C.EXPECTED_N_NODES
    Gd = B.to_igraph(e, N, directed=True)
    Gu = B.to_igraph(e, N, directed=False)

    real = {**dir_metrics(Gd), **und_metrics(Gu)}
    sig = Z.mixing_significance(Gu, seed=null_seed(date, backbone, "perm100", 0), n_perm=100)
    real["r_und_jackknife_sigma"] = sig["r_und"]["jackknife_sigma"]
    real["R_avg_jackknife_sigma"] = sig["R_avg"]["jackknife_sigma"]
    real["R_max_jackknife_sigma"] = sig["R_max"]["jackknife_sigma"]

    Nk, Ek, phi = rich_club_phi(Gu, K_GRID)
    TOPN = np.array(C.RICH_CLUB_TOP_N)
    En_real, phit_real = rich_club_top(Gu, TOPN)
    phit_rew = []

    nulls = {k: [] for k in NULLS}
    phi_rew, phi_er = [], []
    for rep in range(n_null):
        H = rewire(Gu, null_seed(date, backbone, "rew_und", rep), swaps_per_edge=swaps_per_edge)
        nulls["rew_und"].append(und_metrics(H))
        phi_rew.append(rich_club_phi(H, K_GRID)[2])
        phit_rew.append(rich_club_top(H, TOPN)[1])
        H = rewire(Gd, null_seed(date, backbone, "rew_dir", rep), swaps_per_edge=swaps_per_edge)
        d = dir_metrics(H)
        d["E_und"] = B.to_igraph(pd.DataFrame(H.get_edgelist(), columns=["i", "j"]), N, directed=False).ecount()
        nulls["rew_dir"].append(d)
        H = erdos_renyi(N, Gu.ecount(), False, null_seed(date, backbone, "er_und", rep))
        nulls["er_und"].append(und_metrics(H))
        phi_er.append(rich_club_phi(H, K_GRID)[2])
        H = erdos_renyi(N, Gd.ecount(), True, null_seed(date, backbone, "er_dir", rep))
        nulls["er_dir"].append(dir_metrics(H))

    phi_rew = np.array(phi_rew)
    phi_er = np.array(phi_er)
    with np.errstate(invalid="ignore", divide="ignore"):
        rew_mean = np.nanmean(phi_rew, axis=0) if n_null else np.full(len(K_GRID), np.nan)
        rew_sd = np.nanstd(phi_rew, axis=0, ddof=1) if n_null > 1 else np.full(len(K_GRID), np.nan)
        er_mean = np.nanmean(phi_er, axis=0) if n_null else np.full(len(K_GRID), np.nan)
        rho = np.where(rew_mean > 0, phi / rew_mean, np.nan)
        zrc = np.where(rew_sd > 1e-12 * np.maximum(1.0, np.abs(rew_mean)), (phi - rew_mean) / rew_sd, np.nan)
    rc = pd.DataFrame(dict(date=date, backbone=backbone, k=K_GRID, N_k=Nk, E_k=Ek, phi_real=phi,
                           phi_rew_mean=rew_mean, phi_rew_sd=rew_sd, rho=rho, z=zrc, phi_er_mean=er_mean))

    phit_rew = np.array(phit_rew).reshape(-1, len(TOPN))
    pt_mu = phit_rew.mean(axis=0)
    pt_sd = phit_rew.std(axis=0, ddof=1) if n_null > 1 else np.full(len(TOPN), np.nan)
    rank_rc = pd.DataFrame(dict(
        date=date, backbone=backbone, top_n=TOPN, deg_min_in_club=np.sort(Gu.degree())[::-1][TOPN - 1],
        E_n=En_real, phi_real=phit_real, phi_rew_mean=pt_mu, phi_rew_sd=pt_sd,
        phi_rew_min=phit_rew.min(axis=0), phi_rew_max=phit_rew.max(axis=0),
        rho=phit_real / pt_mu, rho_ceiling=1.0 / pt_mu,
        rho_rew_sd=pt_sd / pt_mu,
        z=np.where(pt_sd > 1e-12, (phit_real - pt_mu) / pt_sd, np.nan),
        null="rew_und"))

    rows = []
    for met, val in real.items():
        rows.append(dict(date=date, backbone=backbone, metric=met, null="none", real=float(val),
                         null_mean=np.nan, null_sd=np.nan, z=np.nan, n_null=0))
    for nm, lst in nulls.items():
        if not lst:
            continue
        df = pd.DataFrame(lst)
        for met in df.columns:
            if met in _PRESERVED.get(nm, set()) or met == "N" or met not in real:
                continue
            mu, sd, z = zscore(real[met], df[met].values)
            rows.append(dict(date=date, backbone=backbone, metric=met, null=nm, real=float(real[met]),
                             null_mean=mu, null_sd=sd, z=z, n_null=int(df[met].notna().sum())))
    # Zhou, Cox & Hansen (0903.0687 Sec. 2.2.2) 100-permutation sequence null for r_und, R_avg, R_max.
    for met in ("r_und", "R_avg", "R_max"):
        pm = sig[met]["perm"]
        rows.append(dict(date=date, backbone=backbone, metric=met, null="perm100", real=float(pm["real"]),
                         null_mean=float(pm["mean"]), null_sd=float(pm["sd"]), z=float(pm["z"]),
                         n_null=int(pm["n_perm"])))

    for k in RHO_REPORT_K:
        i = int(np.where(K_GRID == k)[0][0])
        rows.append(dict(date=date, backbone=backbone, metric=f"phi_k{k}", null="rew_und", real=float(phi[i]),
                         null_mean=float(rew_mean[i]), null_sd=float(rew_sd[i]), z=float(zrc[i]), n_null=n_null))
        rows.append(dict(date=date, backbone=backbone, metric=f"rho_k{k}", null="none", real=float(rho[i]),
                         null_mean=np.nan, null_sd=np.nan, z=np.nan, n_null=0))

    for _, q in rank_rc.iterrows():
        n_ = int(q.top_n)
        rows.append(dict(date=date, backbone=backbone, metric=f"phi_top{n_}", null="rew_und", real=q.phi_real,
                         null_mean=q.phi_rew_mean, null_sd=q.phi_rew_sd, z=q.z, n_null=n_null))
        rows.append(dict(date=date, backbone=backbone, metric=f"rho_top{n_}", null="none", real=q.rho,
                         null_mean=np.nan, null_sd=np.nan, z=np.nan, n_null=0))
        rows.append(dict(date=date, backbone=backbone, metric=f"rho_ceiling_top{n_}", null="none", real=q.rho_ceiling,
                         null_mean=np.nan, null_sd=np.nan, z=np.nan, n_null=0))

    cc = []
    for kind, vals in (("in", Gd.indegree()), ("out", Gd.outdegree()), ("und", Gu.degree())):
        k, p = ccdf(vals)
        cc.append(pd.DataFrame(dict(date=date, backbone=backbone, kind=kind, k=k, ccdf=p)))
    return dict(rows=rows, rc=rc, rank_rc=rank_rc, ccdf=pd.concat(cc, ignore_index=True),
                coreness=(f"{date}__{backbone}", np.array(Gu.coreness(), dtype=np.int32)))


def _analyse_star(args):
    return analyse(*args)


# ---------------------------------------------------------------------------
# Zhou-style structural-constraint checks, BOSAM, Zhou rich club, Table 1, power-law (main nights)
# ---------------------------------------------------------------------------
def zhou_task(
    date: str,
    group: str,
    backbone: str,
    n_null: int = C.N_REWIRE,
    swaps_per_edge: int = C.REWIRE_SWAPS_PER_EDGE,
) -> dict:
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    e = pd.read_parquet(NET_DIR / f"edges_{date}__{backbone}.parquet")
    N = C.EXPECTED_N_NODES
    Gd = B.to_igraph(e, N, directed=True)
    Gu = B.to_igraph(e, N, directed=False)
    r_real = r_undirected(Gu)
    R_avg_real, R_max_real = Z.second_order_mixing(Gu)
    R_incl_real = second_order_r_incl(Gu)
    label = f"{group}_{date}__{backbone}"

    # Zhou 2007 Case I: connected max-assortative, max-disassortative, and maximal-random rewiring.
    Ha, ia = Z.extreme_rewire(Gu, "assortative", null_seed(date, backbone, "max_assortative", 0), keep_connected=True)
    Hd, idd = Z.extreme_rewire(Gu, "disassortative", null_seed(date, backbone, "max_disassortative", 0), keep_connected=True)
    assert Ha.degree() == Gu.degree() and Hd.degree() == Gu.degree() and Ha.is_simple() and Hd.is_simple()
    assert Z._n_nontrivial_components(Ha) == 1 and Z._n_nontrivial_components(Hd) == 1
    r_max, r_min = r_undirected(Ha), r_undirected(Hd)

    R_GRID = Z.rich_club_rank_grid(N, extra_ranks=C.RICH_CLUB_TOP_N)
    Er_real, phir_real = rich_club_top(Gu, R_GRID)
    phir_rew = []
    mr_r, mr_disc = [], []
    for rep in range(n_null):
        H, info = Z.random_rewire(Gu, null_seed(date, backbone, "maxrand_conn", rep),
                                  swaps_per_edge=swaps_per_edge, keep_connected=True)
        assert H.degree() == Gu.degree() and H.is_simple() and Z._n_nontrivial_components(H) == 1
        mr_r.append(r_undirected(H))
        mr_disc.append(info["disconnected_rejected"])
        phir_rew.append(rich_club_top(H, R_GRID)[1])
    mr_mu, mr_sd, mr_z = zscore(r_real, np.array(mr_r))

    # Zhou (2004, 2007) rich-club table: raw phi(r) first, normalised rho(r) secondary.
    phir_rew = np.array(phir_rew).reshape(-1, len(R_GRID))
    pr_mu = phir_rew.mean(axis=0)
    pr_sd = phir_rew.std(axis=0, ddof=1) if n_null > 1 else np.full(len(R_GRID), np.nan)
    nld = Z.node_node_link_distribution(Gu, frac=0.05)
    n_clique = Z.rich_club_clique_size(Gu)
    deg_desc = np.sort(Gu.degree())[::-1]
    rc_zhou = pd.DataFrame(dict(
        date=date, group=group, backbone=backbone,
        N=N, L=Gu.ecount(),
        r=R_GRID, r_norm=R_GRID / float(N),
        deg_min_in_club=deg_desc[R_GRID - 1],
        E_r=Er_real,
        phi_real=phir_real,
        phi_rew_mean=pr_mu, phi_rew_sd=pr_sd,
        rho=np.where(pr_mu > 0, phir_real / pr_mu, np.nan),
        rho_ceiling=np.where(pr_mu > 0, 1.0 / pr_mu, np.nan),
        z=np.where(pr_sd > 1e-12, (phir_real - pr_mu) / pr_sd, np.nan),
        top5_n=nld["top_n"],
        l_top5_any=nld["l_top_any"], share_l_top5_any=nld["share_l_top_any"],
        l_top5_both=nld["l_top_both"], share_l_top5_both=nld["share_l_top_both"],
        n_clique=n_clique,
        null="rew_und_connected",
    ))

    # Zhou 2007 Case II: connected P(k) + phi(k)-preserving rewiring ("algorithm 2").
    ks_all = np.arange(0, max(Gu.degree()) + 1)
    Ek0 = rich_club_phi(Gu, ks_all)[1]
    r_rc, R_rc, Rmax_rc, Rincl_rc, acc_rc = [], [], [], [], []
    for rep in range(n_null):
        H, acc = Z.rich_club_preserving_rewire(
            Gu, null_seed(date, backbone, "rc_preserving", rep), swaps_per_edge, keep_connected=True
        )
        assert H.degree() == Gu.degree() and H.is_simple() and Z._n_nontrivial_components(H) == 1
        assert np.array_equal(rich_club_phi(H, ks_all)[1], Ek0), "phi(k) not preserved"
        r_rc.append(r_undirected(H))
        ra, rm = Z.second_order_mixing(H)
        R_rc.append(ra)
        Rmax_rc.append(rm)
        Rincl_rc.append(second_order_r_incl(H))
        acc_rc.append(acc)
    mu, sd, z = zscore(r_real, np.array(r_rc))
    Rmu, Rsd, Rz = zscore(R_avg_real, np.array(R_rc))
    Rmmu, Rmsd, Rmz = zscore(R_max_real, np.array(Rmax_rc))
    Rimu, Risd, Riz = zscore(R_incl_real, np.array(Rincl_rc))

    # Exact iterative BOSAM (Guo, Chen & Zhou 2008 Sec. 2).
    key, binfo = Z.bosam_order(Gu, return_info=True)
    nodes = pd.read_parquet(NET_DIR / "nodes.parquet")
    order = pd.DataFrame(key, columns=["idx", "degree", "max_nb_degree", "max_nb_index"])
    order["max_nb_label"] = binfo["max_nb_label"]
    order.insert(0, "position", np.arange(len(order)))
    order["id"] = nodes.id.values[order.idx.values]
    pos = np.empty(N, dtype=np.int64)
    pos[order.idx.values] = order.position.values
    el = np.array(Gu.get_edgelist(), dtype=np.int64).reshape(-1, 2)
    bitmap = pd.DataFrame({"row": np.concatenate([pos[el[:, 0]], pos[el[:, 1]]]),
                           "col": np.concatenate([pos[el[:, 1]], pos[el[:, 0]]])})

    sc = dict(label=label, date=date, group=group, backbone=backbone, graph="undirected projection",
              keep_connected=True, swaps_per_edge=swaps_per_edge,
              r_real=r_real, r_min_maxdisassortative=r_min, r_max_maxassortative=r_max,
              position_in_range=(r_real - r_min) / (r_max - r_min),
              assort_sweeps=ia["sweeps"], assort_stop=ia["stop"], assort_accepted_last_sweep=ia["accepted_last_sweep"],
              assort_r_change_last_window=ia["r_change_last_window"],
              assort_disconnected_rejected=ia["disconnected_rejected"],
              disassort_sweeps=idd["sweeps"], disassort_stop=idd["stop"],
              disassort_accepted_last_sweep=idd["accepted_last_sweep"],
              disassort_r_change_last_window=idd["r_change_last_window"],
              disassort_disconnected_rejected=idd["disconnected_rejected"],
              maxrand_conn_n=n_null, maxrand_conn_r_mean=mr_mu, maxrand_conn_r_sd=mr_sd, maxrand_conn_r_z=mr_z,
              maxrand_conn_disc_rej_mean=float(np.mean(mr_disc)) if mr_disc else np.nan,
              rcpres_n=n_null, rcpres_r_mean=mu, rcpres_r_sd=sd, rcpres_r_z=z,
              rcpres_R_avg_real=R_avg_real, rcpres_R_avg_mean=Rmu, rcpres_R_avg_sd=Rsd, rcpres_R_avg_z=Rz,
              rcpres_R_max_real=R_max_real, rcpres_R_max_mean=Rmmu, rcpres_R_max_sd=Rmsd, rcpres_R_max_z=Rmz,
              rcpres_R_avg_incl_real=R_incl_real, rcpres_R_avg_incl_mean=Rimu,
              rcpres_R_avg_incl_sd=Risd, rcpres_R_avg_incl_z=Riz,
              rcpres_acceptance_mean=float(np.mean(acc_rc)) if acc_rc else np.nan,
              rcpres_phi_preserved_all_k=True,
              bosam_iterations=binfo["iterations"], bosam_converged=binfo["converged"])

    # Zhou 2007 Table 1 row.
    zrow = Z.zhou_table_row(Gu, date, group, backbone, seed=null_seed(date, backbone, "perm100", 1))

    pl = []
    for kind, vals in (("in", Gd.indegree()), ("out", Gd.outdegree())):
        pl.append(dict(date=date, group=group, backbone=backbone, kind=kind, **Z.powerlaw_vs_lognormal(vals)))
    return dict(sc=sc, label=label, order=order, bitmap=bitmap, pl=pl, rc_zhou=rc_zhou, zrow=zrow)


def _zhou_star(args):
    return zhou_task(*args)


# ---------------------------------------------------------------------------
# Rewiring convergence check (Zhou 2007: up to 1000 x L on ordinary main night)
# ---------------------------------------------------------------------------
def _convergence_rep(date: str, group: str, backbone: str, keep_connected: bool, swaps_per_edge: int, rep: int) -> dict:
    e = pd.read_parquet(NET_DIR / f"edges_{date}__{backbone}.parquet")
    Gu = B.to_igraph(e, C.EXPECTED_N_NODES, directed=False)
    seed = null_seed(date, backbone, "convergence", rep + 10000 * swaps_per_edge + (5000 if keep_connected else 0))
    if keep_connected:
        H, info = Z.random_rewire(Gu, seed=seed, swaps_per_edge=swaps_per_edge, keep_connected=True)
        disc_rej = int(info["disconnected_rejected"])
    else:
        H = rewire(Gu, seed=seed, swaps_per_edge=swaps_per_edge, keep_connected=False)
        disc_rej = 0
    ra, rm = Z.second_order_mixing(H)
    return dict(
        date=date, group=group, backbone=backbone,
        method="zhou_connected" if keep_connected else "igraph_simple",
        keep_connected=keep_connected,
        swaps_per_edge=swaps_per_edge,
        rep=rep,
        E_und=Gu.ecount(),
        r_und=r_undirected(H),
        phi_top50=float(rich_club_top(H, [50])[1][0]),
        R_avg=ra,
        R_max=rm,
        C_avg=float(H.transitivity_avglocal_undirected(mode="zero")),
        disconnected_rejected=disc_rej,
        n_edge_components=Z._n_nontrivial_components(H),
    )


def _convergence_star(args):
    return _convergence_rep(*args)


def run_convergence_check(workers: int = 12) -> pd.DataFrame:
    """Measure r_und and phi(top-50) vs swaps_per_edge on the three main nights (w10).

    Ordinary main night runs up to 1000 x L (Zhou 2007) for both igraph and connected
    rewiring; All Saints and Summer run up to 100 x L (and 300 x L igraph on All Saints).
    """
    plan = [
        # (group, date, backbone, keep_connected, swaps_list, n_reps)
        ("ordinary",  C.MAIN_NIGHTS["ordinary"],  "w10", False, [1, 3, 10, 30, 100, 300, 1000], 20),
        ("ordinary",  C.MAIN_NIGHTS["ordinary"],  "w10", True,  [10, 30, 100, 300, 1000],       10),
        ("allsaints", C.MAIN_NIGHTS["allsaints"], "w10", False, [10, 30, 100, 300],             10),
        ("allsaints", C.MAIN_NIGHTS["allsaints"], "w10", True,  [10, 30, 100],                  10),
        ("summer",    C.MAIN_NIGHTS["summer"],    "w10", False, [10, 30, 100],                  10),
        ("summer",    C.MAIN_NIGHTS["summer"],    "w10", True,  [10, 30, 100],                   5),
    ]
    tasks = []
    for g, d, bb, kc, sw_list, n_reps in plan:
        for sw in sw_list:
            for rep in range(n_reps):
                tasks.append((d, g, bb, kc, sw, rep))
    # Schedule heaviest jobs first for clean parallel load balance.
    tasks.sort(key=lambda t: -(t[4] * (3 if t[3] else 1) * (4 if t[1] == "summer" else 2 if t[1] == "allsaints" else 1)))
    with ProcessPoolExecutor(max_workers=workers) as ex:
        raw = list(ex.map(_convergence_star, tasks))
    df = pd.DataFrame(raw)
    agg = (
        df.groupby(["date", "group", "backbone", "method", "keep_connected", "swaps_per_edge", "E_und"], as_index=False)
        .agg(
            n_reps=("rep", "count"),
            r_mean=("r_und", "mean"),
            r_sd=("r_und", lambda s: float(s.std(ddof=1)) if len(s) > 1 else np.nan),
            phi_top50_mean=("phi_top50", "mean"),
            phi_top50_sd=("phi_top50", lambda s: float(s.std(ddof=1)) if len(s) > 1 else np.nan),
            R_avg_mean=("R_avg", "mean"),
            R_avg_sd=("R_avg", lambda s: float(s.std(ddof=1)) if len(s) > 1 else np.nan),
            R_max_mean=("R_max", "mean"),
            R_max_sd=("R_max", lambda s: float(s.std(ddof=1)) if len(s) > 1 else np.nan),
            C_avg_mean=("C_avg", "mean"),
            C_avg_sd=("C_avg", lambda s: float(s.std(ddof=1)) if len(s) > 1 else np.nan),
            disconnected_rejected_mean=("disconnected_rejected", "mean"),
            max_edge_components=("n_edge_components", "max"),
        )
        .sort_values(["group", "method", "swaps_per_edge"])
        .reset_index(drop=True)
    )
    # Compare every swap count against the largest swap count for that (group, method).
    diff_r, z_r, diff_p, z_p = [], [], [], []
    for _, row in agg.iterrows():
        sub = agg[(agg.group == row.group) & (agg.method == row.method)]
        ref = sub.loc[sub.swaps_per_edge.idxmax()]
        se_r = float(np.sqrt((row.r_sd ** 2) / row.n_reps + (ref.r_sd ** 2) / ref.n_reps))
        se_p = float(np.sqrt((row.phi_top50_sd ** 2) / row.n_reps + (ref.phi_top50_sd ** 2) / ref.n_reps))
        dr = float(row.r_mean - ref.r_mean)
        dp = float(row.phi_top50_mean - ref.phi_top50_mean)
        diff_r.append(dr)
        z_r.append(dr / se_r if se_r > 1e-12 else 0.0)
        diff_p.append(dp)
        z_p.append(dp / se_p if se_p > 1e-12 else 0.0)
    agg["diff_r_vs_max_swaps"] = diff_r
    agg["z_r_vs_max_swaps"] = z_r
    agg["diff_phi50_vs_max_swaps"] = diff_p
    agg["z_phi50_vs_max_swaps"] = z_p
    return agg


# ---------------------------------------------------------------------------
# Representativeness
# ---------------------------------------------------------------------------
def unique_extreme(main: float, values) -> bool:
    """True if ``main`` is the strict, unique minimum or maximum of ``values`` (which include it)."""
    v = np.asarray(values, float)
    v = v[~np.isnan(v)]
    if len(v) < 2 or np.isnan(main):
        return False
    if np.isclose(v.min(), v.max()):
        return False
    n_eq = int(np.isclose(v, main).sum())
    return n_eq == 1 and (np.isclose(main, v.min()) or np.isclose(main, v.max()))


def pairwise_sign_agreement(holiday, ordinary, direction: float) -> tuple[int, int]:
    """Count pairs (h, o) with sign(h - o) == direction. Returns (n_agree, n_pairs)."""
    h = np.asarray(holiday, float)[:, None]
    o = np.asarray(ordinary, float)[None, :]
    s = np.sign(h - o)
    ok = ~np.isnan(s)
    return int((s[ok] == direction).sum()), int(ok.sum())


def representativeness(long: pd.DataFrame, nets: pd.DataFrame) -> pd.DataFrame:
    """One row per (backbone, metric, night type). Metrics: every real value plus z vs rew_und / rew_dir.

    main_is_extreme: the main night is the unique strict min or max over its sister nights.
    direction_main : sign(main holiday night - ordinary main night).
    sister_pairs_agree / sister_pairs: pairs (holiday sister h, ordinary sister o) with
        sign(h - o) == direction_main; all_sisters_agree = every pair agrees.
    years_agree    : sign(holiday_y - ordinary_y) == direction_main for each replication year y present.
    robust         : all_sisters_agree and not main_is_extreme (ordinary rows: not main_is_extreme).
    """
    warnings.filterwarnings("ignore", category=RuntimeWarning)  # all-NaN z rows (null preserves the metric)
    real = long.drop_duplicates(["date", "backbone", "metric"])[["date", "backbone", "metric", "real"]]
    z = long[long.null.isin(["rew_und", "rew_dir"])].copy()
    z["metric"] = z["metric"] + "__z_" + z["null"]
    z = z[["date", "backbone", "metric", "z"]].rename(columns={"z": "real"})
    vals = pd.concat([real, z], ignore_index=True)
    V = vals.pivot_table(index=["backbone", "metric"], columns="date", values="real", aggfunc="first")

    sis = {g: sorted(nets[(nets.group == g) & (nets.role == "sister")].date) for g in C.MAIN_NIGHTS}
    yrs = {g: sorted(nets[(nets.group == g) & (nets.role == "year")].date) for g in C.MAIN_NIGHTS}
    ord_main = C.MAIN_NIGHTS["ordinary"]
    ord_year = {d[:4]: d for d in yrs["ordinary"]}
    out = []
    for (bb, met), row in V.iterrows():
        for g, dmain in C.MAIN_NIGHTS.items():
            sv = row.reindex(sis[g]).values.astype(float)
            mv = float(row.get(dmain, np.nan))
            rec = dict(backbone=bb, metric=met, group=g, main_date=dmain, main=mv,
                       sister_n=int((~np.isnan(sv)).sum()), sister_median=float(np.nanmedian(sv)),
                       sister_min=float(np.nanmin(sv)), sister_max=float(np.nanmax(sv)),
                       main_is_extreme=unique_extreme(mv, sv))
            for y in ("2022", "2023"):
                dy = [d for d in yrs[g] if d.startswith(y)]
                rec[f"year_{y}"] = float(row.get(dy[0], np.nan)) if dy else np.nan
            if g == "ordinary":
                rec.update(direction_main=np.nan, sister_pairs_agree=np.nan, sister_pairs=np.nan,
                           all_sisters_agree=np.nan, years_agree=np.nan,
                           robust=not rec["main_is_extreme"])
            else:
                o_main = float(row.get(ord_main, np.nan))
                direction = float(np.sign(mv - o_main))
                ov = row.reindex(sis["ordinary"]).values.astype(float)
                na, npairs = pairwise_sign_agreement(sv, ov, direction)
                ya = []
                for y in ("2022", "2023"):
                    dy = [d for d in yrs[g] if d.startswith(y)]
                    if dy and y in ord_year:
                        ya.append(np.sign(float(row.get(dy[0], np.nan)) - float(row.get(ord_year[y], np.nan))) == direction)
                rec.update(direction_main=direction, sister_pairs_agree=na, sister_pairs=npairs,
                           all_sisters_agree=bool(npairs > 0 and na == npairs),
                           years_agree=bool(all(ya)) if ya else np.nan)
                rec["robust"] = bool(rec["all_sisters_agree"] and not rec["main_is_extreme"])
            out.append(rec)
    return pd.DataFrame(out)


# ---------------------------------------------------------------------------
# Reconciliation with legacy table and pre-Zhou-fix archive
# ---------------------------------------------------------------------------
LEGACY_JSON = C.ROOT / "explore" / "out_p92_B_props.json"
ARCHIVE_DIR = C.RESULTS / "_archive_pre_rerun"
LEGACY_LABELS = {"Feb_wd_1d": "ordinary", "Aug_15": "summer", "All_Saints": "allsaints"}
# legacy key -> (source, metric, column) in the new outputs
LEGACY_MAP = {
    "total_away_M": ("manifest", "total_weight_full", 1e-6),
    "thr_flow_share_pct": ("manifest", "share_weight_retained", 100.0),
    "E_dir": ("long", "E", "real"),
    "E_und": ("long", "E_und", "real"),
    "density": ("long", "density_und", "real"),
    "k_mean": ("long", "k_mean", "real"),
    "k_med": ("long", "k_median", "real"),
    "k_max": ("long", "k_max", "real"),
    "gcc_size": ("long", "giant_weak", "real"),
    "recip_bin": ("long", "reciprocity", "real"),
    "C_avg": ("long", "C_avg", "real"),
    "C_rew_mean": ("long", "C_avg", "null_mean:rew_und"),
    "C_ER": ("long", "C_avg", "null_mean:er_und"),
    "L_avg": ("long", "L", "real"),
    "diam": ("long", "diameter", "real"),
    "L_rew_mean": ("long", "L", "null_mean:rew_und"),
    "r_deg": ("long", "r_und", "real"),
    "r_deg_rew_mean": ("long", "r_und", "null_mean:rew_und"),
    "r_dir_out_in": ("long", "r_out_in", "real"),
    "R_avg": ("long", "R_avg", "real"),
    "R_avg_rew_mean": ("long", "R_avg", "null_mean:rew_und"),
    "kcore_max": ("long", "kcore_max", "real"),
    "kcore_rew_mean": ("long", "kcore_max", "null_mean:rew_und"),
    "rho_k50": ("long", "rho_k50", "real"),
    "rho_k100": ("long", "rho_k100", "real"),
    "rho_k200": ("long", "rho_k200", "real"),
}
NOT_RECOMPUTED = {
    "k_in_small_med / k_in_small_mean / in_str_small_med / kcore_small_mean":
        "village-level statistics; they belong to R1 (k-core promotion, wave 2). The legacy values used the "
        "buggy buffer coastal flag of p73/p92, so they must be recomputed with coastal_true there.",
    "r_dir_in_out": "not in the b02 brief (only directed out-in is required).",
    "r_pop_weighted": "weighted population mixing on the full matrix; not a binary basic property, not in the brief.",
    "mean_km / share_100km / share_xprov / share_small_pct": "weighted reach statistics; not in the b02 brief.",
}


# Explanations for > 2% differences, keyed by legacy key or (legacy key, night type).
RECON_EXPLAIN: dict = {
    "C_ER": "Legacy C_ER is the analytic ER expectation p = density (N = 2,484); new is the mean C_avg of 50 "
            "G(N, E_und) draws with N = 2,435. Both the N change (density x 1.0407) and analytic-vs-sampled "
            "differ; the N-adjusted legacy value is within 0.5% of the new one.",
    "k_med": "N change: the median is now over 2,435 nodes; the 49 removed nodes all had degree 0, so the "
             "median moves up by one rank position (not a rescaling).",
    "R_avg": "Definition fixed to match Zhou, Cox & Hansen (arXiv:0903.0687): excess degrees (k-1) and "
             "neighbours of a EXCLUDING b. The legacy inclusive plain-degree version is retained as R_avg_incl "
             "and matches legacy p92 to 0.00%.",
    "R_avg_rew_mean": "Definition fixed to match Zhou, Cox & Hansen (arXiv:0903.0687): excess degrees (k-1) and "
                      "neighbours of a EXCLUDING b. Retained inclusive version R_avg_incl matches legacy rew_und.",
    ("kcore_rew_mean", "ordinary"):
        "Null-model sampling, not a method change: the legacy value is the mean of 25 rewirings, the new one of "
        "50 rewirings with different seeds (config.SEED). The difference (0.58) is below one null sd (about 0.7, "
        "column null_sd in b02_basic_properties.csv); the real max k-core (37) is identical and z stays about 10.",
}


N_LEGACY = 2484
# Legacy -> N = 2435 adjustment for metrics that depend on N only through the 49 isolated phantom nodes.
def _n_adjust(key: str, old: float):
    N, n0 = C.EXPECTED_N_NODES, N_LEGACY
    if key in ("density", "C_ER"):
        return old * (n0 * (n0 - 1)) / (N * (N - 1)), "density divides by N(N-1)"
    if key in ("k_mean", "C_avg", "C_rew_mean"):
        return old * n0 / N, "mean over nodes, the 49 removed nodes contributed 0"
    return None, ""


AUDIT_RANK_RC = {  # research_notes/audit_B_existing_code.md section (g), N = 2,435, w >= 10
    25: (1.18, 1.01, 0.97), 50: (1.23, 1.05, 0.92), 100: (1.19, 1.05, 0.81), 200: (1.10, 1.02, 0.59),
}


def _zhou_changes_section(long: pd.DataFrame) -> list[str]:
    """Compare pre-change archived b02 outputs (results/_archive_pre_rerun/) vs new Zhou-aligned outputs."""
    from scipy.stats import kendalltau, spearmanr

    lines = [
        "## Changes to match Zhou's definitions",
        "",
        "Pre-change archive: `results/_archive_pre_rerun/` (N = 2,435, pre-Zhou-fix definitions). "
        "New outputs: `results/b02_*.csv` following Zhou & Mondragon (2004, 2007), Guo, Chen & Zhou (2008), "
        "and Zhou, Cox & Hansen (2009).",
        "",
    ]
    old_bp_path = ARCHIVE_DIR / "b02_basic_properties.csv"
    old_sc_path = ARCHIVE_DIR / "b02_structural_constraints.csv"
    new_sc_path = C.RESULTS / "b02_structural_constraints.csv"
    new_zt_path = C.RESULTS / "b02_zhou_table.csv"
    new_conv_path = C.RESULTS / "b02_rewiring_convergence.csv"

    if old_bp_path.exists():
        old_bp = pd.read_csv(old_bp_path, dtype={"date": str})
        lines += [
            "### 1. Second-order assortative mixing (Zhou, Cox & Hansen, arXiv:0903.0687)",
            "",
            "Old `R_avg` used plain node degrees and included the other endpoint `b` in `k_nn(a)` "
            "(now retained as `R_avg_incl`). New `R_avg` and `R_max` use excess degrees `(k - 1)` "
            "and neighbours of `a` excluding `b`, with jackknife sigma and 100-permutation null.",
            "",
            "| Night | Backbone | Old `R_avg` (= new `R_avg_incl`) | New `R_avg` | Jackknife sigma | `R_avg` rew_und mean | `R_avg` perm100 sd | New `R_max` | Jackknife sigma | `R_max` rew_und mean |",
            "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
        for g, d in C.MAIN_NIGHTS.items():
            for bb in BACKBONES:
                o_ra = float(old_bp[(old_bp.date == d) & (old_bp.backbone == bb) & (old_bp.metric == "R_avg") & (old_bp.null == "none")].real.iloc[0])
                n_ra = float(long[(long.date == d) & (long.backbone == bb) & (long.metric == "R_avg") & (long.null == "none")].real.iloc[0])
                n_ras = float(long[(long.date == d) & (long.backbone == bb) & (long.metric == "R_avg_jackknife_sigma") & (long.null == "none")].real.iloc[0])
                n_rarew = float(long[(long.date == d) & (long.backbone == bb) & (long.metric == "R_avg") & (long.null == "rew_und")].null_mean.iloc[0])
                n_raperm = float(long[(long.date == d) & (long.backbone == bb) & (long.metric == "R_avg") & (long.null == "perm100")].null_sd.iloc[0])
                n_rm = float(long[(long.date == d) & (long.backbone == bb) & (long.metric == "R_max") & (long.null == "none")].real.iloc[0])
                n_rms = float(long[(long.date == d) & (long.backbone == bb) & (long.metric == "R_max_jackknife_sigma") & (long.null == "none")].real.iloc[0])
                n_rmrew = float(long[(long.date == d) & (long.backbone == bb) & (long.metric == "R_max") & (long.null == "rew_und")].null_mean.iloc[0])
                lines.append(
                    f"| {g} ({d}) | {bb} | {o_ra:.6f} | {n_ra:.6f} | {n_ras:.6f} | {n_rarew:.6f} | {n_raperm:.6f} | {n_rm:.6f} | {n_rms:.6f} | {n_rmrew:.6f} |"
                )
        lines += [""]

    if old_sc_path.exists() and new_sc_path.exists():
        old_sc = pd.read_csv(old_sc_path).set_index("label")
        new_sc = pd.read_csv(new_sc_path).set_index("label")
        lines += [
            "### 2. Structural-constraint rewiring with single-component constraint (Zhou & Mondragon 2007)",
            "",
            "Added `keep_connected=True` (rewire only if all non-isolated nodes remain in one connected component) "
            "to max-assortative, max-disassortative, maximal-random Case I, and `P(k)+phi(k)`-preserving Case II.",
            "",
            "| Label | Old `r_min` | New `r_min` | Old `r_max` | New `r_max` | Old `rcpres_r_mean` | New `rcpres_r_mean` | Old `rcpres_R_avg_mean` | New `rcpres_R_avg_mean` | Assort disc. rej. |",
            "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
        for lab in new_sc.index:
            o = old_sc.loc[lab]
            n = new_sc.loc[lab]
            lines.append(
                f"| {lab} | {o.r_min_maxdisassortative:.6f} | {n.r_min_maxdisassortative:.6f} | "
                f"{o.r_max_maxassortative:.6f} | {n.r_max_maxassortative:.6f} | "
                f"{o.rcpres_r_mean:.6f} | {n.rcpres_r_mean:.6f} | "
                f"{o.rcpres_R_avg_mean:.6f} | {n.rcpres_R_avg_mean:.6f} | {int(n.assort_disconnected_rejected)} |"
            )
        lines += [""]

    # BOSAM order comparison vs archive.
    lines += [
        "### 3. BOSAM iterative ordering (Guo, Chen & Zhou 2008 Sec. 2)",
        "",
        "Third sort key changed from raw neighbour label (`max_nb_label`) to largest neighbour index in the "
        "new ordering (`max_nb_index`), iterated until stable. Column `max_nb_label` is kept for traceability.",
        "",
        "| Label | Converged | Iterations | Exact position match | Spearman rho | Kendall tau |",
        "|---|---|---:|---:|---:|---:|",
    ]
    sc_df = pd.read_csv(new_sc_path).set_index("label") if new_sc_path.exists() else None
    for g, d in C.MAIN_NIGHTS.items():
        for bb in BACKBONES:
            lab = f"{g}_{d}__{bb}"
            op = ARCHIVE_DIR / f"b02_bosam_order_{lab}.csv"
            np_ = C.RESULTS / f"b02_bosam_order_{lab}.csv"
            if op.exists() and np_.exists():
                oo = pd.read_csv(op).sort_values("position")
                no = pd.read_csv(np_).sort_values("position")
                exact = float((oo["idx"].values == no["idx"].values).mean())
                opos = np.empty(len(oo), dtype=np.int64); opos[oo["idx"].values] = oo["position"].values
                npos = np.empty(len(no), dtype=np.int64); npos[no["idx"].values] = no["position"].values
                sp = float(spearmanr(opos, npos).statistic)
                kt = float(kendalltau(opos, npos).statistic)
                conv = bool(sc_df.loc[lab, "bosam_converged"]) if sc_df is not None else True
                iters = int(sc_df.loc[lab, "bosam_iterations"]) if sc_df is not None else -1
                lines.append(f"| {lab} | {conv} | {iters} | {100.0 * exact:.2f}% | {sp:.6f} | {kt:.6f} |")
    lines += [""]

    if new_zt_path.exists():
        zt = pd.read_csv(new_zt_path, dtype={"date": str})
        lines += [
            "### 4. Rich club (Zhou 2004, 2007) and Zhou 2007 Table 1",
            "",
            "Raw `phi(r)` vs normalised rank `r/N` is now reported first in `results/b02_rich_club_zhou.csv`, "
            "with normalised `rho(r)` secondary, alongside top-5% node-node link distribution and clique size. "
            "Zhou 2007 Table 1 properties are in `results/b02_zhou_table.csv`.",
            "",
            "| Night | Backbone | N (gc) | L | mean k | raw phi(top-50) | rho(top-50) | phi(k>40) | n_clique | l(top5%, any)/L | l(top5%, top5%)/L | l | l* |",
            "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
        for _, row in zt.iterrows():
            p50 = float(long[(long.date == row.date) & (long.backbone == row.backbone) & (long.metric == "phi_top50") & (long.null == "rew_und")].real.iloc[0])
            r50 = float(long[(long.date == row.date) & (long.backbone == row.backbone) & (long.metric == "rho_top50") & (long.null == "none")].real.iloc[0])
            lines.append(
                f"| {row.group} ({row.date}) | {row.backbone} | {int(row.N)} | {int(row.L)} | {row.k_mean:.2f} | "
                f"{p50:.4f} | {r50:.4f} | {row.phi_k_gt_40:.4f} | {int(row.n_clique)} | "
                f"{100*row.share_l_top5_any:.1f}% | {100*row.share_l_top5_both:.1f}% | {row.l:.3f} | {row.l_star:.3f} |"
            )
        lines += [""]

    if new_conv_path.exists():
        cv = pd.read_csv(new_conv_path, dtype={"date": str})
        lines += [
            "### 5. Rewiring convergence check (`results/b02_rewiring_convergence.csv`)",
            "",
            "| Night | Method | Swaps/edge | Reps | r mean | r sd | z(r vs max swaps) | phi(top-50) mean | phi(top-50) sd | z(phi50 vs max swaps) |",
            "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
        for _, row in cv.iterrows():
            lines.append(
                f"| {row.group} | {row.method} | {int(row.swaps_per_edge)} | {int(row.n_reps)} | "
                f"{row.r_mean:.6f} | {row.r_sd:.6f} | {row.z_r_vs_max_swaps:+.2f} | "
                f"{row.phi_top50_mean:.6f} | {row.phi_top50_sd:.6f} | {row.z_phi50_vs_max_swaps:+.2f} |"
            )
        lines += [""]

    return lines


def reconciliation(long: pd.DataFrame, man: pd.DataFrame, explain: dict | None = None) -> str:
    import json

    legacy = json.loads(LEGACY_JSON.read_text())["overnight"]
    explain = explain or {}
    lines = [
        "# b02 reconciliation with the legacy table",
        "",
        "Legacy: `research_notes/Preguntas alternativas desde cero/network_layer_A_B.md` section 2.1, columns "
        "Feb 14, Aug 15, Nov 1. The table was produced by `explore/p92_B_networks_topology.py`; the legacy values "
        "below are read from its output `explore/out_p92_B_props.json` (full precision; they round to the "
        "values printed in the table). New: `pipeline/b02_basic_properties.py`, main backbone (w >= 10), "
        "outputs `results/b02_basic_properties.csv` and `results/b01_networks/manifest.csv`.",
        "",
        "Generated by `pipeline/b02_basic_properties.py` (function `reconciliation`).",
        "",
        "**Node set.** The legacy table uses N = 2,484, which includes 49 uninhabited shared-territory zones "
        "with 11-digit ids that never appear in the MITMA data (audit M3). They are isolated nodes. The new "
        f"pipeline drops them (config.DROP_LONG_IDS), so N = {C.EXPECTED_N_NODES}. Isolated nodes change only "
        "metrics that divide by N or average over all nodes: density (x 2484*2483 / (2435*2434) = 1.0407), "
        "mean degree and C_avg with k<2 nodes counted as 0 (x 2484/2435 = 1.0201), and the median degree. "
        "They do not change edges, components, L, diameter, r, R_avg, phi(k), rho(k) or k-cores. The column "
        "'Legacy adj. N' applies that factor; where it matches the new value, the N change is the whole explanation.",
        "",
        "| Metric (legacy key) | Night | Legacy (N=2484) | Legacy adj. N | New (N=2435) | Rel. diff | Explanation (if > 2%) |",
        "|---|---|---:|---:|---:|---:|---|",
    ]
    n_big = 0
    for key, (src, met, col) in LEGACY_MAP.items():
        for lab, g in LEGACY_LABELS.items():
            date = C.MAIN_NIGHTS[g]
            old = float(legacy[lab][key])
            if src == "manifest":
                new = float(man[(man.date == date) & (man.backbone == "w10")][met].iloc[0]) * col
            else:
                if col == "real":
                    sel = long[(long.date == date) & (long.backbone == "w10") & (long.metric == met)]
                    new = float(sel.real.iloc[0])
                else:
                    stat, nl = col.split(":")
                    sel = long[(long.date == date) & (long.backbone == "w10") & (long.metric == met) & (long.null == nl)]
                    new = float(sel[stat].iloc[0])
            rel = abs(new - old) / abs(old) if old != 0 else float("nan")
            big = rel > 0.02
            n_big += big
            adj, why = _n_adjust(key, old)
            ex = explain.get((key, g), explain.get(key, "")) if big else ""
            is_null = key.endswith("_rew_mean") or key == "C_ER"
            tol = 0.01 if is_null else 0.001   # null means also carry sampling noise (25 vs 50 draws)
            if big and not ex and adj is not None and abs(new - adj) / abs(adj) < tol:
                ex = (f"N change only ({why}); the N-adjusted legacy value matches to < {100 * tol:.1f}%"
                      + (" (null mean, so sampling noise is allowed)." if is_null else "."))
            if big and not ex:
                ex = "**UNEXPLAINED**"
            adj_s = f"{adj:.6g}" if adj is not None else ""
            lines.append(f"| {key} | {g} | {old:.6g} | {adj_s} | {new:.6g} | {100 * rel:.2f}% | {ex} |")
    lines += ["", f"Differences larger than 2%: {n_big}.", ""]
    rk = C.RESULTS / "b02_rich_club_rank.csv"
    if rk.exists():
        R = pd.read_csv(rk, dtype={"date": str})
        R = R[R.backbone == "w10"].set_index(["date", "top_n"])
        lines += ["## Rich club by rank vs the audit re-run", "",
                  "Audit values: `research_notes/audit_B_existing_code.md` section (g) (N = 2,435, w >= 10, "
                  "undirected rewiring, its own seeds). New: `results/b02_rich_club_rank.csv`.", "",
                  "| Top-n | Feb 14 rho audit | new | Aug 15 rho audit | new | Aug 15 phi_rew audit | new | Aug 15 ceiling 1/phi_rew |",
                  "|---:|---:|---:|---:|---:|---:|---:|---:|"]
        for n_, (fa, aa, pa) in AUDIT_RANK_RC.items():
            f = R.loc[(C.MAIN_NIGHTS["ordinary"], n_)]
            a_ = R.loc[(C.MAIN_NIGHTS["summer"], n_)]
            lines.append(f"| {n_} | {fa:.2f} | {f.rho:.3f} | {aa:.2f} | {a_.rho:.3f} | {pa:.2f} | {a_.phi_rew_mean:.3f} | {a_.rho_ceiling:.3f} |")
        lines += [""]
    lines += _zhou_changes_section(long)
    lines += ["## Legacy rows not recomputed in b02", ""]
    lines += [f"- `{k}`: {v}" for k, v in NOT_RECOMPUTED.items()]
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--only-representativeness", action="store_true")
    ap.add_argument("--skip-convergence", action="store_true")
    ap.add_argument("--swaps-per-edge", type=int, default=C.REWIRE_SWAPS_PER_EDGE)
    a = ap.parse_args()

    nets = pd.read_csv(NET_DIR / "networks.csv", dtype=str)
    man = pd.read_csv(NET_DIR / "manifest.csv", dtype={"date": str})
    nodes = pd.read_parquet(NET_DIR / "nodes.parquet")
    assert len(nodes) == C.EXPECTED_N_NODES
    out_long = C.RESULTS / "b02_basic_properties.csv"

    if not a.only_representativeness:
        tasks = [(d, bb, C.N_REWIRE, a.swaps_per_edge) for d, bb in zip(man.date, man.backbone)]
        # Biggest graphs first for better load balance.
        tasks.sort(key=lambda t: -int(man[(man.date == t[0]) & (man.backbone == t[1])].E.iloc[0]))
        ztasks = [(d, g, bb, C.N_REWIRE, a.swaps_per_edge) for g, d in C.MAIN_NIGHTS.items() for bb in BACKBONES]
        res, zres = [], []
        with ProcessPoolExecutor(max_workers=a.workers) as ex:
            zfut = [ex.submit(_zhou_star, t) for t in ztasks]  # slowest jobs first
            for r in ex.map(_analyse_star, tasks):
                res.append(r)
                print("done", r["coreness"][0], flush=True)
            for f in zfut:
                zres.append(f.result())
                print("done zhou", zres[-1]["label"], flush=True)
        long = pd.DataFrame([x for r in res for x in r["rows"]]).sort_values(["date", "backbone", "metric", "null"])
        long["null_description"] = long["null"].map(NULL_DESCRIPTION)
        long.to_csv(out_long, index=False)
        pd.concat([r["rc"] for r in res]).sort_values(["date", "backbone", "k"]).to_csv(
            C.RESULTS / "b02_rich_club_curves.csv", index=False)
        pd.concat([r["rank_rc"] for r in res]).sort_values(["date", "backbone", "top_n"]).to_csv(
            C.RESULTS / "b02_rich_club_rank.csv", index=False)
        pd.concat([r["ccdf"] for r in res]).sort_values(["date", "backbone", "kind", "k"]).to_parquet(
            C.RESULTS / "b02_degree_ccdf.parquet", index=False)
        node_cols = ["idx", "id", "name", "prov", "ccaa", "population", "n_municipalities", "area_km2", "is_am",
                     "coastal_true", "loss_peak_2021"]
        core = nodes[node_cols].copy()
        for name, arr in sorted(r["coreness"] for r in res):
            core[name] = arr
        core.to_parquet(C.RESULTS / "b02_coreness.parquet", index=False)
        pd.DataFrame([z["sc"] for z in zres]).to_csv(C.RESULTS / "b02_structural_constraints.csv", index=False)
        pd.DataFrame([x for z in zres for x in z["pl"]]).to_csv(C.RESULTS / "b02_powerlaw.csv", index=False)
        pd.concat([z["rc_zhou"] for z in zres]).sort_values(["group", "backbone", "r"]).to_csv(
            C.RESULTS / "b02_rich_club_zhou.csv", index=False)
        pd.DataFrame([z["zrow"] for z in zres]).to_csv(C.RESULTS / "b02_zhou_table.csv", index=False)
        for z in zres:
            o = z["order"].merge(nodes[["idx", "population", "n_municipalities", "area_km2", "is_am"]], on="idx")
            o.sort_values("position").to_csv(C.RESULTS / f"b02_bosam_order_{z['label']}.csv", index=False)
            z["bitmap"].to_parquet(C.RESULTS / f"b02_bosam_edges_{z['label']}.parquet", index=False)

        if not a.skip_convergence:
            conv = run_convergence_check(workers=a.workers)
            conv.to_csv(C.RESULTS / "b02_rewiring_convergence.csv", index=False)
            print("done convergence check:", len(conv), "rows", flush=True)

    long = pd.read_csv(out_long, dtype={"date": str})
    rep = representativeness(long, nets)
    rep.to_csv(C.RESULTS / "b02_representativeness.csv", index=False)
    print(rep[rep.backbone == "w10"].groupby("group").robust.agg(["sum", "count"]))
    (C.RESULTS / "b02_reconciliation.md").write_text(reconciliation(long, man, RECON_EXPLAIN))


if __name__ == "__main__":
    main()
