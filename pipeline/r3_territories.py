"""r3_territories: mesoscale result R3 ("summer territories") for question B.

Question: does the city x destination holiday network have community structure that groups
cities with the distant places their residents come from (non-contiguous "diaspora"
territories), beyond what geography alone produces, and is it stronger on holiday than on
an ordinary night?

Network (one per single night, main + sister nights, results/b01_networks edges = b00 loaders)
  rows    : the 137 mainland Census cities (INE 68545 cities that are zones of the N=2,435 set)
  columns : every mainland zone with inflow > 0 from those cities on that night
  weight  : O[c, d] = persons resident in city c sleeping in zone d (d != c), single night.
            No window means are used anywhere.

Communities: Barber bipartite modularity Q = (1/W) sum_{c,d} (O - E)[c, d] delta(g_c, g_d),
maximised with Leiden (leidenalg, CPM with resolution 0 on a signed two-layer multiplex:
layer weights +1 for B > 0 and -1 for |B| for B < 0, so the quality equals sum B delta exactly).
Same optimiser effort (N_RESTARTS seeds, n_iterations=-1) for every observed and null network;
the best-Q restart is kept and the mean pairwise NMI / ARI over restarts is reported.

Null families
  grav : E = doubly-constrained gravity with a distance-bin deterrence function, fitted by IPF
         (row sums, column sums and distance-bin totals matched; relative error < 1e-6 asserted).
         Fair null: O* ~ negative binomial with mean E and variance V(mu) = a mu^b fitted on the
         observed (O - E)^2 by binned log-log regression (keeps the overdispersion), the gravity
         model is REFITTED on O*, and O* - E* is clustered identically.
  bicm : binary links A = 1[O >= 10] (the w10 rule), BiCM (Saracco et al. 2015) with the
         city's own zone masked. Fair null: A* ~ Bernoulli(P), BiCM refitted on A*, clustered
         identically. Q normalised by the number of links.

Statistics per partition (territory = community with >= 1 city and >= MIN_DEST destinations)
  Q, n_territories, n_noncontig (territories with >= 2 polygon pieces of >= 3 zones),
  diaspora_share (share of within-territory flow, persons for grav / links for bicm, going to
  polygon pieces that neither contain nor touch a city of the territory), diaspora_km,
  roots_enrich = mean roots share r[c, prov(d)] over cross-province (c, d) pairs in the same
  community / the same mean over pairs in different communities, and its within-province
  permutation z (cities permuted only among cities of the same province, 1,000 permutations),
  nmi_prov (destination partition vs province), restart stability (NMI, ARI).

Alternatives (in --summarise): NMI of holiday city communities with roots-network communities,
and the per-city two-axis decomposition (roots axis, beach axis) against the same refitted NB null.

Run:   .venv/bin/python -I pipeline/r3_territories.py                       (all nights, ~1 h on 16 cores)
       .venv/bin/python -I pipeline/r3_territories.py --dates D1,D2 --tag T  (subset run, extra files *_T)
       .venv/bin/python -I pipeline/r3_territories.py --summarise            (tables, alternatives, plots)
The 2024-11-03, 2024-12-25 and 2024-03-29 nights were added with --tag addnights after the
orchestrator's date-convention update; 2024-10-31 (old All Saints sister) is run but dropped.
Outputs: results/r3_*.csv, results/r3_partitions*.parquet, results/r3_diag/*.png, r3_*.log;
per-night cache in results/r3_diag/cache (170 MB, git-ignored, needed by --summarise).
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
import zlib
from concurrent.futures import ProcessPoolExecutor
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C  # noqa: E402
import b00_common as B  # noqa: E402

NET = C.RESULTS / "b01_networks"
DIAG = C.RESULTS / "r3_diag"
CACHE = DIAG / "cache"
ROOTS_CSV = C.ROOT / "_local" / "external" / "ine_birth_residence" / "mun_residence_x_prov_birth_censo_2021_2025.csv"
ROOTS_YEAR = 2024

N_RESTARTS = 10           # Leiden seeds per network, identical for observed and null networks
N_NULL_MAIN = C.N_SPATIAL_DRAWS   # 100 draws per family for the three main nights
N_NULL_SISTER = 20        # draws per family for every sister night
N_PERM_ROOTS = 1000       # within-province permutations of city roots profiles
MIN_DEST = 5              # territory = >= 1 city and >= MIN_DEST destination zones
MIN_PIECE = 3             # a "major" polygon piece has >= MIN_PIECE zones
BIN_EDGES_KM = np.r_[0:30:2, 30:100:5, 100:300:10, 300:1200:25]   # as p93 (fixed bins, km)
BICM_W_MIN = C.BACKBONE_MAIN["w_min"]
FAMILIES = ("grav", "bicm")


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------
def load_inputs():
    nodes = B.load_nodes()
    ref = pd.read_parquet(NET / "nodes.parquet")
    assert (ref["id"].values == nodes["id"].values).all(), "node table differs from b01 output"
    rb = pd.read_csv(ROOTS_CSV, dtype={"ine_code": str, "prov_birth": str})
    pos = dict(zip(nodes["id"], nodes["idx"]))
    cities = sorted(c for c in rb["ine_code"].unique() if c in pos)
    assert len(cities) == 137, len(cities)
    city_idx = np.array([pos[c] for c in cities])
    assert not nodes.loc[city_idx, "is_am"].any()
    provs = sorted(nodes["prov"].unique())
    assert len(provs) == 47
    p2i = {p: k for k, p in enumerate(provs)}
    r = rb[(rb.year == ROOTS_YEAR) & rb.ine_code.isin(cities) & rb.prov_birth.isin(p2i)]
    R = np.zeros((len(cities), len(provs)))
    c2i = {c: k for k, c in enumerate(cities)}
    for c, p, v in zip(r.ine_code, r.prov_birth, r.persons):
        if c[:2] != p:
            R[c2i[c], p2i[p]] = v
    roots_share = R / R.sum(1, keepdims=True)
    return nodes, cities, city_idx, provs, roots_share


def polygon_adjacency(nodes: pd.DataFrame) -> np.ndarray:
    """Undirected pairs (a, b), a < b, of node indices whose polygons intersect (exact, as p93)."""
    import geopandas as gpd

    g = gpd.read_file(C.ZONES_GPKG).set_index("id").loc[nodes["id"]].reset_index()
    a, b = g.sindex.query(g.geometry, predicate="intersects")
    keep = a < b
    return np.unique(np.stack([a[keep], b[keep]], 1), axis=0)


def night_matrix(date: str, city_idx: np.ndarray, n: int, nodes: pd.DataFrame | None = None) -> np.ndarray:
    """City rows of the full single-night edge list: the b01 file if it exists, else b00.build_edges."""
    f = NET / f"edges_{date}__full.parquet"
    if f.exists():
        e = pd.read_parquet(f)
    else:
        e, _ = B.build_edges(date, nodes)
    rowpos = np.full(n, -1)
    rowpos[city_idx] = np.arange(len(city_idx))
    e = e[rowpos[e["i"].values] >= 0]
    O = np.zeros((len(city_idx), n))
    O[rowpos[e["i"].values], e["j"].values] = e["w"].values
    return O


# ---------------------------------------------------------------------------
# Null models
# ---------------------------------------------------------------------------
def fit_gravity(O: np.ndarray, ci, di, bi, tol: float = 1e-6, max_iter: int = 2000):
    """Doubly-constrained gravity with distance-bin deterrence by IPF on the masked cells.

    ci, di, bi: row, column and bin index of every allowed cell. Matches row sums, column sums
    and bin totals of O. Returns (E as a full matrix, relative margin error, iterations).
    """
    obs = O[ci, di]
    keys = (ci, di, bi)
    sizes = (O.shape[0], O.shape[1], int(bi.max()) + 1)
    targets = [np.bincount(k, weights=obs, minlength=s) for k, s in zip(keys, sizes)]
    ex = np.ones(len(obs))
    tot = obs.sum()
    err = np.inf
    for it in range(1, max_iter + 1):
        for k, t, s in zip(keys, targets, sizes):
            cur = np.bincount(k, weights=ex, minlength=s)
            ex *= np.divide(t, cur, out=np.zeros_like(t), where=cur > 0)[k]
        err = max(np.abs(np.bincount(k, weights=ex, minlength=s) - t).sum() / tot
                  for k, t, s in zip(keys, targets, sizes))
        if err < tol:
            break
    assert err < tol, f"gravity IPF did not converge: err={err:.2e} after {it} iterations"
    E = np.zeros_like(O)
    E[ci, di] = ex
    return E, float(err), it


def variance_function(O, E, ci, di, width: float = 0.25, min_cells: int = 50):
    """Empirical mean-variance relation of O around the gravity mean E.

    Cells with E > 1e-4 are binned by ln E in bins of equal width ``width``; bins with fewer
    than ``min_cells`` cells are dropped. Returns the table (bin mean E, bin mean (O - E)^2),
    used by ``var_of_mu`` (linear interpolation in log-log, linear extrapolation with the end
    slopes), and the Pearson dispersion sum((O-E)^2/E)/n. The residual variance includes any
    systematic deviation from gravity, so a null calibrated on it is conservative.
    """
    mu = E[ci, di]
    y = O[ci, di]
    ok = mu > 1e-4
    mu, y = mu[ok], y[ok]
    lm = np.log(mu)
    b = np.floor((lm - lm.min()) / width).astype(int)
    nb = b.max() + 1
    cnt = np.bincount(b, minlength=nb)
    m_mu = np.bincount(b, weights=mu, minlength=nb) / np.maximum(cnt, 1)
    m_v = np.bincount(b, weights=(y - mu) ** 2, minlength=nb) / np.maximum(cnt, 1)
    good = (cnt >= min_cells) & (m_v > 0)
    tab = pd.DataFrame({"mean_mu": m_mu[good], "mean_sq_resid": m_v[good], "n_cells": cnt[good]})
    pearson = float(((y - mu) ** 2 / mu).sum() / len(mu))
    return tab, pearson


def _end_slope(lx, ly, w):
    """WLS slope of log V on log mu over the given end bins, clipped to [1, 2]."""
    X = np.c_[np.ones(len(lx)), lx] * np.sqrt(w)[:, None]
    c, *_ = np.linalg.lstsq(X, ly * np.sqrt(w), rcond=None)
    return float(np.clip(c[1], 1.0, 2.0))


def var_of_mu(mu, var_mu, var_v, var_n=None, k_end: int = 6):
    """V(mu): log-log interpolation of the binned table; beyond its range, power-law extrapolation
    from the end bin with the WLS slope of the k_end end bins clipped to [1, 2]
    (between scaled-Poisson and NB2 behaviour)."""
    lx, ly = np.log(var_mu), np.log(var_v)
    w = np.ones(len(lx)) if var_n is None else np.asarray(var_n, float)
    lm = np.log(mu)
    out = np.interp(lm, lx, ly)
    s0 = _end_slope(lx[:k_end], ly[:k_end], w[:k_end])
    s1 = _end_slope(lx[-k_end:], ly[-k_end:], w[-k_end:])
    lo, hi = lm < lx[0], lm > lx[-1]
    out[lo] = ly[0] + s0 * (lm[lo] - lx[0])
    out[hi] = ly[-1] + s1 * (lm[hi] - lx[-1])
    return np.exp(out)


def draw_nb(E, var_mu, var_v, var_n, rng):
    """Negative binomial with mean E and variance max(V(E), E), V from the empirical table."""
    out = np.zeros_like(E)
    pos = E > 0
    mu = E[pos]
    V = np.maximum(var_of_mu(mu, var_mu, var_v, var_n), mu)
    extra = V - mu
    nb = extra > 1e-9 * mu
    x = np.empty_like(mu)
    r = mu[nb] ** 2 / extra[nb]
    x[nb] = rng.negative_binomial(r, r / (r + mu[nb]))
    x[~nb] = rng.poisson(mu[~nb])
    out[pos] = x
    return out


def fit_bicm(A: np.ndarray, mask: np.ndarray, tol: float = 1e-6, max_iter: int = 20000):
    """Bipartite configuration model with structural zeros (mask False): p = xy/(1+xy)."""
    kr, kc = A.sum(1), A.sum(0)
    # Saturated nodes (degree == number of allowed partners) have p = 1 on every allowed cell
    # (their Lagrange multiplier is infinite); fix them and fit the rest on the residual degrees.
    sat_r = np.zeros(A.shape[0], bool)
    sat_c = np.zeros(A.shape[1], bool)
    while True:
        free = mask & ~sat_r[:, None] & ~sat_c[None, :]
        res_r = kr - (mask & sat_c[None, :]).sum(1)
        res_c = kc - (mask & sat_r[:, None]).sum(0)
        new_r = ~sat_r & (res_r == free.sum(1)) & (free.sum(1) > 0)
        new_c = ~sat_c & (res_c == free.sum(0)) & (free.sum(0) > 0)
        if not new_r.any() and not new_c.any():
            break
        sat_r |= new_r
        sat_c |= new_c
    M = free.astype(float)
    kr_f, kc_f = np.where(sat_r, 0, res_r).astype(float), np.where(sat_c, 0, res_c).astype(float)
    x = np.where(kr_f > 0, kr_f / np.maximum(M.sum(1) - kr_f, 0.5), 0.0)
    y = np.where(kc_f > 0, kc_f / np.maximum(M.sum(0) - kc_f, 0.5), 0.0)
    x, y = np.sqrt(x), np.sqrt(y)
    err = np.inf
    for it in range(1, max_iter + 1):
        x = np.where(kr_f > 0, kr_f / np.maximum((M * y[None, :] / (1 + np.outer(x, y))).sum(1), 1e-300), 0.0)
        y = np.where(kc_f > 0, kc_f / np.maximum((M * x[:, None] / (1 + np.outer(x, y))).sum(0), 1e-300), 0.0)
        if it % 5 == 0:
            xy = np.outer(x, y)
            P = M * xy / (1 + xy)
            err = max(np.abs(P.sum(1) - kr_f).max(), np.abs(P.sum(0) - kc_f).max())
            if err < tol:
                break
    P = P + (mask & (sat_r[:, None] | sat_c[None, :]))
    err = max(np.abs(P.sum(1) - kr).max(), np.abs(P.sum(0) - kc).max())
    assert err < 10 * tol, f"BiCM did not converge: err={err:.2e}"
    assert P.max() <= 1 + 1e-12 and (P[~mask] == 0).all()
    return P, float(err), it


# ---------------------------------------------------------------------------
# Community detection
# ---------------------------------------------------------------------------
def leiden_bipartite(Bm: np.ndarray, n_restarts: int, seed: int):
    """Maximise sum_{c,d} Bm[c,d] delta(g_c, g_d) over partitions of rows+columns.

    Returns (memberships array [n_restarts, R+D], quality per restart).
    """
    import igraph as ig
    import leidenalg as la

    R, D = Bm.shape
    n = R + D
    ii, jj = np.nonzero(Bm)
    w = Bm[ii, jj]
    pos = w > 0
    gp = ig.Graph(n=n, edges=np.c_[ii[pos], R + jj[pos]].tolist())
    gn = ig.Graph(n=n, edges=np.c_[ii[~pos], R + jj[~pos]].tolist())
    wp, wn = w[pos].tolist(), (-w[~pos]).tolist()
    mems, quals = [], []
    ss = np.random.SeedSequence(seed).generate_state(n_restarts)
    for s in ss:
        p1 = la.CPMVertexPartition(gp, weights=wp, resolution_parameter=0.0)
        p2 = la.CPMVertexPartition(gn, weights=wn, resolution_parameter=0.0)
        opt = la.Optimiser()
        opt.set_rng_seed(int(s) % (2 ** 31 - 1))
        opt.optimise_partition_multiplex([p1, p2], layer_weights=[1, -1], n_iterations=-1)
        m = np.array(p1.membership)
        mems.append(m)
        quals.append(modularity_sum(Bm, m[:R], m[R:]))
    return np.array(mems), np.array(quals)


def modularity_sum(Bm, gr, gc) -> float:
    return float(Bm[gr[:, None] == gc[None, :]].sum())


def stability(mems: np.ndarray) -> tuple[float, float]:
    from sklearn.metrics import adjusted_rand_score as ari, normalized_mutual_info_score as nmi

    k = len(mems)
    a = [nmi(mems[i], mems[j]) for i in range(k) for j in range(i + 1, k)]
    b = [ari(mems[i], mems[j]) for i in range(k) for j in range(i + 1, k)]
    return float(np.mean(a)), float(np.mean(b))


# ---------------------------------------------------------------------------
# Partition statistics
# ---------------------------------------------------------------------------
def components(nodes_sub: np.ndarray, adj_list) -> list[np.ndarray]:
    s = set(int(x) for x in nodes_sub)
    seen, out = set(), []
    for u in s:
        if u in seen:
            continue
        st, comp = [u], [u]
        seen.add(u)
        while st:
            v = st.pop()
            for w in adj_list[v]:
                if w in s and w not in seen:
                    seen.add(w)
                    st.append(w)
                    comp.append(w)
        out.append(np.array(comp))
    return out


def partition_stats(F, gr, gc, ctx) -> dict:
    """F: flow matrix used for shares (persons or links), gr/gc: city / destination labels."""
    from sklearn.metrics import normalized_mutual_info_score as nmi

    cols = ctx["cols"]                      # node index of each destination column
    adj = ctx["adj_list"]
    city_idx = ctx["city_idx"]
    n_terr = n_nc = 0
    tot_in = dia = dia_km = 0.0
    for g in np.unique(gr):
        cc = np.where(gr == g)[0]
        dd = np.where(gc == g)[0]
        if len(dd) < MIN_DEST:
            continue
        n_terr += 1
        comps = components(cols[dd], adj)
        if sum(len(p) >= MIN_PIECE for p in comps) >= 2:
            n_nc += 1
        home_zone = set(int(z) for z in city_idx[cc])
        touch = set(home_zone)
        for z in home_zone:
            touch.update(adj[z])
        col_of = {int(z): k for k, z in zip(dd, cols[dd])}
        dia_cols = [col_of[int(z)] for p in comps if not (set(p.tolist()) & touch) for z in p]
        sub = F[np.ix_(cc, dd)]
        tot_in += sub.sum()
        if dia_cols:
            fd = F[np.ix_(cc, dia_cols)]
            dia += fd.sum()
            dia_km += (fd * ctx["D_km"][np.ix_(cc, dia_cols)]).sum()
    same = gr[:, None] == gc[None, :]
    cross = ctx["cross"]
    rs = ctx["rs_mat"]
    enr = roots_ratio(rs, same, cross)
    perm = ctx["perm"]
    null = np.array([roots_ratio(rs[p], same, cross) for p in perm])
    return dict(
        n_comm_with_city=int(len(np.unique(gr))),
        n_territories=n_terr, n_noncontig=n_nc,
        frac_noncontig=n_nc / n_terr if n_terr else np.nan,
        flow_inside=float(F[same].sum() / F.sum()),
        diaspora_share=dia / tot_in if tot_in else np.nan,
        diaspora_km=dia_km / dia if dia else np.nan,
        roots_enrich=enr,
        roots_perm_mean=float(null.mean()), roots_perm_sd=float(null.std()),
        roots_perm_z=float((enr - null.mean()) / null.std()) if null.std() > 0 else np.nan,
        roots_perm_p=float(((null >= enr).sum() + 1) / (len(null) + 1)),
        nmi_prov=float(nmi(ctx["dest_prov"], gc)),
    )


def roots_ratio(rs, same, cross) -> float:
    a = same & cross
    b = ~same & cross
    return float(rs[a].mean() / rs[b].mean()) if a.any() and b.any() else np.nan


# ---------------------------------------------------------------------------
# Per-night context (cached on disk so workers load it once)
# ---------------------------------------------------------------------------
def build_context(date: str, base: dict) -> dict:
    n = len(base["nodes"])
    O_full = night_matrix(date, base["city_idx"], n, base["nodes"])
    keep = O_full.sum(0) > 0
    cols = np.where(keep)[0]
    O = O_full[:, cols]
    mask = cols[None, :] != base["city_idx"][:, None]
    assert (O[~mask] == 0).all()
    D_km = base["D_all"][:, cols]
    bins = np.digitize(D_km, BIN_EDGES_KM)
    ci, di = np.nonzero(mask)
    bi = bins[ci, di]
    ctx = dict(date=date, O=O, cols=cols, mask=mask, D_km=D_km, ci=ci, di=di, bi=bi)
    return ctx


def save_ctx(ctx: dict, extra: dict) -> Path:
    CACHE.mkdir(parents=True, exist_ok=True)
    p = CACHE / f"ctx_{ctx['date']}.npz"
    np.savez_compressed(p, **{k: v for k, v in {**ctx, **extra}.items() if isinstance(v, np.ndarray)})
    return p


@lru_cache(maxsize=4)
def load_ctx(date: str) -> dict:
    z = np.load(CACHE / f"ctx_{date}.npz", allow_pickle=False)
    ctx = {k: z[k] for k in z.files}
    g = np.load(CACHE / "base.npz", allow_pickle=False)
    adj = g["adj"]
    n = int(g["n"])
    adj_list = [[] for _ in range(n)]
    for a, b in adj:
        adj_list[a].append(int(b))
        adj_list[b].append(int(a))
    ctx["adj_list"] = adj_list
    ctx["city_idx"] = g["city_idx"]
    prov_of = g["prov_code_idx"]          # province index per node
    city_prov = prov_of[g["city_idx"]]
    dest_prov = prov_of[ctx["cols"]]
    ctx["dest_prov"] = dest_prov
    ctx["cross"] = city_prov[:, None] != dest_prov[None, :]
    ctx["rs_mat"] = g["roots_share"][:, dest_prov]
    ctx["perm"] = g["perm"]
    return ctx


def task_seed(date: str, fam: str, draw: int) -> int:
    ss = np.random.SeedSequence([C.SEED, zlib.crc32(date.encode()), FAMILIES.index(fam), draw + 1])
    return int(ss.generate_state(1)[0])


def run_task(args):
    """One (date, family, draw) unit. draw = -1 is the observed network."""
    date, fam, draw = args
    t0 = time.time()
    ctx = load_ctx(date)
    O = ctx["O"]
    rng = np.random.default_rng(task_seed(date, fam, draw))
    if fam == "grav":
        if draw < 0:
            X, E = O, ctx["E_grav"]
        else:
            X = draw_nb(ctx["E_grav"], ctx["var_mu"], ctx["var_v"], ctx["var_n"], rng)
            E, _, _ = fit_gravity(X, ctx["ci"], ctx["di"], ctx["bi"])
        Bm = X - E
        norm = X.sum()
        F = X
    else:
        if draw < 0:
            X, P = ctx["A_bicm"], ctx["P_bicm"]
        else:
            X = (rng.random(ctx["P_bicm"].shape) < ctx["P_bicm"]).astype(float)
            P, _, _ = fit_bicm(X, ctx["mask"])
        Bm = X - P
        norm = X.sum()
        F = X
    mems, quals = leiden_bipartite(Bm, N_RESTARTS, task_seed(date, fam, draw) + 7)
    best = int(np.argmax(quals))
    R = O.shape[0]
    gr, gc = mems[best][:R], mems[best][R:]
    st = partition_stats(F, gr, gc, ctx)
    s_nmi, s_ari = stability(mems)
    out = dict(date=date, family=fam, draw=draw, Q=quals[best] / norm,
               Q_restart_sd=float(np.std(quals / norm)), stab_nmi=s_nmi, stab_ari=s_ari,
               n_links_or_W=float(norm), secs=time.time() - t0, **st)
    memb = dict(date=date, family=fam, draw=draw, city_comm=gr.astype(np.int32),
                dest_comm=gc.astype(np.int32) if draw < 0 else None)
    return out, memb


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def night_table() -> pd.DataFrame:
    """Main + sister nights from config (date convention: file d = night ending on d), plus the
    appendix nights (Christmas Eve, Good Friday) as their own groups. Ordinary sisters = the 13
    Tue/Wed/Thu files of Feb 2024 resolved by b01 (networks.csv)."""
    nets = pd.read_csv(NET / "networks.csv", dtype=str)
    ordn = sorted(nets[(nets.group == "ordinary") & nets.role.isin(["main", "sister"])].date.unique())
    sis = dict(C.SISTER_NIGHTS)
    sis["ordinary"] = ordn
    rows = [(d, g, d == C.MAIN_NIGHTS[g]) for g, ds in sis.items() for d in ds]
    rows += [(d, g, True) for g, d in C.APPENDIX_NIGHTS.items()]
    t = pd.DataFrame(rows, columns=["date", "group", "is_main"])
    assert t.groupby("group").is_main.sum().eq(1).all() and t.date.is_unique
    return t


def prepare(base: dict, dates, sfx: str = "") -> pd.DataFrame:
    rows, vtabs = [], []
    for d in dates:
        ctx = build_context(d, base)
        E, err, it = fit_gravity(ctx["O"], ctx["ci"], ctx["di"], ctx["bi"])
        vt, pearson = variance_function(ctx["O"], E, ctx["ci"], ctx["di"])
        vmu, vv, vn = vt["mean_mu"].values, vt["mean_sq_resid"].values, vt["n_cells"].values
        A = (ctx["O"] >= BICM_W_MIN).astype(float)
        P, berr, bit = fit_bicm(A, ctx["mask"])
        # check the variance model on one draw: dispersion, zero share and the binned variance
        # of the draw around its own refitted gravity, against the observed values
        rng = np.random.default_rng(task_seed(d, "grav", 10 ** 6))
        X = draw_nb(E, vmu, vv, vn, rng)
        Ex, _, _ = fit_gravity(X, ctx["ci"], ctx["di"], ctx["bi"])
        vtx, pear_draw = variance_function(X, Ex, ctx["ci"], ctx["di"])
        inr = (vtx["mean_mu"].values >= vmu.min()) & (vtx["mean_mu"].values <= vmu.max())
        ratio = var_of_mu(vtx["mean_mu"].values[inr], vmu, vv, vn) / vtx["mean_sq_resid"].values[inr]
        m = ctx["mask"]
        rows.append(dict(date=d, n_dest=len(ctx["cols"]), W=ctx["O"].sum(), grav_err=err, grav_iter=it,
                         pearson_disp_obs=pearson, pearson_disp_nbdraw=pear_draw,
                         var_ratio_obs_over_draw_median=float(np.median(ratio)),
                         var_ratio_p25=float(np.quantile(ratio, 0.25)), var_ratio_p75=float(np.quantile(ratio, 0.75)),
                         zero_share_obs=float((ctx["O"][m] == 0).mean()), zero_share_nbdraw=float((X[m] == 0).mean()),
                         W_draw=float(X.sum()), n_links_w10=int(A.sum()), bicm_err=berr, bicm_iter=bit))
        vt["date"] = d
        vtabs.append(vt)
        save_ctx(ctx, dict(E_grav=E, A_bicm=A, P_bicm=P, var_mu=vmu, var_v=vv, var_n=vn))
        r_ = rows[-1]
        print(f"  prepared {d}: dest={r_['n_dest']} W={r_['W']:,.0f} pearson obs {pearson:.1f} / draw {pear_draw:.1f}; "
              f"V_obs/V_draw median {r_['var_ratio_obs_over_draw_median']:.2f} [{r_["var_ratio_p25"]:.2f}, {r_["var_ratio_p75"]:.2f}]; "
              f"zeros {r_['zero_share_obs']:.3f}/{r_['zero_share_nbdraw']:.3f}; links={r_['n_links_w10']}", flush=True)
    pd.concat(vtabs).to_csv(C.RESULTS / f"r3_variance_function_bins{sfx}.csv", index=False)
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true", help="main summer night, 3 draws (smoke test)")
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 2))
    ap.add_argument("--dates", default="", help="comma-separated subset of night_table() dates")
    ap.add_argument("--tag", default="", help="suffix of the output files of a subset run")
    args = ap.parse_args()
    t0 = time.time()
    DIAG.mkdir(parents=True, exist_ok=True)
    nodes, cities, city_idx, provs, roots_share = load_inputs()
    p2i = {p: k for k, p in enumerate(provs)}
    prov_idx = nodes["prov"].map(p2i).values
    X, Y = nodes["x"].values, nodes["y"].values
    D_all = np.hypot(X[city_idx][:, None] - X[None, :], Y[city_idx][:, None] - Y[None, :]) / 1000.0
    adj = polygon_adjacency(nodes)
    # within-province permutations of city rows (cities in single-city provinces stay fixed)
    rng = np.random.default_rng(np.random.SeedSequence([C.SEED, 3]))
    cprov = prov_idx[city_idx]
    groups = [np.where(cprov == p)[0] for p in np.unique(cprov) if (cprov == p).sum() >= 2]
    perm = np.tile(np.arange(len(cities)), (N_PERM_ROOTS, 1))
    for k in range(N_PERM_ROOTS):
        for g in groups:
            perm[k, g] = rng.permutation(g)
    CACHE.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(CACHE / "base.npz", adj=adj, n=np.array(len(nodes)), city_idx=city_idx,
                        prov_code_idx=prov_idx, roots_share=roots_share, perm=perm)
    print(f"[r3] {len(cities)} cities, {len(provs)} provinces, {len(adj)} adjacent zone pairs, "
          f"{len(groups)} multi-city provinces with {sum(len(g) for g in groups)} cities", flush=True)
    base = dict(nodes=nodes, city_idx=city_idx, D_all=D_all)

    nights = night_table()
    if args.quick:
        nights = nights[nights.date == C.MAIN_NIGHTS["summer"]]
    if args.dates:
        nights = nights[nights.date.isin(args.dates.split(","))]
    sfx = "_quick" if args.quick else (f"_{args.tag}" if args.tag else "")
    prep = prepare(base, nights.date.tolist(), sfx)
    prep = nights.merge(prep, on="date")
    prep.to_csv(C.RESULTS / f"r3_null_diagnostics{sfx}.csv", index=False)

    tasks = []
    for _, r in nights.iterrows():
        nd = 3 if args.quick else (N_NULL_MAIN if r.is_main else N_NULL_SISTER)
        for fam in FAMILIES:
            tasks += [(r.date, fam, k) for k in range(-1, nd)]
    print(f"[r3] {len(tasks)} clustering tasks on {args.workers} workers", flush=True)
    res, mems = [], []
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        for k, (o, m) in enumerate(ex.map(run_task, tasks, chunksize=1)):
            res.append(o)
            mems.append(m)
            if k % 50 == 0 or o["draw"] < 0:
                print(f"  [{k + 1}/{len(tasks)}] {o['date']} {o['family']} draw {o['draw']}: Q={o['Q']:.4f} "
                      f"terr={o['n_territories']} nc={o['n_noncontig']} dia={o['diaspora_share']:.3f} "
                      f"enr={o['roots_enrich']:.3f} (perm z {o['roots_perm_z']:.1f}) stab={o['stab_nmi']:.2f} "
                      f"{o['secs']:.0f}s  t={time.time() - t0:.0f}s", flush=True)
    df = pd.DataFrame(res).merge(nights[["date", "group", "is_main"]], on="date")
    df.to_csv(C.RESULTS / f"r3_draws{sfx}.csv", index=False)
    pm = pd.DataFrame([dict(date=m["date"], family=m["family"], draw=m["draw"],
                            city_comm=m["city_comm"].tolist(),
                            dest_comm=None if m["dest_comm"] is None else m["dest_comm"].tolist())
                       for m in mems])
    pm.to_parquet(C.RESULTS / f"r3_partitions{sfx}.parquet", index=False)
    print(f"[r3] done in {time.time() - t0:.0f}s", flush=True)


# ---------------------------------------------------------------------------
# Summaries (run with --summarise after the main run)
# ---------------------------------------------------------------------------
STATS = ["Q", "stab_nmi", "n_territories", "n_noncontig", "frac_noncontig", "diaspora_share",
         "diaspora_km", "flow_inside", "roots_enrich", "nmi_prov"]
GROUP_ORDER = ["ordinary", "summer", "allsaints"]
GROUP_ALL = GROUP_ORDER + list(C.APPENDIX_NIGHTS)


def null_comparison(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (d, fam), g in df.groupby(["date", "family"]):
        o = g[g.draw < 0].iloc[0]
        nl = g[g.draw >= 0]
        for s in STATS:
            v = nl[s].dropna().values
            sd = v.std(ddof=1) if len(v) > 1 else np.nan
            rows.append(dict(date=d, group=o.group, is_main=o.is_main, family=fam, stat=s, obs=o[s],
                             null_mean=v.mean(), null_sd=sd, n_null=len(v),
                             z=(o[s] - v.mean()) / sd if sd and sd > 0 else np.nan,
                             excess=o[s] - v.mean(), ratio=o[s] / v.mean() if v.mean() else np.nan,
                             p_upper=((v >= o[s]).sum() + 1) / (len(v) + 1),
                             p_lower=((v <= o[s]).sum() + 1) / (len(v) + 1)))
        for s in ["roots_perm_z", "roots_perm_p", "roots_perm_mean", "roots_perm_sd"]:
            rows.append(dict(date=d, group=o.group, is_main=o.is_main, family=fam, stat=s, obs=o[s],
                             null_mean=nl[s].mean(), null_sd=nl[s].std(ddof=1), n_null=len(nl)))
    return pd.DataFrame(rows)


def sister_summary(nc: pd.DataFrame) -> pd.DataFrame:
    """Median / range over sister nights of obs, excess and z; main-night extremeness; and
    holiday-vs-ordinary sign agreement of the excess over all (holiday sister, ordinary sister) pairs."""
    rows = []
    for (fam, s), g in nc.groupby(["family", "stat"]):
        if s.startswith("roots_perm_") and s != "roots_perm_z":
            continue
        val = "obs" if s == "roots_perm_z" else "excess"
        ordn = g[g.group == "ordinary"][val].values
        ord_main = g[(g.group == "ordinary") & g.is_main][val].values[0]
        for grp in GROUP_ALL:
            h = g[g.group == grp]
            main = h[h.is_main].iloc[0]
            vals = h[val].values
            ext = (main[val] == vals.max() and (vals == vals.max()).sum() == 1) or \
                  (main[val] == vals.min() and (vals == vals.min()).sum() == 1)
            r = dict(family=fam, stat=s, group=grp, n_nights=len(h), main_obs=main.obs,
                     main_null_mean=main.null_mean, main_null_sd=main.null_sd, main_z=main.z,
                     main_excess=main.excess,
                     obs_median=h.obs.median(), obs_min=h.obs.min(), obs_max=h.obs.max(),
                     excess_median=h.excess.median(), excess_min=h.excess.min(), excess_max=h.excess.max(),
                     z_median=h.z.median(), z_min=h.z.min(), z_max=h.z.max(),
                     main_is_extreme=bool(ext))
            if grp != "ordinary":
                sgn = np.sign(main[val] - ord_main)
                agree = int(sum(np.sign(a - b) == sgn for a in vals for b in ordn))
                gt = int(sum(a > b for a in vals for b in ordn))
                r.update(holiday_minus_ordinary_main=main[val] - ord_main, pairs_same_sign=agree,
                         pairs_holiday_gt_ordinary=gt, pairs_total=len(vals) * len(ordn),
                         holiday_minus_ordinary_median=float(np.median(vals) - np.median(ordn)))
            rows.append(r)
    return pd.DataFrame(rows)


def cross_night_stability(pm: pd.DataFrame, nights: pd.DataFrame) -> pd.DataFrame:
    from sklearn.metrics import adjusted_rand_score as ari, normalized_mutual_info_score as nmi

    obs = pm[pm.draw < 0].merge(nights[["date", "group", "is_main"]], on="date")
    cols = {d: np.load(CACHE / f"ctx_{d}.npz")["cols"] for d in obs.date.unique()}
    rows = []
    for fam, g in obs.groupby("family"):
        mains = {r.group: r for r in g[g.is_main].itertuples()}
        for r in g.itertuples():
            for ref_grp, m in mains.items():
                if r.date == m.date:
                    continue
                a, b = np.array(r.city_comm), np.array(m.city_comm)
                ca, cb = cols[r.date], cols[m.date]
                common, ia, ib = np.intersect1d(ca, cb, return_indices=True)
                da, db = np.array(r.dest_comm)[ia], np.array(m.dest_comm)[ib]
                rows.append(dict(family=fam, date=r.date, group=r.group, ref_main=m.date, ref_group=ref_grp,
                                 nmi_cities=nmi(a, b), ari_cities=ari(a, b),
                                 nmi_dests=nmi(da, db), ari_dests=ari(da, db), n_common_dests=len(common)))
    return pd.DataFrame(rows)


# ---- Alternative third result: holiday city communities vs roots-based city communities ----
def roots_partition(nodes, city_idx, provs, n_bins: int = 20):
    """Communities of the city x birth-province roots network (persons, own province masked)
    under the same doubly-constrained gravity null (20 distance quantile bins, as p93), clustered
    with the same Leiden routine and effort."""
    rb = pd.read_csv(ROOTS_CSV, dtype={"ine_code": str, "prov_birth": str})
    ids = nodes["id"].values
    p2i = {p: k for k, p in enumerate(provs)}
    c2i = {ids[c]: k for k, c in enumerate(city_idx)}
    r = rb[(rb.year == ROOTS_YEAR) & rb.ine_code.isin(c2i) & rb.prov_birth.isin(p2i)]
    R = np.zeros((len(city_idx), len(provs)))
    for c, p, v in zip(r.ine_code, r.prov_birth, r.persons):
        if c[:2] != p:
            R[c2i[c], p2i[p]] = v
    w = nodes["population"].fillna(0).values + 1.0
    px = np.array([np.average(nodes.x[nodes.prov == p], weights=w[nodes.prov == p]) for p in provs])
    py = np.array([np.average(nodes.y[nodes.prov == p], weights=w[nodes.prov == p]) for p in provs])
    D = np.hypot(nodes.x.values[city_idx][:, None] - px[None, :], nodes.y.values[city_idx][:, None] - py[None, :]) / 1e3
    cprov = np.array([ids[c][:2] for c in city_idx])
    mask = cprov[:, None] != np.array(provs)[None, :]
    q = np.quantile(D[mask], np.linspace(0, 1, n_bins + 1))
    bins = np.clip(np.digitize(D, q[1:-1]), 0, n_bins - 1)
    ci, di = np.nonzero(mask)
    E, err, it = fit_gravity(R, ci, di, bins[ci, di])
    mems, quals = leiden_bipartite(R - E, N_RESTARTS, C.SEED + 11)
    best = int(np.argmax(quals))
    s_nmi, _ = stability(mems)
    info = dict(Q_roots=quals[best] / R.sum(), grav_err=err, stab_nmi=s_nmi,
                n_city_comms=len(np.unique(mems[best][:len(city_idx)])))
    return mems[best][:len(city_idx)], mems[best][len(city_idx):], info


def alternative_nmi(pm, nights, g_roots, perm, cprov) -> pd.DataFrame:
    from sklearn.metrics import normalized_mutual_info_score as nmi

    rows = []
    for (d, fam), g in pm.groupby(["date", "family"]):
        o = np.array(g[g.draw < 0].city_comm.iloc[0])
        nl = [np.array(x) for x in g[g.draw >= 0].city_comm]
        obs = nmi(o, g_roots)
        null = np.array([nmi(x, g_roots) for x in nl])
        pr = np.array([nmi(o, g_roots[p]) for p in perm[:200]])
        rows.append(dict(date=d, family=fam, nmi_obs=obs, null_mean=null.mean(), null_sd=null.std(ddof=1),
                         z_null=(obs - null.mean()) / null.std(ddof=1),
                         perm_mean=pr.mean(), perm_sd=pr.std(ddof=1),
                         z_perm=(obs - pr.mean()) / pr.std(ddof=1) if pr.std() > 0 else np.nan,
                         nmi_city_province=nmi(o, cprov), nmi_roots_province=nmi(g_roots, cprov)))
    return pd.DataFrame(rows).merge(nights[["date", "group", "is_main"]], on="date")


# ---- Alternative third result (orchestrator suggestion): two-axis decomposition per city ----
def _axes(O, E, ctx_cols, city_prov_k, dest_prov_k, coastal_cols, logR, n_prov):
    """roots axis: Pearson over the other provinces p of log((O_cp+1)/(E_cp+1)) and log(1+born_cp);
    beach axis: log(sum_coastal O / sum_coastal E)."""
    nc = O.shape[0]
    onehot = np.zeros((O.shape[1], n_prov))
    onehot[np.arange(O.shape[1]), dest_prov_k] = 1.0
    Op, Ep = O @ onehot, E @ onehot
    lr = np.log((Op + 1) / (Ep + 1))
    roots = np.full(nc, np.nan)
    for c in range(nc):
        m = np.arange(n_prov) != city_prov_k[c]
        roots[c] = np.corrcoef(lr[c, m], logR[c, m])[0, 1]
    beach = np.log(O[:, coastal_cols].sum(1) / E[:, coastal_cols].sum(1))
    return roots, beach, lr


def _two_axis_night(args):
    date, n_draw = args
    ctx = load_ctx(date)
    g = np.load(CACHE / "base.npz")
    prov_of = g["prov_code_idx"]
    cpk = prov_of[g["city_idx"]]
    dpk = prov_of[ctx["cols"]]
    coastal = np.load(CACHE / "coastal.npy")[ctx["cols"]]
    logR = np.log1p(np.load(CACHE / "roots_persons.npy"))
    n_prov = int(prov_of.max()) + 1
    O, E = ctx["O"], ctx["E_grav"]
    r_obs, b_obs, lr = _axes(O, E, ctx["cols"], cpk, dpk, coastal, logR, n_prov)
    # neighbour-city check: correlation with the roots profile of the other cities of the same province
    nb = np.full(len(cpk), np.nan)
    for c in range(len(cpk)):
        mates = np.where((cpk == cpk[c]) & (np.arange(len(cpk)) != c))[0]
        if len(mates):
            m = np.arange(n_prov) != cpk[c]
            nb[c] = np.mean([np.corrcoef(lr[c, m], logR[k, m])[0, 1] for k in mates])
    rn, bn = [], []
    for k in range(n_draw):
        rng = np.random.default_rng(task_seed(date, "grav", 50_000 + k))
        X = draw_nb(E, ctx["var_mu"], ctx["var_v"], ctx["var_n"], rng)
        Ex, _, _ = fit_gravity(X, ctx["ci"], ctx["di"], ctx["bi"])
        a, b, _ = _axes(X, Ex, ctx["cols"], cpk, dpk, coastal, logR, n_prov)
        rn.append(a)
        bn.append(b)
    rn, bn = np.array(rn), np.array(bn)
    return pd.DataFrame(dict(date=date, city_row=np.arange(len(cpk)), roots_axis=r_obs,
                             roots_null_mean=rn.mean(0), roots_null_sd=rn.std(0, ddof=1),
                             roots_neighbour_mean=nb, roots_own_minus_neighbours=r_obs - nb,
                             beach_axis=b_obs, beach_null_mean=bn.mean(0), beach_null_sd=bn.std(0, ddof=1)))


def two_axis(nights, nodes, cities, city_idx, provs, workers: int):
    from scipy.stats import spearmanr

    rb = pd.read_csv(ROOTS_CSV, dtype={"ine_code": str, "prov_birth": str})
    p2i = {p: k for k, p in enumerate(provs)}
    c2i = {c: k for k, c in enumerate(cities)}
    r = rb[(rb.year == ROOTS_YEAR) & rb.ine_code.isin(c2i) & rb.prov_birth.isin(p2i)]
    R = np.zeros((len(cities), len(provs)))
    for c, p, v in zip(r.ine_code, r.prov_birth, r.persons):
        if c[:2] != p:
            R[c2i[c], p2i[p]] = v
    np.save(CACHE / "roots_persons.npy", R)
    np.save(CACHE / "coastal.npy", nodes["coastal_true"].values.astype(bool))
    tasks = [(d, N_NULL_MAIN if m else N_NULL_SISTER) for d, m in zip(nights.date, nights.is_main)]
    with ProcessPoolExecutor(max_workers=workers) as ex:
        out = pd.concat(list(ex.map(_two_axis_night, tasks)), ignore_index=True)
    out["city"] = np.array(cities)[out.city_row]
    out["name"] = nodes["name"].values[city_idx][out.city_row]
    out["roots_z"] = (out.roots_axis - out.roots_null_mean) / out.roots_null_sd
    out["beach_z"] = (out.beach_axis - out.beach_null_mean) / out.beach_null_sd
    out = out.merge(nights, on="date")
    out.to_csv(C.RESULTS / "r3_alt_two_axis_cities.csv", index=False)
    rows = []
    for d, g in out.groupby("date"):
        rho, p = spearmanr(g.roots_axis, g.beach_axis)
        nbv = g.roots_own_minus_neighbours.dropna()
        rows.append(dict(date=d, group=g.group.iloc[0], is_main=g.is_main.iloc[0],
                         roots_axis_median=g.roots_axis.median(), roots_null_median=g.roots_null_mean.median(),
                         roots_z_median=g.roots_z.median(), share_roots_z_gt2=(g.roots_z > 2).mean(),
                         own_minus_nb_median=nbv.median(), share_own_gt_nb=(nbv > 0).mean(), n_nb=len(nbv),
                         beach_axis_median=g.beach_axis.median(), beach_z_median=g.beach_z.median(),
                         share_beach_z_gt2=(g.beach_z > 2).mean(),
                         spearman_roots_beach=rho, spearman_p=p))
    sm = pd.DataFrame(rows)
    sm.to_csv(C.RESULTS / "r3_alt_two_axis_summary.csv", index=False)
    return sm


def two_axis_contrast(nights, n_thin: int = 20) -> pd.DataFrame:
    """Holiday vs ordinary on the per-city roots axis.

    (1) Paired over cities: d_c = roots_axis(holiday night) - roots_axis(ordinary night), for every
        (holiday sister, ordinary sister) pair: median d, share of cities with d > 0, Wilcoxon p.
    (2) Volume check (attenuation): the ordinary night carries fewer persons, so its residuals are
        noisier and the correlation is attenuated. Each holiday MAIN night is thinned to the
        ordinary main-night total (binomial thinning of round(O) with keep-probability
        W_ordinary / W_holiday), gravity refitted, roots axis recomputed (n_thin replicates).
    """
    from scipy.stats import wilcoxon

    ta = pd.read_csv(C.RESULTS / "r3_alt_two_axis_cities.csv", dtype={"date": str})
    piv = ta.pivot(index="city_row", columns="date", values="roots_axis")
    ordn = nights[nights.group == "ordinary"].date.tolist()
    rows = []
    for grp in [g for g in GROUP_ALL if g != "ordinary"]:
        for h in nights[nights.group == grp].date:
            for o in ordn:
                d = (piv[h] - piv[o]).dropna()
                rows.append(dict(group=grp, holiday=h, ordinary=o, median_diff=d.median(),
                                 share_cities_pos=(d > 0).mean(), wilcoxon_p=wilcoxon(d).pvalue))
    pairs = pd.DataFrame(rows)
    pairs.to_csv(C.RESULTS / "r3_alt_two_axis_pairs.csv", index=False)
    g = np.load(CACHE / "base.npz")
    prov_of = g["prov_code_idx"]
    cpk = prov_of[g["city_idx"]]
    logR = np.log1p(np.load(CACHE / "roots_persons.npy"))
    coastal_all = np.load(CACHE / "coastal.npy")
    n_prov = int(prov_of.max()) + 1
    o_main = C.MAIN_NIGHTS["ordinary"]
    W_o = load_ctx(o_main)["O"].sum()
    thin = []
    for grp in [x for x in GROUP_ALL if x != "ordinary"]:
        d = C.MAIN_NIGHTS.get(grp, C.APPENDIX_NIGHTS.get(grp))
        ctx = load_ctx(d)
        O = ctx["O"]
        keep = W_o / O.sum()
        for k in range(n_thin):
            rng = np.random.default_rng(task_seed(d, "grav", 90_000 + k))
            X = rng.binomial(np.round(O).astype(np.int64), keep).astype(float)
            Ex, _, _ = fit_gravity(X, ctx["ci"], ctx["di"], ctx["bi"])
            r, _, _ = _axes(X, Ex, ctx["cols"], cpk, prov_of[ctx["cols"]], coastal_all[ctx["cols"]], logR, n_prov)
            thin.append(dict(group=grp, date=d, rep=k, keep_prob=keep, W_thinned=X.sum(),
                             roots_axis_median=np.nanmedian(r),
                             share_cities_above_ordinary=float(np.nanmean(r > piv[o_main].values))))
    th = pd.DataFrame(thin)
    th.to_csv(C.RESULTS / "r3_alt_two_axis_thinned.csv", index=False)
    return pairs, th


# ---- Diagnostic maps and plots ----
def diag_plots(df, pm, nights, nodes, city_idx, nc):
    import geopandas as gpd
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    g = gpd.read_file(C.ZONES_GPKG).set_index("id").loc[nodes["id"]].reset_index()
    cmap = matplotlib.colormaps["tab20"]
    for fam in FAMILIES:
        fig, axes = plt.subplots(1, 3, figsize=(18, 6.5))
        for ax, grp in zip(axes, GROUP_ORDER):
            d = C.MAIN_NIGHTS[grp]
            r = pm[(pm.date == d) & (pm.family == fam) & (pm.draw < 0)].iloc[0]
            cols = np.load(CACHE / f"ctx_{d}.npz")["cols"]
            gc, gr = np.array(r.dest_comm), np.array(r.city_comm)
            lab = np.full(len(nodes), -1)
            lab[cols] = gc
            sizes = pd.Series(gc).value_counts()
            top = [c for c in sizes.index if c in set(gr)][:20]
            colour = {c: cmap(k % 20) for k, c in enumerate(top)}
            fc = [colour.get(l, (0.85, 0.85, 0.85, 1)) if l >= 0 else (1, 1, 1, 1) for l in lab]
            g.plot(ax=ax, color=fc, edgecolor="none")
            ax.scatter(nodes.x.values[city_idx], nodes.y.values[city_idx], s=14,
                       c=[colour.get(c, (0.5, 0.5, 0.5, 1)) for c in gr], edgecolors="k", linewidths=0.6, zorder=3)
            o = nc[(nc.date == d) & (nc.family == fam)].set_index("stat")
            ax.set_title(f"{grp} {d} [{fam}]\nQ={o.loc['Q','obs']:.3f} (null {o.loc['Q','null_mean']:.3f}±{o.loc['Q','null_sd']:.3f}); "
                         f"non-contig {int(o.loc['n_noncontig','obs'])}/{int(o.loc['n_territories','obs'])}; "
                         f"roots enr {o.loc['roots_enrich','obs']:.2f} (null {o.loc['roots_enrich','null_mean']:.2f})",
                         fontsize=8.5)
            ax.set_axis_off()
        fig.suptitle(f"R3 bipartite communities, {fam} null: destination zones coloured by community, cities = dots "
                     f"(grey = communities without a city or outside the 20 largest)", fontsize=10)
        fig.tight_layout()
        fig.savefig(DIAG / f"map_main_nights_{fam}.png", dpi=130)
        plt.close(fig)
    show = ["Q", "n_noncontig", "frac_noncontig", "diaspora_share", "roots_enrich", "stab_nmi"]
    for fam in FAMILIES:
        fig, axes = plt.subplots(3, len(show), figsize=(3.1 * len(show), 8))
        for i, grp in enumerate(GROUP_ORDER):
            d = C.MAIN_NIGHTS[grp]
            sub = df[(df.date == d) & (df.family == fam)]
            for j, s in enumerate(show):
                ax = axes[i, j]
                ax.hist(sub[sub.draw >= 0][s].dropna(), bins=20, color="#94a3b8")
                ax.axvline(sub[sub.draw < 0][s].iloc[0], color="#dc2626", lw=2)
                ax.set_title(f"{grp}: {s}", fontsize=8)
                ax.tick_params(labelsize=7)
        fig.suptitle(f"R3 {fam}: null distribution (grey, {N_NULL_MAIN} draws) vs observed (red), main nights", fontsize=10)
        fig.tight_layout()
        fig.savefig(DIAG / f"null_hist_{fam}.png", dpi=120)
        plt.close(fig)
    vt = _concat_parts("r3_variance_function_bins", ".csv")
    fig, ax = plt.subplots(figsize=(6, 5))
    for k, grp in enumerate(GROUP_ORDER):
        t = vt[vt.date == C.MAIN_NIGHTS[grp]]
        ax.loglog(t.mean_mu, t.mean_sq_resid, "o-", ms=3, label=f"{grp} {C.MAIN_NIGHTS[grp]}")
    xx = np.logspace(-3, 4)
    ax.loglog(xx, xx, "k--", lw=1, label="Poisson V = mu")
    ax.set_xlabel("gravity mean mu (persons)")
    ax.set_ylabel("mean (O - mu)^2")
    ax.legend(fontsize=8)
    ax.set_title("Empirical mean-variance relation used by the NB null", fontsize=9)
    fig.tight_layout()
    fig.savefig(DIAG / "variance_function.png", dpi=120)
    plt.close(fig)


def _concat_parts(stem: str, ext: str) -> pd.DataFrame:
    """Main run file + subset-run files (stem_<tag>ext); quick-test files excluded."""
    fs = [f for f in sorted(C.RESULTS.glob(f"{stem}*{ext}")) if "_quick" not in f.name and "_all" not in f.name]
    rd = pd.read_parquet if ext == ".parquet" else pd.read_csv
    return pd.concat([rd(f) for f in fs], ignore_index=True)


def summarise():
    nights = night_table()
    keep = set(nights.date)
    df = _concat_parts("r3_draws", ".csv")
    dropped = sorted(set(df.date) - keep)
    print("dates run but no longer in the night table (dropped):", dropped)
    df = df[df.date.isin(keep)].drop(columns=["group", "is_main"]).merge(nights, on="date")
    assert not df.duplicated(["date", "family", "draw"]).any()
    assert set(df.date) == keep, keep - set(df.date)
    pm = _concat_parts("r3_partitions", ".parquet")
    pm = pm[pm.date.isin(keep)]
    diag = _concat_parts("r3_null_diagnostics", ".csv")
    diag = diag[diag.date.isin(keep)].drop(columns=["group", "is_main", "roles"], errors="ignore").merge(nights, on="date")
    diag.to_csv(C.RESULTS / "r3_null_diagnostics_all.csv", index=False)
    nc = null_comparison(df)
    nc.to_csv(C.RESULTS / "r3_null_comparison.csv", index=False)
    ss = sister_summary(nc)
    ss.to_csv(C.RESULTS / "r3_sister_summary.csv", index=False)
    st = cross_night_stability(pm, nights)
    st.to_csv(C.RESULTS / "r3_partition_stability.csv", index=False)
    nodes, cities, city_idx, provs, _ = load_inputs()
    g_roots, g_roots_prov, info = roots_partition(nodes, city_idx, provs)
    base = np.load(CACHE / "base.npz")
    cprov = np.array([c[:2] for c in cities])
    alt = alternative_nmi(pm, nights, g_roots, base["perm"], cprov)
    alt.to_csv(C.RESULTS / "r3_alt_roots_partition_nmi.csv", index=False)
    pd.DataFrame({"city": cities, "name": nodes["name"].values[city_idx], "roots_comm": g_roots}).to_csv(
        C.RESULTS / "r3_alt_roots_partition.csv", index=False)
    pd.DataFrame({"prov": provs, "roots_comm": g_roots_prov}).to_csv(C.RESULTS / "r3_alt_roots_partition_provinces.csv", index=False)
    print("roots partition:", info)
    ta = two_axis(nights, nodes, cities, city_idx, provs, max(1, (os.cpu_count() or 2) - 2))
    pairs, th = two_axis_contrast(nights)
    print(pairs.groupby("group")[["median_diff", "share_cities_pos", "wilcoxon_p"]].median().to_string())
    print(th.groupby("group")[["keep_prob", "roots_axis_median", "share_cities_above_ordinary"]].agg(["mean", "std"]).to_string())
    diag_plots(df, pm, nights, nodes, city_idx, nc)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 30)
    m = nc[nc.is_main & nc.stat.isin(STATS)][["family", "stat", "group", "obs", "null_mean", "null_sd", "z", "p_upper", "p_lower"]]
    print(m.sort_values(["family", "stat", "group"]).to_string(index=False, float_format=lambda v: f"{v:.4g}"))
    print(ss.to_string(index=False, float_format=lambda v: f"{v:.4g}"))
    print(st.groupby(["family", "group", "ref_group"])[["nmi_cities", "nmi_dests"]].agg(["median", "min", "max"]).to_string())
    print(ta.to_string(index=False, float_format=lambda v: f"{v:.3g}"))
    print(alt.groupby(["family", "group"])[["nmi_obs", "null_mean", "z_null", "perm_mean", "z_perm", "nmi_city_province"]].median().to_string())


if __name__ == "__main__":
    if "--summarise" in sys.argv:
        summarise()
    else:
        main()
