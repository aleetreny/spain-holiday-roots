"""b04_bosam_envelopes: self-similar prediction of BOSAM leaf envelopes (Guo, Chen & Zhou 2008, Sec. 3, App. III).

For each main night (w10 backbone, undirected projection) and for N_ENV_REWIRE degree-preserving rewirings:
  * BOSAM order from b02_zhou.bosam_order (real networks: the saved b02 order is reused and checked).
  * Envelope of the leaf of degree k: E_k = {(i, e_i) | i in I_k}, e_i = largest neighbour index.
  * Root envelope = leaf of the characteristic degree k* (largest N_k), as in the paper.
  * Predicted envelope for degree k (Eq. 14-15): i~ = N_<k + N_k ((i - N_<l) / N_l)^(k/l), keeping e_i.
  * Error of leaf k = mean over predicted points of |e_i - e_real(i~)| / N, where e_real is the real
    envelope of leaf k interpolated at i~. Only leaves with N_k >= MIN_LEAF are scored.
The theorem holds when neighbour degrees are i.i.d. (same for all nodes); a larger error than in the
rewired networks means that nodes of equal degree connect to systematically different neighbours.
Writes results/b04_bosam_envelopes.csv (per leaf) and results/b04_bosam_envelopes_summary.csv.
Run: .venv/bin/python -I pipeline/b04_bosam_envelopes.py
"""
from pathlib import Path
import random
import sys
import numpy as np
import pandas as pd
import igraph as ig

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C  # noqa: E402
import b00_common as B  # noqa: E402
import b02_zhou as Z  # noqa: E402

NET = C.RESULTS / "b01_networks"
N_ENV_REWIRE = 5
MIN_LEAF = 10


def envelopes(rows: np.ndarray):
    """rows: BOSAM key array (node, degree, omega, e) in BOSAM order. Returns dict k -> (i array, e array)."""
    deg, e = rows[:, 1], rows[:, 3]
    out = {}
    for k in np.unique(deg[deg > 0]):
        i = np.flatnonzero(deg == k)
        out[int(k)] = (i, e[i])
    return out


def score(rows: np.ndarray):
    n = len(rows)
    env = envelopes(rows)
    sizes = {k: len(v[0]) for k, v in env.items()}
    l = max(sizes, key=lambda k: (sizes[k], -k))  # characteristic degree k*
    il, el = env[l]
    n_lt_l, n_l = il.min(), len(il)
    res = []
    for k, (ik, ek) in env.items():
        if k == l or len(ik) < MIN_LEAF:
            continue
        n_lt_k, n_k = ik.min(), len(ik)
        it = n_lt_k + n_k * ((il - n_lt_l + 1) / n_l) ** (k / l) - 1  # positions inside leaf k (0-based)
        e_real = np.interp(it, ik, np.maximum.accumulate(ek))
        res.append(dict(k=k, N_k=n_k, err=float(np.mean(np.abs(el - e_real)) / n)))
    return l, res


rows_out = []
nodes = pd.read_parquet(NET / "nodes.parquet")
n = len(nodes)
for group, d in C.MAIN_NIGHTS.items():
    e = pd.read_parquet(NET / f"edges_{d}__w10.parquet")
    G = B.to_igraph(e, n, directed=False)
    rows = Z.bosam_order(G)
    saved = pd.read_csv(C.RESULTS / f"b02_bosam_order_{group}_{d}__w10.csv")
    assert (saved.idx.values == rows[:, 0]).all(), "BOSAM order differs from b02"
    for s in range(-1, N_ENV_REWIRE):
        if s >= 0:
            H = G.copy()
            ig.set_random_number_generator(random.Random(C.SEED + s))
            H.rewire(n=C.REWIRE_SWAPS_PER_EDGE * H.ecount(), allowed_edge_types="simple")
            rows = Z.bosam_order(H)
        l, res = score(rows)
        for r in res:
            rows_out.append(dict(date=d, group=group, series="real" if s < 0 else "rewired", draw=s, k_star=l, **r))
    print(group, "done", flush=True)
ig.set_random_number_generator(random)

df = pd.DataFrame(rows_out)
df.to_csv(C.RESULTS / "b04_bosam_envelopes.csv", index=False)
# Summary: leaf-size-weighted mean error per network; rewired = mean and sd over draws.
w = df.assign(we=df.err * df.N_k).groupby(["group", "date", "series", "draw", "k_star"]).agg(
    we=("we", "sum"), nk=("N_k", "sum"), leaves=("k", "size")).reset_index()
w["err_weighted"] = w.we / w.nk
summ = w.groupby(["group", "date", "series"]).agg(k_star=("k_star", "first"), leaves=("leaves", "mean"),
                                                 err_mean=("err_weighted", "mean"), err_sd=("err_weighted", "std")).reset_index()
summ.to_csv(C.RESULTS / "b04_bosam_envelopes_summary.csv", index=False)
print(summ.to_string())
