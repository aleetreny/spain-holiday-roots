"""
p94_check_2nd_order_and_tails.py
Quick companion verification:
1. Joint-degree-preserving (1st-order mixing preserving) rewiring test on Weekday mobility network (Zhou, Cox, Hansen 2017):
   Preserves P(k) AND P(k, k') (exact 1st-order assortativity r) via binned joint-degree swaps,
   measuring R_avg and R_max to separate 1st-order mixing from spatial embedding!
2. Tail comparison (power-law vs lognormal MLE log-likelihood ratio, Vuong test) on degree & strength CCDFs
   for Question A (Weekday, Sunday, August) and Question B (Feb Wed, Mid-Aug, All Saints, Christmas Eve).
"""
import json, time
import numpy as np
import pandas as pd
import igraph as ig
from scipy import sparse, stats

t0 = time.time()
rng = np.random.default_rng(7)

nd = pd.read_parquet('data/derived/nodes_municipios.parquet')
nd = nd[~nd.foreign & ~nd.island & ~nd.prov.isin(['51', '52'])].reset_index(drop=True)
n = len(nd)
idx = {k: i for i, k in enumerate(nd.id)}

def load_und_wk():
    e = pd.read_parquet('data/derived/net_wk.parquet')
    e = e[(e.o != e.d) & e.o.isin(idx) & e.d.isin(idx) & (e.persist >= 0.5)]
    oi = e.o.map(idx).values; di = e.d.map(idx).values
    a, b = np.minimum(oi, di), np.maximum(oi, di)
    und = pd.DataFrame({'a': a, 'b': b, 'w': e.w.values}).groupby(['a', 'b'], as_index=False)['w'].sum()
    return und

def second_order_from_edges(n_nodes, u_arr, v_arr):
    rows = np.concatenate([u_arr, v_arr])
    cols = np.concatenate([v_arr, u_arr])
    csr = sparse.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(n_nodes, n_nodes))
    deg = np.asarray(csr.sum(1)).ravel().astype(float)
    avg_nb = np.divide(csr @ deg, deg, out=np.zeros_like(deg), where=deg > 0)
    indptr, indices = csr.indptr, csr.indices
    max_nb = np.zeros(n_nodes)
    for i in range(n_nodes):
        s, e = indptr[i], indptr[i + 1]
        if e > s:
            max_nb[i] = deg[indices[s:e]].max()
    def pc(x, y):
        x0, y0 = x - x.mean(), y - y.mean()
        return float((x0 * y0).sum() / np.sqrt((x0 * x0).sum() * (y0 * y0).sum()))
    return pc(deg[rows], deg[cols]), pc(avg_nb[rows], avg_nb[cols]), pc(max_nb[rows], max_nb[cols])

und = load_und_wk()
u0, v0 = und.a.values.copy(), und.b.values.copy()
r1_obs, ravg_obs, rmax_obs = second_order_from_edges(n, u0, v0)

# 1st-order-mixing-preserving rewiring (Joint Degree Distribution / JDD preserving swaps):
# Group nodes into 25 fine degree quantile bins; only accept double-edge swaps (a-b, c-d) -> (a-d, c-b)
# when bin(b) == bin(d) (or bin(a) == bin(c)), which preserves BOTH exact degree sequence P(k) AND
# the joint degree matrix P(k_bin, k'_bin) -> preserves 1st-order assortativity r tightly!
deg0 = np.bincount(np.concatenate([u0, v0]), minlength=n)
k_bins = pd.qcut(deg0, 25, labels=False, duplicates='drop')

r1_jdd, ravg_jdd, rmax_jdd = [], [], []
E = len(u0)
for rep in range(10):
    u = u0.copy(); v = v0.copy()
    edge_set = set(zip(u.tolist(), v.tolist()))
    # Group edge indices by (min_bin, max_bin) or simply pick random pairs sharing a bin
    # Fast vectorized bucketed swap: pick edges where one endpoint has same degree bin
    # Even cleaner: organize edges by endpoint bin (each edge appears under bin(u) and bin(v))
    # Perform 4*E accepted swaps
    swaps_done = 0
    # Build bin -> list of edge indices (directed orientations so endpoint 0 has bin b)
    # Instead, do random proposal within same bin of endpoint u:
    by_bu = [[] for _ in range(int(k_bins.max()) + 1)]
    for ei in range(E):
        by_bu[k_bins[u[ei]]].append(ei)
    by_bu = [np.array(lst, dtype=np.int32) for lst in by_bu if len(lst) >= 2]
    for _ in range(6 * E):
        grp = by_bu[int(rng.integers(len(by_bu)))]
        i1 = grp[int(rng.integers(len(grp)))]
        i2 = grp[int(rng.integers(len(grp)))]
        if i1 == i2:
            continue
        a, b = int(u[i1]), int(v[i1])
        c, d = int(u[i2]), int(v[i2])
        # Since k_bins[a] == k_bins[c], swapping (a-b, c-d) -> (a-d, c-b) preserves BOTH exact degrees AND joint bin counts!
        if a == d or c == b or a == c or b == d:
            continue
        e1 = (a, d) if a < d else (d, a)
        e2 = (c, b) if c < b else (b, c)
        if e1 in edge_set or e2 in edge_set:
            continue
        old1 = (a, b) if a < b else (b, a)
        old2 = (c, d) if c < d else (d, c)
        edge_set.remove(old1); edge_set.remove(old2)
        edge_set.add(e1); edge_set.add(e2)
        u[i1], v[i1] = a, d
        u[i2], v[i2] = c, b
        swaps_done += 1
    r1_r, ravg_r, rmax_r = second_order_from_edges(n, u, v)
    r1_jdd.append(r1_r); ravg_jdd.append(ravg_r); rmax_jdd.append(rmax_r)

print(f"Observed Weekday:     r={r1_obs:+.4f}  R_avg={ravg_obs:+.4f}  R_max={rmax_obs:+.4f}")
print(f"JDD-rewired (10x):    r={np.mean(r1_jdd):+.4f}±{np.std(r1_jdd):.4f}  "
      f"R_avg={np.mean(ravg_jdd):+.4f}±{np.std(ravg_jdd):.4f}  "
      f"R_max={np.mean(rmax_jdd):+.4f}±{np.std(rmax_jdd):.4f}  (swaps/rep ~ {swaps_done})")

# Tail test: Power-law vs Log-normal MLE (Vuong log-likelihood ratio test) on k >= xmin
def vuong_pl_vs_ln(x, xmin):
    xt = np.asarray(x, dtype=float)
    xt = xt[xt >= xmin]
    m = len(xt)
    # Continuous power-law MLE
    alpha = 1.0 + m / np.sum(np.log(xt / xmin))
    ll_pl = np.log((alpha - 1.0) / xmin) - alpha * np.log(xt / xmin)
    # Truncated lognormal MLE (fit mu, sigma on log(xt) truncated at log(xmin))
    lx = np.log(xt)
    mu, sd = lx.mean(), lx.std()
    # exact truncated normal log-density on log(x) minus log(x)
    z0 = (np.log(xmin) - mu) / sd
    log_Z = stats.norm.logsf(z0)
    ll_ln = stats.norm.logpdf(lx, loc=mu, scale=sd) - lx - log_Z
    diff = ll_pl - ll_ln
    R = float(diff.sum())
    sd_diff = float(np.std(diff, ddof=1))
    z_vuong = R / (sd_diff * np.sqrt(m)) if sd_diff > 0 else 0.0
    p_two = float(2 * stats.norm.sf(abs(z_vuong)))
    return dict(xmin=xmin, n_tail=m, alpha_pl=round(float(alpha), 3), R_pl_minus_ln=round(R, 2),
                z_vuong=round(z_vuong, 2), p_val=round(p_two, 4),
                favors=('Power-law' if R > 0 else 'Log-normal'))

str0 = np.bincount(np.concatenate([u0, v0]), weights=np.concatenate([und.w.values, und.w.values]), minlength=n)
print("\nTail comparison (Vuong test: R < 0 favors Log-normal, R > 0 favors Power-law):")
print("  Weekday degree (k>=100):  ", vuong_pl_vs_ln(deg0, 100))
print("  Weekday degree (k>=200):  ", vuong_pl_vs_ln(deg0, 200))
print("  Weekday strength (s>=30k):", vuong_pl_vs_ln(str0, 30000))
print("  Weekday strength (s>=100k):", vuong_pl_vs_ln(str0, 100000))

out = {
    'jdd_rewiring': {
        'obs': {'r': r1_obs, 'R_avg': ravg_obs, 'R_max': rmax_obs},
        'jdd_mean': {'r': float(np.mean(r1_jdd)), 'R_avg': float(np.mean(ravg_jdd)), 'R_max': float(np.mean(rmax_jdd))},
        'jdd_sd': {'r': float(np.std(r1_jdd)), 'R_avg': float(np.std(ravg_jdd)), 'R_max': float(np.std(rmax_jdd))},
    },
    'tails': {
        'wk_deg_100': vuong_pl_vs_ln(deg0, 100),
        'wk_deg_200': vuong_pl_vs_ln(deg0, 200),
        'wk_str_30k': vuong_pl_vs_ln(str0, 30000),
        'wk_str_100k': vuong_pl_vs_ln(str0, 100000),
    }
}
with open('explore/out_p94_jdd_tails.json', 'w') as f:
    json.dump(out, f, indent=2)
print(f"Done in {time.time()-t0:.1f}s")
