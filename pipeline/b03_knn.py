"""b03_knn: average nearest-neighbour degree knn(k) of the undirected projection (Zhou & Mondragon 2007, Fig. 2b).

For the three main nights (w10 backbone): real knn(k), and the mean over N_KNN_REWIRE degree-preserving
rewirings (igraph edge swaps, REWIRE_SWAPS_PER_EDGE x L trials, the rew_und null of b02).
knn(k) is the mean over nodes of degree k of the mean degree of their neighbours; it is also reported
in logarithmic degree bins (mean over the nodes in the bin) for plotting.
Writes results/b03_knn.csv. Run: .venv/bin/python -I pipeline/b03_knn.py
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

NET = C.RESULTS / "b01_networks"
N_KNN_REWIRE = 10
BINS = np.unique(np.round(np.geomspace(1, 2500, 25)).astype(int))


def knn_per_node(G: ig.Graph) -> tuple[np.ndarray, np.ndarray]:
    deg = np.array(G.degree())
    knn = np.array(G.knn()[0], dtype=float)  # mean neighbour degree, NaN for isolated nodes
    ok = deg > 0
    return deg[ok], knn[ok]


def binned(deg, knn):
    idx = np.digitize(deg, BINS) - 1
    out = []
    for b in range(len(BINS) - 1):
        m = idx == b
        if m.any():
            out.append((BINS[b], BINS[b + 1], float(np.exp(np.log(deg[m]).mean())), float(knn[m].mean()), int(m.sum())))
    return out


rows = []
nodes = pd.read_parquet(NET / "nodes.parquet")
n = len(nodes)
for group, d in C.MAIN_NIGHTS.items():
    e = pd.read_parquet(NET / f"edges_{d}__w10.parquet")
    G = B.to_igraph(e, n, directed=False)
    deg, knn = knn_per_node(G)
    for lo, hi, kc, kv, nn in binned(deg, knn):
        rows.append(dict(date=d, group=group, series="real", draw=-1, k_lo=lo, k_hi=hi, k_center=kc, knn=kv, n_nodes=nn))
    for s in range(N_KNN_REWIRE):
        H = G.copy()
        ig.set_random_number_generator(random.Random(C.SEED + s))  # same seeding as b02.rewire
        H.rewire(n=C.REWIRE_SWAPS_PER_EDGE * H.ecount(), allowed_edge_types="simple")
        assert H.degree() == G.degree()
        dh, kh = knn_per_node(H)
        for lo, hi, kc, kv, nn in binned(dh, kh):
            rows.append(dict(date=d, group=group, series="rewired", draw=s, k_lo=lo, k_hi=hi, k_center=kc, knn=kv, n_nodes=nn))
    print(group, "done", flush=True)

ig.set_random_number_generator(random)
pd.DataFrame(rows).to_csv(C.RESULTS / "b03_knn.csv", index=False)
print("wrote results/b03_knn.csv")
