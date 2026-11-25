"""b02_zhou: exact implementations of Prof. Shi Zhou's network methods.

Covers four papers:
1. Zhou & Mondragon (2004, IEEE Comm. Lett., arXiv:cs/0308036):
   Rich-club connectivity phi(r) vs normalised rank r/N, node-node link
   distribution l(r_i <= 5%, r_j) and l(r_i <= 5%, r_j <= 5%), rich-club clique.
2. Zhou & Mondragon (2007, New J. Phys. 9:173, arXiv:physics/0702096):
   Structural constraints via link rewiring preserving P(k) (Case I) and
   P(k) + phi(k) (Case II), with the single-connected-component constraint
   ("a pair of links can be rewired only if the resulted graph remains as a
   single connected component"), plus Table 1 topological properties.
3. Guo, Chen & Zhou (2008, arXiv:0804.0499):
   BOSAM ordering (ascending degree k, then ascending largest neighbour degree
   omega, then ascending largest neighbour index e in the new ordering;
   iterated until the order is stable).
4. Zhou, Cox & Hansen (arXiv:0903.0687):
   Second-order assortative mixing R_avg and R_max using EXCESS degrees (k - 1)
   and neighbours of a EXCLUDING b (N_a \\ {b}), Pearson over both directions of
   each link, jackknife standard deviation sigma, and 100-permutation sequence
   null.

Also includes Clauset-Shalizi-Newman discrete power-law vs log-normal tests.
Imported by pipeline/b02_basic_properties.py.
"""
from __future__ import annotations

import math
import numpy as np

MAX_SWEEPS = 3000


# ---------------------------------------------------------------------------
# Graph helpers & connectivity check for edge swaps
# ---------------------------------------------------------------------------
def _prep(Gu):
    n = Gu.vcount()
    el = [list(e) for e in Gu.get_edgelist()]
    deg = Gu.degree()
    adj = [set() for _ in range(n)]
    for a, b in el:
        adj[a].add(b)
        adj[b].add(a)
    return el, deg, adj


def _to_graph(n, el):
    import igraph as ig
    return ig.Graph(n=n, edges=[tuple(e) for e in el], directed=False)


def _apply(el, adj, i, j, new1, new2):
    a, b = el[i]
    c, d = el[j]
    adj[a].discard(b); adj[b].discard(a)
    adj[c].discard(d); adj[d].discard(c)
    el[i] = list(new1)
    el[j] = list(new2)
    adj[new1[0]].add(new1[1]); adj[new1[1]].add(new1[0])
    adj[new2[0]].add(new2[1]); adj[new2[1]].add(new2[0])


def _n_nontrivial_components(Gu) -> int:
    """Count connected components with >= 2 nodes (isolated nodes have no edges)."""
    return sum(1 for sz in Gu.connected_components().sizes() if sz > 1)


def _assert_single_edge_component(Gu) -> None:
    nc = _n_nontrivial_components(Gu)
    if nc > 1:
        raise ValueError(
            f"keep_connected=True requires all non-isolated nodes to form one connected component "
            f"(found {nc} components with >= 2 nodes)"
        )


def _stays_connected(adj, n1: tuple[int, int], n2: tuple[int, int]) -> bool:
    """Exact connectivity check after replacing (a, b), (c, d) with n1=(u1, v1), n2=(u2, v2).

    When the pre-swap edge-induced subgraph is connected, every vertex in the post-swap
    graph remains connected to either {u1, v1} or {u2, v2}. Thus the graph stays one
    connected component iff there is a path between {u1, v1} and {u2, v2} in adj.
    """
    u1, v1 = n1
    u2, v2 = n2
    # Fast 2-hop common-neighbour checks (succeeds for >99% of swaps in clustered/hub graphs).
    au1, av1, au2, av2 = adj[u1], adj[v1], adj[u2], adj[v2]
    if not au1.isdisjoint(au2) or not av1.isdisjoint(av2):
        return True
    if not au1.isdisjoint(av2) or not av1.isdisjoint(au2):
        return True
    # Bidirectional BFS between {u1, v1} and {u2, v2}.
    vis1 = {u1, v1}
    vis2 = {u2, v2}
    f1 = [u1, v1]
    f2 = [u2, v2]
    while f1 and f2:
        if len(f1) <= len(f2):
            nxt = []
            for u in f1:
                for w in adj[u]:
                    if w in vis2:
                        return True
                    if w not in vis1:
                        vis1.add(w)
                        nxt.append(w)
            f1 = nxt
        else:
            nxt = []
            for u in f2:
                for w in adj[u]:
                    if w in vis1:
                        return True
                    if w not in vis2:
                        vis2.add(w)
                        nxt.append(w)
            f2 = nxt
    return False


# ---------------------------------------------------------------------------
# Rewiring algorithms (Zhou & Mondragon 2007, NJP 9:173)
# ---------------------------------------------------------------------------
def r_from_product_sum(S: float, deg, m: int) -> float:
    """Undirected degree assortativity as an affine function of S = sum over edges of k_u k_v
    (Newman 2002, eq. 4); the other two sums depend only on the degree sequence."""
    k = np.asarray(deg, float)
    a = (k * k).sum() / (2 * m)          # (1/m) sum_e (j + k) / 2
    b = (k ** 3).sum() / (2 * m)         # (1/m) sum_e (j^2 + k^2) / 2
    return (S / m - a * a) / (b - a * a)


def random_rewire(Gu, seed: int, swaps_per_edge: int = 10, keep_connected: bool = True):
    """Degree-preserving maximal-random rewiring (Zhou & Mondragon 2007 Case I).

    When ``keep_connected=True``, a pair of links is rewired only if all non-isolated
    nodes remain in a single connected component (Zhou 2007 Sec. 2.1).
    Processes trials in chunks of ``10 * m`` to keep peak memory small even at 1000 x L.
    Returns ``(graph, info)``.
    """
    if keep_connected:
        _assert_single_edge_component(Gu)
    rng = np.random.default_rng(seed)
    el, deg, adj = _prep(Gu)
    m = len(el)
    T = int(swaps_per_edge * m)
    if m < 2 or T == 0:
        return _to_graph(Gu.vcount(), el), dict(
            swaps_per_edge=swaps_per_edge, trials=T, accepted=0,
            acceptance_rate=0.0, disconnected_rejected=0, keep_connected=keep_connected,
        )
    chunk = min(T, max(10 * m, 10000))
    acc = 0
    disc_rej = 0
    done = 0
    while done < T:
        sz = min(chunk, T - done)
        I = rng.integers(0, m, sz)
        J = rng.integers(0, m, sz)
        R = rng.random(sz) < 0.5
        for i, j, r in zip(I.tolist(), J.tolist(), R.tolist()):
            if i == j:
                continue
            a, b = el[i]
            c, d = el[j]
            if a == c or a == d or b == c or b == d:
                continue
            n1, n2 = ((a, c), (b, d)) if r else ((a, d), (b, c))
            if n1[1] in adj[n1[0]] or n2[1] in adj[n2[0]]:
                continue
            _apply(el, adj, i, j, n1, n2)
            if keep_connected and not _stays_connected(adj, n1, n2):
                _apply(el, adj, i, j, (a, b), (c, d))
                disc_rej += 1
                continue
            acc += 1
        done += sz
    H = _to_graph(Gu.vcount(), el)
    return H, dict(
        swaps_per_edge=swaps_per_edge, trials=T, accepted=acc,
        acceptance_rate=acc / max(T, 1), disconnected_rejected=disc_rej,
        keep_connected=keep_connected,
    )


def extreme_rewire(
    Gu,
    mode: str,
    seed: int,
    max_sweeps: int = MAX_SWEEPS,
    window: int = 25,
    tol: float = 1e-5,
    keep_connected: bool = False,
):
    """Greedy maximum-assortative (mode='assortative') or maximum-disassortative rewiring.

    Sweeps of E random pair trials. When ``keep_connected=True`` (Zhou & Mondragon 2007),
    a swap is accepted only if all non-isolated nodes remain in a single connected component.
    Stops when a sweep accepts no move (converged), or when r has moved less than ``tol``
    over the last ``window`` sweeps (plateau), or at ``max_sweeps``.
    Returns (graph, info).
    """
    assert mode in ("assortative", "disassortative")
    if keep_connected:
        _assert_single_edge_component(Gu)
    rng = np.random.default_rng(seed)
    el, deg, adj = _prep(Gu)
    m = len(el)
    S = float(sum(deg[a] * deg[b] for a, b in el))
    hist = [r_from_product_sum(S, deg, m)]
    acc_tot = 0
    disc_rej = 0
    acc = -1
    sweeps = 0
    stop = "max_sweeps"
    while sweeps < max_sweeps:
        sweeps += 1
        acc = 0
        I = rng.integers(0, m, m)
        J = rng.integers(0, m, m)
        for i, j in zip(I.tolist(), J.tolist()):
            if i == j:
                continue
            a, b = el[i]
            c, d = el[j]
            if a == c or a == d or b == c or b == d:
                continue
            s = sorted((a, b, c, d), key=lambda x: -deg[x])
            if mode == "assortative":
                n1, n2 = (s[0], s[1]), (s[2], s[3])
            else:
                n1, n2 = (s[0], s[3]), (s[1], s[2])
            cur = deg[a] * deg[b] + deg[c] * deg[d]
            new = deg[n1[0]] * deg[n1[1]] + deg[n2[0]] * deg[n2[1]]
            if (mode == "assortative" and new <= cur) or (mode == "disassortative" and new >= cur):
                continue
            if n1[1] in adj[n1[0]] or n2[1] in adj[n2[0]]:
                continue
            _apply(el, adj, i, j, n1, n2)
            if keep_connected and not _stays_connected(adj, n1, n2):
                _apply(el, adj, i, j, (a, b), (c, d))
                disc_rej += 1
                continue
            S += new - cur
            acc += 1
        acc_tot += acc
        hist.append(r_from_product_sum(S, deg, m))
        if acc == 0:
            stop = "no_move_in_sweep"
            break
        if len(hist) > window and abs(hist[-1] - hist[-1 - window]) < tol:
            stop = "plateau"
            break
    return _to_graph(Gu.vcount(), el), dict(
        sweeps=sweeps, accepted_total=acc_tot, accepted_last_sweep=acc,
        disconnected_rejected=disc_rej, keep_connected=keep_connected,
        stop=stop, r_change_last_window=float(hist[-1] - hist[max(0, len(hist) - 1 - window)]),
    )


def rich_club_preserving_rewire(Gu, seed: int, swaps_per_edge: int = 10, keep_connected: bool = False):
    """Random degree- and phi(k)-preserving rewiring ("algorithm 2", Zhou 2007 Case II).

    When ``keep_connected=True``, rejects swaps that disconnect the non-isolated component.
    Returns (graph, acceptance rate).
    """
    if keep_connected:
        _assert_single_edge_component(Gu)
    rng = np.random.default_rng(seed)
    el, deg, adj = _prep(Gu)
    m = len(el)
    T = int(swaps_per_edge * m)
    I = rng.integers(0, m, T)
    J = rng.integers(0, m, T)
    R = rng.random(T) < 0.5
    acc = 0
    for i, j, r in zip(I.tolist(), J.tolist(), R.tolist()):
        if i == j:
            continue
        a, b = el[i]
        c, d = el[j]
        if a == c or a == d or b == c or b == d:
            continue
        # the two alternative wirings of {a, b, c, d}
        n1, n2 = ((a, c), (b, d)) if r else ((a, d), (b, c))
        old = sorted((min(deg[a], deg[b]), min(deg[c], deg[d])))
        new = sorted((min(deg[n1[0]], deg[n1[1]]), min(deg[n2[0]], deg[n2[1]])))
        if old != new:
            continue
        if n1[1] in adj[n1[0]] or n2[1] in adj[n2[0]]:
            continue
        _apply(el, adj, i, j, n1, n2)
        if keep_connected and not _stays_connected(adj, n1, n2):
            _apply(el, adj, i, j, (a, b), (c, d))
            continue
        acc += 1
    return _to_graph(Gu.vcount(), el), acc / max(T, 1)


# ---------------------------------------------------------------------------
# Second-order assortative mixing (Zhou, Cox & Hansen, arXiv:0903.0687)
# ---------------------------------------------------------------------------
def _pearson_sym_pairs(x: np.ndarray, y: np.ndarray) -> float:
    """Eq. (1) of Zhou, Cox & Hansen (0903.0687): Pearson over both directions of each link."""
    L = len(x)
    if L == 0:
        return float("nan")
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    mean = 0.5 * float((x + y).sum()) / L
    num = float((x * y).sum()) / L - mean * mean
    den = 0.5 * float((x * x + y * y).sum()) / L - mean * mean
    return float(num / den) if den > 1e-15 else float("nan")


def second_order_sequences(Gu) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Per-link endpoint second-order excess-degree summaries (Zhou, Cox & Hansen Sec 2.1).

    For each undirected link i = (a, b):
      - Excess degree of any node u is e_u = k_u - 1.
      - Neighbours of a excluding b: N_a \\ {b} (size k_a - 1).
      - K_a^avg = (1 / (k_a - 1)) * sum_{u in N_a \\ {b}} (k_u - 1), or 0 if k_a == 1.
      - K_a^max = max_{u in N_a \\ {b}} (k_u - 1), or 0 if k_a == 1.
      - Symmetrically for b excluding a.

    Computed in O(L) time. Returns (Ka_avg, Kb_avg, Ka_max, Kb_max), each length L.
    """
    el = np.array(Gu.get_edgelist(), dtype=np.int64).reshape(-1, 2)
    L = len(el)
    if L == 0:
        z = np.empty(0, dtype=float)
        return z, z, z, z
    n = Gu.vcount()
    deg = np.array(Gu.degree(), dtype=np.int64)
    ex = np.maximum(deg - 1, 0).astype(float)

    # Sum of neighbour excess degrees per node (for R_avg).
    s = np.zeros(n, dtype=float)
    np.add.at(s, el[:, 0], ex[el[:, 1]])
    np.add.at(s, el[:, 1], ex[el[:, 0]])
    u, v = el[:, 0], el[:, 1]
    du, dv = deg[u], deg[v]
    Ka_avg = np.where(du > 1, (s[u] - ex[v]) / np.maximum(du - 1, 1), 0.0)
    Kb_avg = np.where(dv > 1, (s[v] - ex[u]) / np.maximum(dv - 1, 1), 0.0)

    # Top-2 neighbour excess degrees per node (for R_max in O(L)).
    m1 = np.full(n, -1.0, dtype=float)
    m2 = np.full(n, -1.0, dtype=float)
    c1 = np.zeros(n, dtype=np.int64)
    for a, b in el.tolist():
        eb = ex[b]
        if eb > m1[a]:
            m2[a] = m1[a]; m1[a] = eb; c1[a] = 1
        elif eb == m1[a]:
            c1[a] += 1
        elif eb > m2[a]:
            m2[a] = eb
        ea = ex[a]
        if ea > m1[b]:
            m2[b] = m1[b]; m1[b] = ea; c1[b] = 1
        elif ea == m1[b]:
            c1[b] += 1
        elif ea > m2[b]:
            m2[b] = ea

    Ka_max = np.where(du <= 1, 0.0, np.where((ex[v] == m1[u]) & (c1[u] == 1), m2[u], m1[u]))
    Kb_max = np.where(dv <= 1, 0.0, np.where((ex[u] == m1[v]) & (c1[v] == 1), m2[v], m1[v]))
    return Ka_avg, Kb_avg, Ka_max, Kb_max


def second_order_r_avg(Gu) -> float:
    """Exact Zhou, Cox & Hansen (0903.0687) R_avg (excess degrees, excluding other endpoint)."""
    Ka_avg, Kb_avg, _, _ = second_order_sequences(Gu)
    return _pearson_sym_pairs(Ka_avg, Kb_avg)


def second_order_r_max(Gu) -> float:
    """Exact Zhou, Cox & Hansen (0903.0687) R_max (excess degrees, excluding other endpoint)."""
    _, _, Ka_max, Kb_max = second_order_sequences(Gu)
    return _pearson_sym_pairs(Ka_max, Kb_max)


def second_order_mixing(Gu) -> tuple[float, float]:
    """Compute (R_avg, R_max) in a single pass over edges."""
    Ka_avg, Kb_avg, Ka_max, Kb_max = second_order_sequences(Gu)
    return _pearson_sym_pairs(Ka_avg, Kb_avg), _pearson_sym_pairs(Ka_max, Kb_max)


def jackknife_sigma_pairs(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    """Exact O(L) leave-one-link-out jackknife for symmetrised link-pair correlation.

    Following Newman (2003, PRE 67:026126 eq. 27) and Zhou, Cox & Hansen (0903.0687
    Sec. 2.2.1): sigma^2 = sum_{i=1}^L (r_i - r)^2, where r_i is the correlation after
    removing the i-th link pair (x_i, y_i) (both directions) from the link sequence.
    Returns (r, sigma).
    """
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    L = len(x)
    r = _pearson_sym_pairs(x, y)
    if L < 3 or math.isnan(r):
        return r, float("nan")
    Sxy = float((x * y).sum())
    S1 = float((x + y).sum())
    S2 = float((x * x + y * y).sum())
    Lp = L - 1.0
    Sxy_i = Sxy - x * y
    S1_i = S1 - (x + y)
    S2_i = S2 - (x * x + y * y)
    mean_i = 0.5 * S1_i / Lp
    num_i = Sxy_i / Lp - mean_i * mean_i
    den_i = 0.5 * S2_i / Lp - mean_i * mean_i
    with np.errstate(invalid="ignore", divide="ignore"):
        ri = np.where(den_i > 1e-15, num_i / den_i, np.nan)
    if np.isnan(ri).any():
        return r, float("nan")
    sigma = float(np.sqrt(np.sum((ri - r) ** 2)))
    return r, sigma


def permutation_null_pairs(x: np.ndarray, y: np.ndarray, seed: int, n_perm: int = 100) -> dict:
    """Zhou, Cox & Hansen (0903.0687 Sec. 2.2.2) sequence-permutation null.

    Builds the two directed-link sequences A = [x, y] and B = [y, x] (length 2L each),
    randomly permutes B ``n_perm`` times (default 100), and computes Pearson(A, B_perm).
    Returns dict(real, mean, sd, z, n_perm).
    """
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    obs = _pearson_sym_pairs(x, y)
    if len(x) < 2 or math.isnan(obs):
        return dict(real=obs, mean=float("nan"), sd=float("nan"), z=float("nan"), n_perm=0)
    A = np.concatenate([x, y])
    B = np.concatenate([y, x])
    A0 = A - A.mean()
    B0 = B - B.mean()
    denom = float(np.sqrt((A0 * A0).sum() * (B0 * B0).sum()))
    if denom <= 1e-15:
        return dict(real=obs, mean=float("nan"), sd=float("nan"), z=float("nan"), n_perm=0)
    rng = np.random.default_rng(seed)
    vals = np.empty(n_perm, dtype=float)
    for k in range(n_perm):
        Bp = rng.permutation(B0)
        vals[k] = float((A0 * Bp).sum() / denom)
    mu = float(vals.mean())
    sd = float(vals.std(ddof=1)) if n_perm > 1 else float("nan")
    z = float((obs - mu) / sd) if (not math.isnan(sd) and sd > 1e-12) else float("nan")
    return dict(real=obs, mean=mu, sd=sd, z=z, n_perm=n_perm)


def mixing_significance(Gu, seed: int, n_perm: int = 100) -> dict:
    """Jackknife sigma and 100-permutation null for r_und, R_avg, and R_max."""
    el = np.array(Gu.get_edgelist(), dtype=np.int64).reshape(-1, 2)
    deg = np.array(Gu.degree(), dtype=float)
    ex = np.maximum(deg - 1.0, 0.0)
    r_x = ex[el[:, 0]] if len(el) else np.empty(0, float)
    r_y = ex[el[:, 1]] if len(el) else np.empty(0, float)
    Ka_avg, Kb_avg, Ka_max, Kb_max = second_order_sequences(Gu)

    ss = np.random.SeedSequence(seed)
    s_r, s_avg, s_max = [int(c.generate_state(1)[0]) for c in ss.spawn(3)]

    r_val, r_sig = jackknife_sigma_pairs(r_x, r_y)
    avg_val, avg_sig = jackknife_sigma_pairs(Ka_avg, Kb_avg)
    max_val, max_sig = jackknife_sigma_pairs(Ka_max, Kb_max)
    return {
        "r_und": {
            "real": r_val,
            "jackknife_sigma": r_sig,
            "perm": permutation_null_pairs(r_x, r_y, s_r, n_perm=n_perm),
        },
        "R_avg": {
            "real": avg_val,
            "jackknife_sigma": avg_sig,
            "perm": permutation_null_pairs(Ka_avg, Kb_avg, s_avg, n_perm=n_perm),
        },
        "R_max": {
            "real": max_val,
            "jackknife_sigma": max_sig,
            "perm": permutation_null_pairs(Ka_max, Kb_max, s_max, n_perm=n_perm),
        },
    }


# ---------------------------------------------------------------------------
# BOSAM ordering (Guo, Chen & Zhou 2008, arXiv:0804.0499 Sec. 2)
# ---------------------------------------------------------------------------
def bosam_order(Gu, max_iter: int = 100, return_info: bool = False):
    """Exact BOSAM node order (Guo, Chen & Zhou 2008 Sec. 2).

    Sort keys (all ascending):
      1. node degree k
      2. largest neighbour degree omega = max_{u in N(v)} k_u (-1 if isolated)
      3. largest neighbour index e = max_{u in N(v)} pos[u] IN THE NEW ORDERING (-1 if isolated)

    Iterated until the ordering is stable (capped at ``max_iter``). Ties in (k, omega, e)
    preserve their previous relative order (stable sort, initialised by original node index).
    Returns array of shape (N, 4): (node_idx, degree, max_nb_degree, max_nb_index), in BOSAM order.
    If ``return_info=True``, returns ``(key_array, info_dict)``.
    """
    n = Gu.vcount()
    deg = [int(d) for d in Gu.degree()]
    nbs = [Gu.neighbors(v) for v in range(n)]
    omega = [max((deg[u] for u in nbs[v]), default=-1) for v in range(n)]

    # Initial pass: ascending (k, omega, original node index).
    order = sorted(range(n), key=lambda v: (deg[v], omega[v], v))
    pos = [0] * n
    for p, v in enumerate(order):
        pos[v] = p

    converged = False
    it = 0
    for it in range(1, max_iter + 1):
        e_idx = [max((pos[u] for u in nbs[v]), default=-1) for v in range(n)]
        # Stable sort preserves previous tie order when (k, omega, e) are equal.
        new_order = sorted(order, key=lambda v: (deg[v], omega[v], e_idx[v]))
        if new_order == order:
            converged = True
            break
        order = new_order
        for p, v in enumerate(order):
            pos[v] = p

    # Recompute e_idx and the corresponding neighbour's original label under final pos.
    e_idx = [max((pos[u] for u in nbs[v]), default=-1) for v in range(n)]
    rows = np.array([[v, deg[v], omega[v], e_idx[v]] for v in order], dtype=np.int64)
    if not return_info:
        return rows
    # Original index of the neighbour attaining max_nb_index (for traceability).
    max_nb_label = np.array([order[e_idx[v]] if e_idx[v] >= 0 else -1 for v in order], dtype=np.int64)
    info = dict(iterations=it, converged=converged, max_iter=max_iter, max_nb_label=max_nb_label)
    return rows, info


# ---------------------------------------------------------------------------
# Rich club & node-node link distribution (Zhou & Mondragon 2004, 2007)
# ---------------------------------------------------------------------------
def rich_club_clique_size(Gu) -> int:
    """Largest n >= 2 such that the top-n nodes by degree (ties broken by node index
    ascending) form a clique (phi(n) == 1). Returns 0 if the top-2 nodes are not linked
    (matching Zhou & Mondragon 2007 Table 1).
    """
    n = Gu.vcount()
    if n < 2:
        return 0
    deg = np.array(Gu.degree(), dtype=np.int64)
    order = np.lexsort((np.arange(n), -deg))
    clique_n = 0
    for r in range(2, n + 1):
        v = int(order[r - 1])
        nb_v = set(Gu.neighbors(v))
        if all(int(order[j]) in nb_v for j in range(r - 1)):
            clique_n = r
        else:
            break
    return clique_n


def node_node_link_distribution(Gu, frac: float = 0.05) -> dict:
    """Zhou & Mondragon (2004, IEEE Comm. Lett. Table I) node-node link distribution.

    Nodes are ranked by descending degree (ties broken by node index ascending).
    Top-``frac`` rich nodes = nodes with 1-based normalised rank r/N <= frac
    (i.e. ``top_n = max(1, floor(frac * N))``).
    Returns:
      - top_n, frac
      - L (total undirected links)
      - l_top_any  : l(r_i <= frac, r_j) — links with >= 1 endpoint in top-frac
      - share_l_top_any
      - l_top_both : l(r_i <= frac, r_j <= frac) — links with both endpoints in top-frac
      - share_l_top_both
    """
    n = Gu.vcount()
    L = Gu.ecount()
    top_n = max(1, int(math.floor(frac * n)))
    deg = np.array(Gu.degree(), dtype=np.int64)
    order = np.lexsort((np.arange(n), -deg))
    rank = np.empty(n, dtype=np.int64)
    rank[order] = np.arange(n)
    el = np.array(Gu.get_edgelist(), dtype=np.int64).reshape(-1, 2)
    if L == 0:
        return dict(top_n=top_n, frac=frac, L=0, l_top_any=0, share_l_top_any=float("nan"),
                    l_top_both=0, share_l_top_both=float("nan"))
    r_min = np.minimum(rank[el[:, 0]], rank[el[:, 1]])
    r_max = np.maximum(rank[el[:, 0]], rank[el[:, 1]])
    l_any = int((r_min < top_n).sum())
    l_both = int((r_max < top_n).sum())
    return dict(
        top_n=top_n, frac=frac, L=L,
        l_top_any=l_any, share_l_top_any=l_any / L,
        l_top_both=l_both, share_l_top_both=l_both / L,
    )


def rich_club_rank_grid(n: int, n_log: int = 50, extra_ranks: tuple[int, ...] = (25, 50, 100, 200)) -> np.ndarray:
    """Log-spaced integer ranks r in [2, n] for phi(r) vs normalised rank r/N."""
    base = np.unique(np.round(np.geomspace(2, n, n_log)).astype(int)).tolist()
    pcts = [int(math.floor(f * n)) for f in (0.01, 0.05, 0.10, 0.25, 0.50)]
    grid = sorted({r for r in (base + list(extra_ranks) + pcts + [2, 5, 10, 15, 20, n]) if 2 <= r <= n})
    return np.array(grid, dtype=np.int64)


# ---------------------------------------------------------------------------
# Zhou & Mondragon (2007, NJP 9:173) Table 1 properties
# ---------------------------------------------------------------------------
def _fit_loglog_slope(x: np.ndarray, y: np.ndarray) -> float:
    mask = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
    if int(mask.sum()) < 2:
        return float("nan")
    return float(np.polyfit(np.log(x[mask]), np.log(y[mask]), 1)[0])


def zhou_table_row(Gu, date: str, group: str, backbone: str, seed: int) -> dict:
    """Compute Zhou & Mondragon (2007) Table 1 properties on the giant connected component.

    Also reports Zhou (2004) top-5% link distribution and Zhou, Cox & Hansen (2009)
    second-order mixing R_avg, R_max with jackknife sigmas.
    """
    import b02_basic_properties as P  # local import avoids circularity at module load

    N_total = Gu.vcount()
    gc = Gu.connected_components().giant()
    N = gc.vcount()
    L = gc.ecount()
    deg = np.array(gc.degree(), dtype=np.int64)
    k_mean = float(2.0 * L / N) if N > 0 else float("nan")
    k_max = int(deg.max()) if N > 0 else 0

    # Degree distribution CCDF and PMF for descriptive power-law slope gamma.
    k_u, cc = P.ccdf(deg)
    vals, cnts = np.unique(deg, return_counts=True)
    pmf = cnts / float(N)

    # Literal Zhou 2007 fit range: [k_mean, 40].
    m_cc40 = (k_u >= k_mean) & (k_u <= 40)
    slope_cc40 = _fit_loglog_slope(k_u[m_cc40], cc[m_cc40])
    gamma_ccdf_40 = (1.0 - slope_cc40) if not math.isnan(slope_cc40) else float("nan")

    m_pmf40 = (vals >= k_mean) & (vals <= 40)
    slope_pmf40 = _fit_loglog_slope(vals[m_pmf40], pmf[m_pmf40])
    gamma_pmf_40 = (-slope_pmf40) if not math.isnan(slope_pmf40) else float("nan")

    # Scaled range [k_mean, 6*k_mean] (Zhou's [6, 40] spans ~1x to ~6.3x k_mean) so nights
    # with k_mean >= 40 (All Saints w10, Summer w10) also have a descriptive slope.
    hi_scaled = 6.0 * k_mean
    m_cc_sc = (k_u >= k_mean) & (k_u <= hi_scaled)
    slope_cc_sc = _fit_loglog_slope(k_u[m_cc_sc], cc[m_cc_sc])
    gamma_scaled = (1.0 - slope_cc_sc) if not math.isnan(slope_cc_sc) else float("nan")

    k_cut = float(N ** (1.0 / (gamma_ccdf_40 - 1.0))) if (not math.isnan(gamma_ccdf_40) and gamma_ccdf_40 > 1.0) else float("nan")
    k_cut_pmf = float(N ** (1.0 / (gamma_pmf_40 - 1.0))) if (not math.isnan(gamma_pmf_40) and gamma_pmf_40 > 1.0) else float("nan")
    k_cut_scaled = float(N ** (1.0 / (gamma_scaled - 1.0))) if (not math.isnan(gamma_scaled) and gamma_scaled > 1.0) else float("nan")

    # Rich-club connectivity phi(k > 40) and exponent theta of phi(k) ~ k^theta.
    Nk40, _Ek40, phi40 = P.rich_club_phi(gc, [40])
    ks_all = np.arange(1, max(k_max, 2))
    _, _, phi_all = P.rich_club_phi(gc, ks_all)
    pm40 = (ks_all >= k_mean) & (ks_all <= 40)
    theta_40 = _fit_loglog_slope(ks_all[pm40], phi_all[pm40])
    pm_sc = (ks_all >= k_mean) & (ks_all <= hi_scaled)
    theta_scaled = _fit_loglog_slope(ks_all[pm_sc], phi_all[pm_sc])

    alpha = P.r_undirected(gc)
    n_clique = rich_club_clique_size(gc)
    l_avg, diam = P.path_stats(gc)
    l_star = float(math.log(N) / math.log(k_mean)) if (N > 1 and k_mean > 1.0) else float("nan")
    nld = node_node_link_distribution(gc, frac=0.05)
    sig = mixing_significance(gc, seed=seed, n_perm=100)

    note = "Descriptive slope only (CSN likelihood-ratio tests favour log-normal over power law)"
    if k_mean >= 40.0:
        note += "; mean k >= 40 so [mean k, 40] is empty (see gamma_scaled / theta_scaled on [mean k, 6*mean k])"

    return dict(
        date=date, group=group, backbone=backbone,
        N=N, N_total=N_total, L=L, k_mean=k_mean,
        gamma_kmean_40=gamma_ccdf_40,
        gamma_pmf_kmean_40=gamma_pmf_40,
        gamma_scaled=gamma_scaled,
        k_max=k_max,
        k_cut=k_cut,
        k_cut_pmf=k_cut_pmf,
        k_cut_scaled=k_cut_scaled,
        alpha=alpha,
        alpha_jackknife_sigma=sig["r_und"]["jackknife_sigma"],
        phi_k_gt_40=float(phi40[0]),
        N_k_gt_40=int(Nk40[0]),
        theta_kmean_40=theta_40,
        theta_scaled=theta_scaled,
        n_clique=n_clique,
        l=l_avg,
        l_star=l_star,
        diameter=diam,
        top5_n=nld["top_n"],
        l_top5_any=nld["l_top_any"],
        share_l_top5_any=nld["share_l_top_any"],
        l_top5_both=nld["l_top_both"],
        share_l_top5_both=nld["share_l_top_both"],
        R_avg=sig["R_avg"]["real"],
        R_avg_jackknife_sigma=sig["R_avg"]["jackknife_sigma"],
        R_avg_perm100_mean=sig["R_avg"]["perm"]["mean"],
        R_avg_perm100_sd=sig["R_avg"]["perm"]["sd"],
        R_max=sig["R_max"]["real"],
        R_max_jackknife_sigma=sig["R_max"]["jackknife_sigma"],
        R_max_perm100_mean=sig["R_max"]["perm"]["mean"],
        R_max_perm100_sd=sig["R_max"]["perm"]["sd"],
        gamma_note=note,
    )


# ---------------------------------------------------------------------------
# Power-law vs log-normal (Clauset-Shalizi-Newman)
# ---------------------------------------------------------------------------
def powerlaw_vs_lognormal(values) -> dict:
    """Clauset-Shalizi-Newman fit with the `powerlaw` package: discrete, xmin estimated by KS.

    Zero degrees are excluded (a power law is defined for k >= 1). Likelihood-ratio tests use
    normalized_ratio=True (Vuong); R > 0 favours the power law.
    """
    import powerlaw
    import warnings

    x = np.asarray(values)
    x = x[x > 0]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        f = powerlaw.Fit(x, discrete=True, verbose=False)
        R_ln, p_ln = f.distribution_compare("power_law", "lognormal", normalized_ratio=True)
        R_ex, p_ex = f.distribution_compare("power_law", "exponential", normalized_ratio=True)
        R_tp, p_tp = f.distribution_compare("power_law", "truncated_power_law", normalized_ratio=True)
    return dict(n_pos=int(len(x)), xmin=float(f.xmin), alpha=float(f.power_law.alpha),
                alpha_sigma=float(f.power_law.standard_err), n_tail=int((x >= f.xmin).sum()),
                ks_D=float(f.power_law.D),
                lognormal_mu=float(f.lognormal.mu), lognormal_sigma=float(f.lognormal.sigma),
                R_pl_vs_lognormal=float(R_ln), p_pl_vs_lognormal=float(p_ln),
                R_pl_vs_exponential=float(R_ex), p_pl_vs_exponential=float(p_ex),
                R_pl_vs_truncated=float(R_tp), p_pl_vs_truncated=float(p_tp))
