"""Unit tests for the question-B pipeline core (b00_common, b02_basic_properties).

Run: .venv/bin/python -I -m pytest tests/test_b_metrics.py -q
"""
import math
import sys
from pathlib import Path

import igraph as ig
import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pipeline"))
import config as C  # noqa: E402
import b00_common as B  # noqa: E402
import b02_basic_properties as P  # noqa: E402


# ---------------------------------------------------------------- toy graphs
def triangle_pendant():
    # 0-1-2 triangle, 3 hangs on 0
    return ig.Graph(n=4, edges=[(0, 1), (1, 2), (0, 2), (0, 3)])


def k4_with_leaves():
    # K4 on 0..3, each core node has one leaf 4..7
    e = [(a, b) for a in range(4) for b in range(a + 1, 4)] + [(i, i + 4) for i in range(4)]
    return ig.Graph(n=8, edges=e)


# ---------------------------------------------------------------- edge builder
def test_edge_builder_drops_self_loops_and_outside_nodes():
    df = pd.DataFrame({
        "r": ["a", "a", "b", "b", "b", "x", "c"],
        "s": ["b", "a", "x", "c", "c", "a", "a"],
        "p": [5.0, 7.0, 3.0, 2.0, 1.0, 4.0, 0.5],
    })
    e = B.edges_from_frame(df, ["a", "b", "c"], expected_n=3)
    assert list(zip(e.i, e.j, e.w)) == [(0, 1, 5.0), (1, 2, 3.0), (2, 0, 0.5)]
    assert (e.i != e.j).all()


def test_edge_builder_asserts_n():
    df = pd.DataFrame({"r": ["a"], "s": ["b"], "p": [1.0]})
    with pytest.raises(AssertionError):
        B.edges_from_frame(df, ["a", "b", "c"], expected_n=4)


@pytest.mark.skipif(not C.ZONES_GPKG.exists(), reason="zones GeoPackage not available")
def test_load_nodes_has_expected_n():
    nd = B.load_nodes()
    assert len(nd) == C.EXPECTED_N_NODES
    assert not nd.prov.isin(C.EXCLUDED_PROVINCES).any()
    assert (nd.idx.values == np.arange(len(nd))).all()


# ---------------------------------------------------------------- backbones
def _edges(ws):
    return pd.DataFrame({"i": np.arange(len(ws)), "j": np.arange(len(ws)) + 1, "w": np.array(ws, float)})


def test_backbone_threshold():
    b = B.backbone_threshold(_edges([12, 10, 9.99, 3]), 10)
    assert sorted(b.w.tolist()) == [10.0, 12.0]


@pytest.mark.parametrize("E", [0, 1, 2, 3, 5, 10])
def test_backbone_equal_edges_exact_E_with_ties(E):
    e = _edges([5, 5, 5, 3, 1])
    b = B.backbone_equal_edges(e, E)
    assert len(b) == min(E, len(e))
    if E == 2:  # tie at w=5 broken by (i, j) ascending
        assert b.i.tolist() == [0, 1]


def test_to_igraph_undirected_merges_reciprocal_pairs():
    e = pd.DataFrame({"i": [0, 1, 1], "j": [1, 0, 2], "w": [1.0, 1.0, 1.0]})
    assert B.to_igraph(e, 3, directed=True).ecount() == 3
    Gu = B.to_igraph(e, 3, directed=False)
    assert Gu.ecount() == 2 and Gu.is_simple() and not Gu.is_directed()


# ---------------------------------------------------------------- metrics
def test_density():
    assert P.density_directed(4, 3) == pytest.approx(3 / 12)
    assert P.density_undirected(4, 3) == pytest.approx(6 / 12)


def test_path_stats_path_graph():
    L, d = P.path_stats(ig.Graph(n=4, edges=[(0, 1), (1, 2), (2, 3)]))
    assert L == pytest.approx(10 / 6) and d == 3


def test_path_stats_uses_giant_component():
    L, d = P.path_stats(ig.Graph(n=6, edges=[(0, 1), (1, 2), (4, 5)]))
    assert L == pytest.approx(4 / 3) and d == 2


def test_clustering_triangle_pendant():
    m = P.und_metrics(triangle_pendant())
    assert m["C_avg"] == pytest.approx(7 / 12)
    assert m["C_global"] == pytest.approx(3 / 5)


def test_r_undirected_star_is_minus_one():
    star = ig.Graph.Star(5)
    assert P.r_undirected(star) == pytest.approx(-1.0)
    assert P.r_undirected(star) == pytest.approx(star.assortativity_degree(directed=False))


def test_r_undirected_matches_igraph_on_random_graph():
    ig.set_random_number_generator(__import__("random").Random(1))
    G = ig.Graph.Barabasi(200, 3)
    assert P.r_undirected(G) == pytest.approx(G.assortativity_degree(directed=False))


def test_r_out_in_known_value():
    # edges 0->1, 0->2, 3->1: (k_out src, k_in tgt) = (2,2), (2,1), (1,2) -> r = -0.5
    G = ig.Graph(n=4, edges=[(0, 1), (0, 2), (3, 1)], directed=True)
    assert P.r_out_in(G) == pytest.approx(-0.5)


def test_reciprocity():
    G = ig.Graph(n=3, edges=[(0, 1), (1, 0), (1, 2)], directed=True)
    assert P.reciprocity(G) == pytest.approx(2 / 3)


def test_rich_club_phi():
    Nk, Ek, phi = P.rich_club_phi(k4_with_leaves(), [0, 1, 3, 4])
    assert Nk.tolist() == [8, 4, 4, 0]
    assert Ek.tolist() == [10, 6, 6, 0]
    assert phi[0] == pytest.approx(20 / 56)
    assert phi[1] == pytest.approx(1.0) and phi[2] == pytest.approx(1.0)
    assert math.isnan(phi[3])


def test_kcore():
    G = k4_with_leaves()
    core = G.coreness()
    assert max(core) == 3 and core[:4] == [3, 3, 3, 3] and core[4:] == [1, 1, 1, 1]
    assert P.und_metrics(G)["kcore_max"] == 3


def test_giant_components():
    Gd = ig.Graph(n=5, edges=[(0, 1), (1, 2), (2, 0), (2, 3)], directed=True)
    assert P.dir_metrics(Gd)["giant_strong"] == 3
    Gu = B.to_igraph(pd.DataFrame(Gd.get_edgelist(), columns=["i", "j"]), 5, directed=False)
    assert P.und_metrics(Gu)["giant_weak"] == 4


def test_degree_stats():
    m = P.und_metrics(ig.Graph.Star(5))
    assert (m["k_mean"], m["k_median"], m["k_max"]) == (pytest.approx(8 / 5), 1.0, 4)
    d = P.dir_metrics(ig.Graph(n=3, edges=[(0, 1), (0, 2)], directed=True))
    assert d["kout_max"] == 2 and d["kin_max"] == 1 and d["kin_mean"] == pytest.approx(2 / 3)


def test_ccdf():
    k, p = P.ccdf([0, 1, 1, 3])
    assert k.tolist() == [0, 1, 3]
    assert p.tolist() == pytest.approx([1.0, 0.75, 0.25])


def test_zscore():
    mu, sd, z = P.zscore(5.0, np.array([1.0, 2.0, 3.0]))
    assert (mu, sd, z) == (pytest.approx(2.0), pytest.approx(1.0), pytest.approx(3.0))
    assert math.isnan(P.zscore(1.0, np.array([2.0, 2.0]))[2])


# ---------------------------------------------------------------- null models
def _hetero_directed(seed=7, n=150):
    rng = np.random.default_rng(seed)
    w = rng.pareto(1.5, n) + 1
    p = np.outer(w, w)
    p = np.clip(0.02 * p, 0, 0.9)
    A = rng.random((n, n)) < p
    np.fill_diagonal(A, False)
    i, j = np.nonzero(A)
    return ig.Graph(n=n, edges=list(zip(i.tolist(), j.tolist())), directed=True)


def test_rewire_directed_preserves_in_out_degrees_exactly():
    G = _hetero_directed()
    assert G.ecount() > 300
    H = P.rewire(G, seed=C.SEED)
    assert H.indegree() == G.indegree()
    assert H.outdegree() == G.outdegree()
    assert H.is_simple() and not any(H.is_loop())
    assert H.ecount() == G.ecount()
    assert set(H.get_edgelist()) != set(G.get_edgelist())


def test_rewire_undirected_preserves_degrees_exactly():
    Gu = _hetero_directed().as_undirected(mode="collapse")
    H = P.rewire(Gu, seed=C.SEED)
    assert H.degree() == Gu.degree()
    assert H.is_simple() and not any(H.is_loop()) and not any(H.is_multiple())
    assert set(H.get_edgelist()) != set(Gu.get_edgelist())


def test_rewire_is_reproducible_with_seed():
    G = _hetero_directed()
    assert P.rewire(G, 123).get_edgelist() == P.rewire(G, 123).get_edgelist()
    assert P.rewire(G, 123).get_edgelist() != P.rewire(G, 124).get_edgelist()


def test_null_seeds_distinct():
    s = {P.null_seed(d, b, nl, r) for d in ("2024-02-14", "2024-08-15") for b in P.BACKBONES
         for nl in P.NULLS for r in range(5)}
    assert len(s) == 2 * 2 * 4 * 5


def test_erdos_renyi_same_n_e():
    for directed in (False, True):
        H = P.erdos_renyi(100, 400, directed, seed=1)
        assert H.vcount() == 100 and H.ecount() == 400 and H.is_simple()
        assert H.is_directed() == directed


# ---------------------------------------------------------------- representativeness logic
def test_unique_extreme():
    assert P.unique_extreme(5, [1, 2, 5])
    assert not P.unique_extreme(2, [1, 2, 5])
    assert not P.unique_extreme(5, [5, 5, 1])      # tied extreme is not "the" extreme
    assert not P.unique_extreme(3, [3, 3, 3])


def test_pairwise_sign_agreement():
    assert P.pairwise_sign_agreement([5, 6], [1, 2, 3], 1.0) == (6, 6)
    assert P.pairwise_sign_agreement([2.5, 6], [1, 2, 3], 1.0) == (5, 6)


# ---------------------------------------------------------------- audit C3: window means
def test_window_mean_divides_by_nights_in_window_not_nights_present():
    df = pd.DataFrame({
        "date": ["d1", "d2", "d3", "d1"],
        "r": ["a", "a", "a", "b"],
        "s": ["b", "b", "b", "c"],
        "p": [3.0, 3.0, 3.0, 12.0],   # (b, c) present on one night only
    })
    m = B.window_mean_from_frame(df, n_nights=3).set_index(["r", "s"]).p
    assert m[("a", "b")] == pytest.approx(3.0)
    assert m[("b", "c")] == pytest.approx(4.0)       # 12 / 3, not 12 / 1
    m4 = B.window_mean_from_frame(df, n_nights=4).set_index(["r", "s"]).p   # one night missing for everyone
    assert m4[("b", "c")] == pytest.approx(3.0)
    with pytest.raises(AssertionError):
        B.window_mean_from_frame(df, n_nights=2)


@pytest.mark.skipif(not C.ZONES_GPKG.exists(), reason="zones GeoPackage not available")
def test_node_set_drops_long_ids_and_carries_aggregation():
    nd = B.load_nodes()
    assert len(nd) == C.EXPECTED_N_NODES == 2435
    assert (nd.id.str.len() <= 8).all()
    assert (nd.n_municipalities >= 1).all() and (nd.area_km2 > 0).all()
    # Only non-'_AM' exception in the gpkg: 03015 Almoradi also lists 03903 (a municipality split across zones).
    multi = nd[~nd.is_am & (nd.n_municipalities != 1)]
    assert multi.id.tolist() == ["03015"]


# ---------------------------------------------------------------- rich club by rank, R_avg
def test_rich_club_top_rank_with_tie_break():
    En, phi = P.rich_club_top(k4_with_leaves(), [4, 5])
    assert En.tolist() == [6, 7]           # top-5 adds leaf 4 (lowest index among degree-1 ties)
    assert phi.tolist() == pytest.approx([1.0, 0.7])


def test_second_order_r_path():
    # Path 0-1-2-3: degrees [1,2,2,1], excess [0,1,1,0].
    # Excluding other endpoint: (0,1)->(0,1), (1,2)->(0,0), (2,3)->(1,0) -> R_avg = R_max = -0.5.
    G = ig.Graph(n=4, edges=[(0, 1), (1, 2), (2, 3)])
    assert P.second_order_r(G) == pytest.approx(-0.5)
    assert P.second_order_r_max(G) == pytest.approx(-0.5)
    assert P.second_order_r_incl(G) == pytest.approx(-0.5)


def test_second_order_r_avg_and_r_max_hand_computed():
    # 6-node graph: (0,1),(0,2),(0,3),(1,2),(2,4),(4,5).
    # Degrees: [3,2,3,1,2,1], excess degrees: [2,1,2,0,1,0].
    # Hand-computed (excluding other endpoint):
    #   R_avg pairs: (1,2),(0.5,1),(1.5,0),(2,1.5),(1.5,0),(2,0) -> R_avg = -37/83
    #   R_max pairs: (2,2),(1,1),(2,0),(2,2),(2,0),(2,0)         -> R_max = -5/13
    G = ig.Graph(n=6, edges=[(0, 1), (0, 2), (0, 3), (1, 2), (2, 4), (4, 5)])
    Ka, Kb, Ma, Mb = Z.second_order_sequences(G)
    assert P.second_order_r(G) == pytest.approx(-37.0 / 83.0)
    assert P.second_order_r_max(G) == pytest.approx(-5.0 / 13.0)
    assert Z.second_order_mixing(G) == (pytest.approx(-37.0 / 83.0), pytest.approx(-5.0 / 13.0))
    # Distinct from pre-change inclusive plain-degree R_avg_incl:
    assert abs(P.second_order_r(G) - P.second_order_r_incl(G)) > 0.05

    # Exact O(L) jackknife matches brute-force leave-one-link-out loop.
    for x, y in ((Ka, Kb), (Ma, Mb)):
        r_full, sig_fast = Z.jackknife_sigma_pairs(x, y)
        ri_brute = []
        for idx in range(len(x)):
            m = np.ones(len(x), dtype=bool)
            m[idx] = False
            ri_brute.append(Z._pearson_sym_pairs(x[m], y[m]))
        sig_brute = math.sqrt(sum((ri - r_full) ** 2 for ri in ri_brute))
        assert sig_fast == pytest.approx(sig_brute, rel=1e-12, abs=1e-12)

    # Sequence-permutation null has mean near 0 and 100 draws.
    sig = Z.mixing_significance(G, seed=C.SEED, n_perm=100)
    for met in ("r_und", "R_avg", "R_max"):
        assert sig[met]["perm"]["n_perm"] == 100
        assert abs(sig[met]["perm"]["mean"]) < 0.2
        assert sig[met]["jackknife_sigma"] > 0


# ---------------------------------------------------------------- Zhou-style rewirings and BOSAM
import b02_zhou as Z  # noqa: E402


def _hetero_undirected():
    return _hetero_directed(seed=3, n=120).as_undirected(mode="collapse")


def test_r_from_product_sum_matches_pearson():
    G = _hetero_undirected()
    d = G.degree()
    S = sum(d[a] * d[b] for a, b in G.get_edgelist())
    assert Z.r_from_product_sum(S, d, G.ecount()) == pytest.approx(P.r_undirected(G))


def test_extreme_rewiring_bounds():
    G = _hetero_undirected()
    Ha, ia = Z.extreme_rewire(G, "assortative", seed=1)
    Hd, idd = Z.extreme_rewire(G, "disassortative", seed=1)
    for H in (Ha, Hd):
        assert H.degree() == G.degree() and H.is_simple()
    r = P.r_undirected(G)
    assert P.r_undirected(Hd) < r < P.r_undirected(Ha)


def test_extreme_rewiring_two_stars():
    # two disjoint 3-stars + one edge between leaves; max-assortative must join the two hubs
    G = ig.Graph(n=8, edges=[(0, 1), (0, 2), (0, 3), (4, 5), (4, 6), (4, 7)])
    Ha, _ = Z.extreme_rewire(G, "assortative", seed=2)
    assert Ha.degree() == G.degree() and Ha.is_simple()
    assert P.r_undirected(Ha) >= P.r_undirected(G)


def test_keep_connected_rewiring_preserves_single_component():
    # Connected graph with two K4 cliques linked by a bridge plus two leaves (would disconnect
    # into a separate leaf-leaf edge under unconstrained max-assortative rewiring).
    edges = (
        [(a, b) for a in range(4) for b in range(a + 1, 4)]
        + [(a, b) for a in range(4, 8) for b in range(a + 1, 8)]
        + [(3, 4), (0, 8), (5, 9)]
    )
    G = ig.Graph(n=10, edges=edges, directed=False)
    assert Z._n_nontrivial_components(G) == 1
    for mode in ("assortative", "disassortative"):
        H, info = Z.extreme_rewire(G, mode, seed=11, keep_connected=True)
        assert H.degree() == G.degree() and H.is_simple()
        assert Z._n_nontrivial_components(H) == 1
        assert info["keep_connected"] is True
    Hr, rinfo = Z.random_rewire(G, seed=11, swaps_per_edge=20, keep_connected=True)
    assert Hr.degree() == G.degree() and Hr.is_simple()
    assert Z._n_nontrivial_components(Hr) == 1
    Hrc, _ = Z.rich_club_preserving_rewire(G, seed=11, swaps_per_edge=20, keep_connected=True)
    assert Hrc.degree() == G.degree() and Hrc.is_simple()
    assert Z._n_nontrivial_components(Hrc) == 1


def test_rich_club_preserving_rewiring_keeps_degrees_and_phi():
    G = _hetero_undirected()
    H, acc = Z.rich_club_preserving_rewire(G, seed=5)
    ks = np.arange(0, max(G.degree()) + 1)
    assert H.degree() == G.degree() and H.is_simple()
    assert np.array_equal(P.rich_club_phi(H, ks)[1], P.rich_club_phi(G, ks)[1])
    assert acc > 0 and set(H.get_edgelist()) != set(G.get_edgelist())


def test_bosam_order():
    G = ig.Graph(n=5, edges=[(0, 1), (0, 2), (0, 3), (3, 4)])
    key = Z.bosam_order(G)
    assert key[:, 0].tolist() == [4, 1, 2, 3, 0]
    assert key[0].tolist() == [4, 1, 2, 3]


def test_bosam_order_iterative_third_key_hand_worked():
    # 10-node graph where old max_nb_label and Zhou's iterative max_nb_index strictly disagree:
    # node 5 (k=5): nbs {0,4,6,7,8}; node 4 (k=3): nbs {1,5,9};
    # node 0 (k=2, omega=5): nbs {2,5}; node 1 (k=2, omega=3): nbs {3,4};
    # leaves 2 (nb 0) and 3 (nb 1) tie on k=1 and omega=2.
    # Because pos[1]=6 < pos[0]=7 in BOSAM order, e_3=6 < e_2=7 -> leaf 3 MUST precede leaf 2
    # (even though original label 0 < 1).
    edges = [(5, 0), (5, 4), (5, 6), (5, 7), (5, 8), (4, 1), (4, 9), (0, 2), (1, 3)]
    G = ig.Graph(n=10, edges=edges, directed=False)
    key, info = Z.bosam_order(G, return_info=True)
    assert info["converged"] is True
    assert info["iterations"] >= 2  # required iterative pass to swap leaves 2 and 3
    order = key[:, 0].tolist()
    assert order == [3, 2, 9, 6, 7, 8, 1, 0, 4, 5]
    # Verify exact BOSAM invariants: for all i < j, (k_i, omega_i, e_i) <= (k_j, omega_j, e_j).
    tuples = [tuple(row[1:].tolist()) for row in key]
    assert tuples == sorted(tuples)
    # Check specific third-key values e_3 = pos[1] = 6 and e_2 = pos[0] = 7.
    assert key[0].tolist() == [3, 1, 2, 6]
    assert key[1].tolist() == [2, 1, 2, 7]


def test_rich_club_zhou_clique_and_link_distribution():
    # K4 on 0..3 (deg 4) with leaves 4..7 (deg 1) -> clique size = 4.
    G = k4_with_leaves()
    assert Z.rich_club_clique_size(G) == 4
    # Remove edge (0,1): degs [3,3,4,4,1,1,1,1]; top-3 = {2,3,0} is a clique, top-4 misses (0,1) -> 3.
    G3 = G.copy()
    G3.delete_edges([(0, 1)])
    assert Z.rich_club_clique_size(G3) == 3
    # Top-2 nodes (0 and 1, deg 3 each) NOT linked -> clique size = 0 (as Zhou 2007 Table 1).
    G0 = ig.Graph(n=8, edges=[(0, 2), (0, 3), (0, 4), (1, 5), (1, 6), (1, 7), (2, 5)])
    assert Z.rich_club_clique_size(G0) == 0

    # Top-25% of K4+leaves = top 2 nodes {0, 1}: 1 link between them, 7 links attached to {0, 1}.
    nld = Z.node_node_link_distribution(G, frac=0.25)
    assert nld["top_n"] == 2
    assert nld["L"] == 10
    assert nld["l_top_both"] == 1
    assert nld["l_top_any"] == 7


def test_zhou_table_row_toy():
    G = k4_with_leaves()
    row = Z.zhou_table_row(G, "2024-02-14", "ordinary", "w10", seed=C.SEED)
    assert row["N"] == 8 and row["L"] == 10
    assert row["k_mean"] == pytest.approx(2.5)
    assert row["k_max"] == 4
    assert row["n_clique"] == 4
    assert row["l_star"] == pytest.approx(math.log(8) / math.log(2.5))
    assert "log-normal" in row["gamma_note"].lower()


def test_powerlaw_fit_recovers_exponent():
    x = np.random.default_rng(C.SEED).zipf(2.5, 5000)
    f = Z.powerlaw_vs_lognormal(x)
    assert abs(f["alpha"] - 2.5) < 0.15
    assert f["xmin"] >= 1
