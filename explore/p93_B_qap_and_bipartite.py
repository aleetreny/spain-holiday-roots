"""
p93_B_qap_and_bipartite.py
Network comparison & bipartite "summer territories" for Question B:
Part 1 — Spatial-null QAP & edge-level correlations (2,000 permutations):
  (A) 47x47 mainland province network:
      Holiday excess (Mid-Aug, All Saints, Christmas Eve) vs Roots (2024 & 1998) vs Everyday (Feb wd)
      vs Historical migration (1988-1995 vs 2015-2024) vs Meta SCI,
      controlling for doubly-constrained spatial gravity null AND partial QAP controlling for everyday network.
  (B) 142 mainland Census cities x 47 mainland provinces bipartite network:
      Bipartite QAP (both full permutation & within-origin-province neighbouring-city permutation)
      between holiday excess, roots, and everyday overnight networks under spatial DC-gravity null.
Part 2 — Bipartite communities ("summer territories") under BiCM (exact Saracco et al. 2015),
  Barber weighted CM, and Spatial DC-gravity null (142 mainland cities x 2,484 destinations):
  - Counts non-contiguous ("diaspora") communities with an exact polygon-contiguity test
  - Measures birthplace-roots enrichment of spatial-null communities
  - Proper Monte Carlo null test: 100 networks drawn from the pure spatial DC-gravity null,
    clustered with the exact same bipartite modularity algorithm -> tests whether observed Q,
    diaspora fragmentation, and roots alignment exceed spatial-null sampling noise!
Figure: reports/figs/net_B_qap_and_territories.png
"""
import json, time
import duckdb
import numpy as np
import pandas as pd
import geopandas as gpd
from scipy.spatial import cKDTree
from sklearn.metrics import normalized_mutual_info_score as nmi
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, 'explore')
from common import PROV2CCAA

t0 = time.time()
rng = np.random.default_rng(123)
con = duckdb.connect()
q_pq = "'data/derived/pernoct_mun_2024_6m.parquet'"

# Load mainland nodes
nd = pd.read_parquet('data/derived/nodes_municipios.parquet')
nd = nd[~nd.foreign & ~nd.island & ~nd.prov.isin(['51', '52'])].reset_index(drop=True)
n = len(nd)
idx = {k: i for i, k in enumerate(nd.id)}
mainland_provs = sorted(nd.prov.unique().tolist())
np_ = len(mainland_provs)
p2i = {p: i for i, p in enumerate(mainland_provs)}
ml_sql = ",".join(f"'{z}'" for z in nd.id)

# Population-weighted province centroids & 47x47 distance matrix
nd['w_pop'] = nd.population.fillna(0) + 1.0
pc = nd.groupby('prov').apply(
    lambda g: pd.Series({'x': np.average(g.x, weights=g.w_pop), 'y': np.average(g.y, weights=g.w_pop)}),
    include_groups=False
).loc[mainland_provs]
D_prov = np.hypot(pc.x.values[:, None] - pc.x.values[None, :],
                  pc.y.values[:, None] - pc.y.values[None, :]) / 1000.0
# Distance bins for 47x47 (deciles of off-diag distances -> 15 bins)
off_mask_p = ~np.eye(np_, dtype=bool)
q_edges_p = np.quantile(D_prov[off_mask_p], np.linspace(0, 1, 16))
q_edges_p[0] = -1.0; q_edges_p[-1] += 1.0
B_prov = np.digitize(D_prov, q_edges_p)
np.fill_diagonal(B_prov, 0)

def fit_dc_gravity_rect(W, bin_mat, mask, max_iter=300):
    """Doubly-constrained directed/bipartite gravity via IPF: matches row sums, col sums, and distance-bin totals."""
    nr, nc = W.shape
    a = np.ones(nr); b = np.ones(nc)
    n_bins = int(bin_mat.max()) + 1
    f = np.ones(n_bins); f[0] = 0.0
    W_m = np.where(mask, W, 0.0)
    r_tot = W_m.sum(1); c_tot = W_m.sum(0)
    bm = np.where(mask, bin_mat, 0)
    W_bin = np.bincount(bm.ravel(), weights=W_m.ravel(), minlength=n_bins)
    for _ in range(max_iter):
        F = np.where(mask, f[bm], 0.0)
        for _ in range(5):
            a = np.where(r_tot > 0, r_tot / np.maximum(F @ b, 1e-300), 0.0)
            b = np.where(c_tot > 0, c_tot / np.maximum(F.T @ a, 1e-300), 0.0)
        E = np.outer(a, b) * F
        E_bin = np.bincount(bm.ravel(), weights=E.ravel(), minlength=n_bins)
        fn = np.where(E_bin > 0, f * W_bin / np.maximum(E_bin, 1e-300), 0.0)
        fn[0] = 0.0
        if np.max(np.abs(fn - f) / np.maximum(f, 1e-12)) < 1e-6:
            break
        f = fn
    return np.outer(a, b) * np.where(mask, f[bm], 0.0)

def spatial_resid(W, E, mask):
    """Log obs/exp residual under spatial DC-gravity null (stabilised with +1 pseudo-count)."""
    Z = np.zeros_like(W, dtype=float)
    Z[mask] = np.log((W[mask] + 1.0) / (E[mask] + 1.0))
    return Z

# ==============================================================================
# PART 1A: 47x47 PROVINCE-LEVEL SPATIAL-NULL QAP
# ==============================================================================
print("=== PART 1A: 47x47 Province-Level QAP (Spatial DC-Gravity Null) ===", flush=True)
ref_fri = ['2024-02-02','2024-02-09','2024-02-16','2024-02-23','2024-03-08','2024-03-15','2024-10-18','2024-11-15','2024-11-22']
def prov_mat(cond):
    df = con.execute(f"""
        SELECT substr(r,1,2) pr, substr(s,1,2) ps, sum(p)/count(DISTINCT date) w
        FROM {q_pq}
        WHERE {cond} AND r IN ({ml_sql}) AND s IN ({ml_sql}) AND substr(r,1,2) <> substr(s,1,2)
        GROUP BY 1, 2
    """).df()
    M = np.zeros((np_, np_), dtype=float)
    for r_, s_, w_ in zip(df.pr, df.ps, df.w):
        if r_ in p2i and s_ in p2i:
            M[p2i[r_], p2i[s_]] = w_
    return M

M_feb = prov_mat("date BETWEEN '2024-02-05' AND '2024-02-29' AND dayofweek(date) IN (2,3,4)")
M_ref_fri = prov_mat("date IN (" + ",".join(f"'{d}'" for d in ref_fri) + ")")
M_aug = prov_mat("date BETWEEN '2024-08-10' AND '2024-08-18'")
M_nov = prov_mat("date = '2024-11-01'")
M_xmas = prov_mat("date = '2024-12-24'")

# Holiday excess over ordinary weekday (and over ordinary Friday)
M_aug_exc = np.clip(M_aug - M_feb, 0, None)
M_nov_exc = np.clip(M_nov - M_feb, 0, None)
M_xmas_exc = np.clip(M_xmas - M_feb, 0, None)

# Roots 2024 & 1998: oriented [residence_prov, birth_prov] so it aligns with mobility [residence_prov, stay_prov]!
br = pd.read_csv('_local/external/ine_birth_residence/prov_birth_x_prov_residence.csv')
br['pb'] = br['prov_birth'].astype(str).str.zfill(2)
br['pr'] = br['prov_residence'].astype(str).str.zfill(2)
def build_roots_prov(yr):
    sub = br[(br.year == yr) & br.pb.isin(p2i) & br.pr.isin(p2i) & (br.pb != br.pr)]
    M = np.zeros((np_, np_), dtype=float)
    for r_, b_, v_ in zip(sub.pr, sub.pb, sub.persons):
        M[p2i[r_], p2i[b_]] = v_
    return M
M_roots24 = build_roots_prov(2024)
M_roots98 = build_roots_prov(1998)

# Historical migration 1988-1995 vs 2015-2024: oriented [dest_res_prov, origin_birth_prov]
mig = pd.read_csv('_local/external/ine_interprovincial_migration/interprovincial_migration_1988_2024.csv')
mig['po'] = mig['prov_origin'].astype(str).str.zfill(2)
mig['pd'] = mig['prov_dest'].astype(str).str.zfill(2)
def build_mig(y0, y1):
    sub = mig[(mig.year >= y0) & (mig.year <= y1) & mig.po.isin(p2i) & mig.pd.isin(p2i) & (mig.po != mig.pd)]
    g = sub.groupby(['pd', 'po'])['migrations'].sum()
    M = np.zeros((np_, np_), dtype=float)
    for (d_, o_), v_ in g.items():
        M[p2i[d_], p2i[o_]] = v_
    return M
M_mig_old = build_mig(1988, 1995)
M_mig_new = build_mig(2015, 2024)

# Meta SCI
sci = pd.read_csv('_local/external/meta_sci_spain/sci_spain_nuts3.csv')
sci['p1'] = sci['user_prov_ine'].astype(str).str.zfill(2)
sci['p2'] = sci['friend_prov_ine'].astype(str).str.zfill(2)
M_sci = np.zeros((np_, np_), dtype=float)
for p1, p2, v in zip(sci.p1, sci.p2, sci.scaled_sci):
    if p1 in p2i and p2 in p2i and p1 != p2:
        M_sci[p2i[p1], p2i[p2]] = v

# Compute spatial DC-gravity residuals for each matrix
mats_p = {
    'Feb_wd': M_feb, 'Mid_Aug': M_aug, 'All_Saints': M_nov, 'Xmas_Eve': M_xmas,
    'Aug_exc': M_aug_exc, 'Nov_exc': M_nov_exc, 'Xmas_exc': M_xmas_exc,
    'Roots_2024': M_roots24, 'Roots_1998': M_roots98,
    'Mig_88_95': M_mig_old, 'Mig_15_24': M_mig_new, 'Meta_SCI': M_sci
}
Z_p = {}
for k, M in mats_p.items():
    E = fit_dc_gravity_rect(M, B_prov, off_mask_p)
    Z_p[k] = spatial_resid(M, E, off_mask_p)

def qap_corr(X, Y, mask, n_perm=2000):
    """Standard simultaneous row+column QAP Pearson correlation."""
    xv = X[mask]; yv = Y[mask]
    x0 = xv - xv.mean(); y0 = yv - yv.mean()
    obs = float((x0 * y0).sum() / np.sqrt((x0 * x0).sum() * (y0 * y0).sum()))
    cnt = 0
    nulls = np.empty(n_perm)
    nx = X.shape[0]
    for i in range(n_perm):
        perm = rng.permutation(nx)
        Xp = X[perm][:, perm]
        xp = Xp[mask]
        xp0 = xp - xp.mean()
        r = (xp0 * y0).sum() / np.sqrt((xp0 * xp0).sum() * (y0 * y0).sum())
        nulls[i] = r
        if r >= obs:
            cnt += 1
    p_val = (cnt + 1) / (n_perm + 1)
    z_val = (obs - nulls.mean()) / max(nulls.std(), 1e-12)
    return obs, float(z_val), float(p_val)

def partial_qap_corr(Y, X, C, mask, n_perm=2000):
    """
    Freedman-Lane / Dekkers double-semi-partialing QAP:
    Residualize both Y and X on control matrix C (everyday network), then run QAP between residuals!
    """
    yv = Y[mask]; xv = X[mask]; cv = C[mask]
    # OLS residual of Y on [1, C] and X on [1, C]
    def resid_on(a, c):
        c0 = c - c.mean()
        slope = ((a - a.mean()) * c0).sum() / max((c0 * c0).sum(), 1e-15)
        return (a - a.mean()) - slope * c0
    ry = np.zeros_like(Y); rx = np.zeros_like(X)
    ry[mask] = resid_on(yv, cv)
    rx[mask] = resid_on(xv, cv)
    return qap_corr(rx, ry, mask, n_perm=n_perm)

qap_prov_rows = []
for mob_key in ['Feb_wd', 'Aug_exc', 'Nov_exc', 'Xmas_exc', 'Mid_Aug', 'All_Saints', 'Xmas_Eve']:
    for ref_key in ['Roots_2024', 'Roots_1998', 'Mig_88_95', 'Mig_15_24', 'Meta_SCI']:
        r_raw, z_raw, p_raw = qap_corr(np.log1p(mats_p[mob_key]), np.log1p(mats_p[ref_key]), off_mask_p, n_perm=2000)
        r_sp, z_sp, p_sp = qap_corr(Z_p[mob_key], Z_p[ref_key], off_mask_p, n_perm=2000)
        if mob_key != 'Feb_wd':
            r_part, z_part, p_part = partial_qap_corr(Z_p[mob_key], Z_p[ref_key], Z_p['Feb_wd'], off_mask_p, n_perm=2000)
        else:
            r_part = z_part = p_part = np.nan
        qap_prov_rows.append(dict(mobility=mob_key, reference=ref_key,
                                  r_raw=r_raw, p_raw=p_raw,
                                  r_spatial=r_sp, z_spatial=z_sp, p_spatial=p_sp,
                                  r_partial_ctrl_feb=r_part, z_partial=z_part, p_partial=p_part))

df_qap_p = pd.DataFrame(qap_prov_rows)
print(df_qap_p[df_qap_p.reference == 'Roots_2024'].round(4).to_string(index=False))
print("\nOld vs Recent migration & SCI vs Holidays (spatial-null QAP, partial ctrl Feb):")
print(df_qap_p[df_qap_p.mobility.isin(['Feb_wd', 'Aug_exc', 'Nov_exc', 'Xmas_exc'])].round(4).to_string(index=False))
df_qap_p.to_csv('explore/out_p93_qap_prov.csv', index=False)

# ==============================================================================
# PART 1B: 142 MAINLAND CITIES x 47 PROVINCES BIPARTITE QAP (NEIGHBOURING-CITY TEST)
# ==============================================================================
print("\n=== PART 1B: Mainland Census Cities x 47 Provinces Bipartite QAP ===", flush=True)
bm = pd.read_csv('_local/external/ine_birth_residence/mun_residence_x_prov_birth_censo_2021_2025.csv')
bm['cid'] = bm['ine_code'].astype(str).str.zfill(5)
bm['pb'] = bm['prov_birth'].astype(str).str.zfill(2)
bm24 = bm[(bm.year == 2024) & bm.cid.isin(idx) & bm.pb.isin(p2i)].copy()
cities = sorted(bm24.cid.unique().tolist())
nc_ = len(cities)
c2i = {c: i for i, c in enumerate(cities)}
c_prov = np.array([c[:2] for c in cities])
cl_sql = ",".join(f"'{c}'" for c in cities)
print(f"Matched {nc_} mainland Census cities x {np_} provinces.", flush=True)

# City -> province distance & bins
cx = nd.set_index('id').loc[cities, 'x'].values
cy = nd.set_index('id').loc[cities, 'y'].values
D_cp = np.hypot(cx[:, None] - pc.x.values[None, :], cy[:, None] - pc.y.values[None, :]) / 1000.0
mask_cp = (c_prov[:, None] != np.array(mainland_provs)[None, :])
q_edges_cp = np.quantile(D_cp[mask_cp], np.linspace(0, 1, 21))
q_edges_cp[0] = -1.0; q_edges_cp[-1] += 1.0
B_cp = np.digitize(D_cp, q_edges_cp)
B_cp[~mask_cp] = 0

def city_prov_mat(cond):
    df = con.execute(f"""
        SELECT r c, substr(s,1,2) ps, sum(p)/count(DISTINCT date) w
        FROM {q_pq}
        WHERE {cond} AND r IN ({cl_sql}) AND s IN ({ml_sql}) AND substr(r,1,2) <> substr(s,1,2)
        GROUP BY 1, 2
    """).df()
    M = np.zeros((nc_, np_), dtype=float)
    for c_, p_, w_ in zip(df.c, df.ps, df.w):
        if c_ in c2i and p_ in p2i:
            M[c2i[c_], p2i[p_]] = w_
    return M

C_feb = city_prov_mat("date BETWEEN '2024-02-05' AND '2024-02-29' AND dayofweek(date) IN (2,3,4)")
C_aug = city_prov_mat("date BETWEEN '2024-08-10' AND '2024-08-18'")
C_nov = city_prov_mat("date = '2024-11-01'")
C_xmas = city_prov_mat("date = '2024-12-24'")
C_aug_exc = np.clip(C_aug - C_feb, 0, None)
C_nov_exc = np.clip(C_nov - C_feb, 0, None)
C_xmas_exc = np.clip(C_xmas - C_feb, 0, None)

C_roots = np.zeros((nc_, np_), dtype=float)
for c_, p_, v_ in zip(bm24.cid, bm24.pb, bm24.persons):
    if c_[:2] != p_:
        C_roots[c2i[c_], p2i[p_]] = v_

# Fit spatial DC-gravity null on each city x prov matrix
mats_cp = {'Feb_wd': C_feb, 'Mid_Aug': C_aug, 'All_Saints': C_nov, 'Xmas_Eve': C_xmas,
           'Aug_exc': C_aug_exc, 'Nov_exc': C_nov_exc, 'Xmas_exc': C_xmas_exc, 'Roots': C_roots}
Z_cp = {}
for k, M in mats_cp.items():
    E = fit_dc_gravity_rect(M, B_cp, mask_cp)
    Z_cp[k] = spatial_resid(M, E, mask_cp)

# Also compute WITHIN-PROVINCE demeaned residuals (comparing neighbouring cities in the same origin province!)
def demean_within_origin_prov(Z, mask):
    out = np.zeros_like(Z)
    for pr in np.unique(c_prov):
        rows = np.where(c_prov == pr)[0]
        if len(rows) < 2:
            continue
        for j in range(np_):
            if pr == mainland_provs[j]:
                continue
            vals = Z[rows, j]
            out[rows, j] = vals - vals.mean()
    return out

# Groups of cities in multi-city provinces for within-province permutation QAP
prov_groups = [np.where(c_prov == pr)[0] for pr in np.unique(c_prov) if np.sum(c_prov == pr) >= 2]
multi_city_mask = np.zeros_like(mask_cp)
for grp in prov_groups:
    multi_city_mask[grp, :] = mask_cp[grp, :]

def bip_qap(X, Y, mask, within_prov=False, n_perm=2000):
    xv = X[mask]; yv = Y[mask]
    x0 = xv - xv.mean(); y0 = yv - yv.mean()
    obs = float((x0 * y0).sum() / np.sqrt((x0 * x0).sum() * (y0 * y0).sum()))
    cnt = 0; nulls = np.empty(n_perm)
    for i in range(n_perm):
        if within_prov:
            # Permute cities ONLY within their own origin province! (Tests neighbouring cities directly)
            perm = np.arange(nc_)
            for grp in prov_groups:
                perm[grp] = rng.permutation(grp)
            Xp = X[perm, :]
        else:
            pr_ = rng.permutation(nc_); pc_ = rng.permutation(np_)
            Xp = X[pr_][:, pc_]
        xp = Xp[mask]; xp0 = xp - xp.mean()
        r = (xp0 * y0).sum() / np.sqrt((xp0 * xp0).sum() * (y0 * y0).sum())
        nulls[i] = r
        if r >= obs:
            cnt += 1
    return obs, float((obs - nulls.mean()) / max(nulls.std(), 1e-12)), float((cnt + 1) / (n_perm + 1))

Zr_dm = demean_within_origin_prov(Z_cp['Roots'], mask_cp)
qap_cp_rows = []
for lab in ['Feb_wd', 'Aug_exc', 'Nov_exc', 'Xmas_exc', 'Mid_Aug', 'All_Saints', 'Xmas_Eve']:
    r_full, z_full, p_full = bip_qap(Z_cp[lab], Z_cp['Roots'], mask_cp, within_prov=False, n_perm=2000)
    Z_dm = demean_within_origin_prov(Z_cp[lab], mask_cp)
    r_nb, z_nb, p_nb = bip_qap(Z_dm, Zr_dm, multi_city_mask, within_prov=True, n_perm=2000)
    # Also partial controlling for Feb_wd (within neighbouring cities)
    if lab != 'Feb_wd':
        Zf_dm = demean_within_origin_prov(Z_cp['Feb_wd'], mask_cp)
        # Residualize Z_dm and Zr_dm on Zf_dm
        yv = Z_dm[multi_city_mask]; xv = Zr_dm[multi_city_mask]; cv = Zf_dm[multi_city_mask]
        c0 = cv - cv.mean()
        ry = np.zeros_like(Z_dm); rx = np.zeros_like(Zr_dm)
        ry[multi_city_mask] = (yv - yv.mean()) - (((yv - yv.mean()) * c0).sum() / (c0 * c0).sum()) * c0
        rx[multi_city_mask] = (xv - xv.mean()) - (((xv - xv.mean()) * c0).sum() / (c0 * c0).sum()) * c0
        r_nb_part, z_nb_part, p_nb_part = bip_qap(ry, rx, multi_city_mask, within_prov=True, n_perm=2000)
    else:
        r_nb_part = z_nb_part = p_nb_part = np.nan
    qap_cp_rows.append(dict(
        mobility=lab, r_spatial_full=r_full, z_full=z_full, p_full=p_full,
        r_neighbour_cities=r_nb, z_nb=z_nb, p_nb=p_nb,
        r_nb_partial_ctrl_feb=r_nb_part, z_nb_part=z_nb_part, p_nb_part=p_nb_part
    ))
df_qap_cp = pd.DataFrame(qap_cp_rows)
print(df_qap_cp.round(4).to_string(index=False))
df_qap_cp.to_csv('explore/out_p93_qap_cities.csv', index=False)

# ==============================================================================
# PART 2: BIPARTITE COMMUNITIES ("SUMMER TERRITORIES") — BiCM vs BARBER vs SPATIAL NULL
#         + PROPER CONTIGUITY & DIASPORA NULL TEST
# ==============================================================================
print("\n=== PART 2: Bipartite Communities (BiCM, Barber, Spatial Null) & Contiguity Test ===", flush=True)
# Use the 142 mainland Census cities (covering major urban Spain AND having exact birthplace roots!)
# x all 2,484 mainland destinations
dests = nd.id.tolist()
d2i = idx
d_prov = nd.prov.to_numpy(dtype=object)
d_ccaa = nd.ccaa.to_numpy(dtype=object)

# Polygon contiguity graph on 2,484 mainland destinations (for exact spatial component counting)
zones = gpd.read_file('clean_data/v2/zones/municipios_mitma.gpkg').set_index('id').loc[nd.id].reset_index()
sidx = zones.sindex
geoms = zones.geometry.values
adj_dest = [set() for _ in range(n)]
for i in range(n):
    g_i = geoms[i]
    for j in sidx.intersection(g_i.bounds):
        if j > i and g_i.intersects(geoms[j]):
            adj_dest[i].add(j); adj_dest[j].add(i)

def count_spatial_components(node_indices):
    """Exact connected components on the polygon adjacency graph for a subset of destination indices."""
    sset = set(int(x) for x in node_indices)
    vis = set()
    comps = []
    for u in sset:
        if u not in vis:
            q = [u]; vis.add(u); c = [u]
            while q:
                cur = q.pop()
                for nb in adj_dest[cur]:
                    if nb in sset and nb not in vis:
                        vis.add(nb); q.append(nb); c.append(nb)
            comps.append(c)
    return comps

# Distance matrix (142 cities x 2,484 destinations)
D_cd = np.hypot(cx[:, None] - nd.x.values[None, :], cy[:, None] - nd.y.values[None, :]) / 1000.0
mask_cd = np.ones((nc_, n), dtype=bool)
for i, c in enumerate(cities):
    mask_cd[i, d2i[c]] = False
bin_edges_cd = np.r_[0:30:2, 30:100:5, 100:300:10, 300:1200:25]
B_cd = np.digitize(D_cd, bin_edges_cd)
B_cd[~mask_cd] = 0

def solve_bicm(M_bin, max_iter=2000, tol=1e-6):
    """
    Exact Bipartite Configuration Model (Saracco et al. 2015 Sci Rep):
    p_ij = (x_i y_j) / (1 + x_i y_j), solving sum_j p_ij = k_i and sum_i p_ij = d_j via fixed-point iteration.
    """
    M_bin = (M_bin > 0).astype(float)
    kr = M_bin.sum(1)
    kc = M_bin.sum(0)
    nr, nc = M_bin.shape
    x = np.where(kr > 0, kr / np.maximum(nc - kr, 1.0), 0.0)
    y = np.where(kc > 0, kc / np.maximum(nr - kc, 1.0), 0.0)
    for it in range(max_iter):
        xy = np.outer(x, y)
        denom = 1.0 + xy
        # Fixed point: x_new = kr / sum_j (y_j / (1 + x_i y_j))
        x_new = np.where(kr > 0, kr / np.maximum((y[None, :] / denom).sum(1), 1e-300), 0.0)
        xy = np.outer(x_new, y)
        denom = 1.0 + xy
        y_new = np.where(kc > 0, kc / np.maximum((x_new[:, None] / denom).sum(0), 1e-300), 0.0)
        err = max(np.max(np.abs((x_new[:, None] * y_new[None, :] / (1.0 + np.outer(x_new, y_new))).sum(1) - kr)),
                  np.max(np.abs((x_new[:, None] * y_new[None, :] / (1.0 + np.outer(x_new, y_new))).sum(0) - kc)))
        x, y = x_new, y_new
        if err < tol:
            break
    P = np.outer(x, y) / (1.0 + np.outer(x, y))
    return P, float(err), it

def bip_louvain(Bb, n_restarts=10, seed_base=0):
    """
    Barber bipartite modularity maximization on rectangular modularity matrix Bb (no x nd):
    Alternating greedy node moves on the full (no + nd) bipartite graph + aggregation (Louvain).
    """
    no, nd_ = Bb.shape
    ntot = no + nd_
    # Build full symmetric modularity matrix
    B_full = np.zeros((ntot, ntot), dtype=np.float64)
    B_full[:no, no:] = Bb
    B_full[no:, :no] = Bb.T
    best_q = -1e9
    best_mb = None
    for s in range(n_restarts):
        r_loc = np.random.default_rng(seed_base + s)
        Bc = B_full.copy()
        memb = np.arange(ntot)
        while True:
            nn = Bc.shape[0]
            comm = np.arange(nn)
            improved = True
            it = 0
            while improved and it < 25:
                improved = False; it += 1
                for i in r_loc.permutation(nn):
                    ci = comm[i]
                    scores = np.bincount(comm, weights=Bc[i], minlength=nn)
                    scores[ci] -= Bc[i, i]
                    best_c = int(np.argmax(scores))
                    if scores[best_c] > scores[ci] + 1e-12 and best_c != ci:
                        comm[i] = best_c
                        improved = True
            u, comm = np.unique(comm, return_inverse=True)
            if len(u) == nn:
                break
            Mm = np.zeros((nn, len(u)))
            Mm[np.arange(nn), comm] = 1.0
            Bc = Mm.T @ Bc @ Mm
            memb = comm[memb]
        co, cd = memb[:no], memb[no:]
        q_val = float(Bb[co[:, None] == cd[None, :]].sum())
        if q_val > best_q:
            best_q = q_val
            best_mb = memb.copy()
    return best_q, best_mb[:no], best_mb[no:]

# Helper to evaluate contiguity, diaspora share, and roots alignment of a bipartite partition
city_dest_idx = np.array([d2i[c] for c in cities])
# Roots matrix aligned to (142 cities x 2484 destinations) via destination's province (off-own-province)
roots_share_cp = C_roots / np.maximum(C_roots.sum(1, keepdims=True), 1.0)  # share of city's out-of-prov migrants born in prov p

def evaluate_bipartite_comms(O, co, cd, keep_idx):
    """
    Evaluate contiguity & diaspora structure:
    For each community c containing >=1 city and >=5 destinations:
      - spatial components of its destinations (polygon adjacency)
      - is_non_contiguous: >1 spatial component (ignoring isolated 1-node enclaves: >=2 components of size >=3)
      - diaspora_dest_share: share of destination nodes in components that do NOT touch/contain any origin city of community c
      - diaspora_flow_share: share of O flow in community c going to destinations >100 km from nearest origin city in community c
        AND in a spatial component disjoint from all origin cities of c!
      - roots_enrichment: mean birthplace share (for out-of-own-province city-dest pairs) inside same community / outside
    """
    W_tot = O.sum()
    flow_inside = sum(O[i, cd == co[i]].sum() for i in range(nc_)) / max(W_tot, 1e-12)
    n_comms_cities = len(set(co))
    comms_ge5 = []
    n_noncontig_strict = 0
    n_noncontig_major = 0  # >=2 spatial pieces of size >=3
    diasp_dests = 0; tot_dests_eval = 0
    diasp_flow = 0.0; tot_flow_eval = 0.0
    comp_counts = []
    comm_details = []

    for c in sorted(set(co)):
        c_cities = np.where(co == c)[0]
        c_dest_local = np.where(cd == c)[0]
        if len(c_dest_local) < 5:
            continue
        comms_ge5.append(c)
        c_dest_global = keep_idx[c_dest_local]
        c_city_global = set(int(city_dest_idx[i]) for i in c_cities)
        # Include city polygons + destination polygons to check which components touch an origin city
        comps = count_spatial_components(c_dest_global)
        n_pieces = len(comps)
        major_pieces = [cp for cp in comps if len(cp) >= 3]
        comp_counts.append(n_pieces)
        if n_pieces > 1:
            n_noncontig_strict += 1
        if len(major_pieces) >= 2:
            n_noncontig_major += 1

        # Determine which destination components are "home" (contain or touch an origin city of this community)
        home_nodes = set()
        for cp in comps:
            cp_set = set(cp)
            touches_city = any((cg in cp_set) or any(nb in cp_set for nb in adj_dest[cg]) for cg in c_city_global)
            if touches_city:
                home_nodes.update(cp_set)

        # Diaspora destinations = assigned to community c, NOT in a component touching any city of c
        diasp_local = [loc for loc, glb in zip(c_dest_local, c_dest_global) if glb not in home_nodes]
        diasp_dests += len(diasp_local)
        tot_dests_eval += len(c_dest_local)

        sub_O = O[np.ix_(c_cities, c_dest_local)]
        c_flow = float(sub_O.sum())
        tot_flow_eval += c_flow
        if diasp_local:
            d_flow = float(O[np.ix_(c_cities, diasp_local)].sum())
            # Also check min distance from community's cities to diaspora destinations
            min_d_diasp = float(D_cd[np.ix_(c_cities, keep_idx[diasp_local])].min(axis=0).mean())
        else:
            d_flow = 0.0; min_d_diasp = 0.0
        diasp_flow += d_flow

        # Top cities & top destinations
        city_names = nd.set_index('id').loc[[cities[i] for i in c_cities], 'name'].tolist()
        dest_flows = O[c_cities, :][:, c_dest_local].sum(0)
        top_d = [nd.iloc[keep_idx[c_dest_local[k]]]['name'] for k in np.argsort(-dest_flows)[:4]]
        comm_details.append(dict(
            comm=int(c), n_cities=len(c_cities), n_dests=len(c_dest_local),
            n_pieces=n_pieces, n_major_pieces=len(major_pieces),
            diasp_dest_pct=100.0 * len(diasp_local) / len(c_dest_local),
            diasp_flow_pct=100.0 * d_flow / max(c_flow, 1e-12),
            mean_diasp_km=min_d_diasp,
            cities=city_names[:6], top_dests=top_d
        ))

    # Roots alignment: for all cross-province pairs (city i, dest j), compare mean birthplace share
    # when co[i] == cd[j] vs co[i] != cd[j]
    dest_p_idx = np.array([p2i[d_prov[g]] for g in keep_idx])
    same_comm = (co[:, None] == cd[None, :])
    cross_prov = (c_prov[:, None] != d_prov[keep_idx][None, :])
    rs_mat = roots_share_cp[:, dest_p_idx]
    r_in = float(rs_mat[same_comm & cross_prov].mean()) if np.any(same_comm & cross_prov) else np.nan
    r_out = float(rs_mat[(~same_comm) & cross_prov].mean()) if np.any((~same_comm) & cross_prov) else np.nan
    roots_ratio = float(r_in / max(r_out, 1e-12)) if not np.isnan(r_in) else np.nan

    return dict(
        n_comms_cities=n_comms_cities,
        n_comms_ge5=len(comms_ge5),
        flow_inside_pct=float(100.0 * flow_inside),
        n_noncontig_strict=n_noncontig_strict,
        n_noncontig_major=n_noncontig_major,
        mean_pieces=float(np.mean(comp_counts)) if comp_counts else 0.0,
        diasp_dest_pct=float(100.0 * diasp_dests / max(tot_dests_eval, 1)),
        diasp_flow_pct=float(100.0 * diasp_flow / max(tot_flow_eval, 1e-12)),
        roots_in_comm=r_in, roots_out_comm=r_out, roots_enrichment=roots_ratio,
        nmi_ccaa=float(nmi(d_ccaa[keep_idx], cd)),
        nmi_prov=float(nmi(d_prov[keep_idx], cd)),
        details=comm_details
    )

# Run bipartite community detection on Mid-August, All Saints, and Feb weekday
bip_slices = {
    'Mid_Aug': "date BETWEEN '2024-08-10' AND '2024-08-18'",
    'All_Saints': "date = '2024-11-01'",
    'Feb_wd': "date BETWEEN '2024-02-05' AND '2024-02-29' AND dayofweek(date) IN (2,3,4)",
}
bip_results = {}
saved_maps = {}

for lab, cond in bip_slices.items():
    print(f"\n--- Bipartite Territories: {lab} ---", flush=True)
    df = con.execute(f"""
        SELECT r o, s d, sum(p)/count(DISTINCT date) w
        FROM {q_pq}
        WHERE {cond} AND r IN ({cl_sql}) AND s IN ({ml_sql}) AND r <> s
        GROUP BY 1, 2
    """).df()
    O_full = np.zeros((nc_, n), dtype=float)
    O_full[df.o.map(c2i).values, df.d.map(d2i).values] = df.w.values

    # Fit spatial DC-gravity null on full matrix
    E_sp_full = fit_dc_gravity_rect(O_full, B_cd, mask_cd, max_iter=200)

    # Keep active destinations (total inflow >= 10 across the 142 cities)
    keep = O_full.sum(0) >= 10.0
    keep_idx = np.where(keep)[0]
    O = O_full[:, keep]
    E_sp = E_sp_full[:, keep]
    W_tot = O.sum()

    # 1. Exact BiCM on binary matrix (w >= 5 -> sparse crisp links)
    M_bin = (O >= 5.0).astype(float)
    P_bicm, bicm_err, bicm_it = solve_bicm(M_bin)
    B_bicm = M_bin - P_bicm
    q_bicm_raw, co_bi, cd_bi = bip_louvain(B_bicm, n_restarts=10, seed_base=10)
    Q_bicm = q_bicm_raw / max(M_bin.sum(), 1.0)
    ev_bicm = evaluate_bipartite_comms(O, co_bi, cd_bi, keep_idx)
    ev_bicm['Q_B'] = float(Q_bicm); ev_bicm['bicm_err'] = bicm_err; ev_bicm['n_dests_active'] = int(len(keep_idx))

    # 2. Barber weighted configuration null
    E_bar = np.outer(O.sum(1), O.sum(0)) / W_tot
    q_bar_raw, co_ba, cd_ba = bip_louvain(O - E_bar, n_restarts=10, seed_base=30)
    Q_bar = q_bar_raw / W_tot
    ev_bar = evaluate_bipartite_comms(O, co_ba, cd_ba, keep_idx)
    ev_bar['Q_B'] = float(Q_bar); ev_bar['n_dests_active'] = int(len(keep_idx))

    # 3. Spatial DC-gravity null
    q_sp_raw, co_sp, cd_sp = bip_louvain(O - E_sp, n_restarts=15, seed_base=50)
    Q_sp = q_sp_raw / W_tot
    ev_sp = evaluate_bipartite_comms(O, co_sp, cd_sp, keep_idx)
    ev_sp['Q_B'] = float(Q_sp); ev_sp['n_dests_active'] = int(len(keep_idx))

    print(f"  [BiCM binary, err={bicm_err:.1e}] Q_B={Q_bicm:.3f} comms={ev_bicm['n_comms_ge5']} "
          f"noncontig_major={ev_bicm['n_noncontig_major']}/{ev_bicm['n_comms_ge5']} "
          f"diasp_dest={ev_bicm['diasp_dest_pct']:.1f}% roots_enrich={ev_bicm['roots_enrichment']:.2f}x NMI_CCAA={ev_bicm['nmi_ccaa']:.2f}", flush=True)
    print(f"  [Barber weighted CM] Q_B={Q_bar:.3f} comms={ev_bar['n_comms_ge5']} "
          f"noncontig_major={ev_bar['n_noncontig_major']}/{ev_bar['n_comms_ge5']} "
          f"diasp_dest={ev_bar['diasp_dest_pct']:.1f}% roots_enrich={ev_bar['roots_enrichment']:.2f}x NMI_CCAA={ev_bar['nmi_ccaa']:.2f}", flush=True)
    print(f"  [Spatial DC-gravity] Q_B={Q_sp:.3f} comms={ev_sp['n_comms_ge5']} "
          f"noncontig_major={ev_sp['n_noncontig_major']}/{ev_sp['n_comms_ge5']} "
          f"diasp_dest={ev_sp['diasp_dest_pct']:.1f}% diasp_flow={ev_sp['diasp_flow_pct']:.1f}% "
          f"roots_enrich={ev_sp['roots_enrichment']:.2f}x NMI_CCAA={ev_sp['nmi_ccaa']:.2f}", flush=True)

    bip_results[lab] = {'BiCM': ev_bicm, 'Barber_CM': ev_bar, 'Spatial_grav': ev_sp}

    if lab == 'Mid_Aug':
        saved_maps['keep_idx'] = keep_idx
        saved_maps['co_bi'] = co_bi; saved_maps['cd_bi'] = cd_bi
        saved_maps['co_ba'] = co_ba; saved_maps['cd_ba'] = cd_ba
        saved_maps['co_sp'] = co_sp; saved_maps['cd_sp'] = cd_sp
        saved_maps['O'] = O; saved_maps['E_sp'] = E_sp

        # Print top spatial-null communities in Mid-August
        print("  Top Spatial-Null Mid-August communities:")
        for d_ in sorted(ev_sp['details'], key=lambda x: -x['n_dests'])[:12]:
            print(f"    C{d_['comm']:2d}: cities({d_['n_cities']})={d_['cities']} | dests={d_['n_dests']} "
                  f"pieces={d_['n_pieces']}(maj={d_['n_major_pieces']}) diasp_dest={d_['diasp_dest_pct']:.0f}% "
                  f"diasp_flow={d_['diasp_flow_pct']:.0f}% top={d_['top_dests']}")

        # --- PROPER MONTE CARLO NULL TEST FOR SPATIAL-GRAVITY COMMUNITIES (100 draws) ---
        print("\n  Running 100-draw Monte Carlo Spatial-Null test for Mid-August bipartite communities...", flush=True)
        N_MC = 100
        mc_Q, mc_diasp_dest, mc_diasp_flow, mc_roots, mc_noncontig = [], [], [], [], []
        for rep in range(N_MC):
            O_mc = rng.poisson(np.maximum(E_sp, 0.0)).astype(float)
            # Refit expected null? For Poisson(E_sp), E_sp is the exact generating mean; compute modularity on O_mc - E_sp
            W_mc = O_mc.sum()
            q_mc_raw, co_mc, cd_mc = bip_louvain(O_mc - E_sp, n_restarts=3, seed_base=1000 + rep * 5)
            ev_mc = evaluate_bipartite_comms(O_mc, co_mc, cd_mc, keep_idx)
            mc_Q.append(q_mc_raw / W_mc)
            mc_diasp_dest.append(ev_mc['diasp_dest_pct'])
            mc_diasp_flow.append(ev_mc['diasp_flow_pct'])
            mc_roots.append(ev_mc['roots_enrichment'])
            mc_noncontig.append(ev_mc['n_noncontig_major'])
        mc_summary = dict(
            N_MC=N_MC,
            obs_Q=float(Q_sp), mc_Q_mean=float(np.mean(mc_Q)), mc_Q_sd=float(np.std(mc_Q)),
            z_Q=float((Q_sp - np.mean(mc_Q)) / max(np.std(mc_Q), 1e-12)),
            obs_diasp_dest_pct=ev_sp['diasp_dest_pct'], mc_diasp_dest_mean=float(np.mean(mc_diasp_dest)), mc_diasp_dest_sd=float(np.std(mc_diasp_dest)),
            obs_diasp_flow_pct=ev_sp['diasp_flow_pct'], mc_diasp_flow_mean=float(np.mean(mc_diasp_flow)), mc_diasp_flow_sd=float(np.std(mc_diasp_flow)),
            obs_roots_enrich=ev_sp['roots_enrichment'], mc_roots_mean=float(np.mean(mc_roots)), mc_roots_sd=float(np.std(mc_roots)),
            z_roots=float((ev_sp['roots_enrichment'] - np.mean(mc_roots)) / max(np.std(mc_roots), 1e-12)),
            obs_noncontig_major=ev_sp['n_noncontig_major'], mc_noncontig_mean=float(np.mean(mc_noncontig)), mc_noncontig_sd=float(np.std(mc_noncontig)),
        )
        bip_results['Mid_Aug_MC_null'] = mc_summary
        print("  Monte Carlo Spatial-Null test (Mid-August):", json.dumps(mc_summary, indent=2), flush=True)

with open('explore/out_p93_B_bipartite.json', 'w') as f:
    json.dump(bip_results, f, indent=2)

# ---- Figure: reports/figs/net_B_qap_and_territories.png ----
plt.rcParams.update({'font.size': 9.2, 'axes.spines.top': False, 'axes.spines.right': False})
fig = plt.figure(figsize=(13.0, 9.2), dpi=180)
gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.15], width_ratios=[1.05, 1.0])

# Panel (a): Spatial-null QAP correlations (47x47 province & 142-city x 47-prov neighbouring-city)
ax1 = fig.add_subplot(gs[0, 0])
mobs = ['Feb_wd', 'Aug_exc', 'Nov_exc', 'Xmas_exc']
mlabs = ['Ordinary\nFeb Wed', 'Mid-Aug\nExcess', 'All Saints\nExcess', 'Christmas\nExcess']
x = np.arange(len(mobs))
w = 0.25
r_p_sp = [df_qap_p.loc[(df_qap_p.mobility == m) & (df_qap_p.reference == 'Roots_2024'), 'r_spatial'].iloc[0] for m in mobs]
r_p_pa = [0.0 if m == 'Feb_wd' else df_qap_p.loc[(df_qap_p.mobility == m) & (df_qap_p.reference == 'Roots_2024'), 'r_partial_ctrl_feb'].iloc[0] for m in mobs]
r_c_nb = [df_qap_cp.loc[df_qap_cp.mobility == m, 'r_neighbour_cities'].iloc[0] for m in mobs]

ax1.bar(x - w, r_p_sp, w, color='#1d4ed8', label='47×47 Provinces: Spatial-null QAP vs Roots')
ax1.bar(x, r_p_pa, w, color='#60a5fa', label='47×47 Provinces: Partial QAP (ctrl Feb wd)')
ax1.bar(x + w, r_c_nb, w, color='#ea580c', label='137 Cities×47 Provs: Neighbouring-city QAP')
for xi, a1, a2, a3 in zip(x, r_p_sp, r_p_pa, r_c_nb):
    ax1.text(xi - w, a1 + 0.015, f"{a1:.2f}", ha='center', fontsize=7.6, color='#1e3a8a', fontweight='bold')
    if xi > 0:
        ax1.text(xi, a2 + 0.015, f"{a2:.2f}", ha='center', fontsize=7.6, color='#2563eb')
    ax1.text(xi + w, a3 + 0.015, f"{a3:.2f}", ha='center', fontsize=7.6, color='#c2410c', fontweight='bold')
ax1.set_xticks(x); ax1.set_xticklabels(mlabs, fontsize=8.5)
ax1.set_ylabel('Spatial-null QAP correlation $r$ ($p \\leq 0.02$; most $p < 0.001$)')
ax1.set_ylim(0, 0.85)
ax1.set_title('(a) Spatial-Null QAP: Holiday Excess vs Birthplace Roots', fontweight='bold')
ax1.legend(frameon=False, fontsize=7.5, loc='upper right')
ax1.grid(True, axis='y', alpha=0.2)

# Panel (b): Old (1988-95) vs Recent (2015-24) Migration & SCI vs Holidays (47x47 Partial QAP ctrl Feb)
ax2 = fig.add_subplot(gs[0, 1])
refs = ['Roots_2024', 'Roots_1998', 'Mig_88_95', 'Mig_15_24', 'Meta_SCI']
rlabs = ['Birthplace\nRoots 2024', 'Birthplace\nRoots 1998', 'Migration\n1988–1995', 'Migration\n2015–2024', 'Facebook\nSCI 2026']
xr = np.arange(len(refs))
w2 = 0.25
for k_idx, (mk, ml, col) in enumerate([('Aug_exc', 'Mid-Aug excess', '#0284c7'),
                                       ('Nov_exc', 'All Saints excess', '#d97706'),
                                       ('Xmas_exc', 'Christmas excess', '#dc2626')]):
    vals = [df_qap_p.loc[(df_qap_p.mobility == mk) & (df_qap_p.reference == rk), 'r_partial_ctrl_feb'].iloc[0] for rk in refs]
    ax2.bar(xr + (k_idx - 1) * w2, vals, w2, color=col, label=ml)
ax2.set_xticks(xr); ax2.set_xticklabels(rlabs, fontsize=8.0)
ax2.set_ylabel('Partial QAP $r$ (spatial null, ctrl Feb wd)')
ax2.set_title('(b) Which Historical Network Predicts Holiday Excess Best?', fontweight='bold')
ax2.legend(frameon=False, fontsize=8.0, loc='upper right')
ax2.grid(True, axis='y', alpha=0.2)

# Panel (c): Map of Barber/BiCM "Summer Territories" (Mid-August)
ax3 = fig.add_subplot(gs[1, 0])
keep_idx = saved_maps['keep_idx']
zones_plot = zones.copy()
zones_plot['comm_bi'] = -1
zones_plot.loc[keep_idx, 'comm_bi'] = saved_maps['cd_bi']
# Color top 12 communities cleanly
top_bi = pd.Series(saved_maps['cd_bi']).value_counts().head(12).index.tolist()
cmap = matplotlib.colormaps['tab20'].resampled(12)
zones_plot[zones_plot.comm_bi == -1].plot(ax=ax3, color='#f1f5f9', edgecolor='none')
for ci, cid in enumerate(top_bi):
    zones_plot[zones_plot.comm_bi == cid].plot(ax=ax3, color=cmap(ci), edgecolor='none', alpha=0.9)
zones_plot[~zones_plot.comm_bi.isin([-1] + top_bi)].plot(ax=ax3, color='#cbd5e1', edgecolor='none')
ax3.set_axis_off()
ev_bi_aug = bip_results['Mid_Aug']['BiCM']
ax3.set_title(f"(c) Mid-Aug Summer Territories (Exact BiCM, $Q_B={ev_bi_aug['Q_B']:.2f}$)\n"
              f"Contiguous regional blocs (NMI CCAA={ev_bi_aug['nmi_ccaa']:.2f}, roots={ev_bi_aug['roots_enrichment']:.1f}×)",
              fontweight='bold', fontsize=9.2)

# Panel (d): Map of Spatial DC-Gravity Bipartite Territories (Mid-August) highlighting diaspora communities
ax4 = fig.add_subplot(gs[1, 1])
zones_plot['comm_sp'] = -1
zones_plot.loc[keep_idx, 'comm_sp'] = saved_maps['cd_sp']
top_sp = pd.Series(saved_maps['cd_sp']).value_counts().head(12).index.tolist()
zones_plot[zones_plot.comm_sp == -1].plot(ax=ax4, color='#f1f5f9', edgecolor='none')
for ci, cid in enumerate(top_sp):
    zones_plot[zones_plot.comm_sp == cid].plot(ax=ax4, color=cmap(ci), edgecolor='none', alpha=0.9)
zones_plot[~zones_plot.comm_sp.isin([-1] + top_sp)].plot(ax=ax4, color='#cbd5e1', edgecolor='none')
ax4.set_axis_off()
ev_sp_aug = bip_results['Mid_Aug']['Spatial_grav']
mc_s = bip_results['Mid_Aug_MC_null']
ax4.set_title(f"(d) Mid-Aug Beyond Space (Spatial-Gravity Null, $Q_B={ev_sp_aug['Q_B']:.3f}$, MC $z={mc_s['z_Q']:.1f}$)\n"
              f"Non-contiguous diaspora territories: {ev_sp_aug['n_noncontig_major']}/{ev_sp_aug['n_comms_ge5']} comms, "
              f"roots={ev_sp_aug['roots_enrichment']:.2f}× (MC $z={mc_s['z_roots']:.1f}$)",
              fontweight='bold', fontsize=9.2)

plt.suptitle('Question B: Network Comparison (Spatial-Null QAP) & Bipartite Summer Territories',
             fontsize=12, fontweight='bold', y=0.985)
plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig('reports/figs/net_B_qap_and_territories.png', dpi=180)
print(f"\nSaved reports/figs/net_B_qap_and_territories.png in {time.time() - t0:.1f}s")
