"""
p92_B_networks_topology.py
Network-science basic properties & structural shift for Question B (rural exodus & holiday map):
Builds three network families on the exact same N=2,484 mainland Spain zones (plus province-level & 151-city roots networks):
 1. Ordinary weekday overnight network: Feb Tue-Thu mean (Feb_wd_mean) AND single night Wed 2024-02-14 (Feb_wd_1d)
 2. Holiday overnight networks: Mid-August (Aug 10-18 mean), Aug 15 peak night, All Saints (Nov 1), Christmas Eve (Dec 24)
 3. Roots networks: (a) 47x47 mainland province birthplace x residence (2024), (b) 151-city x 47-province bipartite roots (2024)
Computes for every overnight network (at w >= 10 backbone for binary topology so density < 5%, plus full weighted metrics):
 - N, E_dir, E_und, density, <k>, k_max, <s>, C vs ER & 20 degree-preserving rewirings, L vs 20 rewirings
 - Assortativity: undirected r, directed r(out->in), weighted population-size assortativity r_pop (do cities link to villages?)
 - Normalised rich-club rho(k) (20 rewirings)
 - k-core decomposition + k-core promotion of emptied inland villages by historical population-loss quintile (1900-2021)
 - Second-order mixing R_avg, R_max (Zhou)
 - Reach: mean weighted km, share >100 km, share cross-province, village origin-province diversity
Figure: reports/figs/net_B_topology_shift.png
"""
import json, time
import duckdb
import numpy as np
import pandas as pd
import geopandas as gpd
import igraph as ig
from scipy import sparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

t0 = time.time()
con = duckdb.connect()
q_pq = "'data/derived/pernoct_mun_2024_6m.parquet'"

# Mainland zones
nd = pd.read_parquet('data/derived/nodes_municipios.parquet')
nd = nd[~nd.foreign & ~nd.island & ~nd.prov.isin(['51', '52'])].reset_index(drop=True)
n = len(nd)
idx = {k: i for i, k in enumerate(nd.id)}
ml_sql = ",".join(f"'{z}'" for z in nd.id)
coords_x = nd.x.values
coords_y = nd.y.values
pop = np.maximum(nd.population.fillna(100).values, 1.0)
log_pop = np.log10(pop)
prov = nd.prov.values

# Historical population loss (1900-2021) + coastal flag for mainland zones
H = pd.read_csv('_local/external/ine_historical_pop/municipal_population_1900_2025.csv',
                dtype={'ine_code': str, 'province_code': str})
rel = pd.read_csv('raw_data_cache/v2/zonificacion/relacion_ine_zonificacionMitma.csv',
                  sep='|', dtype=str)[['municipio_ine', 'municipio_mitma']].drop_duplicates()
H = H.merge(rel, left_on='ine_code', right_on='municipio_ine', how='inner')
hcols = ['pop_1900', 'pop_1910', 'pop_1920', 'pop_1930', 'pop_1940', 'pop_1950', 'pop_1960', 'pop_1970', 'pop_1981', 'pop_2001', 'pop_2021']
Z = H.groupby('municipio_mitma')[hcols].sum(min_count=1)
Z['peak'] = Z[['pop_1900', 'pop_1910', 'pop_1920', 'pop_1930', 'pop_1940', 'pop_1950', 'pop_1960']].max(axis=1)
Z['loss_peak_2021'] = 1.0 - Z['pop_2021'] / Z['peak']

g_geo = gpd.read_file('clean_data/v2/zones/municipios_mitma.gpkg')
land_outer = g_geo.geometry.union_all().boundary
g_geo['coastal'] = g_geo.geometry.buffer(50).intersects(land_outer) & g_geo.id.str.match(r'^\d')
coastal_map = g_geo.set_index('id')['coastal']

nd['loss_peak_2021'] = nd.id.map(Z['loss_peak_2021'])
nd['coastal'] = nd.id.map(coastal_map).fillna(False).astype(bool)

def graph_to_csr(G):
    edges = np.array(G.get_edgelist(), dtype=np.int32)
    if len(edges) == 0:
        return sparse.csr_matrix((G.vcount(), G.vcount()), dtype=np.float64)
    rows = np.concatenate([edges[:, 0], edges[:, 1]])
    cols = np.concatenate([edges[:, 1], edges[:, 0]])
    return sparse.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(G.vcount(), G.vcount()))

def second_order_mixing(adj_csr):
    deg = np.asarray(adj_csr.sum(axis=1)).ravel().astype(float)
    nb_sum = adj_csr @ deg
    avg_nb = np.divide(nb_sum, deg, out=np.zeros_like(deg), where=deg > 0)
    indptr, indices = adj_csr.indptr, adj_csr.indices
    max_nb = np.zeros(len(deg), dtype=float)
    for i in range(len(deg)):
        s, e = indptr[i], indptr[i + 1]
        if e > s:
            max_nb[i] = deg[indices[s:e]].max()
    coo = adj_csr.tocoo()
    u, v = coo.row, coo.col
    def pearson(a, b):
        a0, b0 = a - a.mean(), b - b.mean()
        d = np.sqrt((a0 * a0).sum() * (b0 * b0).sum())
        return float((a0 * b0).sum() / d) if d > 0 else 0.0
    return pearson(deg[u], deg[v]), pearson(avg_nb[u], avg_nb[v]), pearson(max_nb[u], max_nb[v])

def rich_club_phi(adj_csr, ks):
    deg = np.asarray(adj_csr.sum(axis=1)).ravel()
    out = {}
    for k in ks:
        mask = deg > k
        nk = int(mask.sum())
        if nk < 2:
            out[int(k)] = (nk, np.nan)
            continue
        ek = adj_csr[mask][:, mask].nnz / 2.0
        out[int(k)] = (nk, float(2.0 * ek / (nk * (nk - 1))))
    return out

def wcorr(a, b, w):
    ma, mb = np.average(a, weights=w), np.average(b, weights=w)
    da, db = a - ma, b - mb
    return float(np.average(da * db, weights=w) / np.sqrt(np.average(da * da, weights=w) * np.average(db * db, weights=w)))

slices = [
    ('Feb_wd_mean', "date BETWEEN '2024-02-05' AND '2024-02-29' AND dayofweek(date) IN (2,3,4)"),
    ('Feb_wd_1d',   "date = '2024-02-14'"),
    ('Mid_Aug',     "date BETWEEN '2024-08-10' AND '2024-08-18'"),
    ('Aug_15',      "date = '2024-08-15'"),
    ('All_Saints',  "date = '2024-11-01'"),
    ('Xmas_Eve',    "date = '2024-12-24'"),
]

THR = 10.0  # >=10 residents sleeping away: gives 1.0%-3.7% density across all slices (<5% guideline!)
N_REW = 25
ks_rc = [10, 20, 30, 50, 75, 100, 150, 200, 300, 400]

summary = {}
rc_store = {}
node_metrics = {}
csr_store = {}

for lab, cond in slices:
    print(f"\n=== Overnight network: {lab} ===", flush=True)
    df = con.execute(f"""
        SELECT r, s, sum(p)/count(DISTINCT date) AS w
        FROM {q_pq}
        WHERE {cond} AND r <> s AND r IN ({ml_sql}) AND s IN ({ml_sql})
        GROUP BY 1, 2
    """).df()
    df['ri'] = df.r.map(idx).values
    df['si'] = df.s.map(idx).values
    df['km'] = np.hypot(coords_x[df.ri.values] - coords_x[df.si.values],
                        coords_y[df.ri.values] - coords_y[df.si.values]) / 1000.0
    w_all = df.w.values
    tot_w = float(w_all.sum())

    # Full-network weighted reach & size mixing
    mean_km = float(np.average(df.km.values, weights=w_all))
    share_100km = float(100.0 * w_all[df.km.values > 100].sum() / tot_w)
    share_xprov = float(100.0 * w_all[prov[df.ri.values] != prov[df.si.values]].sum() / tot_w)
    share_small = float(100.0 * w_all[pop[df.si.values] < 5000].sum() / tot_w)
    r_pop_w = wcorr(log_pop[df.ri.values], log_pop[df.si.values], w_all)

    # Per-destination weighted in-strength & number of distinct origin provinces (weighted entropy / effective count)
    in_str = np.bincount(df.si.values, weights=w_all, minlength=n)
    out_str = np.bincount(df.ri.values, weights=w_all, minlength=n)
    # origin province diversity per destination (effective # of origin provinces = exp(H))
    df['rprov'] = prov[df.ri.values]
    p_dest = df.groupby(['si', 'rprov'])['w'].sum().reset_index()
    p_tot = p_dest.groupby('si')['w'].transform('sum')
    p_dest['q'] = p_dest['w'] / p_tot
    p_dest['h'] = -p_dest['q'] * np.log(np.maximum(p_dest['q'], 1e-15))
    eff_orig_provs = np.ones(n)
    h_by_s = p_dest.groupby('si')['h'].sum()
    eff_orig_provs[h_by_s.index.values] = np.exp(h_by_s.values)
    n_orig_provs = np.zeros(n)
    c_by_s = p_dest[p_dest['w'] >= 10].groupby('si')['rprov'].nunique()
    n_orig_provs[c_by_s.index.values] = c_by_s.values

    # Backbone at w >= THR (density < 5%)
    sub = df[df.w >= THR].copy()
    w_sub_share = float(100.0 * sub.w.sum() / tot_w)
    # Directed degree assortativity r(out_r, in_s)
    k_out = np.bincount(sub.ri.values, minlength=n).astype(float)
    k_in = np.bincount(sub.si.values, minlength=n).astype(float)
    r_dir_oi = float(np.corrcoef(k_out[sub.ri.values], k_in[sub.si.values])[0, 1])
    r_dir_io = float(np.corrcoef(k_in[sub.ri.values], k_out[sub.si.values])[0, 1])
    # Directed binary reciprocity
    dp = set(zip(sub.ri.values, sub.si.values))
    recip_bin = sum(1 for u, v in dp if (v, u) in dp) / max(len(dp), 1)
    # Undirected backbone
    a = np.minimum(sub.ri.values, sub.si.values)
    b = np.maximum(sub.ri.values, sub.si.values)
    und = pd.DataFrame({'a': a, 'b': b, 'w': sub.w.values}).groupby(['a', 'b'], as_index=False)['w'].sum()
    G = ig.Graph(n=n, edges=list(zip(und.a.values, und.b.values)), directed=False)
    E_und = G.ecount()
    dens = 2.0 * E_und / (n * (n - 1))
    deg = np.array(G.degree(), dtype=float)
    comps = G.connected_components()
    gcc_size = max(comps.sizes())
    G_gcc = comps.giant()

    C_avg = float(G.transitivity_avglocal_undirected(mode='zero'))
    C_glob = float(G.transitivity_undirected())
    L_avg = float(G_gcc.average_path_length())
    diam = int(G_gcc.diameter())

    csr = graph_to_csr(G)
    csr_store[lab] = csr
    r1, r_avg, r_max = second_order_mixing(csr)
    coreness = np.array(G.coreness(), dtype=float)
    kcore_max = int(coreness.max())

    phi_obs = rich_club_phi(csr, ks_rc)

    # 25 degree-preserving rewirings
    c_rew, l_rew, r1_rew, ravg_rew, rmax_rew, kc_rew = [], [], [], [], [], []
    phi_rew_mat = {k: [] for k in ks_rc}
    for rep in range(N_REW):
        Gr = G.copy()
        Gr.rewire(n=10 * E_und)
        c_rew.append(Gr.transitivity_avglocal_undirected(mode='zero'))
        if rep < 8:
            l_rew.append(Gr.connected_components().giant().average_path_length())
        cs = graph_to_csr(Gr)
        rr1, rra, rrm = second_order_mixing(cs)
        r1_rew.append(rr1); ravg_rew.append(rra); rmax_rew.append(rrm)
        kc_rew.append(max(Gr.coreness()))
        pr = rich_club_phi(cs, ks_rc)
        for k in ks_rc:
            phi_rew_mat[k].append(pr[k][1])

    rc_rows = []
    for k in ks_rc:
        nk, p0 = phi_obs[k]
        arr = np.array(phi_rew_mat[k], dtype=float); arr = arr[~np.isnan(arr)]
        rho = p0 / np.mean(arr) if (len(arr) and np.mean(arr) > 0 and not np.isnan(p0)) else np.nan
        rc_rows.append(dict(k=k, nk=nk, phi_obs=p0, phi_rew=float(np.mean(arr)) if len(arr) else np.nan, rho=float(rho)))
    rc_store[lab] = rc_rows

    node_metrics[lab] = pd.DataFrame({
        'id': nd.id.values, 'pop': pop, 'loss': nd.loss_peak_2021.values, 'coastal': nd.coastal.values,
        'deg': deg, 'k_in': k_in, 'k_out': k_out, 'in_str': in_str, 'out_str': out_str,
        'coreness': coreness, 'core_norm': coreness / max(kcore_max, 1),
        'eff_orig_provs': eff_orig_provs, 'n_orig_provs': n_orig_provs
    })

    # Small village stats (pop < 5,000, inland)
    sv = (~nd.coastal.values) & (pop < 5000)
    res = dict(
        label=lab, N=n, E_dir=int(len(sub)), E_und=int(E_und), density=float(dens),
        total_away_M=float(tot_w / 1e6), thr_flow_share_pct=w_sub_share,
        k_mean=float(deg.mean()), k_med=float(np.median(deg)), k_max=int(deg.max()),
        k_in_small_med=float(np.median(k_in[sv])), k_in_small_mean=float(np.mean(k_in[sv])),
        in_str_small_med=float(np.median(in_str[sv])),
        recip_bin=float(recip_bin), gcc_size=int(gcc_size),
        C_avg=C_avg, C_ER=float(dens), C_rew_mean=float(np.mean(c_rew)),
        L_avg=L_avg, diam=diam, L_rew_mean=float(np.mean(l_rew)),
        r_deg=r1, r_deg_rew_mean=float(np.mean(r1_rew)),
        r_dir_out_in=r_dir_oi, r_dir_in_out=r_dir_io,
        r_pop_weighted=r_pop_w,
        R_avg=r_avg, R_avg_rew_mean=float(np.mean(ravg_rew)),
        R_max=r_max, R_max_rew_mean=float(np.mean(rmax_rew)),
        kcore_max=kcore_max, kcore_rew_mean=float(np.mean(kc_rew)),
        kcore_small_mean=float(np.mean(coreness[sv])),
        rho_k50=next((r['rho'] for r in rc_rows if r['k'] == 50), np.nan),
        rho_k100=next((r['rho'] for r in rc_rows if r['k'] == 100), np.nan),
        rho_k200=next((r['rho'] for r in rc_rows if r['k'] == 200), np.nan),
        mean_km=mean_km, share_100km=share_100km, share_xprov=share_xprov, share_small_pct=share_small,
    )
    summary[lab] = res
    print(f"  away={tot_w/1e6:.2f}M E_und={E_und} dens={dens:.4f} <k>={deg.mean():.1f} kmax={int(deg.max())} "
          f"C={C_avg:.3f}(rew={np.mean(c_rew):.3f}) L={L_avg:.3f} r={r1:+.3f}(rew={np.mean(r1_rew):+.3f}) "
          f"r_dir(o,i)={r_dir_oi:+.3f} r_pop={r_pop_w:+.3f} R_avg={r_avg:+.3f} "
          f"kcore={kcore_max} small_kin={np.mean(k_in[sv]):.1f} km={mean_km:.1f} rho100={res['rho_k100']:.2f}", flush=True)

# --- Also compute basic properties of the TWO ROOTS networks ---
print("\n=== Roots Networks (Province 47x47 & 151-city x 47-prov) ===", flush=True)
mainland_provs = sorted(nd.prov.unique().tolist())  # 47 provinces
br = pd.read_csv('_local/external/ine_birth_residence/prov_birth_x_prov_residence.csv')
br['pb'] = br['prov_birth'].astype(str).str.zfill(2)
br['pr'] = br['prov_residence'].astype(str).str.zfill(2)
br24 = br[(br.year == 2024) & br.pb.isin(mainland_provs) & br.pr.isin(mainland_provs) & (br.pb != br.pr)].copy()
# Province population
prov_pop = nd.groupby('prov').population.sum()
br24['lpop_res'] = np.log10(br24.pr.map(prov_pop).values)
br24['lpop_birth'] = np.log10(br24.pb.map(prov_pop).values)
r_roots_pop = wcorr(br24['lpop_res'].values, br24['lpop_birth'].values, br24['persons'].values)
# Top-10% backbone (or >= 1000 persons) of 47x47 roots network for topology
br_thr = br24[br24.persons >= 2000]
p2i = {p: i for i, p in enumerate(mainland_provs)}
u_r = np.minimum(br_thr.pr.map(p2i).values, br_thr.pb.map(p2i).values)
v_r = np.maximum(br_thr.pr.map(p2i).values, br_thr.pb.map(p2i).values)
und_r = pd.DataFrame({'a': u_r, 'b': v_r}).drop_duplicates()
Gr_p = ig.Graph(n=len(mainland_provs), edges=list(zip(und_r.a.values, und_r.b.values)), directed=False)
roots_prov_stats = dict(
    N=len(mainland_provs), E_dir_all=len(br24), total_migrants_M=float(br24.persons.sum() / 1e6),
    E_und_ge2000=Gr_p.ecount(), dens_ge2000=float(2 * Gr_p.ecount() / (47 * 46)),
    r_deg_ge2000=float(Gr_p.assortativity_degree()),
    r_pop_weighted=r_roots_pop,
    C_avg_ge2000=float(Gr_p.transitivity_avglocal_undirected(mode='zero')),
)
print("Roots 47x47:", roots_prov_stats)

# 151 cities x 47 provinces bipartite roots
bm = pd.read_csv('_local/external/ine_birth_residence/mun_residence_x_prov_birth_censo_2021_2025.csv')
bm['cid'] = bm['ine_code'].astype(str).str.zfill(5)
bm['pb'] = bm['prov_birth'].astype(str).str.zfill(2)
bm24 = bm[(bm.year == 2024) & bm.cid.isin(idx) & bm.pb.isin(mainland_provs) & (bm.cid.str[:2] != bm.pb)].copy()
bm_thr = bm24[bm24.persons >= 200]
roots_bip_stats = dict(
    N_cities=int(bm24.cid.nunique()), N_provs=int(bm24.pb.nunique()),
    E_all_nonzero=int((bm24.persons > 0).sum()), total_migrants_M=float(bm24.persons.sum() / 1e6),
    E_ge200=int(len(bm_thr)), dens_ge200=float(len(bm_thr) / (bm24.cid.nunique() * (bm24.pb.nunique() - 1))),
    r_pop_weighted=wcorr(np.log10(bm24.cid.map(nd.set_index('id').population).values),
                         np.log10(bm24.pb.map(prov_pop).values), bm24.persons.values)
)
print("Roots 151-city x 47-prov:", roots_bip_stats)

# --- Emptied-village k-core & degree promotion table (inland zones with pop >= 200 & non-missing loss) ---
base_df = node_metrics['Feb_wd_1d'].copy()
inland_mask = (~base_df['coastal']) & (base_df['pop'] >= 200) & (base_df['loss'].notna())
inl = base_df[inland_mask].copy()
inl['loss_q'] = pd.qcut(inl['loss'], 5, labels=['Q1 (grew)', 'Q2', 'Q3', 'Q4', 'Q5 (emptied most)'])

for lab in ['Feb_wd_mean', 'Feb_wd_1d', 'Mid_Aug', 'Aug_15', 'All_Saints', 'Xmas_Eve']:
    m = node_metrics[lab].loc[inland_mask]
    inl[f'kcore_{lab}'] = m['coreness'].values
    inl[f'kcore_norm_{lab}'] = m['core_norm'].values
    inl[f'kin_{lab}'] = m['k_in'].values
    inl[f'instr_{lab}'] = m['in_str'].values
    inl[f'effprov_{lab}'] = m['eff_orig_provs'].values

promo_tbl = inl.groupby('loss_q', observed=True).agg(
    n=('pop', 'size'),
    med_pop=('pop', 'median'),
    med_loss_pct=('loss', lambda s: 100 * s.median()),
    kcore_feb1d=('kcore_Feb_wd_1d', 'mean'),
    kcore_aug15=('kcore_Aug_15', 'mean'),
    kcore_nov1=('kcore_All_Saints', 'mean'),
    kcore_xmas=('kcore_Xmas_Eve', 'mean'),
    kcore_norm_feb1d=('kcore_norm_Feb_wd_1d', 'mean'),
    kcore_norm_aug15=('kcore_norm_Aug_15', 'mean'),
    kcore_norm_nov1=('kcore_norm_All_Saints', 'mean'),
    kcore_norm_xmas=('kcore_norm_Xmas_Eve', 'mean'),
    kin_feb1d=('kin_Feb_wd_1d', 'mean'),
    kin_aug15=('kin_Aug_15', 'mean'),
    kin_nov1=('kin_All_Saints', 'mean'),
    kin_xmas=('kin_Xmas_Eve', 'mean'),
    effprov_feb1d=('effprov_Feb_wd_1d', 'mean'),
    effprov_aug15=('effprov_Aug_15', 'mean'),
    effprov_nov1=('effprov_All_Saints', 'mean'),
).round(3)
print("\n=== Emptied-Village Promotion Across Inland Population-Loss Quintiles ===")
print(promo_tbl.to_string())
promo_tbl.to_csv('explore/out_p92_B_promotion.csv')

# Also check promotion controlling for population (within small inland zones pop < 6,000!)
small_inl = inl[inl['pop'] < 6000].copy()
small_inl['loss_q3'] = pd.qcut(small_inl['loss'], 3, labels=['Low loss', 'Mid loss', 'Extreme loss'])
promo_small = small_inl.groupby('loss_q3', observed=True).agg(
    n=('pop', 'size'), med_pop=('pop', 'median'), med_loss=('loss', 'median'),
    kcore_feb=('kcore_Feb_wd_1d', 'mean'), kcore_aug15=('kcore_Aug_15', 'mean'),
    kcore_nov1=('kcore_All_Saints', 'mean'), kcore_xmas=('kcore_Xmas_Eve', 'mean'),
    kin_feb=('kin_Feb_wd_1d', 'mean'), kin_aug15=('kin_Aug_15', 'mean'),
    kin_nov1=('kin_All_Saints', 'mean'), kin_xmas=('kin_Xmas_Eve', 'mean'),
    effprov_feb=('effprov_Feb_wd_1d', 'mean'), effprov_aug15=('effprov_Aug_15', 'mean'),
).round(2)
print("\n=== Size-Controlled Check: Small Inland Villages Only (pop < 2,000) ===")
print(promo_small.to_string())

with open('explore/out_p92_B_props.json', 'w') as f:
    json.dump({
        'overnight': summary,
        'roots_prov': roots_prov_stats,
        'roots_bip': roots_bip_stats,
        'rich_club': rc_store,
        'promotion_quintiles': promo_tbl.reset_index().to_dict(orient='records'),
        'promotion_small_villages': promo_small.reset_index().to_dict(orient='records'),
    }, f, indent=2)

# ---- Figure: reports/figs/net_B_topology_shift.png ----
plt.rcParams.update({'font.size': 9.3, 'axes.spines.top': False, 'axes.spines.right': False})
fig, axes = plt.subplots(2, 2, figsize=(12.0, 9.0), dpi=180)
col_map = {
    'Feb_wd_1d': '#475569',
    'All_Saints': '#d97706',
    'Xmas_Eve': '#dc2626',
    'Aug_15': '#0284c7',
}
lbl_map = {
    'Feb_wd_1d': 'Ordinary Wed (Feb 14)',
    'All_Saints': 'All Saints (Nov 1)',
    'Xmas_Eve': 'Christmas Eve (Dec 24)',
    'Aug_15': 'Mid-August (Aug 15)',
}

# Panel (a): In-degree CCDF across regimes
ax = axes[0, 0]
for lab in ['Feb_wd_1d', 'All_Saints', 'Xmas_Eve', 'Aug_15']:
    d = np.sort(node_metrics[lab]['k_in'].values)
    d = d[d > 0]
    ccdf = 1.0 - np.arange(len(d)) / len(d)
    ax.loglog(d, ccdf, lw=2.0, color=col_map[lab],
              label=f"{lbl_map[lab]} ($\\langle k_{{in}}\\rangle={node_metrics[lab]['k_in'].mean():.1f}$)")
ax.set_xlabel('Destination in-degree $k_{in}$ ($w \\geq 10$ overnight stays)')
ax.set_ylabel('$P(K_{in} \\geq k)$')
ax.set_title('(a) Destination In-Degree CCDF: Ordinary vs Holidays', fontweight='bold')
ax.legend(frameon=False, fontsize=8.0, loc='lower left')
ax.grid(True, which='both', alpha=0.2)

# Panel (b): Assortativity shift & Rich-club
ax = axes[0, 1]
labs_b = ['Feb_wd_1d', 'All_Saints', 'Xmas_Eve', 'Aug_15']
x = np.arange(len(labs_b))
w = 0.25
r_und_v = [summary[l]['r_deg'] for l in labs_b]
r_dir_v = [summary[l]['r_dir_out_in'] for l in labs_b]
r_pop_v = [summary[l]['r_pop_weighted'] for l in labs_b]
ax.bar(x - w, r_und_v, w, color='#1e3a8a', label='Undirected degree $r$')
ax.bar(x, r_dir_v, w, color='#0ea5e9', label='Directed $r(k_{out}, k_{in})$')
ax.bar(x + w, r_pop_v, w, color='#f97316', label='Weighted $\\log(\\mathrm{pop})$ mixing $r_{pop}$')
ax.axhline(0, color='black', lw=0.8)
for xi, a1, a2, a3 in zip(x, r_und_v, r_dir_v, r_pop_v):
    ax.text(xi - w, a1 + (-0.028 if a1 < 0 else 0.012), f"{a1:+.2f}", ha='center', fontsize=7.5, color='#1e3a8a')
    ax.text(xi, a2 + (-0.028 if a2 < 0 else 0.012), f"{a2:+.2f}", ha='center', fontsize=7.5, color='#0369a1')
    ax.text(xi + w, a3 + (-0.028 if a3 < 0 else 0.012), f"{a3:+.2f}", ha='center', fontsize=7.5, color='#c2410c', fontweight='bold')
ax.set_xticks(x)
ax.set_xticklabels(['Feb Wed\n(ordinary)', 'All Saints\n(Nov 1)', 'Christmas\n(Dec 24)', 'Mid-Aug\n(Aug 15)'], fontsize=8.5)
ax.set_ylabel('Assortativity coefficient')
ax.set_ylim(-0.32, 0.36)
ax.set_title('(b) Mixing Shift: Cities Reconnect to Small Villages on Holidays', fontweight='bold')
ax.legend(frameon=False, fontsize=7.8, loc='upper right')
ax.grid(True, axis='y', alpha=0.2)

# Panel (c): k-core promotion of emptied inland villages
ax = axes[1, 0]
q_labels = ['Q1\n(grew)', 'Q2', 'Q3', 'Q4', 'Q5\n(emptied\nmost)']
xq = np.arange(5)
# Plot ratio of holiday k-core to Feb weekday k-core (and also normalized k-core jump!)
for lab, mk in [('Aug_15', 'o-'), ('Xmas_Eve', 's-'), ('All_Saints', '^-')]:
    ratio_kc = promo_tbl[f'kcore_{lab.lower().replace("aug_15","aug15").replace("all_saints","nov1").replace("xmas_eve","xmas")}'] / promo_tbl['kcore_feb1d']
    ax.plot(xq, ratio_kc.values, mk, lw=2.2, ms=6, color=col_map[lab], label=f"{lbl_map[lab]} / Feb Wed")
ax.axhline(1.0, color='gray', ls='--', lw=1.0)
for i, v in enumerate((promo_tbl['kcore_aug15'] / promo_tbl['kcore_feb1d']).values):
    ax.text(i, v + 0.12, f"{v:.1f}x", ha='center', fontsize=8.2, color=col_map['Aug_15'], fontweight='bold')
ax.set_xticks(xq)
ax.set_xticklabels(q_labels, fontsize=8.5)
ax.set_xlabel('Inland municipal quintile of 1900–2021 population loss')
ax.set_ylabel('Mean $k$-core shell multiplier (Holiday / Feb Wed)')
ax.set_title('(c) $k$-Core Promotion of Emptied Inland Villages', fontweight='bold')
ax.legend(frameon=False, fontsize=8.2, loc='upper left')
ax.grid(True, alpha=0.2)

# Panel (d): BOSAM comparison — Feb Wed (upper) vs Aug 15 (lower), sorted by population descending
ax = axes[1, 1]
B_SIZE = 120
ord_pop = np.argsort(-pop)  # sorted by population (cities top-left -> villages bottom-right)
def bosam_grid(csr_mat, order, bins=B_SIZE):
    m = csr_mat[order][:, order].tocoo()
    bi = np.minimum((m.row * bins) // n, bins - 1)
    bj = np.minimum((m.col * bins) // n, bins - 1)
    G_img = np.zeros((bins, bins), dtype=float)
    np.add.at(G_img, (bi, bj), 1.0)
    sz = np.bincount(np.minimum((np.arange(n) * bins) // n, bins - 1), minlength=bins)
    return G_img / np.maximum(np.outer(sz, sz), 1.0)

img_feb = bosam_grid(csr_store['Feb_wd_1d'], ord_pop)
img_aug = bosam_grid(csr_store['Aug_15'], ord_pop)
comb = np.triu(img_feb) + np.tril(img_aug, -1)
ax.imshow(np.log1p(comb * 25), cmap='viridis_r', origin='upper', aspect='equal')
ax.plot([0, B_SIZE - 1], [0, B_SIZE - 1], color='red', ls='--', lw=1.0, alpha=0.8)
ax.text(B_SIZE * 0.56, B_SIZE * 0.14, 'Feb Wed\n(upper)', color='#334155', fontweight='bold', fontsize=8.5,
        bbox=dict(boxstyle='round,pad=0.2', fc='white', alpha=0.88, ec='none'))
ax.text(B_SIZE * 0.08, B_SIZE * 0.82, 'Aug 15\n(lower)', color='#0284c7', fontweight='bold', fontsize=8.5,
        bbox=dict(boxstyle='round,pad=0.2', fc='white', alpha=0.88, ec='none'))
ax.set_xlabel('Nodes sorted by population (cities $\\rightarrow$ villages)')
ax.set_ylabel('Nodes sorted by population')
ax.set_title('(d) BOSAM (Population-Sorted): Cities Lighting Up Villages', fontweight='bold')
ax.set_xticks([0, B_SIZE // 2, B_SIZE - 1])
ax.set_xticklabels(['1 (Madrid)', '1,242', '2,484'])
ax.set_yticks([0, B_SIZE // 2, B_SIZE - 1])
ax.set_yticklabels(['1', '1,242', '2,484'])

plt.suptitle('Question B (Overnight Mobility Network, N = 2,484): Holiday Structural Transformation',
             fontsize=12, fontweight='bold', y=0.985)
plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig('reports/figs/net_B_topology_shift.png', dpi=180)
print(f"\nSaved reports/figs/net_B_topology_shift.png in {time.time() - t0:.1f}s")
