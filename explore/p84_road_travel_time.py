# Check 4: Long-distance road travel time for 151 cities x all 2,484 mainland destination zones
# Method:
#   (1) Shortest-path Dijkstra on the 2,484-node, 677,116-edge exact OSM <=150km road graph (road_tt_150km.parquet),
#       with empirically calibrated intermediate-hop urban-access penalty (calibrated directly on 100-150km pairs
#       where both exact 1-hop OSM tt_min and best 2-hop chained tt_min are observed!).
#   (2) Empirical circuity + speed approximation fitted on road_tt_150km.parquet (for comparison).
#   (3) Re-estimate the city roots PPML model at BOTH:
#       - City x destination-province (population-weighted road travel time from city to province),
#       - City x all 2,484 mainland destination zones (exact zone-to-zone road travel time bins).
import sys, time
sys.path.insert(0, 'explore')
import duckdb, numpy as np, pandas as pd, pyfixest as pf
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra
from common import nodes

t0 = time.time()
tt = pd.read_parquet('data/derived/road_tt_150km.parquet')
B = pd.read_csv('_local/external/ine_birth_residence/mun_residence_x_prov_birth_censo_2021_2025.csv',
                dtype={'ine_code': str, 'prov_birth': str})
B24 = B[(B.year == 2024) & B.prov_birth.str.match(r'^\d\d$')].copy()

N = nodes().set_index('id')
M = N[(~N.foreign) & (~N.island) & (~N.prov.isin(['51', '52']))].copy()
main_zones = sorted(set(tt.o) | set(tt.d))
assert len(main_zones) == 2484
z2i = {z: i for i, z in enumerate(main_zones)}
n_z = len(main_zones)

cities_all = sorted(set(B24.ine_code) & set(N.index))
cities_main = sorted(set(cities_all) & set(main_zones))
print(f"Cities matched: {len(cities_all)} total, {len(cities_main)} mainland (excluded islands/Ceuta/Melilla: {sorted(set(cities_all)-set(cities_main))})")

# 1. Calibrate per-hop urban-access penalty on 100-150km pairs (exact 1-hop vs best 2-hop via midpoint)
g_tt = csr_matrix((tt.tt_min.values.astype(np.float64), (tt.o.map(z2i).values, tt.d.map(z2i).values)), shape=(n_z, n_z))
g_km = csr_matrix((tt.road_km.values.astype(np.float64), (tt.o.map(z2i).values, tt.d.map(z2i).values)), shape=(n_z, n_z))

# Sample 200 pairs at 120-150km to compare direct 1-hop vs forced 2-hop (removing direct edge > 85km)
tt_short = tt[tt.eucl_km <= 90.0]
g_tt_90 = csr_matrix((tt_short.tt_min.values.astype(np.float64), (tt_short.o.map(z2i).values, tt_short.d.map(z2i).values)), shape=(n_z, n_z))
long150 = tt[(tt.eucl_km >= 120.0) & (tt.eucl_km <= 150.0) & (tt.o.isin(cities_main))].copy()
src_u = sorted(set(long150.o.map(z2i)))
d90 = dijkstra(g_tt_90, directed=True, indices=src_u)
smap = {s: k for k, s in enumerate(src_u)}
two_hop = np.array([d90[smap[z2i[o]], z2i[d]] for o, d in zip(long150.o, long150.d)])
diff = two_hop - long150.tt_min.values
valid = np.isfinite(diff) & (diff >= 0) & (diff < 60)
hop_pen_med = float(np.median(diff[valid]))
hop_pen_mean = float(np.mean(diff[valid]))
print(f"Empirical 1-extra-hop urban-access penalty (n={valid.sum()} pairs at 120-150km): median={hop_pen_med:.2f} min, mean={hop_pen_mean:.2f} min, p25={np.percentile(diff[valid],25):.2f}, p75={np.percentile(diff[valid],75):.2f}")

# 2. Compute Dijkstra shortest paths from ALL 138 mainland cities to ALL 2,484 mainland zones!
# Also compute hop count and road_km along fastest-time route.
city_idx = [z2i[c] for c in cities_main]
# Subtract hop_pen_med from every edge (floored at 0.5*tt_min) when routing multi-hop? Better: route on raw g_tt,
# count hops h, and define tt_adj = tt_raw - (h - 1) * hop_pen_med (exact for h=1, adjusts off-motorway dips for h>=2).
D_tt, pred = dijkstra(g_tt, directed=True, indices=city_idx, return_predecessors=True)

# Vectorized hop count and road_km extraction
hops = np.zeros_like(D_tt, dtype=np.int16)
D_rkm = np.zeros_like(D_tt, dtype=np.float32)
# Build fast dict/lookup for edge road_km
edge_km = {}
for oi, di, rk in zip(tt.o.map(z2i).values, tt.d.map(z2i).values, tt.road_km.values):
    edge_km[(int(oi), int(di))] = float(rk)

for ci_pos, ci in enumerate(city_idx):
    p_row = pred[ci_pos]
    for j in range(n_z):
        if j == ci:
            continue
        cur = j
        h = 0
        rk = 0.0
        while cur != ci and cur >= 0:
            prv = int(p_row[cur])
            rk += edge_km[(prv, cur)]
            cur = prv
            h += 1
        hops[ci_pos, j] = h
        D_rkm[ci_pos, j] = rk

D_tt_adj = D_tt - np.maximum(hops - 1, 0) * hop_pen_med

# Sanity check major corridors
print("\n=== Long-Distance Road Travel Time Sanity Checks (Chained OSM vs Hop-Adjusted vs Real Car) ===")
checks = [
    ('28079','08019','Madrid-Barcelona (~620km, ~375m)'),
    ('28079','41091','Madrid-Sevilla (~530km, ~315m)'),
    ('28079','46250','Madrid-Valencia (~355km, ~210m)'),
    ('28079','48020','Madrid-Bilbao (~400km, ~240m)'),
    ('28079','15030','Madrid-A Coruna (~595km, ~345m)'),
    ('08019','41091','Barcelona-Sevilla (~995km, ~570m)'),
    ('08019','29067','Barcelona-Malaga (~990km, ~570m)'),
    ('08019','48020','Barcelona-Bilbao (~610km, ~350m)'),
    ('01059','10037','Vitoria-Caceres (~560km, ~330m)'),
    ('48020','11012','Bilbao-Cadiz (~990km, ~570m)'),
]
c2pos = {c: i for i, c in enumerate(cities_main)}
for a, b, ref in checks:
    if a in c2pos and b in z2i:
        p_, j_ = c2pos[a], z2i[b]
        print(f"  {ref:40s} | raw={D_tt[p_,j_]:5.1f}m  adj={D_tt_adj[p_,j_]:5.1f}m  road_km={D_rkm[p_,j_]:6.1f}  hops={hops[p_,j_]}")

# Save 138 cities x 2,484 mainland zones travel-time table for reuse (Check 1, Check 4, Check 6)
coords = M.reindex(main_zones)[['x', 'y']].values
rows_tt = []
for ci_pos, c in enumerate(cities_main):
    cx, cy = coords[z2i[c]]
    eucl = np.hypot(coords[:, 0] - cx, coords[:, 1] - cy) / 1000.0
    rows_tt.append(pd.DataFrame({
        'c': c,
        'z': main_zones,
        'eucl_km': eucl.astype(np.float32),
        'tt_chained': D_tt[ci_pos].astype(np.float32),
        'tt_adj': D_tt_adj[ci_pos].astype(np.float32),
        'road_km': D_rkm[ci_pos].astype(np.float32),
        'hops': hops[ci_pos],
    }))
tt_long = pd.concat(rows_tt, ignore_index=True)
tt_long = tt_long[tt_long.c != tt_long.z].copy()

# Also compute empirical circuity + speed approximation fitted on road_tt_150km (log-log OLS on 50-150km)
fit_df = tt[tt.eucl_km >= 30.0].copy()
# Fit tt_min = a * eucl_km^b (power law accounts for faster average speed at longer motorway distances)
b_pow, log_a = np.polyfit(np.log(fit_df.eucl_km.values), np.log(fit_df.tt_min.values), 1)
tt_long['tt_circuity_approx'] = np.where(
    tt_long.hops == 1,
    tt_long.tt_chained,
    np.exp(log_a) * (tt_long.eucl_km.values ** b_pow)
).astype(np.float32)
print(f"\nEmpirical power-law fit on 30-150km road_tt: tt_min = {np.exp(log_a):.3f} * eucl_km^{b_pow:.3f}")
print(f"Correlation(tt_adj, eucl_km) = {tt_long.tt_adj.corr(tt_long.eucl_km):.4f}; "
      f"Correlation(tt_adj, tt_chained) = {tt_long.tt_adj.corr(tt_long.tt_chained):.4f}")
tt_long.to_parquet('explore/out_p84_city_zone_tt.parquet', index=False)
print(f"Saved explore/out_p84_city_zone_tt.parquet ({len(tt_long):,} rows) in {time.time()-t0:.1f}s")

# 3. Re-run City Roots PPML model replacing straight-line distance bins with ROAD TRAVEL TIME bins!
con = duckdb.connect()
q = "'data/derived/pernoct_mun_2024_6m.parquet'"
cl = ",".join(f"'{c}'" for c in cities_main)
slices = {
    'Feb weekday': "date BETWEEN '2024-02-05' AND '2024-02-29' AND dayofweek(date) IN (2,3,4)",
    'ordinary Sat': "date IN ('2024-02-10','2024-02-17','2024-02-24','2024-10-19','2024-11-16','2024-11-23')",
    'Good Friday': "date='2024-03-29'",
    'mid-Aug': "date BETWEEN '2024-08-10' AND '2024-08-18'",
    'All Saints': "date='2024-11-01'",
    'Christmas Eve': "date='2024-12-24'",
}

# (A) City x Destination-Province level (mainland 138 cities x 47 mainland provinces, population-weighted road TT)
tt_long['p'] = tt_long.z.str[:2]
tt_long['zpop'] = M.loc[tt_long.z, 'population'].fillna(0).values + 1.0
prov_tt = tt_long.groupby(['c', 'p']).apply(
    lambda g_: pd.Series({
        'eucl_km': np.average(g_.eucl_km, weights=g_.zpop),
        'tt_chained': np.average(g_.tt_chained, weights=g_.zpop),
        'tt_adj': np.average(g_.tt_adj, weights=g_.zpop),
        'road_km': np.average(g_.road_km, weights=g_.zpop),
        'tt_approx': np.average(g_.tt_circuity_approx, weights=g_.zpop),
    }), include_groups=False
).reset_index()
prov_tt['cprov'] = prov_tt.c.str[:2]
prov_tt = prov_tt[prov_tt.cprov != prov_tt.p].copy()

Bm = B24.set_index(['ine_code', 'prov_birth']).persons
prov_tt['lborn'] = np.log1p(Bm.reindex(list(zip(prov_tt.c, prov_tt.p))).fillna(0).values)
prov_tt['pairFE'] = prov_tt.cprov + '_' + prov_tt.p

# Fine bins (~30-40 bins) for each distance/time metric
prov_tt['dbin_eucl'] = np.digitize(prov_tt.eucl_km, np.r_[0:100:10, 100:400:20, 400:1200:50])
prov_tt['ttbin_adj'] = np.digitize(prov_tt.tt_adj, np.r_[0:120:10, 120:360:15, 360:900:30])
prov_tt['ttbin_raw'] = np.digitize(prov_tt.tt_chained, np.r_[0:120:10, 120:360:15, 360:1000:30])
prov_tt['rkmbin'] = np.digitize(prov_tt.road_km, np.r_[0:120:15, 120:500:25, 500:1500:50])
prov_tt['ttbin_approx'] = np.digitize(prov_tt.tt_approx, np.r_[0:120:10, 120:360:15, 360:900:30])

print("\n=== PART A: City x Province PPML (138 mainland cities x 47 mainland provinces, other-province flows) ===")
rows_a = []
for lab, cond in slices.items():
    f = con.execute(f"""
        SELECT r c, substr(s,1,2) p, sum(p)/count(DISTINCT date) w
        FROM {q} WHERE ({cond}) AND r IN ({cl}) AND r<>s
        GROUP BY 1,2
    """).df().set_index(['c', 'p']).w
    g_ = prov_tt.copy()
    g_['w'] = f.reindex(list(zip(g_.c, g_.p))).fillna(0).values
    res_s = {'slice': lab}
    for key, fml in [
        ('eucl_bin', "w ~ lborn | c + p + dbin_eucl"),
        ('tt_adj_bin', "w ~ lborn | c + p + ttbin_adj"),
        ('tt_raw_bin', "w ~ lborn | c + p + ttbin_raw"),
        ('road_km_bin', "w ~ lborn | c + p + rkmbin"),
        ('tt_approx', "w ~ lborn | c + p + ttbin_approx"),
        ('pairFE', "w ~ lborn | c + pairFE"),
        ('pairFE_plus_tt', "w ~ lborn | c + pairFE + ttbin_adj"),
    ]:
        t = pf.fepois(fml, data=g_, vcov={'CRV1': 'c'}).tidy()
        res_s[key] = t.loc['lborn', 'Estimate']
        res_s[key + '_se'] = t.loc['lborn', 'Std. Error']
    rows_a.append(res_s)
df_a = pd.DataFrame(rows_a)
df_a.to_csv('explore/out_p84_prov_tt_roots.csv', index=False)
pd.set_option('display.width', 240)
print(df_a[['slice','eucl_bin','eucl_bin_se','tt_adj_bin','tt_adj_bin_se','tt_raw_bin','road_km_bin','tt_approx','pairFE','pairFE_plus_tt']].round(3).to_string(index=False))

# (B) City x ALL 2,484 Mainland Destination Zones (342,654 pairs! Other-province pairs only)
print("\n=== PART B: City x ALL 2,484 Mainland Destination Zones PPML (138 x 2,484, other-province pairs) ===")
cz = tt_long[['c', 'z', 'p', 'eucl_km', 'tt_adj', 'tt_chained', 'road_km']].copy()
cz['cprov'] = cz.c.str[:2]
cz = cz[cz.cprov != cz.p].copy()
cz['lborn'] = np.log1p(Bm.reindex(list(zip(cz.c, cz.p))).fillna(0).values)
cz['pairFE'] = cz.cprov + '_' + cz.p
cz['dbin_eucl'] = np.digitize(cz.eucl_km, np.r_[0:100:10, 100:400:20, 400:1200:50])
cz['ttbin_adj'] = np.digitize(cz.tt_adj, np.r_[0:120:10, 120:360:15, 360:1000:30])
cz['rkmbin'] = np.digitize(cz.road_km, np.r_[0:120:15, 120:500:25, 500:1500:50])

rows_b = []
for lab in ['Feb weekday', 'mid-Aug', 'All Saints', 'Christmas Eve']:
    cond = slices[lab]
    f = con.execute(f"""
        SELECT r c, s z, sum(p)/count(DISTINCT date) w
        FROM {q} WHERE ({cond}) AND r IN ({cl}) AND r<>s
        GROUP BY 1,2
    """).df().set_index(['c', 'z']).w
    g_ = cz.copy()
    g_['w'] = f.reindex(list(zip(g_.c, g_.z))).fillna(0).values
    res_s = {'slice': lab, 'n_pairs': len(g_), 'nz_pairs': int((g_.w > 0).sum())}
    for key, fml in [
        ('z_eucl', "w ~ lborn | c + z + dbin_eucl"),
        ('z_tt_adj', "w ~ lborn | c + z + ttbin_adj"),
        ('z_rkm', "w ~ lborn | c + z + rkmbin"),
        ('z_pairFE_tt', "w ~ lborn | c + z + pairFE + ttbin_adj"),
    ]:
        t = pf.fepois(fml, data=g_, vcov={'CRV1': 'c'}).tidy()
        res_s[key] = t.loc['lborn', 'Estimate']
        res_s[key + '_se'] = t.loc['lborn', 'Std. Error']
    rows_b.append(res_s)
    print(f"  {lab:15s}: z_eucl={res_s['z_eucl']:.3f} ({res_s['z_eucl_se']:.3f}) | "
          f"z_tt_adj={res_s['z_tt_adj']:.3f} ({res_s['z_tt_adj_se']:.3f}) | "
          f"z_rkm={res_s['z_rkm']:.3f} | z_pairFE+tt={res_s['z_pairFE_tt']:.3f} ({res_s['z_pairFE_tt_se']:.3f})", flush=True)

df_b = pd.DataFrame(rows_b)
df_b.to_csv('explore/out_p84_zone_tt_roots.csv', index=False)
print(f"All Check 4 done in {time.time()-t0:.1f}s")
