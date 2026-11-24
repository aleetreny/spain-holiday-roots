# Check 6 (Part 1, baseline before ADRH income): Tourism confound check
# Restrict city-level roots model to (a) inland destination zones (true maritime coast removed, AND p73 flag),
# (b) strictly inland destination provinces (excluding all 19 coastal provinces), and (c) All Saints' Day.
import sys
sys.path.insert(0, 'explore')
import duckdb, geopandas as gpd, numpy as np, pandas as pd, pyfixest as pf
from common import nodes

# True mainland coastal provinces (19 provinces with sea coastline)
COASTAL_PROVS = {'03','04','08','11','12','15','17','18','20','21','27','29','30','33','36','39','43','46','48'}

# Build both p73 coastal flag AND true maritime coastal zone flag
g = gpd.read_file('clean_data/v2/zones/municipios_mitma.gpkg')
land = g.geometry.union_all()
outer_p73 = land.boundary
# Exterior ring of Iberian+FR+PT mainland polygon (drops internal polygon slivers/holes)
from shapely.geometry import MultiPolygon
main_poly = max(land.geoms, key=lambda p: p.area) if isinstance(land, MultiPolygon) else land
main_ext = main_poly.exterior

g['prov'] = g.id.str[:2]
g['coastal_p73'] = (g.geometry.buffer(50).intersects(outer_p73) & g.id.str.match(r'^\d')).astype(int)
# True maritime coast: touches exterior of Iberian+PT+FR land mass AND is in a maritime province (drops PT/FR land border & inland slivers)
# Also include coastal zones whose polygon is an island/peninsula piece in a maritime province (e.g. Cadiz bay)
spain_main = g[g.id.str.match(r'^\d') & ~g.prov.isin(['07','35','38','51','52'])].copy()
# Check exterior of Spain+PT+FR vs Spain-only? Actually if zone touches main_ext and prov in COASTAL_PROVS:
# Check also small coastal peninsulas/islets in COASTAL_PROVS disconnected from main_poly
g['coastal_true'] = (
    g.id.str.match(r'^\d')
    & g.prov.isin(COASTAL_PROVS)
    & (g.geometry.buffer(50).intersects(main_ext) | ~g.geometry.intersects(main_poly))
).astype(int)
# Also mark 10km coastal buffer (coastal strip / beach hinterland like Benidorm-adjacent or Costa del Sol second row)
sea_edge = main_ext  # wait, restrict to sea portion? We can buffer true coastal zones by 10km within coastal provinces!
coastal_union = g.loc[g.coastal_true == 1, 'geometry'].union_all()
g['coastal_10km'] = (
    g.id.str.match(r'^\d')
    & g.prov.isin(COASTAL_PROVS)
    & g.geometry.intersects(coastal_union.buffer(10000))
).astype(int)

cf = g.set_index('id')[['coastal_p73', 'coastal_true', 'coastal_10km']]
cf.to_csv('explore/out_p86_coastal_flags.csv')
print("Coastal zone counts (mainland Spanish zones):")
m_ids = spain_main.id.values
print(cf.loc[m_ids].sum().to_dict(), "out of", len(m_ids))

# 1. Re-check p73 village-refill gradient with TRUE coastal flag vs p73 flag!
H = pd.read_csv('_local/external/ine_historical_pop/municipal_population_1900_2025.csv',
                dtype={'ine_code': str, 'province_code': str})
rel = pd.read_csv('raw_data_cache/v2/zonificacion/relacion_ine_zonificacionMitma.csv', sep='|', dtype=str)[
    ['municipio_ine', 'municipio_mitma']
].drop_duplicates()
H = H.merge(rel, left_on='ine_code', right_on='municipio_ine', how='inner')
cols = ['pop_1900','pop_1910','pop_1920','pop_1930','pop_1940','pop_1950','pop_1960','pop_1970','pop_1981','pop_2001','pop_2021']
Z = H.groupby('municipio_mitma')[cols].sum(min_count=1)
Z['peak'] = Z[['pop_1900','pop_1910','pop_1920','pop_1930','pop_1940','pop_1950','pop_1960']].max(axis=1)
Z['loss_peak_2021'] = 1 - Z.pop_2021 / Z.peak

R = pd.read_csv('explore/out_refill.csv', dtype={'z': str, 'prov': str}).set_index('z')
D = R.join(Z, how='inner').join(cf, how='inner')
D['prov'] = D.index.str[:2]
D = D[(D.resid >= 200) & ~D.prov.isin(['07','35','38','51','52'])].copy()
D['lmult_aug'] = np.log(D.mult_aug)
D['lmult_nov1'] = np.log(D.mult_nov1)
D['lres'] = np.log(D.resid)
D['loss_rank'] = D.loss_peak_2021.rank(pct=True)

print("\n=== Village Refill Gradient: p73 flag vs True Maritime Coast vs 10km Coastal Buffer vs Inland Provinces Only ===")
for label, mask in [
    ('p73 inland (original)', D.coastal_p73 == 0),
    ('TRUE inland (no sea coast)', D.coastal_true == 0),
    ('STRICT inland (>10km from coast)', D.coastal_10km == 0),
    ('INLAND PROVINCES ONLY (28 provs)', ~D.prov.isin(COASTAL_PROVS)),
]:
    sub = D[mask].copy()
    sub['loss_q'] = pd.qcut(sub.loss_peak_2021, 5, labels=['Q1','Q2','Q3','Q4','Q5'])
    qtab = sub.groupby('loss_q', observed=True).agg(
        med_loss=('loss_peak_2021','median'), aug=('mult_aug','median'), nov1=('mult_nov1','median'), feb=('mult_feb','median')
    )
    m_aug = pf.feols("lmult_aug ~ loss_rank + lres | prov", data=sub, vcov={'CRV1':'prov'}).tidy()
    m_nov = pf.feols("lmult_nov1 ~ loss_rank + lres | prov", data=sub, vcov={'CRV1':'prov'}).tidy()
    print(f"\n  [{label}] n={len(sub)}: Q1->Q5 Aug {qtab.loc['Q1','aug']:.2f} -> {qtab.loc['Q5','aug']:.2f} "
          f"(Nov1 {qtab.loc['Q1','nov1']:.2f} -> {qtab.loc['Q5','nov1']:.2f}, Feb {qtab.loc['Q1','feb']:.2f} -> {qtab.loc['Q5','feb']:.2f})")
    print(f"    OLS loss_rank (prov FE + lres): Aug={m_aug.loc['loss_rank','Estimate']:+.3f} (se {m_aug.loc['loss_rank','Std. Error']:.3f})  "
          f"Nov1={m_nov.loc['loss_rank','Estimate']:+.3f} (se {m_nov.loc['loss_rank','Std. Error']:.3f})")

# 2. City roots elasticity (p72) restricted to inland destinations!
B = pd.read_csv('_local/external/ine_birth_residence/mun_residence_x_prov_birth_censo_2021_2025.csv',
                dtype={'ine_code': str, 'prov_birth': str})
B = B[(B.year == 2024) & B.prov_birth.str.match(r'^\d\d$')]
N = nodes().set_index('id')
cities = sorted(set(B.ine_code) & set(N.index))
M = N[(~N.foreign) & (~N.island) & (~N.prov.isin(['51','52']))].copy()
M['w'] = M.population.fillna(0) + 1
pc = M.groupby('prov').apply(lambda g_: pd.Series({'x': np.average(g_.x, weights=g_.w), 'y': np.average(g_.y, weights=g_.w)}))

con = duckdb.connect()
q = "'data/derived/pernoct_mun_2024_6m.parquet'"
cl = ",".join(f"'{c}'" for c in cities)
# Register coastal flags in duckdb
con.register('cf_tbl', cf.reset_index().rename(columns={'id': 'z'}))

slices = {
    'Feb weekday': "date BETWEEN '2024-02-05' AND '2024-02-29' AND dayofweek(date) IN (2,3,4)",
    'ordinary Sat': "date IN ('2024-02-10','2024-02-17','2024-02-24','2024-10-19','2024-11-16','2024-11-23')",
    'Good Friday': "date='2024-03-29'",
    'mid-Aug': "date BETWEEN '2024-08-10' AND '2024-08-18'",
    'All Saints': "date='2024-11-01'",
    'Christmas Eve': "date='2024-12-24'",
}

grid = pd.MultiIndex.from_product([cities, sorted(pc.index)], names=['c', 'p']).to_frame(index=False)
grid['cprov'] = grid.c.str[:2]
grid = grid[grid.cprov != grid.p].copy()
grid['km'] = np.hypot(N.loc[grid.c, 'x'].values - pc.loc[grid.p, 'x'].values,
                      N.loc[grid.c, 'y'].values - pc.loc[grid.p, 'y'].values) / 1000
grid['dbin'] = np.digitize(grid.km, np.r_[0:100:10, 100:400:20, 400:1200:50])
Bm = B.set_index(['ine_code', 'prov_birth']).persons
grid['lborn'] = np.log1p(Bm.reindex(list(zip(grid.c, grid.p))).fillna(0).values)
grid['pairFE'] = grid.cprov + '_' + grid.p

print("\n=== City Roots Elasticity (151 cities x mainland dest provinces, neighbouring-city pairFE): All vs Inland-Zone Flows vs Inland-Province Flows ===")
rows = []
for mode_name, where_s, prov_filter in [
    ('All mainland zones (p72 ref)', "1=1", None),
    ('Inland zones only (p73 flag=0)', "cf_tbl.coastal_p73 = 0", None),
    ('Inland zones only (true sea coast=0)', "cf_tbl.coastal_true = 0", None),
    ('Strict inland zones (>10km from coast)', "cf_tbl.coastal_10km = 0", None),
    ('Inland provinces only (28 non-coastal provs)', "1=1", lambda df: ~df.p.isin(COASTAL_PROVS)),
    ('Coastal zones ONLY (placebo/contrast)', "cf_tbl.coastal_true = 1", None),
]:
    for lab, cond in slices.items():
        f = con.execute(f"""
            SELECT x.r c, substr(x.s,1,2) p, sum(x.p)/count(DISTINCT x.date) w
            FROM {q} x
            JOIN cf_tbl ON x.s = cf_tbl.z
            WHERE ({cond}) AND x.r IN ({cl}) AND x.r <> x.s AND ({where_s})
            GROUP BY 1,2
        """).df().set_index(['c','p']).w
        g_ = grid.copy()
        if prov_filter is not None:
            g_ = g_[prov_filter(g_)].copy()
        g_['w'] = f.reindex(list(zip(g_.c, g_.p))).fillna(0).values
        m1 = pf.fepois("w ~ lborn | c + p + dbin", data=g_, vcov={'CRV1': 'c'}).tidy()
        m2 = pf.fepois("w ~ lborn | c + pairFE", data=g_, vcov={'CRV1': 'c'}).tidy()
        rows.append(dict(
            sample=mode_name, slice=lab, flow_M=g_.w.sum()/1e6,
            elast_dist=m1.loc['lborn','Estimate'], se_dist=m1.loc['lborn','Std. Error'],
            elast_pair=m2.loc['lborn','Estimate'], se_pair=m2.loc['lborn','Std. Error'],
        ))

res = pd.DataFrame(rows)
res.to_csv('explore/out_p86_tourism_confound.csv', index=False)
pd.set_option('display.width', 220)
print(res.round(3).to_string(index=False))
