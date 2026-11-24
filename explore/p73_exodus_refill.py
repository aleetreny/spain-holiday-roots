# Do the places that emptied most in the rural exodus refill most on holidays (beyond size, coast, province)?
import pandas as pd, numpy as np, geopandas as gpd, sys, pyfixest as pf
sys.path.insert(0,'explore'); from common import nodes
H = pd.read_csv('_local/external/ine_historical_pop/municipal_population_1900_2025.csv', dtype={'ine_code':str,'province_code':str})
rel = pd.read_csv('raw_data_cache/v2/zonificacion/relacion_ine_zonificacionMitma.csv', sep='|', dtype=str)[['municipio_ine','municipio_mitma']].drop_duplicates()
H = H.merge(rel, left_on='ine_code', right_on='municipio_ine', how='inner')
cols = ['pop_1900','pop_1910','pop_1920','pop_1930','pop_1940','pop_1950','pop_1960','pop_1970','pop_1981','pop_2001','pop_2021']
Z = H.groupby('municipio_mitma')[cols].sum(min_count=1)
Z['peak'] = Z[['pop_1900','pop_1910','pop_1920','pop_1930','pop_1940','pop_1950','pop_1960']].max(axis=1)
Z['loss_peak_2021'] = 1 - Z.pop_2021/Z.peak            # share of peak population lost
Z['l_1950_1981'] = np.log(Z.pop_1981/Z.pop_1950)       # exodus-era change
R = pd.read_csv('explore/out_refill.csv', dtype={'z':str}).set_index('z')
D = R.join(Z, how='inner')
# coastal flag: zone touches the outer boundary of Spain+PT+FR land union
g = gpd.read_file('clean_data/v2/zones/municipios_mitma.gpkg')
land = g.geometry.union_all()
outer = land.boundary
g['coastal'] = g.geometry.buffer(50).intersects(outer) & g.id.str.match(r'^\d')
D = D.join(g.set_index('id')['coastal'])
D = D[(D.resid>=200) & D.coastal.notna() & ~D.prov.isin(['07','35','38','51','52'])].copy()
D['coastal'] = D.coastal.astype(int)
D['lmult_aug'] = np.log(D.mult_aug); D['lmult_nov1'] = np.log(D.mult_nov1); D['lres'] = np.log(D.resid)
print('zones', len(D), 'coastal', D.coastal.sum(), '; median loss from peak', round(D.loss_peak_2021.median(),2))
for y in ['lmult_aug','lmult_nov1']:
    m = pf.feols(f"{y} ~ loss_peak_2021 + lres + coastal | prov", data=D, vcov={'CRV1':'prov'}); t = m.tidy()
    print(f"{y}: loss_from_peak {t.loc['loss_peak_2021','Estimate']:+.3f} (se {t.loc['loss_peak_2021','Std. Error']:.3f})  coastal {t.loc['coastal','Estimate']:+.3f}  lres {t.loc['lres','Estimate']:+.3f}  R2 {m._r2:.2f}")
inl = D[D.coastal==0].copy()
inl['loss_q'] = pd.qcut(inl.loss_peak_2021, 5, labels=['Q1 grew','Q2','Q3','Q4','Q5 emptied most'])
print(inl.groupby('loss_q', observed=True).agg(n=('resid','size'), median_loss=('loss_peak_2021','median'), mult_aug=('mult_aug','median'), mult_nov1=('mult_nov1','median'), mult_feb=('mult_feb','median')).round(2))
D['loss_c'] = D.loss_peak_2021.clip(-1, 1)
D['loss_rank'] = D.loss_peak_2021.rank(pct=True)
inl = D[D.coastal==0]
for y in ['lmult_aug','lmult_nov1']:
    for x in ['loss_c','loss_rank']:
        m = pf.feols(f"{y} ~ {x} + lres | prov", data=inl, vcov={'CRV1':'prov'}); t = m.tidy()
        print(f"INLAND {y} ~ {x}: {t.loc[x,'Estimate']:+.3f} (se {t.loc[x,'Std. Error']:.3f}), lres {t.loc['lres','Estimate']:+.3f}")
# size-matched: within resident-size terciles
inl = inl.copy(); inl['size_t'] = pd.qcut(inl.resid, 3, labels=['small','mid','large'])
inl['loss_q'] = pd.qcut(inl.loss_peak_2021, 5, labels=False)
print(inl.pivot_table(index='loss_q', columns='size_t', values='mult_aug', aggfunc='median', observed=True).round(2))
print(inl.pivot_table(index='loss_q', columns='size_t', values='mult_nov1', aggfunc='median', observed=True).round(2))
