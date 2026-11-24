# Do holiday overnight flows follow birthplace ("roots") ties beyond distance? Province-pair and zone-level PPML.
import pandas as pd, numpy as np, duckdb, sys, pyfixest as pf, unicodedata, re
sys.path.insert(0,'explore'); from common import nodes
d = pd.read_csv('_local/external/ine_birth_prov/0ccaa001.csv', sep=';', encoding='utf-8-sig', dtype=str)
d = d[d.Sexo=='Ambos sexos'].copy()
d['n'] = d.Total.str.replace(',','').astype(float)
# residence province code
uni = {'03':'33','04':'07','06':'39','13':'28','14':'30','15':'31','17':'26','18':'51','19':'52'}
def rescode(row):
    if isinstance(row.Provincias,str): return row.Provincias[:2]
    c = row['Comunidades y Ciudades Autónomas']
    return uni.get(c[:2]) if isinstance(c,str) else None
d['res'] = d.apply(rescode, axis=1)
# birth province code from names
def norm(s): return unicodedata.normalize('NFKD', s).encode('ascii','ignore').decode().lower()
provnames = {}
for p in d.Provincias.dropna().unique(): provnames[norm(p[3:])] = p[:2]
extra = {'asturias, principado de':'33','balears, illes':'07','cantabria':'39','madrid, comunidad de':'28','murcia, region de':'30',
         'navarra, comunidad foral de':'31','rioja, la':'26','ceuta':'51','melilla':'52'}
provnames.update(extra)
d['birth'] = d['Lugar de nacimiento'].map(lambda s: provnames.get(norm(s)))
R = d.dropna(subset=['res','birth'])
# keep province-level rows only (residence rows of multi-province CCAA have Provincias filled; CCAA total rows have NaN and are dropped unless uniprovincial)
R = R[(R.Provincias.notna()) | (R['Comunidades y Ciudades Autónomas'].str[:2].isin(uni.keys()))]
R = R.groupby(['birth','res']).n.sum().reset_index()
R.to_csv('_local/external/ine_birth_prov/birth_x_residence_2021_long.csv', index=False)
print('roots matrix', R.shape, 'total', R.n.sum()/1e6, 'M;  example born Jaen living Barcelona:', R[(R.birth=='23')&(R.res=='08')].n.values)
Rm = R.set_index(['birth','res']).n

N = nodes().set_index('id'); N = N[(~N.foreign)&(~N.island)&(~N.prov.isin(['51','52']))]
con = duckdb.connect(); q = "'data/derived/pernoct_mun_2024_6m.parquet'"
mainland = tuple(N.index)
def prov_flows(cond):
    df = con.execute(f"SELECT substr(r,1,2) rp, substr(s,1,2) sp, sum(p)/count(DISTINCT date) w FROM {q} WHERE {cond} GROUP BY 1,2").df()
    return df.set_index(['rp','sp']).w
slices = {
 'ordinary weekday (Feb Tue-Thu)': "date BETWEEN '2024-02-05' AND '2024-02-29' AND dayofweek(date) IN (2,3,4)",
 'ordinary Saturday (Feb/Oct/Nov)': "dayofweek(date)=6 AND date IN ('2024-02-10','2024-02-17','2024-02-24','2024-10-19','2024-11-16','2024-11-23')",
 'Good Friday': "date='2024-03-29'",
 'mid-August (10-18 Aug)': "date BETWEEN '2024-08-10' AND '2024-08-18'",
 'All Saints (1 Nov)': "date='2024-11-01'",
 'Christmas Eve (24 Dec)': "date='2024-12-24'",
}
# population-weighted province centroids for distance
N['w'] = N.population.fillna(0)
pc = N.groupby('prov').apply(lambda g: pd.Series({'x':np.average(g.x,weights=g.w+1),'y':np.average(g.y,weights=g.w+1),'pop':g.w.sum()}))
provs = pc.index.tolist()
rows=[]
pairs = pd.DataFrame([(a,b) for a in provs for b in provs if a!=b], columns=['rp','sp'])
pairs['km'] = np.hypot(pc.loc[pairs.rp,'x'].values-pc.loc[pairs.sp,'x'].values, pc.loc[pairs.rp,'y'].values-pc.loc[pairs.sp,'y'].values)/1000
pairs['roots'] = Rm.reindex(list(zip(pairs.sp, pairs.rp))).fillna(0).values   # born in DESTINATION, living in RESIDENCE province
pairs['roots_rev'] = Rm.reindex(list(zip(pairs.rp, pairs.sp))).fillna(0).values # born in residence prov, living in destination (family members there)
pairs['lroots'] = np.log1p(pairs.roots); pairs['lroots_rev'] = np.log1p(pairs.roots_rev)
pairs['lkm'] = np.log(pairs.km); pairs['kmbin'] = pd.qcut(pairs.km, 20, labels=False)
ref = prov_flows(slices['ordinary weekday (Feb Tue-Thu)'])
for lab, cond in slices.items():
    f = prov_flows(cond)
    p = pairs.copy(); p['w'] = f.reindex(list(zip(p.rp,p.sp))).fillna(0).values
    p['w_ref'] = ref.reindex(list(zip(p.rp,p.sp))).fillna(0).values
    p['excess'] = (p.w - p.w_ref).clip(lower=0)
    for dep in (['w'] if lab.startswith('ordinary weekday') else ['w','excess']):
        m = pf.fepois(f"{dep} ~ lroots + lroots_rev | rp + sp + kmbin", data=p, vcov='hetero')
        t = m.tidy()
        m0 = pf.fepois(f"{dep} ~ 1 | rp + sp + kmbin", data=p, vcov='hetero')
        # pseudo-R2 gain: deviance reduction share
        dev0, dev1 = float(np.sum(m0.deviance)), float(np.sum(m.deviance))
        rows.append(dict(slice=lab, dep=dep, b_roots=t.loc['lroots','Estimate'], se=t.loc['lroots','Std. Error'],
                         b_rev=t.loc['lroots_rev','Estimate'], dev_explained_by_roots=1-dev1/dev0, total=p[dep].sum()))
out = pd.DataFrame(rows); out.to_csv('explore/out_roots_prov.csv', index=False)
pd.set_option('display.width',200); print(out.round(3).to_string())
