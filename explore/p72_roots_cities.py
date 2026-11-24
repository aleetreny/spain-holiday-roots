# Roots test with municipal birthplace (Census 2021-2025, 151 large cities x province of birth).
# Do residents of a city sleep (on holidays) in the provinces where its residents were born, beyond distance -- and comparing NEIGHBOURING cities?
import pandas as pd, numpy as np, duckdb, sys, pyfixest as pf
sys.path.insert(0,'explore'); from common import nodes
B = pd.read_csv('_local/external/ine_birth_residence/mun_residence_x_prov_birth_censo_2021_2025.csv', dtype={'ine_code':str,'prov_birth':str})
B = B[(B.year==2024) & B.prov_birth.str.match(r'^\d\d$')]
N = nodes().set_index('id')
cities = sorted(set(B.ine_code) & set(N.index))   # big cities are single-municipality MITMA zones
print('cities matched to MITMA zones:', len(cities), 'of', B.ine_code.nunique())
con = duckdb.connect(); q = "'data/derived/pernoct_mun_2024_6m.parquet'"
cl = ",".join(f"'{c}'" for c in cities)
# population-weighted province centroids
M = N[~N.foreign].copy(); M['w'] = M.population.fillna(0)+1
pc = M.groupby('prov').apply(lambda g: pd.Series({'x':np.average(g.x,weights=g.w),'y':np.average(g.y,weights=g.w)}))
slices = {'Feb weekday': "date BETWEEN '2024-02-05' AND '2024-02-29' AND dayofweek(date) IN (2,3,4)",
          'ordinary Sat': "date IN ('2024-02-10','2024-02-17','2024-02-24','2024-10-19','2024-11-16','2024-11-23')",
          'Good Friday': "date='2024-03-29'", 'mid-Aug': "date BETWEEN '2024-08-10' AND '2024-08-18'",
          'All Saints': "date='2024-11-01'", 'Christmas Eve': "date='2024-12-24'"}
grid = pd.MultiIndex.from_product([cities, sorted(pc.index)], names=['c','p']).to_frame(index=False)
grid = grid[grid.p.isin([f'{i:02d}' for i in range(1,51)])]
grid['cprov'] = grid.c.str[:2]
grid['same_prov'] = (grid.cprov==grid.p).astype(int)
grid['km'] = np.hypot(N.loc[grid.c,'x'].values-pc.loc[grid.p,'x'].values, N.loc[grid.c,'y'].values-pc.loc[grid.p,'y'].values)/1000
grid['dbin'] = np.digitize(grid.km, np.r_[0:100:10, 100:400:20, 400:1200:50])
Bm = B.set_index(['ine_code','prov_birth']).persons
grid['born'] = Bm.reindex(list(zip(grid.c, grid.p))).fillna(0).values
grid['lborn'] = np.log1p(grid.born)
grid['pairFE'] = grid.cprov + '_' + grid.p
rows=[]
for lab, cond in slices.items():
    f = con.execute(f"SELECT r c, substr(s,1,2) p, sum(p)/count(DISTINCT date) w FROM {q} WHERE {cond} AND r IN ({cl}) AND r<>s GROUP BY 1,2").df().set_index(['c','p']).w
    g = grid.copy(); g['w'] = f.reindex(list(zip(g.c,g.p))).fillna(0).values
    g = g[g.same_prov==0]   # other provinces only: where do people go when they LEAVE their province?
    m1 = pf.fepois("w ~ lborn | c + p + dbin", data=g, vcov={'CRV1':'c'})
    m2 = pf.fepois("w ~ lborn | c + pairFE", data=g, vcov={'CRV1':'c'})   # neighbouring cities: same origin province x destination province
    t1, t2 = m1.tidy(), m2.tidy()
    rows.append(dict(slice=lab, total_M=g.w.sum()/1e6, elast_distFE=t1.loc['lborn','Estimate'], se1=t1.loc['lborn','Std. Error'],
                     elast_within_province_pair=t2.loc['lborn','Estimate'], se2=t2.loc['lborn','Std. Error']))
out = pd.DataFrame(rows); out.to_csv('explore/out_roots_cities.csv', index=False)
pd.set_option('display.width',200); print(out.round(3).to_string())
# illustration: Barcelona-area cities -> Andalusian provinces, August, vs share born in Andalusia
AND = ['04','11','14','18','21','23','29','41']
f = con.execute(f"SELECT r c, substr(s,1,2) p, sum(p)/count(DISTINCT date) w FROM {q} WHERE date BETWEEN '2024-08-10' AND '2024-08-18' AND r IN ({cl}) AND r<>s GROUP BY 1,2").df()
f = f[f.p!=f.c.str[:2]]
tot = f.groupby('c').w.sum(); toAND = f[f.p.isin(AND)].groupby('c').w.sum()
bornAND = B[B.prov_birth.isin(AND)].groupby('ine_code').persons.sum()/B.groupby('ine_code').persons.sum()
ill = pd.DataFrame({'pct_Aug_nights_out_of_prov_in_Andalusia':100*toAND/tot, 'pct_residents_born_in_Andalusia':100*bornAND}).dropna()
nm = B.drop_duplicates('ine_code').set_index('ine_code').municipality_name
bcn = ill[ill.index.str[:2]=='08'].copy(); bcn['city'] = nm.reindex(bcn.index).values
print(bcn.sort_values('pct_residents_born_in_Andalusia', ascending=False).round(1).to_string())
print('corr within Barcelona province:', round(bcn.iloc[:,0].corr(bcn.iloc[:,1]),3))
madr = ill[ill.index.str[:2]=='28'].copy(); madr['city']=nm.reindex(madr.index).values
print(madr.sort_values('pct_residents_born_in_Andalusia', ascending=False).round(1).to_string())
