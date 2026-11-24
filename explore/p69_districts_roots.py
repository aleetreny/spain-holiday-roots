# Within-metro contrast: same location, different roots. Share of August / All Saints away-nights spent in Andalusia+Extremadura+CLM+CyL+Galicia (exodus origins)
import duckdb, pandas as pd, numpy as np, sys
sys.path.insert(0,'explore'); from common import PROV2CCAA
con = duckdb.connect(); q="'data/derived/pernoct_2024_6m.parquet'"
dn = pd.read_csv('raw_data_cache/v2/zonificacion/zonificacion_distritos/nombres_distritos.csv', sep='|', dtype=str)
dn = dn.set_index(dn.columns[1])[dn.columns[2]] if dn.shape[1]>=3 else dn.set_index(dn.columns[0])[dn.columns[1]]
pop = pd.read_csv('raw_data_cache/v2/zonificacion/zonificacion_distritos/poblacion_distritos.csv', sep='|', dtype=str)
print(pop.head(2).to_string())
EXO = {'AND','EXT','CLM','CYL','GAL','MUR','ARA'}  # main origin regions of 1950-75 migrants to Catalonia
def shares(cond, prefix_list):
    pl = " OR ".join(f"r LIKE '{p}%'" for p in prefix_list)
    df = con.execute(f"""SELECT r, substr(s,1,2) sp, sum(p)/count(DISTINCT date) w FROM {q}
        WHERE ({cond}) AND ({pl}) AND substr(r,1,5)<>substr(s,1,5) GROUP BY 1,2""").df()
    df['ccaa'] = df.sp.map(PROV2CCAA)
    g = df.groupby('r')
    out = pd.DataFrame({'away': g.w.sum(),
                        'pct_exodus_regions': 100*df[df.ccaa.isin(EXO)].groupby('r').w.sum()/g.w.sum(),
                        'pct_AND': 100*df[df.ccaa=='AND'].groupby('r').w.sum()/g.w.sum(),
                        'pct_CAT': 100*df[df.ccaa=='CAT'].groupby('r').w.sum()/g.w.sum()}).fillna(0)
    return out
prefixes = ['08019','08101','08073','08015','08245','08279','08187','08205','08217','08113','08200','08121']  # BCN, Hospitalet, Cornellà, Badalona, Sta Coloma, Terrassa, Sabadell, Sant Cugat, Sant Boi, Mataró, Sant Adrià?, Mollet?
a = shares("date BETWEEN '2024-08-10' AND '2024-08-18'", prefixes)
b = shares("date='2024-11-01'", prefixes)
c = shares("date BETWEEN '2024-02-05' AND '2024-02-29' AND dayofweek(date) IN (2,3,4)", prefixes)
t = a[['away','pct_exodus_regions','pct_AND','pct_CAT']].join(b[['pct_exodus_regions','pct_AND']], rsuffix='_nov1').join(c[['pct_AND']], rsuffix='_feb')
t['name'] = [dn.get(i, i) for i in t.index]
pd.set_option('display.width',220)
print(t.sort_values('pct_AND', ascending=False).round(1).to_string())
