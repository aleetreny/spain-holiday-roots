# Village refill: present population / resident population; where do visitors come from; distance of holiday-excess edges
import duckdb, pandas as pd, numpy as np, sys
sys.path.insert(0,'explore'); from common import nodes
con = duckdb.connect(); q="'data/derived/pernoct_mun_2024_6m.parquet'"
N = nodes().set_index('id')
nm = pd.read_csv('raw_data_cache/v2/zonificacion/zonificacion_municipios/nombres_municipios.csv', sep='|', dtype=str).set_index('ID')['name']
ref = "date BETWEEN '2024-02-01' AND '2024-02-29' AND dayofweek(date) IN (2,3,4)"
def present(cond):
    return con.execute(f"SELECT s z, sum(p)/count(DISTINCT date) present FROM {q} WHERE {cond} GROUP BY 1").df().set_index('z').present
def resid(cond):
    return con.execute(f"SELECT r z, sum(p)/count(DISTINCT date) resid FROM {q} WHERE {cond} GROUP BY 1").df().set_index('z').resid
R = resid(ref); P0 = present(ref)
out = pd.DataFrame({'resid':R,'present_feb':P0})
for lab, cond in [('aug15',"date='2024-08-15'"),('aug_mid',"date BETWEEN '2024-08-10' AND '2024-08-18'"),('nov1',"date='2024-11-01'"),('dec24',"date='2024-12-24'")]:
    out['present_'+lab] = present(cond)
out = out[out.index.str.match(r'^\d')]
out['mult_aug'] = out.present_aug_mid/out.resid; out['mult_nov1']=out.present_nov1/out.resid; out['mult_feb']=out.present_feb/out.resid
out['name'] = nm.reindex(out.index.str.replace('_AM','')).values
out['prov'] = out.index.str[:2]
out.to_csv('explore/out_refill.csv')
print('zones', len(out))
print('multiplier mid-Aug quantiles', out.mult_aug.quantile([.1,.25,.5,.75,.9,.99]).round(2).to_dict())
print('Nov1 quantiles', out.mult_nov1.quantile([.1,.5,.9,.99]).round(2).to_dict())
# by resident size class
out['size'] = pd.cut(out.resid, [0,1000,5000,20000,100000,1e8], labels=['<1k','1-5k','5-20k','20-100k','>100k'])
print(out.groupby('size', observed=True)[['mult_aug','mult_nov1','mult_feb']].median().round(3))
print(out.sort_values('mult_aug', ascending=False).head(15)[['name','prov','resid','mult_aug','mult_nov1']].round(2).to_string())
# national: share of nights in <5k zones
for c in ['present_feb','present_aug_mid','present_nov1']:
    print(c, 'share of population sleeping in zones with <5k residents:', round(100*out.loc[out.resid<5000, c].sum()/out[c].sum(),2), '%')
