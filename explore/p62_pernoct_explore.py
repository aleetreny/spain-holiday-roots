import duckdb, pandas as pd, numpy as np, sys
sys.path.insert(0,'explore'); from common import nodes
con = duckdb.connect(); q="'data/derived/pernoct_mun_2024_6m.parquet'"
N = nodes().set_index('id'); name = N['name'] if 'name' in N else None
nm = pd.read_csv('raw_data_cache/v2/zonificacion/zonificacion_municipios/nombres_municipios.csv', sep='|', dtype=str)
print(nm.head(2)); nm = nm.set_index(nm.columns[0])[nm.columns[1]] if nm.shape[1]>=2 else None
PN = {'01':'Álava','02':'Albacete','03':'Alicante','04':'Almería','05':'Ávila','06':'Badajoz','07':'Baleares','08':'Barcelona','09':'Burgos','10':'Cáceres','11':'Cádiz','12':'Castellón','13':'C.Real','14':'Córdoba','15':'A Coruña','16':'Cuenca','17':'Girona','18':'Granada','19':'Guadalajara','20':'Gipuzkoa','21':'Huelva','22':'Huesca','23':'Jaén','24':'León','25':'Lleida','26':'La Rioja','27':'Lugo','28':'Madrid','29':'Málaga','30':'Murcia','31':'Navarra','32':'Ourense','33':'Asturias','34':'Palencia','35':'Las Palmas','36':'Pontevedra','37':'Salamanca','38':'Tenerife','39':'Cantabria','40':'Segovia','41':'Sevilla','42':'Soria','43':'Tarragona','44':'Teruel','45':'Toledo','46':'Valencia','47':'Valladolid','48':'Bizkaia','49':'Zamora','50':'Zaragoza','51':'Ceuta','52':'Melilla'}
# reference = median over Feb Tue-Thu nights; holiday nights
ref_days = "date BETWEEN '2024-02-01' AND '2024-02-29' AND dayofweek(date) IN (2,3,4)"
def provmat(cond):
    return con.execute(f"""SELECT substr(r,1,2) rp, substr(s,1,2) sp, sum(p)/count(DISTINCT date) p FROM {q} WHERE {cond} GROUP BY 1,2""").df()
ref = provmat(ref_days).set_index(['rp','sp']).p
for lab, cond in [('Aug 15', "date='2024-08-15'"), ('Nov 1 (All Saints)', "date='2024-11-01'"), ('Dec 24', "date='2024-12-24'"), ('Good Fri', "date='2024-03-29'")]:
    h = provmat(cond).set_index(['rp','sp']).p
    ex = (h - ref.reindex(h.index).fillna(0)).rename('excess').reset_index()
    ex = ex[ex.rp!=ex.sp]
    print(f'\n=== {lab}: total excess other-province sleepers {ex.excess.sum()/1e6:.2f}M')
    for city in ['28','08','48','46','41','20']:
        e = ex[ex.rp==city].sort_values('excess', ascending=False).head(8)
        tot = ex[ex.rp==city].excess.sum()
        print(f"  residents of {PN[city]:10s} excess {tot/1e3:7.0f}k -> " + ', '.join(f"{PN.get(s,s)} {100*v/tot:.0f}%" for s,v in zip(e.sp,e.excess)))
