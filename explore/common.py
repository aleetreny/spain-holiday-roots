import duckdb, pandas as pd, numpy as np, geopandas as gpd
P = "'data/derived/od_daily_2024.parquet'"
PROV2CCAA = {**{p:'AND' for p in ['04','11','14','18','21','23','29','41']}, **{p:'ARA' for p in ['22','44','50']},
 '33':'AST','07':'BAL','35':'CAN','38':'CAN','39':'CANT', **{p:'CYL' for p in ['05','09','24','34','37','40','42','47','49']},
 **{p:'CLM' for p in ['02','13','16','19','45']}, **{p:'CAT' for p in ['08','17','25','43']}, **{p:'VAL' for p in ['03','12','46']},
 '06':'EXT','10':'EXT', **{p:'GAL' for p in ['15','27','32','36']}, '28':'MAD','30':'MUR','31':'NAV', **{p:'PVA' for p in ['01','20','48']},
 '26':'RIO','51':'CEU','52':'MEL'}
PROV_CAPITAL_INE = {'01':'01059','02':'02003','03':'03014','04':'04013','05':'05019','06':'06015','07':'07040','08':'08019','09':'09059','10':'10037',
 '11':'11012','12':'12040','13':'13034','14':'14021','15':'15030','16':'16078','17':'17079','18':'18087','19':'19130','20':'20069','21':'21041',
 '22':'22125','23':'23050','24':'24089','25':'25120','26':'26089','27':'27028','28':'28079','29':'29067','30':'30030','31':'31201','32':'32054',
 '33':'33044','34':'34120','35':'35016','36':'36038','37':'37274','38':'38038','39':'39075','40':'40194','41':'41091','42':'42173','43':'43148',
 '44':'44216','45':'45168','46':'46250','47':'47186','48':'48020','49':'49275','50':'50297','51':'51001','52':'52001'}
def nodes():
    g = gpd.read_file('clean_data/v2/zones/municipios_mitma.gpkg')
    g['area_km2'] = g.geometry.area/1e6
    cen = g.geometry.representative_point()
    g['x'] = cen.x; g['y'] = cen.y
    g['foreign'] = ~g.id.str.match(r'^\d')
    g['prov'] = np.where(g.foreign, g.id.str[:2], g.id.str[:2])
    g['ccaa'] = g.prov.map(PROV2CCAA).fillna(g.id.str[:2])
    g['population'] = pd.to_numeric(g.population, errors='coerce')
    g['island'] = g.prov.isin(['07','35','38'])
    g['canary'] = g.prov.isin(['35','38'])
    caps = set(PROV_CAPITAL_INE.values())
    g['prov_capital'] = g.municipalities.fillna('').apply(lambda s: any(m.strip() in caps for m in s.split(';')))
    return pd.DataFrame(g.drop(columns='geometry'))
def con():
    return duckdb.connect()
