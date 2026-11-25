import duckdb, pandas as pd
con=duckdb.connect()
rel=con.execute("""SELECT DISTINCT distrito_mitma dz, municipio_mitma mz FROM read_csv('raw_data_cache/v2/zonificacion/relacion_ine_zonificacionMitma.csv', delim='|', header=true, all_varchar=true)""").df()
print('rel rows',len(rel),'distinct dz',rel.dz.nunique(),'distinct mz',rel.mz.nunique())
m=rel.groupby('dz').mz.nunique(); print('districts mapping to >1 mz:', (m>1).sum()); print(rel[rel.dz.isin(m[m>1].index)].head(20))
print(rel.head())
raw="'data/derived/pernoct_2024_6m.parquet'"
print(con.execute(f"select count(distinct r), count(distinct s), min(r),max(r),min(s),max(s) from {raw}").fetchall())
# r codes not in rel
print('unmapped r', con.execute(f"select r, sum(p) from {raw} where r not in (select dz from rel) group by 1 order by 2 desc limit 20").df())
print('s not municipal-zone codes', con.execute(f"select count(distinct s) from {raw} where s not in (select mz from rel)").fetchall())
print(con.execute(f"select s, sum(p) from {raw} where s not in (select mz from rel) group by 1 order by 2 desc limit 10").df())
# totals before/after mapping on one day
mun="'data/derived/pernoct_mun_2024_6m.parquet'"
for d in ['2024-02-14','2024-08-15']:
    print(d, con.execute(f"select sum(p) from {raw} where date='{d}'").fetchall(), con.execute(f"select sum(p) from {mun} where date='{d}'").fetchall())
print(con.execute(f"select count(distinct date), min(date), max(date) from {mun}").fetchall())
print(con.execute(f"select strftime(date,'%Y-%m') m, count(distinct date) from {mun} group by 1 order by 1").df())
