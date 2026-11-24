# Residence district -> municipal zone; build municipal residence x overnight matrix; first descriptive pass
import duckdb, pandas as pd
con = duckdb.connect()
con.execute("""CREATE TABLE rel AS SELECT DISTINCT distrito_mitma AS dz, municipio_mitma AS mz
  FROM read_csv('raw_data_cache/v2/zonificacion/relacion_ine_zonificacionMitma.csv', delim='|', header=true, all_varchar=true)""")
print(con.execute("SELECT count(*), count(DISTINCT dz) FROM rel").fetchall())
miss = con.execute("""SELECT count(DISTINCT r) FROM 'data/derived/pernoct_2024_6m.parquet' p LEFT JOIN rel ON p.r=rel.dz WHERE rel.mz IS NULL""").fetchall()
print('unmapped residence codes', miss)
con.execute("""COPY (SELECT date, rel.mz AS r, s, sum(p) p FROM 'data/derived/pernoct_2024_6m.parquet' x JOIN rel ON x.r=rel.dz GROUP BY 1,2,3)
  TO 'data/derived/pernoct_mun_2024_6m.parquet' (FORMAT PARQUET)""")
q="'data/derived/pernoct_mun_2024_6m.parquet'"
d = con.execute(f"""SELECT date, dayname(date) dow, sum(p) tot, sum(CASE WHEN r<>s THEN p END) away,
  round(100*sum(CASE WHEN r<>s THEN p END)/sum(p),2) pct_away,
  round(100*sum(CASE WHEN substr(r,1,2)<>substr(s,1,2) THEN p END)/sum(p),2) pct_other_prov
  FROM {q} GROUP BY 1,2 ORDER BY 1""").df()
d.to_csv('explore/out_pernoct_daily.csv', index=False)
print(d.sort_values('pct_other_prov', ascending=False).head(15).to_string())
print(d.sort_values('pct_other_prov').head(5).to_string())
print(d.groupby('dow')[['pct_away','pct_other_prov']].median())
