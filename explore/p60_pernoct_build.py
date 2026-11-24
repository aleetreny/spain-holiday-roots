# Build one parquet with all downloaded daily municipal overnight-stay files (residence zone -> overnight zone, persons)
import duckdb, glob
con = duckdb.connect()
files = sorted(glob.glob('_local/external/mitma_pernoct/2024*_Pernoctaciones_municipios.csv.gz'))
print(len(files))
con.execute(f"""
COPY (SELECT strptime(CAST(fecha AS VARCHAR),'%Y%m%d')::DATE AS date, CAST(zona_residencia AS VARCHAR) AS r,
             CAST(zona_pernoctacion AS VARCHAR) AS s, CAST(personas AS DOUBLE) AS p
      FROM read_csv({files}, delim='|', header=true, columns={{'fecha':'VARCHAR','zona_residencia':'VARCHAR','zona_pernoctacion':'VARCHAR','personas':'DOUBLE'}}))
TO 'data/derived/pernoct_2024_6m.parquet' (FORMAT PARQUET)""")
print(con.execute("SELECT count(*), count(DISTINCT date), min(date), max(date), count(DISTINCT r), count(DISTINCT s) FROM 'data/derived/pernoct_2024_6m.parquet'").fetchall())
print(con.execute("""SELECT date, sum(p) tot, sum(CASE WHEN r<>s THEN p END) away, round(100*sum(CASE WHEN r<>s THEN p END)/sum(p),2) pct_away
 FROM 'data/derived/pernoct_2024_6m.parquet' GROUP BY 1 ORDER BY pct_away DESC LIMIT 12""").df())
print(con.execute("""SELECT date, round(100*sum(CASE WHEN r<>s THEN p END)/sum(p),2) pct_away
 FROM 'data/derived/pernoct_2024_6m.parquet' GROUP BY 1 ORDER BY pct_away ASC LIMIT 6""").df())
