from pathlib import Path
import duckdb
import geopandas as gpd
import numpy as np
import pandas as pd

audit = Path(__file__).resolve().parent
raw = Path("raw_data_cache/v2/zonificacion")
relation = pd.read_csv(raw/"relacion_ine_zonificacionMitma.csv", sep="|", dtype=str)
population = pd.read_csv(raw/"zonificacion_municipios/poblacion_municipios.csv",
                         sep="|", names=["id", "pop"], dtype={"id": str})
boundary = gpd.read_file("clean_data/v2/zones/municipios_mitma.gpkg").set_index("id")
nodes = population.loc[population.id.str.match(r"^\d") &
                       ~population.id.str[:2].isin(["07","35","38","51","52"]) &
                       population.id.str.len().le(8)].sort_values("id").reset_index(drop=True)
nodes["idx"] = nodes.index
nodes["prov"] = nodes.id.str[:2]
geometry = boundary.geometry.reindex(nodes.id)
nodes["name"] = nodes.id.map(boundary.name)
for label, points in [("", geometry.representative_point()), ("c", geometry.centroid)]:
    nodes[label+"x"] = points.x.to_numpy()
    nodes[label+"y"] = points.y.to_numpy()
nodes["area"] = geometry.area.to_numpy()/1e6
mun_map = relation[["municipio_ine","municipio_mitma"]].drop_duplicates()
nodes["n_mun"] = nodes.id.map(mun_map.groupby("municipio_mitma").size())
coast = pd.read_csv("explore/out_p86_coastal_flags.csv",dtype={"id":str}).set_index("id")
nodes["coast"] = nodes.id.map(coast.coastal_true).astype(bool)
nodes["city"] = nodes["pop"] >= 100_000
nodes["kind"] = np.select([nodes.coast,nodes.city],["coast","city"],default="inland")
hist = pd.read_csv("_local/external/ine_historical_pop/municipal_population_1900_2025.csv",
                   dtype={"ine_code":str})
peaks = [f"pop_{year}" for year in range(1900,1970,10)]
hist = (hist.merge(mun_map,left_on="ine_code",right_on="municipio_ine")
        .groupby("municipio_mitma")[peaks+["pop_2021"]].sum(min_count=1))
nodes["loss"] = nodes.id.map(1-hist.pop_2021/hist[peaks].max(axis=1))
nodes.to_parquet(audit/"nodes.parquet",index=False)
ref = pd.read_parquet("results/b01_networks/nodes.parquet")
print("nodes",len(nodes),"coast",nodes.coast.sum(),"city",nodes.city.sum(),
      "id equal",nodes.id.equals(ref.id),"population max error",max(abs(nodes["pop"]-ref.population)),
      "coast equal",nodes.coast.equals(ref.coastal_true),"loss max error",
      (nodes.loss-ref.loss_peak_2021).abs().max(),"municipal count differences",
      (nodes.n_mun!=ref.n_municipalities).sum(),flush=True)

db = duckdb.connect(config={"temp_directory":str(audit/"duckdb_tmp"),"threads":4})
db.register("zones",nodes)
mapping = relation[["distrito_mitma","municipio_mitma"]].drop_duplicates()
print("nonunique mapped source districts", mapping.distrito_mitma.duplicated().sum(),flush=True)
db.register("mapping",mapping)
source = """FROM 'data/derived/pernoct_2024_6m.parquet' f
JOIN mapping m ON f.r=m.distrito_mitma
JOIN zones a ON m.municipio_mitma=a.id
JOIN zones b ON f.s=b.id
WHERE a.idx!=b.idx"""
db.execute(f"""CREATE TEMP TABLE ods AS
    SELECT date, a.idx::INTEGER i, b.idx::INTEGER j, CAST(SUM(CAST(p AS DECIMAL(18,3))) AS DOUBLE) w
    {source} GROUP BY 1,2,3""")
print(db.sql("SELECT count(*) rows,count(distinct date) dates FROM ods").fetchall(),flush=True)
for date in ["2024-02-14","2024-08-15","2024-11-01","2024-12-25","2024-03-29"]:
    frame = db.execute("SELECT i,j,w FROM ods WHERE date=? ORDER BY i,j",[date]).df()
    frame.to_parquet(audit/f"edges_{date}.parquet",index=False)
    for mode in ["full","w10"]:
        path = Path(f"results/b01_networks/edges_{date}__{mode}.parquet")
        if not path.exists():
            continue
        edges = frame if mode=="full" else frame[frame.w>=10]
        reference = pd.read_parquet(path)
        pair = edges.merge(reference,on=["i","j"],how="outer",indicator=True)
        print(date,mode,len(edges),len(reference),edges.w.sum(),reference.w.sum(),
              pair._merge.value_counts().to_dict(),(pair.w_x-pair.w_y).abs().max(),
              len(set(edges.i)|set(edges.j)),len(set(reference.i)|set(reference.j)),flush=True)

def out(query, name):
    db.sql(query).df().to_parquet(audit/f"{name}.parquet",index=False)

out("""SELECT date,i,prov,coast,SUM(w) w FROM ods JOIN zones ON j=idx GROUP BY ALL""","city_prov")
out("""SELECT date,idx,
       COALESCE(v.visitors,0) visitors, COALESCE(t.travelers,0) travelers,
       COALESCE(v.in_degree,0) in_degree, COALESCE(t.out_degree,0) out_degree
       FROM (SELECT DISTINCT date FROM ods) d CROSS JOIN zones
       LEFT JOIN (SELECT date,j,SUM(w) visitors,COUNT_IF(w>=10) in_degree FROM ods GROUP BY ALL) v
          ON v.date=d.date AND v.j=zones.idx
       LEFT JOIN (SELECT date,i,SUM(w) travelers,COUNT_IF(w>=10) out_degree FROM ods GROUP BY ALL) t
          ON t.date=d.date AND t.i=zones.idx""".replace("SELECT date,idx", "SELECT d.date,idx"),"node_daily")
out("""SELECT date,b.kind,SUM(w) total,
    SUM(w) FILTER (WHERE (a.x-b.x)^2+(a.y-b.y)^2>=50000^2) dist50,
    SUM(w) FILTER (WHERE (a.x-b.x)^2+(a.y-b.y)^2>=100000^2) dist100,
    SUM(w) FILTER (WHERE (a.x-b.x)^2+(a.y-b.y)^2>=200000^2) dist200,
    SUM(w) FILTER (WHERE (a.cx-b.cx)^2+(a.cy-b.cy)^2>=100000^2) centroid100
    FROM ods JOIN zones a ON i=a.idx JOIN zones b ON j=b.idx GROUP BY ALL""","calendar")
db.execute(f"""COPY (
    SELECT date,i,j,w FROM ods WHERE w>=10 AND (
        (month(date)=2 AND isodow(date) in (2,3,4)) OR
        date BETWEEN '2024-08-10' AND '2024-08-18' OR
        date BETWEEN '2024-11-01' AND '2024-11-03')
    ) TO '{audit/"sister_w10.parquet"}' (FORMAT PARQUET)""")
db.execute(f"""COPY (
    SELECT date,a.idx::INTEGER i,b.idx::INTEGER j,CAST(SUM(CAST(p AS DECIMAL(18,3))) AS DOUBLE) w
    FROM 'data/derived/pernoct_mun_2022_2023.parquet' f
    JOIN zones a ON f.r=a.id JOIN zones b ON f.s=b.id
    WHERE a.idx!=b.idx
    GROUP BY ALL) TO '{audit/"earlier.parquet"}' (FORMAT PARQUET)""")
