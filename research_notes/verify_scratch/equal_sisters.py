from pathlib import Path
import duckdb
import igraph as ig
import numpy as np
import pandas as pd

audit=Path(__file__).resolve().parent
db=duckdb.connect(config={"threads":4,"temp_directory":str(audit/"duckdb_tmp")})
nodes=pd.read_parquet(audit/"nodes.parquet")
db.register("zone",nodes[["id","idx"]])
mapping=pd.read_csv("raw_data_cache/v2/zonificacion/relacion_ine_zonificacionMitma.csv",sep="|",dtype=str)
db.register("crosswalk",mapping[["distrito_mitma","municipio_mitma"]].drop_duplicates())
data=db.sql("""
    WITH flows AS (
        SELECT f.date,a.idx i,b.idx j,SUM(CAST(p AS DECIMAL(18,3))) w
        FROM 'data/derived/pernoct_2024_6m.parquet' f
        JOIN crosswalk c ON f.r=c.distrito_mitma
        JOIN zone a ON a.id=c.municipio_mitma
        JOIN zone b ON b.id=f.s
        WHERE a.idx!=b.idx AND (
           (month(f.date)=2 AND isodow(f.date) IN (2,3,4)) OR
           f.date BETWEEN '2024-08-10' AND '2024-08-18' OR
           f.date BETWEEN '2024-11-01' AND '2024-11-03'
        )
        GROUP BY ALL HAVING w>0
    ) SELECT date,i,j,CAST(w AS DOUBLE) w
    FROM flows
    QUALIFY ROW_NUMBER() OVER(PARTITION BY date ORDER BY w DESC,i,j)<=35812
""").df()
data.date=pd.to_datetime(data.date).dt.strftime("%Y-%m-%d")
rows=[]
for date, frame in data.groupby("date"):
    test=Path(f"results/b01_networks/edges_{date}__equal_edges.parquet")
    if test.exists():
        ref=pd.read_parquet(test)
        compared=frame.merge(ref,on=["i","j"],how="outer",indicator=True)
        print(date,compared._merge.value_counts().to_dict(),(compared.w_x-compared.w_y).abs().max())
        if compared._merge.ne("both").any():
            print(compared[compared._merge.ne("both")].to_dict(orient="records"))
    graph=ig.Graph(len(nodes),frame[["i","j"]].to_numpy(),directed=True).as_undirected()
    rank=np.lexsort((nodes.idx,-np.asarray(graph.degree())))
    clique=0
    for i in range(1,len(nodes)):
        if not set(rank[:i]).issubset(graph.neighbors(int(rank[i]))):
            break
        clique=i+1
    rows.append(dict(date=date,diameter=graph.connected_components().giant().diameter(),
                     clique=clique,phi50=graph.induced_subgraph(rank[:50]).density(),
                     C=graph.transitivity_avglocal_undirected(mode="zero"),
                     links=graph.ecount()))
rows=pd.DataFrame(rows)
rows.to_csv(audit/"sister_eq_network.csv",index=False)
for x in ["diameter","clique","phi50","C"]:
    f=rows[rows.date.str[5:7]=="02"][x].to_numpy()
    a=rows[rows.date.str[5:7]=="08"][x].to_numpy()
    d=(a[:,None]-f[None,:]).ravel()
    print("AugEq",x,int(sum(d<0)),int(sum(d==0)),int(sum(d>0)))
