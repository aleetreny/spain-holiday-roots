"""Verify manuscript figures and graph measurements from immutable saved inputs."""
from pathlib import Path
import json
import re

import igraph as ig
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

root = Path.cwd()
scratch = root / "research_notes/review_scratch"
snapshot = scratch / "snapshot"
res = root / "results"
net = res / "b01_networks"
fig = root / "reports/figures"
nodes = pd.read_parquet(net / "nodes.parquet").sort_values("idx")
N = len(nodes)
city = nodes.population.to_numpy() >= 100_000
coast = nodes.coastal_true.to_numpy()
x, y = nodes.x.to_numpy(), nodes.y.to_numpy()
madrid = int(nodes.loc[nodes.name == "Madrid", "idx"].iloc[0])
dates = ["2024-02-14", "2024-03-29", "2024-08-15", "2024-11-01", "2024-12-25"]

def edge(day, name):
    return pd.read_parquet(net / f"edges_{day}__{name}.parquet")

eligible = nodes.loss_peak_2021.notna() & ~nodes.id.isin(["31905_AM", "19103_AM", "40149_AM"])
quintile = pd.qcut(nodes.loc[eligible, "loss_peak_2021"].rank(method="first"), 5, labels=False)
last_quintile = quintile.index[quintile == 4]
print("Nodes", N, "municipality memberships", nodes.n_municipalities.sum(), "population", nodes.population.sum(),
      "grouped", (nodes.n_municipalities > 1).sum(), "missing loss", nodes.loss_peak_2021.isna().sum(),
      "n for R1", eligible.sum(), "quintiles", quintile.value_counts().sort_index().tolist())
calendar = pd.read_csv(fig / "fig1_calendar_data.csv", parse_dates=["date"])
for date in [calendar.date, calendar.date - pd.Timedelta(days=1)]:
    print("Aug coast wins", ((date.dt.month == 8) & (calendar.coast > calendar.towns)).sum(), "/", (date.dt.month == 8).sum(),
          "outside Aug inland losses", calendar.loc[(date.dt.month != 8) & (calendar.coast >= calendar.towns), "date"].tolist(),
          "/", (date.dt.month != 8).sum())
print("calendar count", len(calendar), calendar.groupby(calendar.date.dt.month).size().to_dict())
print("Nov missing", sorted(set(range(1, 31)) - set(calendar.loc[calendar.date.dt.month == 11].date.dt.day)))
manifest = pd.read_csv(net / "manifest.csv", dtype={"date": str})
print("sisters by file and starting weekday",
      manifest.loc[(manifest.groups.str.contains("ordinary")) & (manifest.roles.str.contains("sister|main")) & (manifest.backbone == "w10"), "date"]
      .map(lambda day: (day, pd.Timestamp(day).day_name(), (pd.Timestamp(day)-pd.Timedelta(days=1)).day_name())).tolist())
for day in dates:
    path = net / f"edges_{day}__full.parquet"
    if not path.exists():
        continue
    e = edge(day, "full")
    i, j, w = [e[c].to_numpy() for c in ("i", "j", "w")]
    assert (i != j).all() and (w > 0).all()
    inc, out = np.bincount(j, weights=w, minlength=N), np.bincount(i, weights=w, minlength=N)
    far = np.hypot(x[i] - x[j], y[i] - y[j]) >= 100000
    parts = [w[far & ~coast[j] & ~city[j]].sum(), w[far & ~coast[j] & city[j]].sum(), w[far & coast[j]].sum()]
    print(day, "full weight", w.sum(), "far parts towns/cities/coast", parts,
          "net receivers", (inc > out).sum(), "Q5", (inc[last_quintile] > out[last_quintile]).sum(),
          "Madrid balance", (inc[madrid] - out[madrid]) / (inc[madrid] + out[madrid]))

for day, name in [(dates[0],"w10"),(dates[2],"w10"),(dates[3],"w10"),(dates[2],"equal_edges")]:
    b = edge(day,name)
    f = edge(day,"full")
    if name == "w10":
        assert b.sort_values(["i","j"]).reset_index(drop=True).equals(f[f.w >= 10].sort_values(["i","j"]).reset_index(drop=True))
    else:
        assert len(b) == len(edge(dates[0],"w10")) and b.w.min() >= f.sort_values("w",ascending=False).iloc[len(b)-1].w
    i, j = b.i.to_numpy(), b.j.to_numpy()
    D = ig.Graph(n=N, edges=list(zip(map(int,i), map(int,j))), directed=True)
    G = D.as_undirected(mode="collapse")
    deg = np.array(G.degree())
    out, inc = np.array(D.outdegree()), np.array(D.indegree())
    a, c = np.array(G.get_edgelist()).T
    rank = np.lexsort((np.arange(N), -deg))
    R = np.empty(N,dtype=int)
    R[rank] = np.arange(N)
    L = G.ecount()
    Gc = G.connected_components().giant()
    ex = deg - 1
    s = np.bincount(a, weights=ex[c], minlength=N) + np.bincount(c, weights=ex[a], minlength=N)
    xa = (s[a]-ex[c])/np.maximum(ex[a],1)
    xc = (s[c]-ex[a])/np.maximum(ex[c],1)
    ra = np.corrcoef(np.r_[xa,xc], np.r_[xc,xa])[0,1]
    rk = np.maximum(R[a], R[c])
    rich_sizes = np.cumsum(np.bincount(rk,minlength=N))
    prefix = np.where(rich_sizes == np.arange(N)*(np.arange(N)+1)//2)[0].max()+1
    print(day,name,"kept",b.w.sum()/f.w.sum(),"arcs",len(b),"pairs",L,
          "dir/und densities",len(b)/(N*(N-1)),G.density(),"giant",Gc.vcount(),"isolates",(deg==0).sum(),
          "mean median max",deg.mean(),np.median(deg),deg.max(),"degree medians",np.median(out),np.median(inc),
          "Madrid out/in",out[madrid],inc[madrid],"Spearman",spearmanr(out,inc).statistic,
          "coast top 50 in",coast[np.lexsort((np.arange(N),-inc))[:50]].sum(),
          "reciprocity", D.reciprocity(), "path",Gc.average_path_length(),"lstar",np.log(Gc.vcount())/np.log(np.mean(Gc.degree())),
          "diameter",Gc.diameter(),"clustering",G.transitivity_avglocal_undirected(mode="zero"),
          "r",G.assortativity_degree(directed=False),"Ravg",ra,"phi50/121",
          rich_sizes[49]/1225,rich_sizes[120]/7260,"prefix",prefix,
          "touch top121", (np.minimum(R[a],R[c])<121).mean(), "kcore",max(G.coreness()))
    print("top 50 kinds, cities/coast",city[rank[:50]].sum(),coast[rank[:50]].sum())

sources = {}
for path in sorted(snapshot.glob("*.tex")):
    for ln, line in enumerate(path.read_text().splitlines(),1):
        for match in re.finditer(r"\\cite\w*(?:\[[^\]]*\])*\{([^}]+)\}",line):
            for key in match.group(1).split(","):
                sources.setdefault(key.strip(),[]).append(f"{path.name}:{ln}")
bib = {m[1]: ln for ln, line in enumerate((snapshot/"references.bib").read_text().splitlines(),1)
       if (m := re.match(r"@\w+\{([^,]+),",line))}
for key, loc in sorted(sources.items()):
    print(key, bib.get(key,"MISSING"), ", ".join(loc))
print("keys total",len(sources),"missing",sorted(set(sources)-set(bib)))
for f in sorted(scratch.glob("crossref_*.json")):
    if f.stat().st_size:
        m = json.loads(f.read_text())["message"]
        print("Crossref",m["DOI"],m["title"],m.get("volume"),m.get("page"),re.sub("<[^>]*>", "", m.get("abstract","")))
