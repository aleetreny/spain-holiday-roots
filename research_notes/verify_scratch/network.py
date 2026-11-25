from pathlib import Path
import math
import random
import warnings
import igraph as ig
import numpy as np
import pandas as pd
import powerlaw
from scipy.stats import spearmanr

audit = Path(__file__).resolve().parent

def mixing(x,y):
    center = (x.mean()+y.mean())/2
    return ((x-center)*(y-center)).mean()/(((x-center)**2).mean()/2+((y-center)**2).mean()/2)

def second(g):
    u,v = np.asarray(g.get_edgelist()).T
    degrees = np.array(g.degree())
    totals = (np.bincount(u,degrees[v],minlength=g.vcount())+
              np.bincount(v,degrees[u],minlength=g.vcount()))
    valid = (degrees[u]>1)&(degrees[v]>1)
    a,b = u[valid],v[valid]
    return dict(Ravg=mixing((totals[a]-degrees[b])/(degrees[a]-1),
                            (totals[b]-degrees[a])/(degrees[b]-1)),
                Ravg_paper_literal=mixing((totals[u]-degrees[v]-(degrees[u]-1))/degrees[u],
                                           (totals[v]-degrees[u]-(degrees[v]-1))/degrees[v]),
                Ravg_incl=mixing(totals[u]/degrees[u],totals[v]/degrees[v]))

def measures(g, order, extra=True):
    giant = g.connected_components().giant()
    row = dict(N=g.vcount(), giant=giant.vcount(), L=g.ecount(), density=g.density(),
               k_mean=np.mean(g.degree()), k_median=np.median(g.degree()), k_max=g.maxdegree(),
               r=g.assortativity_degree(directed=False), path=giant.average_path_length(),
               path_star=np.log(g.vcount())/np.log(np.mean(g.degree())),
               diam=giant.diameter(), C=g.transitivity_avglocal_undirected(mode="zero"),
               C_drop=g.transitivity_avglocal_undirected(mode="nan"), core=max(g.coreness()))
    if extra:
        row |= {f"phi{n}":g.induced_subgraph(order[:n]).density() for n in (50,121)}
        row["clique50"] = g.induced_subgraph(order[:50]).clique_number()
        rank = np.lexsort((np.arange(g.vcount()), -np.array(g.degree())))
        row["clique_prefix"] = 0
        for end in range(1,len(rank)):
            if not set(rank[:end]).issubset(g.neighbors(int(rank[end]))):
                break
            row["clique_prefix"] = end+1
        top = np.zeros(g.vcount(),bool)
        top[order[:121]]=True
        edges = np.asarray(g.get_edgelist())
        row["touch121"] = np.any(top[edges],axis=1).mean()
    return row

def graph_for(date, condition, n, reference):
    df = pd.read_parquet(audit/f"edges_{date}.parquet")
    cut = df[df.w>=10] if condition=="w10" else df.sort_values(["w","i","j"],ascending=[False,True,True]).head(reference)
    g = ig.Graph(n,cut[["i","j"]].to_numpy(),directed=True)
    return df,cut,g,g.as_undirected()

def run():
    nodes = pd.read_parquet(audit/"nodes.parquet")
    london_seed = 734815
    ig.set_random_number_generator(random.Random(london_seed))
    rand = []
    real = []
    tails = []
    for label, date, mode in [
        ("Feb","2024-02-14","w10"),("Aug","2024-08-15","w10"),
        ("Nov","2024-11-01","w10"),("AugEq","2024-08-15","equal_edges")]:
        n = len(nodes)
        df,cut,di,g = graph_for(date,mode,n,35812)
        indeg,outdeg = np.array(di.indegree()), np.array(di.outdegree())
        degree = np.array(g.degree())
        greatest_neighbor = np.array([max(degree[g.neighbors(i)],default=-1) for i in range(n)])
        order = np.lexsort((nodes.idx,-greatest_neighbor,-degree))
        rank_in = np.lexsort((nodes.idx,-indeg))
        madrid=int(nodes.index[nodes.id=="28079"][0])
        row = dict(label=label,date=date,mode=mode, full=df.w.sum(), retained=cut.w.sum()/df.w.sum(),
                   directed=di.ecount(), med_in=np.median(indeg),med_out=np.median(outdeg),
                   recip=di.reciprocity(), spear=spearmanr(indeg,outdeg).statistic,
                   madrid_in=indeg[madrid],madrid_out=outdeg[madrid],
                   madrid_neigh=degree[madrid]/n,
                   coast50=nodes.coast.iloc[rank_in[:50]].mean(),
                   coast121=nodes.coast.iloc[order[:121]].mean(),
                   city121=nodes.city.iloc[order[:121]].mean(),
                   above_equal_boundary=(df.w>cut.w.min()).sum(),
                   ties_at_boundary=(df.w==cut.w.min()).sum())
        row |= measures(g,order)|second(g)
        real.append(row)
        print(label,row,flush=True)
        pd.DataFrame(real).to_csv(audit/"network_real.csv",index=False)
        np.savetxt(audit/f"und_{label}.txt",np.asarray(g.get_edgelist()),fmt="%d")
        for kind,k in [("in",indeg),("out",outdeg),("und",degree)]:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                fit = powerlaw.Fit(k[k>0], discrete=True,verbose=False)
                lr,p = fit.distribution_compare("power_law","lognormal")
            tails.append(dict(label=label,kind=kind,xmin=fit.xmin,R=lr,p=p))
        pd.DataFrame(tails).to_csv(audit/"tails.csv",index=False)
        for draw in range(16):
            rew = g.copy()
            rew.rewire(n=10*g.ecount(),allowed_edge_types="simple")
            tries=0
            assert sorted(g.degree())==sorted(rew.degree()) and rew.is_simple()
            dr = di.copy()
            dr.rewire(n=10*di.ecount(),allowed_edge_types="simple")
            assert di.indegree()==dr.indegree() and di.outdegree()==dr.outdegree()
            er = ig.Graph.Erdos_Renyi(n=n,m=g.ecount())
            rand.append(dict(label=label,draw=draw,model="rew",tries=tries,recip=dr.reciprocity(),**measures(rew,order)))
            rand.append(dict(label=label,draw=draw,model="er",**measures(er,order,False)))
            pd.DataFrame(rand).to_csv(audit/"network_null.csv",index=False)

if __name__=="__main__":
    run()
