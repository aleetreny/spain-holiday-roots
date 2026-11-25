from pathlib import Path
import igraph as ig
import numpy as np
import pandas as pd
import subprocess

audit=Path(__file__).resolve().parent

def pearson(x,y):
    mean=(x.mean()+y.mean())/2
    return ((x-mean)*(y-mean)).mean()/(((x-mean)**2).mean()/2+((y-mean)**2).mean()/2)

def stats(g):
    degree=np.array(g.degree(), dtype=float)
    edges=np.array(g.get_edgelist())
    u,v=edges.T
    total=(np.bincount(u,degree[v]-1,minlength=g.vcount())+
           np.bincount(v,degree[u]-1,minlength=g.vcount()))
    a=(total[u]-(degree[v]-1))/np.maximum(1,degree[u]-1)
    b=(total[v]-(degree[u]-1))/np.maximum(1,degree[v]-1)
    order=np.lexsort((np.arange(len(degree)),-degree))
    stop=0
    for index in range(1,len(order)):
        if not set(order[:index]).issubset(g.neighbors(int(order[index]))):
            break
        stop=index+1
    return dict(R_avg=pearson(a,b),clique_prefix=stop,r=g.assortativity_degree(),
                giant=g.connected_components().giant().vcount())

def main():
    out=[]
    for label in ("Feb","Aug","Nov","AugEq"):
        path=audit/f"und_{label}.txt"
        edges=np.loadtxt(path,dtype=int)
        g=ig.Graph(2435,edges)
        real=stats(g)
        out.append(dict(label=label,model="real",**real))
        d=np.array(g.degree())
        threshold=np.sort(np.min(d[edges],axis=1))
        rng=np.random.default_rng(75333+len(edges))
        for model,sign,tries,draws in (("max",1,900,2),("min",-1,900,2),("phi",0,10,20)):
            for draw in range(draws):
                perm=edges[rng.permutation(len(edges))]
                flip = rng.random(len(perm)) < 0.5
                perm[flip] = perm[flip, ::-1]
                input_path=audit/"swaps_in.txt"
                output_path=audit/"swaps_out.txt"
                np.savetxt(input_path,perm,fmt="%d")
                subprocess.run([str(audit/"swaps"), str(input_path),str(output_path),"2435",
                                str(tries),str(sign)],check=True,capture_output=True)
                new=np.loadtxt(output_path,dtype=int)
                h=ig.Graph(2435,new)
                assert h.degree()==g.degree() and h.is_simple()
                if model=="phi":
                    assert np.array_equal(np.sort(np.min(d[new],axis=1)),threshold)
                out.append(dict(label=label,model=model,draw=draw,**stats(h)))
            pd.DataFrame(out).to_csv(audit/"constraints.csv",index=False)

if __name__=="__main__":
    main()
