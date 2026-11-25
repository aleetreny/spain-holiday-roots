#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <random>
#include <numeric>
#include <string>
#include <vector>
using namespace std;

struct Edge { int u, v; };
int main(int argc, char **argv) {
    if (argc != 6) return 1;
    ifstream input(argv[1]);
    ofstream output(argv[2]);
    int nodes = stoi(argv[3]), attempts_per_edge = stoi(argv[4]), sign = stoi(argv[5]);
    vector<Edge> edges;
    vector<int> deg(nodes);
    vector<uint8_t> adj(nodes*nodes);
    auto set = [&](Edge edge, bool value) {
        adj[edge.u*nodes+edge.v] = adj[edge.v*nodes+edge.u] = value;
    };
    Edge e;
    while (input >> e.u >> e.v) {
        edges.push_back(e);
        deg[e.u]++; deg[e.v]++;
        set(e,true);
    }
    mt19937_64 rng(0x9fd8e25b + edges.size());
    uniform_int_distribution<size_t> choose(0,edges.size()-1);
    const auto connected = [&]() {
        vector<int> parent(nodes);
        iota(parent.begin(),parent.end(),0);
        auto root = [&](int i) {
            while (i != parent[i]) {
                i=parent[i]=parent[parent[i]];
            }
            return i;
        };
        for (auto [u,v]:edges) parent[root(u)]=root(v);
        int start=edges.front().u;
        for (int i=0;i<nodes;i++)
            if (deg[i] && root(i)!=root(start)) return false;
        return true;
    };
    vector<Edge> checkpoint=edges;
    auto validate = [&]() {
        if (connected()) { checkpoint=edges; return; }
        for (auto edge:edges) set(edge,false);
        edges=checkpoint;
        for (auto edge:edges) set(edge,true);
    };
    long accepted = 0;
    for (long t=0, limit=long(edges.size())*attempts_per_edge; t<limit; t++) {
        if (sign>0 && t%5000==0) validate();
        size_t i=choose(rng), j=choose(rng);
        auto [a,b]=edges[i];
        auto [c,d]=edges[j];
        if (a==c || a==d || b==c || b==d) continue;
        if (rng()&1) swap(c,d);
        if (adj[a*nodes+c] || adj[b*nodes+d]) continue;
        if (sign) {
            // The endpoints keep their degrees, so the change in Pearson r has this sign.
            const int gain=(deg[a]-deg[d])*(deg[c]-deg[b]);
            if (sign*gain<=0) continue;
        } else {
            // Degree-threshold rich-club counts are the survival curve of min endpoint degree.
            int p=min(deg[a],deg[b]), q=min(deg[c],deg[d]);
            int x=min(deg[a],deg[c]), y=min(deg[b],deg[d]);
            if (!((p==x && q==y)||(p==y && q==x))) continue;
        }
        set({a,b},false); set({c,d},false);
        set({a,c},true); set({b,d},true);
        edges[i]={a,c}; edges[j]={b,d};
        accepted++;
    }
    if (sign>0) validate();
    for (auto [a,b]:edges) output << a << ' ' << b << '\n';
    cerr << "accepted " << accepted << '\n';
}
