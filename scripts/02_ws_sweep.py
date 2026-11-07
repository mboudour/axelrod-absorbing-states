import os
from _paths import RES, FIG, TAB
import json, numpy as np, networkx as nx
from workloads import run_set
ps=[0.0,0.005,0.01,0.02,0.05,0.1,0.2,0.5,1.0]
res={}
for p in ps:
    rows=[]
    for g in range(8):
        G=nx.connected_watts_strogatz_graph(200,6,p,seed=100+g)
        C=nx.transitivity(G)
        for d in run_set(G,3,12,10,seed0=5000+100*g):
            d['C']=C; rows.append(d)
    res[str(p)]=rows; print(p,flush=True)
json.dump(res,open(os.path.join(RES,'sweep.json'),'w'),default=float)
