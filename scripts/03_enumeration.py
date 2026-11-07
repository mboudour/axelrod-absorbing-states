import os
from _paths import RES, FIG, TAB
import itertools, numpy as np, networkx as nx, json
from counting import set_partitions, chrom, formula, brute
graphs = {'$K_3$': nx.complete_graph(3), '$P_4$': nx.path_graph(4), '$K_{1,3}$': nx.star_graph(3), '$C_4$': nx.cycle_graph(4),
          'paw': nx.Graph([(0,1),(1,2),(0,2),(2,3)]), 'diamond': nx.Graph([(0,1),(1,2),(0,2),(1,3),(2,3)]), '$K_4$': nx.complete_graph(4),
          '$C_5$': nx.cycle_graph(5), 'bull': nx.bull_graph()}
rows = []
for name, G in graphs.items():
    for F, q in [(2,2),(2,3),(3,2),(3,3)]:
        if (q**F) ** G.number_of_nodes() > 2e7: continue
        b, f = brute(G, F, q), formula(G, F, q)
        rows.append((name, G.number_of_nodes(), F, q, b, f)); print(name, F, q, b, f, b == f, flush=True)
json.dump(rows, open(os.path.join(RES, 'enum.json'), 'w'))
