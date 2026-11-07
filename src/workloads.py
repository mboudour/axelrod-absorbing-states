"""Graph set and batch runner used in Section 7."""
import json, numpy as np, networkx as nx
from axsim import prepare, simulate, diagnose, Q_of

def graphs():
    G = {}
    G['Florentine'] = nx.florentine_families_graph()
    G['Karate'] = nx.karate_club_graph()
    G['Les Mis'] = nx.les_miserables_graph()
    er = nx.gnp_random_graph(100, 0.06, seed=7)
    G['Erdos-Renyi'] = er.subgraph(max(nx.connected_components(er), key=len)).copy()
    G['Ring lattice'] = nx.watts_strogatz_graph(100, 6, 0.0, seed=1)
    G['Random 6-regular'] = nx.random_regular_graph(6, 100, seed=3)
    G['Triangular lattice'] = nx.triangular_lattice_graph(8, 20)
    G['Holme-Kim'] = nx.powerlaw_cluster_graph(100, 3, 0.6, seed=5)
    G['Caveman'] = nx.connected_caveman_graph(10, 10)
    return G

def run_set(G, F, q, runs, seed0=0):
    P = prepare(G); H = P[0]
    out = []
    for r in range(runs):
        s, t, ev, Q0 = simulate(P, F, q, seed0 + r)
        d = diagnose(H, s)
        d.pop('dom'); d.pop('H')
        d.update(t=t, Q0=Q0, Qtau=Q_of(s, q), n=H.number_of_nodes(), m=H.number_of_edges())
        out.append(d)
    return out

