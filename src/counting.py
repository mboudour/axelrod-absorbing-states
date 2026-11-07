"""Absorbing configurations of the Axelrod model on small graphs: brute force and the chromatic-sum formula."""
import itertools
import numpy as np
import networkx as nx

def set_partitions(items):
    if not items: yield []; return
    first, rest = items[0], items[1:]
    for p in set_partitions(rest):
        for i in range(len(p)):
            yield p[:i] + [[first] + p[i]] + p[i+1:]
        yield [[first]] + p
def chrom(H, q):
    nodes = list(H); idx = {v: i for i, v in enumerate(nodes)}
    return sum(all(c[idx[u]] != c[idx[v]] for u, v in H.edges()) for c in itertools.product(range(q), repeat=len(nodes)))
def formula(G, F, q):
    tot = 0
    for P in set_partitions(list(G)):
        if not all(nx.is_connected(G.subgraph(B)) for B in P): continue
        b = {x: i for i, B in enumerate(P) for x in B}
        H = nx.Graph(); H.add_nodes_from(range(len(P)))
        H.add_edges_from((b[u], b[v]) for u, v in G.edges() if b[u] != b[v])
        tot += chrom(H, q) ** F
    return tot
def brute(G, F, q):
    cults = list(itertools.product(range(q), repeat=F)); C = len(cults)
    O = np.array([[sum(a == b for a, b in zip(c1, c2)) for c2 in cults] for c1 in cults])
    ok = (O == 0) | (O == F)
    n = G.number_of_nodes()
    grids = np.indices((C,) * n).reshape(n, -1)
    mask = np.ones(grids.shape[1], bool)
    for u, v in G.edges(): mask &= ok[grids[u], grids[v]]
    return int(mask.sum())

def tree_count(n, F, q):
    """|A(T)| = q^F (1+(q-1)^F)^(n-1) for any tree on n vertices."""
    return q**F * (1 + (q - 1)**F) ** (n - 1)

def cycle_count(n, F, q):
    """|A(C_n)| = tr M^n with M = I + (J_q - I_q)^{(x)F}."""
    from math import comb
    return sum(comb(F, j) * (q - 1)**j * (1 + (-1)**j * (q - 1)**(F - j))**n for j in range(F + 1))
