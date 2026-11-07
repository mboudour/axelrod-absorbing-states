"""Edge-driven Axelrod dynamics (Lanchier's convention) and absorbing-state diagnostics."""
import numpy as np
import networkx as nx
from numba import njit


@njit(cache=True)
def _run(sigma, eu, ev, inc_ptr, inc_edges, F, seed, max_events=-1):
    np.random.seed(seed)
    n = sigma.shape[0]
    m = eu.shape[0]
    o = np.zeros(m, dtype=np.int64)
    for e in range(m):
        c = 0
        for f in range(F):
            if sigma[eu[e], f] == sigma[ev[e], f]:
                c += 1
        o[e] = c
    w = np.zeros(m, dtype=np.int64)
    tot = 0
    for e in range(m):
        if 0 < o[e] < F:
            w[e] = o[e]
            tot += w[e]
    t = 0.0
    events = 0
    diff = np.empty(F, dtype=np.int64)
    while tot > 0 and (max_events < 0 or events < max_events):
        t += -np.log(np.random.random()) / (tot / F)
        r = np.random.randint(tot)
        acc = 0
        e = m - 1
        for k in range(m):
            acc += w[k]
            if acc > r:
                e = k
                break
        if np.random.random() < 0.5:
            x, y = eu[e], ev[e]
        else:
            x, y = ev[e], eu[e]
        nd = 0
        for f in range(F):
            if sigma[x, f] != sigma[y, f]:
                diff[nd] = f
                nd += 1
        f = diff[np.random.randint(nd)]
        sigma[x, f] = sigma[y, f]
        events += 1
        for k in range(inc_ptr[x], inc_ptr[x + 1]):
            e2 = inc_edges[k]
            c = 0
            for g in range(F):
                if sigma[eu[e2], g] == sigma[ev[e2], g]:
                    c += 1
            o[e2] = c
            tot -= w[e2]
            if 0 < c < F:
                w[e2] = c
            else:
                w[e2] = 0
            tot += w[e2]
    return t, events


def prepare(G):
    G = nx.convert_node_labels_to_integers(G)
    n = G.number_of_nodes()
    edges = np.array(list(G.edges()), dtype=np.int64)
    eu, ev = edges[:, 0].copy(), edges[:, 1].copy()
    inc = [[] for _ in range(n)]
    for k, (u, v) in enumerate(edges):
        inc[u].append(k)
        inc[v].append(k)
    ptr = np.zeros(n + 1, dtype=np.int64)
    for i in range(n):
        ptr[i + 1] = ptr[i] + len(inc[i])
    flat = np.array([k for l in inc for k in l], dtype=np.int64)
    return G, eu, ev, ptr, flat


def Q_of(sigma, q):
    F = sigma.shape[1]
    return sum(int((np.bincount(sigma[:, f], minlength=q) ** 2).sum()) for f in range(F))


def simulate(prep, F, q, seed, sigma0=None):
    G, eu, ev, ptr, flat = prep
    n = G.number_of_nodes()
    rng = np.random.default_rng(seed)
    sigma = rng.integers(0, q, size=(n, F)).astype(np.int64) if sigma0 is None else sigma0.copy()
    Q0 = Q_of(sigma, q)
    t, ev_count = _run(sigma, eu, ev, ptr, flat, F, int(seed) % (2**31 - 1))
    return sigma, t, ev_count, Q0


def r_of_partition(G, labels):
    """Newman's categorical assortativity: (e_in - s)/(1 - s)."""
    m = G.number_of_edges()
    ein = sum(1 for u, v in G.edges() if labels[u] == labels[v]) / m
    vol = {}
    for x, d in G.degree():
        vol[labels[x]] = vol.get(labels[x], 0) + d
    s = sum((vv / (2 * m)) ** 2 for vv in vol.values())
    return (ein - s) / (1 - s) if s < 1 else np.nan, ein, s


def open_wedges(G):
    tri = nx.triangles(G)
    m2 = sum(d * (d - 1) // 2 for _, d in G.degree())
    return m2 - sum(tri.values())  # m2 - 3*tau


def diagnose(G, sigma):
    F = sigma.shape[1]
    n = G.number_of_nodes()
    cult = {x: tuple(sigma[x]) for x in G}
    o = {(u, v): sum(sigma[u, f] == sigma[v, f] for f in range(F)) for u, v in G.edges()}
    assert all(val in (0, F) for val in o.values()), "not absorbing"
    I = [e for e, val in o.items() if val == F]
    GI = nx.Graph()
    GI.add_nodes_from(G)
    GI.add_edges_from(I)
    dom = {}
    for k, comp in enumerate(nx.connected_components(GI)):
        for x in comp:
            dom[x] = k
    ndom = max(dom.values()) + 1
    ncul = len(set(cult.values()))
    r_dom, iota, s_dom = r_of_partition(G, dom)
    r_cul, _, s_cul = r_of_partition(G, cult)
    rf = []
    for f in range(F):
        rf.append(r_of_partition(G, {x: sigma[x, f] for x in G})[0])
    # triangle split by number of boundary edges
    split = [0, 0, 0, 0]
    bset = {frozenset(e) for e, val in o.items() if val == 0}
    for c in nx.enumerate_all_cliques(G):
        if len(c) == 3:
            a, b, d = c
            k = sum(frozenset(p) in bset for p in [(a, b), (b, d), (a, d)])
            split[k] += 1
        elif len(c) > 3:
            break
    # quotient graph and colourability witness (feature 1 colours)
    H = nx.Graph()
    H.add_nodes_from(range(ndom))
    for u, v in G.edges():
        if dom[u] != dom[v]:
            H.add_edge(dom[u], dom[v])
    rep = {}
    for x in G:
        rep.setdefault(dom[x], x)
    proper = all(sigma[rep[a], f] != sigma[rep[b], f] for a, b in H.edges() for f in range(F))
    return dict(iota=iota, ndom=ndom, ncul=ncul, r_dom=r_dom, r_cul=r_cul, r_f=float(np.mean(rf)),
                r_f_max=float(np.max(rf)), split=split, omega_G=open_wedges(G), omega_GI=open_wedges(GI),
                quotient_proper=proper, quotient_clique=max((len(c) for c in nx.find_cliques(H)), default=1),
                dom=dom, H=H)
