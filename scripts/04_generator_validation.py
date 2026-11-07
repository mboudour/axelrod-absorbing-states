import os
from _paths import RES, FIG, TAB
import numpy as np, networkx as nx, itertools, json
from axsim import prepare, _run
def exact(G, s, F):
    out = {}
    for x, y in G.edges():
        o = int((s[x] == s[y]).sum())
        if 0 < o < F:
            for u, v in [(x, y), (y, x)]:
                for f in np.flatnonzero(s[u] != s[v]):
                    t = s.copy(); t[u, f] = s[v, f]
                    out[t.tobytes()] = out.get(t.tobytes(), 0) + o / (2 * F * (F - o))
    Z = sum(out.values()); return {k: v / Z for k, v in out.items()}, Z
res = []
cases = {'K2': (nx.complete_graph(2), [[0,0,1],[0,1,0]]), 'P3': (nx.path_graph(3), [[0,0,1],[0,1,1],[1,1,2]]),
         'K3': (nx.complete_graph(3), [[0,0,1],[0,1,1],[0,2,2]])}
for name, (G, s0) in cases.items():
    P = prepare(G); s0 = np.array(s0, dtype=np.int64); F = s0.shape[1]
    p, Z = exact(P[0], s0, F); N = 200000; cnt = {}; ts = []
    for r in range(N):
        s = s0.copy(); t, ev = _run(s, P[1], P[2], P[3], P[4], F, r + 1, 1)
        cnt[s.tobytes()] = cnt.get(s.tobytes(), 0) + 1; ts.append(t)
    tv = 0.5 * sum(abs(cnt.get(k, 0) / N - p.get(k, 0)) for k in set(cnt) | set(p))
    res.append(dict(graph=name, outcomes=len(p), tv=tv, mean_wait=float(np.mean(ts)), exact_wait=1 / Z))
    print(res[-1])
json.dump(res, open(os.path.join(RES, 'genval.json'), 'w'))
