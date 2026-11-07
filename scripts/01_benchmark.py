"""Section 7.2: 40 runs per graph and (F,q) on the nine graphs of Table 3 (seeds 1000+j)."""
import json, os
from _paths import RES
from workloads import graphs, run_set

res = {}
for name, G in graphs().items():
    for F, q in [(2, 2), (3, 6), (3, 12)]:
        res[f'{name}|{F}|{q}'] = run_set(G, F, q, 40, seed0=1000)
        print(name, F, q, flush=True)
json.dump(res, open(os.path.join(RES, 'bench.json'), 'w'), default=float)
