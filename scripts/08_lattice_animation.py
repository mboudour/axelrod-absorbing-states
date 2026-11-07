import os
from _paths import RES, FIG, TAB
import numpy as np, networkx as nx, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.animation import FFMpegWriter
from axsim import prepare, _run

F, q, seed = 3, 10, 0
G0 = nx.triangular_lattice_graph(20, 40)
pos0 = nx.get_node_attributes(G0, 'pos')
P = prepare(G0); G, eu, ev, ptr, flat = P
mapping = {i: v for i, v in enumerate(G0.nodes())}
pos = np.array([pos0[mapping[i]] for i in range(G.number_of_nodes())])
n, m = G.number_of_nodes(), len(eu)
sigma = np.random.default_rng(seed).integers(0, q, size=(n, F)).astype(np.int64)
tris = np.array([c for c in nx.enumerate_all_cliques(G) if len(c) == 3])
eidx = {}
for k, (u, v) in enumerate(zip(eu, ev)):
    eidx[(u, v)] = k; eidx[(v, u)] = k
tri_e = np.array([[eidx[(a, b)], eidx[(b, c)], eidx[(a, c)]] for a, b, c in tris])

def state(s):
    o = (s[eu] == s[ev]).sum(1)
    inter, bnd = o == F, o == 0
    act = ~inter & ~bnd
    nb = bnd[tri_e].sum(1)
    H = nx.Graph(); H.add_nodes_from(range(n)); H.add_edges_from(zip(eu[inter], ev[inter]))
    return o, inter, act, bnd, (nb == 1).mean(), nx.number_connected_components(H)

# record frames at geometrically spaced event counts until absorption
frames, times, t = [], [], 0.0
targets = np.unique(np.round(np.geomspace(1, 400000, 260)).astype(int))
done_events, k = 0, 0
frames.append(sigma.copy()); times.append(0.0)
while True:
    step = int(targets[k] - done_events) if k < len(targets) else 5000
    k += 1
    if step <= 0: continue
    dt, e = _run(sigma, eu, ev, ptr, flat, F, 7919 * k + 1, step)
    t += dt; done_events += e
    frames.append(sigma.copy()); times.append(t)
    if e < step: break
print('absorbed after', done_events, 'events, time', round(t, 1), 'frames', len(frames))

stats = [state(s) for s in frames]
cols = lambda s: 0.12 + 0.8 * s / (q - 1)
fig = plt.figure(figsize=(12.8, 7.2), dpi=100)
axL = fig.add_axes([0.02, 0.10, 0.56, 0.80]); axR = fig.add_axes([0.65, 0.30, 0.32, 0.55])
segs = np.stack([pos[eu], pos[ev]], 1)
T = np.array(times); T[0] = T[1] / 2
series = np.array([[st[1].mean(), st[2].mean(), st[3].mean(), st[4]] for st in stats])
writer = FFMpegWriter(fps=15, bitrate=4000)
hold = 60
txt = fig.text(0.65, 0.17, '', fontsize=12, va='top')
with writer.saving(fig, os.path.join(FIG, 'axelrod_triangular_lattice.mp4'), dpi=100):
    for i in list(range(len(frames))) + [len(frames) - 1] * hold:
        s = frames[i]; o, inter, act, bnd, d1, ndom = stats[i]
        axL.clear(); axR.clear(); axL.set_aspect('equal'); axL.axis('off')
        axL.add_collection(LineCollection(segs[bnd], colors='#c8c8c8', linewidths=0.6, linestyles=(0, (2, 2))))
        axL.add_collection(LineCollection(segs[inter], colors='#333333', linewidths=1.2))
        axL.add_collection(LineCollection(segs[act], colors='#eb6834', linewidths=1.6))
        axL.scatter(pos[:, 0], pos[:, 1], s=34, c=cols(s), edgecolors='white', linewidths=0.4, zorder=3)
        axL.set_xlim(pos[:, 0].min() - .5, pos[:, 0].max() + .5); axL.set_ylim(pos[:, 1].min() - .5, pos[:, 1].max() + .5)
        axL.set_title(f'Axelrod model on a triangular lattice  (n={n}, F={F}, q={q})', fontsize=13, loc='left')
        axL.text(0, -0.06, 'Vertex colour = culture (features → RGB).  Solid dark: interior edge (all features shared).  '
                 'Orange: active edge.  Dashed grey: boundary edge (nothing shared).', transform=axL.transAxes, fontsize=8.5, color='#555555')
        x = T[:i + 1]
        for j, (lab, c) in enumerate([('interior edges', '#333333'), ('active edges', '#eb6834'), ('boundary edges', '#2a78d6'),
                                      ('triangles with exactly one\nboundary edge', '#1baf7a')]):
            axR.plot(x, series[:i + 1, j], color=c, lw=2, label=lab)
        axR.set_xscale('log'); axR.set_xlim(T[0], T[-1] * 1.2); axR.set_ylim(0, 1)
        axR.set_xlabel('time'); axR.set_ylabel('fraction'); axR.spines[['top', 'right']].set_visible(False)
        axR.grid(color='#eeeeee'); axR.legend(frameon=False, fontsize=9, loc='upper left')
        msg = f't = {times[i]:.1f}    domains = {ndom}'
        if i == len(frames) - 1:
            msg += '\nAbsorbed: no active edges, and no triangle\nhas exactly one boundary edge (Theorem 4.2).'
        txt.set_text(msg)
        writer.grab_frame()
print('final', stats[-1][5], 'domains, Δ1 share', stats[-1][4], 'boundary share', round(stats[-1][3].mean(), 3))
