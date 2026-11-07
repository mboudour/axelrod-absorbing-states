import os
from _paths import RES, FIG, TAB
import json, numpy as np, networkx as nx
from workloads import graphs
from axsim import open_wedges
R=json.load(open(os.path.join(RES, 'bench.json'))); Gs=graphs()
L0=[]
for name,G in Gs.items():
    tau=sum(nx.triangles(G).values())//3
    L0.append(f"{name} & {G.number_of_nodes()} & {G.number_of_edges()} & {tau} & {open_wedges(G)} & {nx.transitivity(G):.3f} \\\\")
open(os.path.join(TAB,'tab_graphs.tex'),'w').write('\n'.join(L0)+'\n')
L1=[];L2=[]
prev=None
for k,v in R.items():
    name,F,q=k.split('|'); F=int(F); q=int(q)
    nm=name if name!=prev else ''; prev=name
    tri=np.array([d['split'] for d in v]); T=tri.sum(1)
    cons=np.mean([d['ndom']==1 for d in v])
    nc=[d for d in v if d['ndom']>1]
    f=lambda key,vs: np.mean([d[key] for d in vs]) if vs else float('nan')
    L1.append(f"{nm} & $({F},{q})$ & {cons:.2f} & {f('ndom',v):.1f} & {1-f('iota',v):.2f} & {np.mean(tri[:,2]/T):.2f} & {np.mean(tri[:,3]/T):.3f} & {int(tri[:,1].sum())} & {np.mean([d['omega_GI']/d['omega_G'] for d in v]):.2f} & {max(d['quotient_clique'] for d in v)} \\\\")
    bound=np.mean([F*(F*d['n']**2-d['Q0'])/2 for d in v])
    def fmt(x): return '--' if np.isnan(x) else f"{x:.2f}"
    bstr = f"{bound:,.0f}".replace(',', '{,}')
    L2.append(f"{nm} & $({F},{q})$ & {len(nc)} & {fmt(f('r_f',nc))} & {fmt(f('r_cul',nc))} & {fmt(f('r_dom',nc))} & {f('t',v):.1f} & {bstr} \\\\")
    if q==12 and name!='Caveman':
        L1.append('\\midrule'); L2.append('\\midrule')
open(os.path.join(TAB,'tab_structure.tex'),'w').write('\n'.join(L1)+'\n')
open(os.path.join(TAB,'tab_assort.tex'),'w').write('\n'.join(L2)+'\n')

# typeset graph names as in the manuscript
NAMES = {'Erdos-Renyi': "Erd\\H{o}s--R\\'enyi", 'Holme-Kim': 'Holme--Kim', 'Les Mis ': "Les Mis\\'erables ", 'Random 6-regular': 'Random $6$-regular'}
for fn in ('tab_graphs.tex', 'tab_structure.tex', 'tab_assort.tex'):
    p = os.path.join(TAB, fn); t = open(p).read()
    for k, v in NAMES.items():
        t = t.replace(k, v)
    open(p, 'w').write(t)
