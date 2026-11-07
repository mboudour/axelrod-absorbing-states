import os
from _paths import RES, FIG, TAB
import networkx as nx, numpy as np, pygraphviz as pgv, json
from axsim import *
PAL=['#2a78d6','#eb6834','#1baf7a','#eda100','#e87ba4','#008300']
G0=nx.karate_club_graph(); P=prepare(G0); G=P[0]
F,q,seed=3,6,23
s,t,ev,Q0=simulate(P,F,q,seed); d=diagnose(G,s)
dom=d['dom']; cult={x:tuple(int(c) for c in s[x]) for x in G}
cults=sorted(set(cult.values())); cidx={c:i for i,c in enumerate(cults)}
# order domains by size, label D1..
sizes={}
for x,k in dom.items(): sizes[k]=sizes.get(k,0)+1
order=sorted(sizes,key=lambda k:-sizes[k]); dname={k:f'D{i+1}' for i,k in enumerate(order)}
A=pgv.AGraph(strict=True); A.graph_attr.update(overlap='false',splines='true',outputorder='edgesfirst',pad='0.1')
for x in G:
    A.add_node(x,label=str(x+1),shape='circle',style='filled',fillcolor=PAL[cidx[cult[x]]],fontcolor='white',fontsize='9',width='0.32',fixedsize='true',color='white',penwidth='1.5',fontname='Helvetica')
for u,v in G.edges():
    if dom[u]==dom[v]: A.add_edge(u,v,color='#333333',penwidth='1.4')
    else: A.add_edge(u,v,color='#9a9a9a',style='dashed',penwidth='0.8')
A.layout(prog='neato',args='-Gseed=4')
A.draw(os.path.join(FIG, 'karate_absorbing.pdf'))
# quotient
from collections import Counter
W=Counter()
for u,v in G.edges():
    a,b=sorted((dname[dom[u]],dname[dom[v]])); W[(a,b)]+=1
Q=pgv.AGraph(strict=False); Q.graph_attr.update(overlap='false',pad='0.2',splines='true')
rep={dname[dom[x]]:x for x in G}
for k in order:
    nm=dname[k]; c=cidx[cult[rep[nm]]]; sz=sizes[k]
    Q.add_node(nm,label=f'{nm}\\n{sz}',shape='circle',style='filled',fillcolor=PAL[c],fontcolor='white',fontsize='22',width=str(round(0.75*sz**0.5,3)),fixedsize='shape',color='white',fontname='Helvetica')
for (a,b),w in W.items():
    if a!=b: Q.add_edge(a,b,penwidth=str(1.5+0.6*w),color='#777777',label=str(w),fontsize='22',fontname='Helvetica')
Q.layout(prog='neato',args='-Gseed=2 -Gsep=+12')
Q.draw(os.path.join(FIG, 'karate_quotient.pdf'))
# weighted culture assortativity of quotient (self-loop weights = interior edges)
labs=sorted(set(cidx.values())); L=len(labs); e=np.zeros((L,L))
for (a,b),w in W.items():
    ca,cb=cidx[cult[rep[a]]],cidx[cult[rep[b]]]
    e[ca,cb]+=w; e[cb,ca]+=w
e/=e.sum(); aa=e.sum(1); rq=(np.trace(e)-(aa**2).sum())/(1-(aa**2).sum())
info=dict(F=F,q=q,seed=seed,ndom=d['ndom'],ncul=d['ncul'],split=d['split'],r_f=d['r_f'],r_cul=d['r_cul'],r_dom=d['r_dom'],iota=d['iota'],rq=rq,
 rf_each=[r_of_partition(G,{x:s[x,f] for x in G})[0] for f in range(F)],nbound=sum(1 for u,v in G.edges() if dom[u]!=dom[v]),m=G.number_of_edges(),
 omega_G=d['omega_G'],omega_GI=d['omega_GI'],t=t,sizes=[sizes[k] for k in order],cultures={dname[k]:cult[rep[dname[k]]] for k in order},W={f'{a}-{b}':w for (a,b),w in W.items()})
json.dump(info,open(os.path.join(RES, 'karate_info.json'),'w'),default=float,indent=1); print(info)
