import os
from _paths import RES, FIG, TAB
import json, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'serif','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'axes.grid':True,'grid.color':'#e6e6e6','grid.linewidth':0.6,'axes.edgecolor':'#555555'})
PAL=['#2a78d6','#eb6834','#1baf7a','#eda100']
R=json.load(open(os.path.join(RES, 'sweep.json'))); ps=[float(p) for p in R]; x=np.array(ps)
rng=np.random.default_rng(2026)
def ms(vals):
    # cluster bootstrap over the 8 graph realisations (10 runs each), percentile 95% interval
    v=np.array(vals,float).reshape(8,10)
    est=np.nanmean(v); bs=[]
    for _ in range(2000):
        idx=rng.integers(0,8,8); bs.append(np.nanmean(v[idx]))
    lo,hi=np.percentile(bs,[2.5,97.5]); return est,lo,hi
fig,ax=plt.subplots(1,3,figsize=(7.2,2.5))
def line(a,ys,lab,c,mk):
    m=np.array([y[0] for y in ys]); lo=np.array([y[1] for y in ys]); hi=np.array([y[2] for y in ys])
    a.plot(x,m,color=c,lw=1.6,marker=mk,ms=4,label=lab); a.fill_between(x,lo,hi,color=c,alpha=0.18,lw=0)
for k,lab,c,mk in [(0,r'$0$ boundary edges',PAL[0],'o'),(2,r'$2$ boundary edges',PAL[1],'s'),(3,r'$3$ boundary edges',PAL[2],'^')]:
    line(ax[0],[ms([d['split'][k]/sum(d['split']) for d in R[p]]) for p in R],lab,c,mk)
ax[0].set_ylabel('fraction of triangles'); ax[0].set_title('(a) triangle composition',fontsize=9)
for k,lab,c,mk in [('r_f',r'$\bar r_{\mathrm{feat}}$',PAL[0],'o'),('r_cul',r'$r_{\mathrm{cult}}$',PAL[1],'s'),('r_dom',r'$r_{\mathrm{dom}}$',PAL[2],'^')]:
    line(ax[1],[ms([d[k] for d in R[p]]) for p in R],lab,c,mk)
ax[1].set_ylabel('assortativity'); ax[1].set_title('(b) assortativity chain',fontsize=9)
line(ax[2],[ms([1-d['iota'] for d in R[p]]) for p in R],r'$|\partial|/m$',PAL[0],'o')
line(ax[2],[ms([d['omega_GI']/d['omega_G'] for d in R[p]]) for p in R],r'$\omega(G_I)/\omega(G)$',PAL[1],'s')
ax[2].set_ylabel('fraction'); ax[2].set_title('(c) boundary and openness',fontsize=9)
for a in ax:
    a.set_xscale('symlog',linthresh=0.005); a.set_xlim(-0.0005,1.05); a.set_xlabel(r'rewiring probability $p$')
    a.set_xticks([0,0.01,0.1,1]); a.set_xticklabels(['0','0.01','0.1','1'])
    a.legend(frameon=False,fontsize=7,loc='upper left' if a is ax[0] else 'best')
ax[0].set_ylim(0,1.32); ax[0].set_yticks([0,0.2,0.4,0.6,0.8,1.0]); ax[2].set_ylim(0,1)
fig.tight_layout(); fig.savefig(os.path.join(FIG, 'ws_sweep.pdf'))
