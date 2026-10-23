"""Generate vector figures directly from reconstructed graphs and CSV data."""
import csv
from fractions import Fraction
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from verify import construct

ROOT=Path(__file__).resolve().parent
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,
                     'pdf.fonttype':42,'ps.fonttype':42,'savefig.bbox':'tight'})
colors=['#006B74','#B54132']

def graph_figure():
    fig,axs=plt.subplots(1,3,figsize=(7.1,2.5))
    for ax,(n,st,title) in zip(axs,[(2,False,'(a) Coronene: $H_2$'),(2,True,'(b) Face-stellated: $T_2$'),(3,True,'(c) Coordinate shells: $T_3$')]):
        V,E,A,C=construct(n,st)
        pos=np.array([(x+y/2,np.sqrt(3)*y/2) for x,y in V])
        ax.add_collection(LineCollection([pos[[u,v]] for u,v in E],colors='#92989A',linewidths=.65))
        if n==3:
            rho=[max(abs(x),abs(y),abs(x+y)) for x,y in V]
            ax.scatter(*pos.T,c=rho,cmap='viridis',s=12,zorder=3)
        else:
            ax.scatter(*pos.T,c=[colors[1] if v in C else colors[0] for v in V],s=20,zorder=3)
        ax.autoscale();ax.set_aspect('equal');ax.axis('off');ax.set_title(title,fontsize=10)
    fig.tight_layout();fig.savefig(ROOT/'graphs.pdf');plt.close(fig)

def figures():
    graph_figure()
    rows=list(csv.DictReader((ROOT/'indices.csv').open()))
    fig,axs=plt.subplots(1,3,figsize=(7.1,2.55))
    for family,color,label in [('Cor',colors[0],'$H_n$'),('St',colors[1],'$T_n$')]:
        r=[v for v in rows if v['family']==family]
        n=np.array([int(v['n']) for v in r])
        axs[0].plot(n,[float(Fraction(v['AD']))/int(v['n']) for v in r],color=color,label=label,marker='.',ms=4)
        axs[1].plot(n,[int(v['M1'])/int(v['vertices']) for v in r],color=color,marker='.',ms=4)
        axs[2].plot(n,[int(v['Mo'])/(int(v['edges'])*int(v['vertices'])) for v in r],color=color,marker='.',ms=4)
    for ax,label in zip(axs,['$AD(G)/n$','$M_1(G)/|V(G)|$','$Mo(G)/(|E(G)||V(G)|)$']):
        ax.set_xlabel('$n$');ax.set_ylabel(label);ax.grid(alpha=.2)
    axs[0].legend(frameon=False);fig.tight_layout();fig.savefig(ROOT/'normalized_indices.pdf');plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(6.8,2.5))
    st=[v for v in rows if v['family']=='St'];n=np.array([int(v['n']) for v in st])
    for ax,name,base,amp in [(axs[0],'xi',171*n**3-423*n*n/2+102*n-75/4,.75),
                            (axs[1],'zeta',57*n**3/2-81*n*n/4+13*n/2-11/8,.375)]:
        actual=np.array([int(v[name]) for v in st])-base
        ax.axhline(0,color='#888888',lw=.6);ax.plot(n,actual,'o-',color=colors[0],ms=4,lw=.8)
        ax.set_xlabel('$n$');ax.set_ylabel('$'+('\\xi' if name=='xi' else '\\zeta')+'(T_n)$ minus cubic part')
        ax.set_yticks([-amp,0,amp]);ax.set_ylim(-amp*1.4,amp*1.4);ax.grid(alpha=.2)
    fig.tight_layout();fig.savefig(ROOT/'parity.pdf');plt.close(fig)
    hist=list(csv.DictReader((ROOT/'distance_partitions.csv').open()))
    fig,ax=plt.subplots(figsize=(5.8,2.5))
    for family,color in [('Cor',colors[0]),('St',colors[1])]:
        r=[v for v in hist if v['family']==family and v['n']=='5'];total=sum(int(v['pairs']) for v in r)
        ax.plot([int(v['distance']) for v in r],[int(v['pairs'])/total for v in r],'.-',color=color,label='$H_5$' if family=='Cor' else '$T_5$')
    ax.set_xlabel('Shortest-path distance');ax.set_ylabel('Fraction of unordered pairs');ax.legend(frameon=False);ax.grid(alpha=.2)
    fig.tight_layout();fig.savefig(ROOT/'distance_distribution.pdf');plt.close(fig)
    print('FIGURES GENERATED: 4')

if __name__=='__main__':figures()
