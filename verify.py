"""Rebuild graph data and check the manuscript without its original code.

Run: python verify.py --max-n 15
Dependencies: numpy, sympy, matplotlib. Uses exact integer/rational arithmetic
for the rational indices and double precision only for irrational/log values.
"""
from collections import Counter, deque
from fractions import Fraction
from pathlib import Path
import argparse
import csv
import json
import math
import platform
import numpy as np
import sympy
from symbolic_checks import formulas, n as symbol_n, check as symbolic_check

ROOT=Path(__file__).resolve().parent
OFFSETS=((1,0),(0,1),(-1,1),(-1,0),(0,-1),(1,-1))

def construct(n, st=False):
    centers={(2*i+j,-i+j) for i in range(1-n,n) for j in range(1-n,n)
             if max(abs(i),abs(j),abs(i+j))<n}
    vertices,edges=set(),set()
    for x,y in centers:
        ring=[(x+a,y+b) for a,b in OFFSETS]
        vertices.update(ring)
        edges.update(tuple(sorted((u,v))) for u,v in zip(ring,ring[1:]+ring[:1]))
        if st:
            vertices.add((x,y))
            edges.update(tuple(sorted(((x,y),u))) for u in ring)
    vertices=sorted(vertices); ids={v:i for i,v in enumerate(vertices)}
    edges=sorted((ids[u],ids[v]) for u,v in edges)
    adj=[[] for _ in vertices]
    for u,v in edges:
        adj[u].append(v); adj[v].append(u)
    return vertices,edges,adj,centers

def bfs(adj):
    rows=[]
    for source in range(len(adj)):
        dist=[-1]*len(adj); dist[source]=0; q=deque([source])
        while q:
            u=q.popleft()
            for v in adj[u]:
                if dist[v]<0:
                    dist[v]=dist[u]+1; q.append(v)
        assert min(dist)==0
        rows.append(dist)
    return np.array(rows,dtype=np.int32)

def bounds(n,x):
    if abs(x)>2*n-1:
        return 1,0
    B=3*n-1
    return max(-((B+x)//2),-B-2*x,x-B),min((B-x)//2,B-2*x,x+B)

def coordinate_distance(V,st):
    if st:
        z=np.asarray(V)
        d=z[:,None,:]-z[None,:,:]
        return np.maximum(np.abs(d).max(axis=2),np.abs(d.sum(axis=2)))
    z=np.array([((2*x+y)//3,(x+2*y)//3,(x-y)//3) for x,y in V])
    return np.abs(z[:,None,:]-z[None,:,:]).sum(axis=2)

def expected_partition(n,st):
    if not st:
        d={(2,2):6,(2,3):12*(n-1),(3,3):9*n*n-15*n+6}
    else:
        d={(3,3):6,(3,5):12*(n-1),(3,6):6*n,(5,6):18*(n-1),(6,6):27*n*n-57*n+30}
    return {k:v for k,v in d.items() if v}

def degree_indices(part):
    weights={
        'M1':lambda a,b:a+b,'M2':lambda a,b:a*b,'F':lambda a,b:a*a+b*b,
        'R':lambda a,b:1/math.sqrt(a*b),'SC':lambda a,b:1/math.sqrt(a+b),
        'GA':lambda a,b:2*math.sqrt(a*b)/(a+b),
        'ABC':lambda a,b:math.sqrt((a+b-2)/(a*b))}
    return {name:math.fsum(count*f(*pair) for pair,count in part.items()) for name,f in weights.items()}

def indices(D,edges,deg):
    triu=np.triu_indices(len(D),1); ds=D[triu].astype(np.int64)
    histogram=np.bincount(ds)
    W=int(ds.sum()); WW=int(((ds*ds+ds)//2).sum())
    S=int(((deg[triu[0]]+deg[triu[1]])*ds).sum())
    SS=int((deg[triu[0]]*deg[triu[1]]*ds).sum())
    part=Counter(tuple(sorted((int(deg[u]),int(deg[v])))) for u,v in edges)
    out=degree_indices(part)
    for name in ('M1','M2','F'):out[name]=int(out[name])
    out.update(W=W,WW=WW,S=S,Sstar=SS,xi=int(deg@D.max(axis=1)),zeta=int(D.max(axis=1).sum()),
               H=str(sum((Fraction(int(v),k) for k,v in enumerate(histogram) if k),Fraction(0))),
               AD=str(Fraction(2*W,len(D)*(len(D)-1))),
               log10_pi=math.fsum(int(v)*math.log10(k) for k,v in enumerate(histogram) if k))
    counts=[]
    for u,v in edges:
        nu=int(np.count_nonzero(D[u]<D[v])); nv=int(np.count_nonzero(D[v]<D[u]))
        counts.append((nu,nv))
    out.update(Sz=sum(a*b for a,b in counts),PI=sum(a+b for a,b in counts),Mo=sum(abs(a-b) for a,b in counts))
    return out,histogram,counts

def st_histogram(n):
    R=2*n-1; hist=[0]*(2*R+1)
    for d in range(1,2*R+1):
        for j in range(d):
            a=d-j
            for x in range(-R,R-a+1):
                l,u=bounds(n,x); l2,u2=bounds(n,x+a)
                hist[d]+=3*max(0,min(u,u2-j)-max(l,l2-j)+1)
    return np.asarray(hist)

def st_edge_counts(n,x,y):
    R=2*n-1; U=V=0
    for r in range(-R,x+1):
        l,u=bounds(n,r)
        U+=max(0,min(u,x+y-r)-l+1)
    for r in range(x+1,R+1):
        l,u=bounds(n,r)
        V+=max(0,u-max(l,x+y+1-r)+1)
    return U,V

def run(max_n):
    certs=symbolic_check(); data=[]; verifications=[]; distances=[]; comparisons=[]
    pair_count=0; edge_count=0
    for n in range(1,max_n+1):
        both={}
        for st in (False,True):
            family='St' if st else 'Cor'
            V,E,A,C=construct(n,st); D=bfs(A); deg=np.array(list(map(len,A)),dtype=np.int64)
            N=9*n*n-3*n+1 if st else 6*n*n
            M=27*n*n-21*n+6 if st else 9*n*n-3*n
            assert (len(V),len(E))==(N,M)
            actual_part=Counter(tuple(sorted((int(deg[u]),int(deg[v])))) for u,v in E)
            assert actual_part==expected_partition(n,st)
            CD=coordinate_distance(V,st)
            assert np.array_equal(D,CD),(n,st,'coordinate metric')
            pair_count+=N*(N-1)//2
            out,hist,counts=indices(D,E,deg)
            # Check degree formulas independently against vertex/edge sums.
            closed=degree_indices(expected_partition(n,st))
            direct=dict(M1=int((deg**2).sum()),F=int((deg**3).sum()),
                        M2=sum(int(deg[u]*deg[v]) for u,v in E),
                        R=math.fsum(1/math.sqrt(int(deg[u]*deg[v])) for u,v in E),
                        SC=math.fsum(1/math.sqrt(int(deg[u]+deg[v])) for u,v in E),
                        GA=math.fsum(2*math.sqrt(int(deg[u]*deg[v]))/int(deg[u]+deg[v]) for u,v in E),
                        ABC=math.fsum(math.sqrt(int(deg[u]+deg[v]-2)/int(deg[u]*deg[v])) for u,v in E))
            for name,value in closed.items():
                error=abs(value-direct[name])
                assert math.isclose(value,direct[name],rel_tol=1e-12,abs_tol=1e-10)
                verifications.append(dict(n=n,family=family,index=name,formula=str(value),bfs=str(direct[name]),absolute_error=str(error)))
            alternative,_,_=indices(CD,E,deg)
            assert out==alternative
            expected=formulas(st)
            for name,f in expected.items():
                target=f.subs(symbol_n,n)
                assert target.is_integer and int(target)==out[name],(n,st,name,target,out[name])
                verifications.append(dict(n=n,family=family,index=name,formula=str(target),bfs=str(out[name]),absolute_error='0'))
            if st:
                polygon={(x,y) for x in range(-2*n+1,2*n) for y in range(bounds(n,x)[0],bounds(n,x)[1]+1)}
                assert polygon==set(V)
                sh=st_histogram(n)
                assert np.array_equal(hist,sh)
                for (u,v),(nu,nv) in zip(E,counts):
                    x,y=V[u]; xx,yy=V[v]
                    if (xx-x,yy-y)==(1,0):
                        assert st_edge_counts(n,x,y)==(nu,nv)
                        edge_count+=1
                R=2*n-1; B=3*n-1
                rho=np.array([max(abs(x),abs(y),abs(x+y)) for x,y in V])
                assert np.array_equal(D.max(axis=1),R+rho)
                for r in range(R+1):
                    cap=max(0,2*r-B)
                    assert int(np.count_nonzero(rho<=r))==1+3*r*(r+1)-3*cap*(cap+1)
                for x in range(R+1):
                    mx=(B+int((x-B)%2==0)) if x<n else 2*B-3*x+1
                    wx=(6*B-2+2*int((x-B)%2==0)) if x<n else (6 if x==R else 6*mx-8)
                    assert sum(v[0]==x for v in V)==mx
                    assert sum(len(a) for v,a in zip(V,A) if v[0]==x)==wx
                old=[i for i,v in enumerate(V) if v not in C]
                both['old_st_D']=D[np.ix_(old,old)]
            else:
                both['old_D']=D
                z3=np.array([((2*x+y)//3,(x+2*y)//3,(x-y)//3) for x,y in V])
                a,b,c=z3.T; ss=a-b-c
                eccentric=2*n+np.maximum.reduce([2*a-ss+1,-2*b-ss-1,-2*c-ss-1,
                                                 ss-2*a-2,ss+2*b,ss+2*c])
                assert np.array_equal(D.max(axis=1),eccentric)
                # Recover all cut profiles from the integer embedding.
                z=np.array([(2*x+y)//3 for x,y in V])
                for j in range(1,n+1):
                    left=z<=-n+j-1
                    assert int(left.sum())==j*(2*n+j)
                    assert int(deg[left].sum())==3*j*(2*n+j)-n-2*j
                    assert sum(bool(left[u])!=bool(left[v]) for u,v in E)==n+j
            data.append(dict(n=n,family=family,vertices=N,edges=M,diameter=int(D.max()),**out))
            distances.extend(dict(n=n,family=family,distance=k,pairs=int(v)) for k,v in enumerate(hist) if k)
            print('CHECKED',family,n,'vertices',N,'W',out['W'],flush=True)
        original=both['old_D']; st_original=both['old_st_D']
        assert np.all(st_original<=original)
        tri=np.triu_indices(len(original),1)
        diff=(original-st_original)[tri]
        comparisons.append(dict(n=n,old_pairs=len(diff),shortened_pairs=int(np.count_nonzero(diff)),
                                maximum_shortening=int(diff.max()),
                                W_old_in_St=int(st_original[tri].sum())))
    for name,rows in [('indices.csv',data),('verification.csv',verifications),('distance_partitions.csv',distances),('old_vertex_comparison.csv',comparisons)]:
        with (ROOT/name).open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    summary=dict(max_n=max_n,graphs=2*max_n,indices_per_graph=19,
                 coordinate_pair_checks=pair_count,horizontal_edge_partition_checks=edge_count,
                 closed_formula_checks=len(verifications),symbolic_certificates=len(certs),
                 max_absolute_degree_error=max(float(v['absolute_error']) for v in verifications),
                 python=platform.python_version(),numpy=np.__version__,sympy=sympy.__version__,status='PASS')
    (ROOT/'verification_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))
    return data

if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('--max-n',type=int,default=15)
    args=p.parse_args()
    if args.max_n<1: p.error('--max-n must be positive')
    run(args.max_n)
