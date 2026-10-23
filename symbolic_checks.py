"""Symbolic certificates for the polynomial identities derived in the paper."""
import sympy as s

n, j, h, t, k = s.symbols('n j h t k', integer=True, positive=True)

def formulas(st):
    Q=s.Rational
    if not st:
        return dict(W=Q(164,5)*n**5-6*n**3+n/5,
            S=Q(984,5)*n**5-82*n**4-18*n**3+10*n**2+Q(6,5)*n,
            Sstar=Q(1476,5)*n**5-246*n**4+46*n**3+15*n**2-Q(11,5)*n,
            Sz=54*n**6-Q(3,2)*n**4+Q(3,2)*n**2,
            PI=54*n**4-18*n**3, Mo=27*n**4-18*n**3-9*n**2,
            xi=60*n**3-33*n**2+9*n,
            zeta=20*n**3-3*n**2+n)
    return dict(W=Q(639,10)*n**5-Q(213,4)*n**4+28*n**3-Q(39,4)*n**2+Q(11,10)*n,
        S=Q(3834,5)*n**5-1065*n**4+730*n**3-294*n**2+Q(331,5)*n-6,
        Sstar=Q(11502,5)*n**5-4473*n**4+4012*n**3-1989*n**2+Q(2668,5)*n-60,
        xi=171*n**3-Q(423,2)*n**2+102*n-Q(75,4)+Q(3,4)*(-1)**n,
        zeta=Q(57,2)*n**3-Q(81,4)*n**2+Q(13,2)*n-Q(11,8)+Q(3,8)*(-1)**n)

def check():
    N=6*n*n; T=18*n*n-6*n; L=j*(2*n+j); D=3*L-n-2*j; c=n+j
    terms=dict(W=L*(N-L),S=D*(N-L)+(T-D)*L,Sstar=D*(T-D),
               Sz=c*L*(N-L),PI=c*N,Mo=c*(N-2*L))
    checks=[]
    for name,f in terms.items():
        value=s.expand(3*(2*s.summation(f,(j,1,n-1))+f.subs(j,n)))
        assert s.simplify(value-formulas(False)[name])==0, name
        checks.append('Cor '+name)
    z=s.summation(6*j*(2*(n+j)-1)+6*(j-1)*(2*(n+j)-2),(j,1,n))
    assert s.expand(z-formulas(False)['zeta'])==0
    assert s.expand(3*z-6*n*(4*n-1)-formulas(False)['xi'])==0
    checks.extend(['Cor zeta','Cor xi'])
    for parity in (0,1):
        nh=2*h+parity; B=3*nh-1; R=2*nh-1; N=9*nh**2-3*nh+1; T=54*nh**2-42*nh+12
        p=1-parity
        m0=B+1-p; w0=6*B-2*p
        def prefix(z, parity_z):
            floor=(z+p-((parity_z+p)%2))/2
            return (N+m0)/2+B*z+floor, (T+w0)/2+(6*B-2)*z+2*floor
        ln,dn=prefix(nh-1,(parity-1)%2)
        lo=ln+s.summation(2*B+1-3*j,(j,nh,t))
        do=dn+s.summation(6*(2*B+1-3*j)-8,(j,nh,t))
        def integrands(l,d):
            return dict(W=l*(N-l),S=d*(N-l)+(T-d)*l,Sstar=d*(T-d))
        vals={name:s.summation(f,(t,nh,2*nh-2)) for name,f in integrands(lo,do).items()}
        for tp in (0,1):
            maxk=h-1+(1 if parity==1 and tp==0 else 0)
            l,d=prefix(2*k+tp,tp)
            for name,f in integrands(l,d).items():
                vals[name]+=s.summation(f,(k,0,maxk))
        for name,v in vals.items():
            assert s.simplify(3*v-formulas(True)[name].subs(n,nh))==0,(parity,name)
            checks.append(f'St {name} parity {parity}')
        # The caps start at floor(B/2)+1. Work separately for each parity of n.
        r0=3*h+(1 if parity else -1)
        cap=2*t-B
        total_c=s.summation(1+3*t*(t+1),(t,0,R-1))-3*s.summation(cap*(cap+1),(t,r0+1,R-1))
        zeta=s.expand(2*R*N-total_c)
        midA=3*h+(1 if parity else -1)
        midB=3*h+(0 if parity else -1)
        A=s.summation(R+B-j,(j,nh,midA))+s.summation(R+j,(j,midA+1,2*nh-1))
        BB=s.summation(R+B-1-j,(j,nh,midB))+s.summation(R+j,(j,midB+1,2*nh-2))
        xi=6*zeta-18*A-6*BB
        for name,v in [('zeta',zeta),('xi',xi)]:
            expected=s.simplify(formulas(True)[name].subs(n,nh))
            assert s.simplify(v-expected)==0,(parity,name,s.simplify(v-expected))
            checks.append(f'St {name} parity {parity}')
    # Explicit expanded degree formulas, separately from the partition form.
    for st,part in [(False,{(2,2):6,(2,3):12*(n-1),(3,3):9*n*n-15*n+6}),
                    (True,{(3,3):6,(3,5):12*(n-1),(3,6):6*n,(5,6):18*(n-1),(6,6):27*n*n-57*n+30})]:
        expected=([54*n*n-30*n,81*n*n-63*n+6,162*n*n-114*n] if not st else
                  [324*n*n-336*n+102,972*n*n-1224*n+414,1944*n*n-2328*n+762])
        for label,f,target in zip(['M1','M2','F'],[lambda a,b:a+b,lambda a,b:a*b,lambda a,b:a*a+b*b],expected):
            assert s.expand(sum(v*f(*pair) for pair,v in part.items())-target)==0
            checks.append(('St ' if st else 'Cor ')+label)
    print('SYMBOLIC CHECKS PASSED:',len(checks))
    return checks

if __name__=='__main__':
    check()
