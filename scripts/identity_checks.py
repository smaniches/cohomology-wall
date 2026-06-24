from itertools import product
from math import comb, factorial
# (R1) finite-difference identity: prod(a+1)-3prod(a-1)+3prod(a-3)-prod(a-5) == 48 for ALL a (identity in a)
viol=0
for a in product(range(-3,25),repeat=3):
    f=lambda t:(a[0]+1-2*t)*(a[1]+1-2*t)*(a[2]+1-2*t)
    if f(0)-3*f(1)+3*f(2)-f(3)!=48: viol+=1
print("R1 finite-difference identity, 28^3 lattice incl. negatives:", "PASS" if viol==0 else f"FAIL {viol}")
# (R2) min determination: src<tgt for all a_i>=6  (src=prod(a+1), tgt=3prod(a-1))
viol=[a for a in product(range(6,40),repeat=3) if (a[0]+1)*(a[1]+1)*(a[2]+1) >= 3*(a[0]-1)*(a[1]-1)*(a[2]-1)]
print("R2 min=src for a_i>=6 (34^3 sweep):", "PASS" if not viol else f"FAIL {viol[:3]}")
# (R3) H^1(Y,O(a-4))=0 for a_i>=4: Kunneth criterion = some entry <= -2; check none
viol=[a for a in product(range(4,30),repeat=3) if any(t-4<=-2 for t in a)]
print("R3 H1(a-4)=0 on domain:", "PASS" if not viol else f"FAIL {viol[:3]}")
# (R4) Thom-Porteous determinant: c1^3-2c1c2+c3 with cj=C(b-1,j) ell^j; bracket = (b-1)^2+C(b-1,3); x48
def tp(b):
    n=b-1
    c1,c2,c3=n,comb(n,2),comb(n,3)   # coefficients of ell^j
    det=c1**3-2*c1*c2+c3             # coefficient of ell^3
    return det*48, ((n*n+comb(n,3))*48)
ok=all(tp(b)[0]==tp(b)[1] for b in range(2,12))
print("R4 TP det == 48[(b-1)^2+C(b-1,3)] for b=2..11:", "PASS" if ok else "FAIL")
print("    TP ceilings b=2..6:", [tp(b)[0] for b in range(2,7)], " factorial:", [factorial(b+1)*8 for b in range(2,7)])
# (R5) TP retrodiction: b=2 -> 48, b=3 -> 192; divergence at b=4: 480 vs 960
print("R5 retrodiction:", "PASS" if tp(2)[0]==48 and tp(3)[0]==192 and tp(4)[0]==480 else "FAIL")
# (R6) explicit 3x3 TP determinant symbolically in exterior algebra (h_i^2=0) for b=3: expect 192
# c(F)=(1+L)^2, L=2h1+2h2+2h3; compute det[[c1,c2,c3],[1,c1,c2],[0,1,c1]] coefficient of h1h2h3, integrate=1
from itertools import combinations
def mul(p,q):  # polynomials as dict: frozenset of factors-> but with squares vanish; represent monomial as tuple of exponents (e1,e2,e3) capped
    out={}
    for m1,c1_ in p.items():
        for m2,c2_ in q.items():
            m=tuple(m1[i]+m2[i] for i in range(3))
            if any(e>1 for e in m): continue
            out[m]=out.get(m,0)+c1_*c2_
    return out
one={(0,0,0):1}; L={(1,0,0):2,(0,1,0):2,(0,0,1):2}
def power(p,k):
    r=one
    for _ in range(k): r=mul(r,p)
    return r
b=3; n=b-1
cF=[one,  {k:n*v for k,v in L.items()}, {k: comb(n,2)*v for k,v in power(L,2).items()}, {k: comb(n,3)*v for k,v in power(L,3).items()}]
c1,c2,c3=cF[1],cF[2],cF[3]
det=mul(c1,mul(c1,c1))
for k,v in mul(c1,c2).items(): det[k]=det.get(k,0)-2*v
for k,v in c3.items(): det[k]=det.get(k,0)+v
top=det.get((1,1,1),0)
print("R6 symbolic exterior-algebra TP at b=3:", top, "PASS" if top==192 else "FAIL")
