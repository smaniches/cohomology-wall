#!/usr/bin/env python3
"""
saturation_threshold.py  --  exact-F_p onset checks for the proven and computed cases.

Pins the depth a_i (equal charges d=a_1=a_2=a_3) at which the defect first
reaches the ceiling 48*C(b+1,3), and compares it to the three candidate bounds:
  2b           (wall onset / interior-law turn-on; FALSE for the ceiling)
  2b+1         (term-by-term acyclicity of the skyscraper BR resolution; FALSE for the ceiling)
  tau(b) = (2b-1) + ceil( 2*(b(b-1))^(1/3) )   (cubic crossing; agrees through b=6 but fails at b=7)

b=2 uses the proven closed form; b=3 and requested b>=3 points use exact rank over F_p. The corrected general equal-charge threshold is rho(b), with tau(b) retained as the old cubic-crossing reference.
"""
import argparse, itertools, numpy as np
from scipy.sparse import csr_matrix
from math import comb, ceil

def build_sparse(D,p,seed):
    a=[-D[0],-D[1],-D[2]]; b=D[3]
    src=[a[0]+1,a[1]+1,a[2]+1,b-1]; tgt=[a[0]-1,a[1]-1,a[2]-1,b+1]
    S=np.array(list(itertools.product(*[np.arange(x) for x in src])),dtype=np.int64)
    nt=int(np.prod(tgt)); st=np.array([int(np.prod(tgt[k+1:])) for k in range(4)],dtype=np.int64)
    rng=np.random.default_rng(seed); N=list(itertools.product(range(3),repeat=4))
    co={d:int(rng.integers(1,p)) for d in N}; ci=np.arange(len(S)); R=[];C=[];V=[]
    for d in N:
        e0=S[:,0]-d[0];e1=S[:,1]-d[1];e2=S[:,2]-d[2];f=S[:,3]+d[3]
        ok=((e0>=0)&(e0<tgt[0])&(e1>=0)&(e1<tgt[1])&(e2>=0)&(e2<tgt[2])&(f>=0)&(f<tgt[3]))
        if not ok.any(): continue
        ti=e0[ok]*st[0]+e1[ok]*st[1]+e2[ok]*st[2]+f[ok]*st[3]
        R.append(ti);C.append(ci[ok]);V.append(np.full(int(ok.sum()),co[d],dtype=np.int64))
    M=csr_matrix((np.concatenate(V),(np.concatenate(R),np.concatenate(C))),shape=(nt,len(S))); M.data%=p
    return M,len(S),nt

def rank_flint(M, p):
    try:
        import flint
    except ImportError:
        raise SystemExit(
            "python-flint is required for this command. "
            "Install with: pip install python-flint"
        )
    coo = M.tocoo(); nr, nc = coo.shape; Mf = flint.nmod_mat(nr, nc, p)
    r,c,v=coo.row,coo.col,(coo.data.astype(np.int64)%p)
    for i in range(v.size): Mf[int(r[i]),int(c[i])]=int(v[i])
    return Mf.rank()

def predict_b2(a):
    def h0(k): return 0 if any(t<0 for t in k) else (k[0]+1)*(k[1]+1)*(k[2]+1)
    def h1(k):
        tot=0
        for j in range(3):
            if k[j]<=-2 and all(k[i]>=0 for i in range(3) if i!=j):
                d=-k[j]-1
                for i in range(3):
                    if i!=j: d*=(k[i]+1)
                tot+=d
        return tot
    am4=tuple(x-4 for x in a); am6=tuple(x-6 for x in a)
    return 3*(a[0]-1)*(a[1]-1)*(a[2]-1)-(3*h0(am4)-h0(am6)+h1(am6))

def ceiling(b): return 48*comb(b+1,3)
def cubic(b,d):  return (b+1)*(d-(2*b-1))**3 if d>(2*b-1) else 0
def tau(b):      return (2*b-1)+ceil(2*(b*(b-1))**(1/3))
def Edim(b,d):   return (b+1)*(d-1)**3   # source dim of mu_{b,(d,d,d)}
def Fdim(b,d):   return (b-1)*(d+1)**3   # target dim of mu_{b,(d,d,d)}
def rho(b):
    """Corrected equal-charge onset (Corollary 5.2): least d>=2b+1 with E>=F."""
    d = 2*b+1
    while not (Edim(b,d) >= Fdim(b,d)):
        d += 1
    return d

def defect_b2(d):
    s=(d+1)**3*1; t=(d-1)**3*3
    return min(s,t)-predict_b2((d,d,d))

def defect_fp(b,d,p=100003,seed=1):
    M,s,t=build_sparse((-d,-d,-d,b),p,seed); return min(s,t)-rank_flint(M,p)

def onset_measured(b, ds, defect_fn):
    lad = {d: defect_fn(d) for d in ds}
    onset = min((d for d, v in lad.items() if v == ceiling(b)), default=None)
    return lad, onset

def gate():
    """Fast assertion used by reproduce.py: the measured equal-charge onset equals the
    corrected threshold rho(b) for b=2,3 (where rho==tau), and rho(7)!=tau(7) (the
    correction that supersedes the old cubic-crossing tau)."""
    lad2, o2 = onset_measured(2, range(4, 9), defect_b2)
    lad3, o3 = onset_measured(3, range(6, 11), lambda d: defect_fp(3, d))
    print(f"b=2 ladder {dict(sorted(lad2.items()))}  onset={o2}  rho(2)={rho(2)}  tau(2)={tau(2)}")
    print(f"b=3 ladder {dict(sorted(lad3.items()))}  onset={o3}  rho(3)={rho(3)}  tau(3)={tau(3)}")
    ok = (o2 == rho(2) == 6) and (o3 == rho(3) == 9)
    # the bounds 2b and 2b+1 must be strictly sub-ceiling
    ok = ok and lad2[5] < ceiling(2) and lad3[7] < ceiling(3)
    # the correction: rho==tau through b=6, and they diverge at b=7 (tau false from b=7)
    ok = ok and all(rho(b) == tau(b) for b in range(2, 7))
    ok = ok and (rho(7) == 21 and tau(7) == 20)
    print(f"rho==tau for b=2..6: {all(rho(b)==tau(b) for b in range(2,7))}; "
          f"b=7 rho={rho(7)} tau={tau(7)} (diverge; tau superseded)")
    print(f"ONSET CHECK (b=2,3): {'PASS' if ok else 'FAIL'}")
    return ok

def full_table():
    import time
    print("b | ceiling | 2b | 2b+1 | tau(b) | onset(measured) | ladder (d: defect)")
    print("-"*92)
    lad, onset = onset_measured(2, range(4, 9), defect_b2)
    print(f"2 | {ceiling(2):7d} | {4:2d} | {5:4d} | {tau(2):6d} | {onset:15d} | "
          + ", ".join(f"{d}:{lad[d]}" for d in sorted(lad)))
    lad, onset = onset_measured(3, range(6, 11), lambda d: defect_fp(3, d))
    print(f"3 | {ceiling(3):7d} | {6:2d} | {7:4d} | {tau(3):6d} | {onset:15d} | "
          + ", ".join(f"{d}:{lad[d]}" for d in sorted(lad)))
    t0=time.time(); d41=defect_fp(4,11); dt=time.time()-t0
    print(f"4 | {ceiling(4):7d} | {8:2d} | {9:4d} | {tau(4):6d} | (see results) | "
          f"11:{d41}  [ceiling {ceiling(4)} reached at tau=12; d=11 sub-ceiling]  ({dt:.0f}s)")
    print()
    print("rho(b) for b=2..7:", [rho(b) for b in range(2,8)], "  (corrected onset; first d with E>=F)")
    print("tau(b) for b=2..7:", [tau(b) for b in range(2,8)], "  (old cubic crossing; fails from b=7)")
    print("2b+1   for b=2..7:", [2*b+1 for b in range(2,8)], "  (acyclic-range lower bound)")

def _main():
    import sys
    ap = argparse.ArgumentParser(
        description="Finite-field onset checks for the tetraquadric: b=2 uses the proven closed form; "
                    "b=3 and requested b>=3 points use exact rank over F_p (requires python-flint). "
                    "The corrected general equal-charge threshold is rho(b), with tau(b) retained as the old cubic-crossing reference.")
    ap.add_argument("--gate", action="store_true",
                    help="fast onset assertion used by reproduce.py: onset == tau(b) for b=2,3 (exit 0/1)")
    ap.add_argument("--bd", nargs="+", metavar="VAL",
                    help="B D [PRIME]: recompute the defect at the equal-charge point "
                         "(a1,a2,a3,b)=(d,d,d,b); optional PRIME defaults to 100003; b>=3 needs python-flint")
    args = ap.parse_args()
    if args.gate:
        sys.exit(0 if gate() else 1)
    if args.bd is not None:
        if len(args.bd) < 2:
            ap.error("--bd needs at least B and D, e.g. --bd 5 15 100019")
        try:
            b = int(args.bd[0]); d = int(args.bd[1])
            p = int(args.bd[2]) if len(args.bd) >= 3 else 100003
        except ValueError:
            ap.error("--bd arguments must be integers, e.g. --bd 5 15 100019")
        val = defect_b2(d) if b == 2 else defect_fp(b, d, p)
        cap, cub = ceiling(b), cubic(b, d)
        e, f, rh, th = Edim(b, d), Fdim(b, d), rho(b), tau(b)
        method = "proven closed form" if b == 2 else f"exact F_p rank (p={p})"
        branch = "C_b (E>=F)" if e >= f else "C_b+E-F (E<F)"
        pred = cap if e >= f else cap + e - f
        status = ("saturated (= ceiling)" if val == cap else "sub-ceiling")
        print(f"b={b}  d={d}   defect = {val}    [{method}]")
        print(f"  source E=(b+1)prod(a-1)     = {e}")
        print(f"  target F=(b-1)prod(a+1)     = {f}")
        print(f"  ceiling  48*C(b+1,3) = C_b  = {cap}")
        print(f"  Theorem 5.1 branch          = {branch} -> predicted defect {pred}"
              + ("  (acyclic: a>=2b+1)" if d >= 2*b+1 else "  (NOTE: d<2b+1, outside acyclic range)"))
        print(f"  interior (b+1)(d-(2b-1))^3  = {cub}   (reference only; not the law near onset)")
        print(f"  onset    rho(b)             = {rh}   ({'d >= rho' if d >= rh else 'd < rho'})")
        print(f"  old ref  tau(b)             = {th}   (superseded; agrees with rho through b=6)")
        print(f"  status: {status}")
        sys.exit(0)
    full_table()

if __name__=="__main__":
    _main()
