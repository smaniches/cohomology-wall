# Pin the TRUE source/target dims of the section map for D=(-a1,-a2,-a3,2).
# Source = H^3(A, O(D-N)) ; Target = H^3(A, O(D)).  (q=3 block)
# By Kunneth on (P^1)^4: choose H^1 on the 3 negative factors, H^0 on factor 4.
#   D = (-a1,-a2,-a3, 2):   factor4 charge = 2  -> H^0(P^1,O(2)) dim 3
#       neg factors -ai     -> H^1(P^1,O(-ai)) dim ai-1
#   => H^3(A,O(D)) dim = (a1-1)(a2-1)(a3-1) * 3
#   D-N = (-a1-2,-a2-2,-a3-2, 0): factor4 charge 0 -> H^0(P^1,O(0)) dim 1
#       neg factors -ai-2   -> H^1(P^1,O(-ai-2)) dim ai+1
#   => H^3(A,O(D-N)) dim = (a1+1)(a2+1)(a3+1) * 1
def dims(a):
    src=(a[0]+1)*(a[1]+1)*(a[2]+1)*1
    tgt=(a[0]-1)*(a[1]-1)*(a[2]-1)*3
    return src,tgt
def h0_3(k): return 0 if any(t<0 for t in k) else (k[0]+1)*(k[1]+1)*(k[2]+1)
def h1_3(k):
    tot=0
    for j in range(3):
        if k[j]<=-2 and all(k[i]>=0 for i in range(3) if i!=j):
            d=-k[j]-1
            for i in range(3):
                if i!=j: d*=(k[i]+1)
            tot+=d
    return tot
def rank(a):
    tgt=3*(a[0]-1)*(a[1]-1)*(a[2]-1)
    am4=tuple(t-4 for t in a); am6=tuple(t-6 for t in a)
    return tgt-(3*h0_3(am4)-h0_3(am6)+h1_3(am6))

print("CORRECTED dims and defect:")
print(f"{'a':12s} {'src':>6s} {'tgt':>6s} {'rank':>6s} {'min-rank=DEFECT':>16s}")
meas={(4,4,4):3,(4,5,6):18,(4,6,6):28,(4,6,7):38,(5,5,5):24,(5,5,6):36,
      (6,6,6):48,(7,7,7):48,(4,7,8):48,(5,6,8):48}
allok=True
for a,md in meas.items():
    s,t=dims(a); r=rank(a); d=min(s,t)-r
    ok=(d==md); allok&=ok
    print(f"{str(a):12s} {s:6d} {t:6d} {r:6d} {d:16d}  prev_meas={md} {'OK' if ok else 'STILL OFF'}")
print("\nDefect with CORRECT min(src,tgt):", "matches original measurements" if allok else "mismatch")
print("\nNote: src=(a+1)^3, tgt=3(a-1)^3. For a_i>=4, tgt vs src:")
for a in [(4,4,4),(5,5,5),(6,6,6),(7,7,7),(10,10,10)]:
    s,t=dims(a); print(f"  a={a}: src={s} tgt={t} min={'src' if s<t else 'tgt'}")
