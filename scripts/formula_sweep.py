# REFEREE 2 (Arithmetician) + REFEREE 1 (Logician edge cases): recompute independently.
# 1. Verify the 48 = 3! * 2^3 intersection number claim via a DIFFERENT method:
#    length of base locus of 3 generic (2,2,2) divisors on (P^1)^3 = the coefficient
#    extraction / mixed volume. (2H1+2H2+2H3)^3 expanded, keep H1 H2 H3 term (top class).
from itertools import product
from math import factorial
# On (P^1)^3, H_i^2=0, top class H1 H2 H3 has integral 1. (2H1+2H2+2H3)^3:
# multinomial: only the term with one H1,one H2,one H3 survives -> coeff = 3! * 2*2*2
coeff = factorial(3) * (2**3)
print("Referee 2 — intersection number recomputation:")
print(f"  (2H1+2H2+2H3)^3 surviving term = 3! * 2^3 = {coeff}   [paper says 48]  {'OK' if coeff==48 else 'FAIL'}")

# 2. Verify b=3 ceiling 192 = 4! * 2^3
print(f"  b=3 ceiling 4!*2^3 = {factorial(4)*2**3}   [computed 192]  {'OK' if factorial(4)*8==192 else 'FAIL'}")

# 3. Logician edge case: does the rank formula stay non-negative and <= min(src,tgt)
#    everywhere in the stated domain a_i>=4? (a formula giving negative rank = fatal)
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
def dims(a): return (a[0]+1)*(a[1]+1)*(a[2]+1), 3*(a[0]-1)*(a[1]-1)*(a[2]-1)
bad=[]
for a in product(range(4,16),repeat=3):
    s,t=dims(a); r=rank(a); d=min(s,t)-r
    if r<0 or d<0 or r>min(s,t): bad.append((a,s,t,r,d))
print(f"\nReferee 1 — edge-case sweep over a_i in [4,15]^3 ({12**3} points):")
print(f"  formula violations (neg rank, neg defect, rank>min): {len(bad)}")
for b in bad[:10]: print("   ",b)

# 4. Logician: boundary a_i=4 exactly (smallest allowed). Defect should be 3 at (4,4,4).
print(f"\n  boundary (4,4,4): defect = {min(*dims((4,4,4)))-rank((4,4,4))}  [expect 3]")
# 5. What about a_i in {4,5} only, never reaching 6 -> ceiling never active, pure cubic region
print(f"  pure-cubic check (5,5,5): defect={min(*dims((5,5,5)))-rank((5,5,5))}  3*prod(a-3)={3*2*2*2}")
print(f"  pure-cubic check (4,5,5): defect={min(*dims((4,5,5)))-rank((4,5,5))}  3*prod(a-3)={3*1*2*2}")

# 6. Does defect ever EXCEED 48 anywhere? (ceiling claim: max defect = 48 for b=2)
maxd=max(min(*dims(a))-rank(a) for a in product(range(4,16),repeat=3))
print(f"\n  Referee 1 ceiling stress: max defect over [4,15]^3 = {maxd}  [ceiling claims 48]  {'OK' if maxd==48 else 'EXCEEDS — investigate'}")
