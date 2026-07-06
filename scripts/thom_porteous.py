#!/usr/bin/env python3
"""Thom-Porteous ceiling prediction (Prop. 5.1 / Conj. 5.2 of the note).

Computes deg D_{b-2}(phi) = 48*[(b-1)^2 + C(b-1,3)] two ways:
(1) the closed bracket; (2) the 3x3 Giambelli determinant evaluated symbolically
in the exterior algebra of (P^1)^3 (h_i^2 = 0). Prints both predictions vs the
factorial pattern (b+1)!*8 and marks the b=4 divergence (480 vs 960).
Runtime: instant. Exit 0 iff internal agreement for b=2..11.
"""
import sys
from math import comb, factorial

def mul(p, q):
    out = {}
    for m1, c1 in p.items():
        for m2, c2 in q.items():
            m = tuple(m1[i] + m2[i] for i in range(3))
            if any(e > 1 for e in m): continue
            out[m] = out.get(m, 0) + c1 * c2
    return out

ONE = {(0, 0, 0): 1}
L = {(1, 0, 0): 2, (0, 1, 0): 2, (0, 0, 1): 2}   # c1(O(2,2,2))

def Lpow(k):
    r = ONE
    for _ in range(k): r = mul(r, L)
    return r

def tp_symbolic(b):
    n = b - 1
    c1 = {k: n * v for k, v in L.items()}
    c2 = {k: comb(n, 2) * v for k, v in Lpow(2).items()}
    c3 = {k: comb(n, 3) * v for k, v in Lpow(3).items()}
    det = mul(c1, mul(c1, c1))
    for k, v in mul(c1, c2).items(): det[k] = det.get(k, 0) - 2 * v
    for k, v in c3.items(): det[k] = det.get(k, 0) + v
    return det.get((1, 1, 1), 0)   # integral of top class

def tp_bracket(b):
    n = b - 1
    return 48 * (n * n + comb(n, 3))

ok = True
print(f"{'b':>3s} {'TP (symbolic)':>14s} {'TP (bracket)':>13s} {'factorial (b+1)!*8':>19s}")
for b in range(2, 12):
    s, br, fa = tp_symbolic(b), tp_bracket(b), factorial(b + 1) * 8
    ok &= (s == br)
    tag = '   <- measured' if b in (2, 3) else ('   <- DECISIVE' if b == 4 else '')
    print(f"{b:3d} {s:14d} {br:13d} {fa:19d}{tag}")
print("\nInternal agreement (symbolic == bracket):", "PASS" if ok else "FAIL")
print("b=2,3 retrodicted by BOTH laws; first divergence at b=4: 480 (TP) vs 960 (factorial).")
sys.exit(0 if ok else 1)
