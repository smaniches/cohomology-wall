#!/usr/bin/env python3
"""Builds sparse exact finite-field Koszul section-map matrices and validates
the b=4 construction against the proven b=2 closed-form formula.

The archived b=4 finite-field ranks are reproduced by scripts/run_decisive_flint.py
using python-flint.

Pipeline:
  build_sparse(D, p, seed)  -> scipy CSR matrix of the Koszul section map mod p
  rank_gfp(M, p)            -> exact rank over GF(p) via python-flint when
                               available, otherwise sparse Gaussian elimination

Ceiling candidates for b=4:
    Thom-Porteous : 480
    factorial     : 960
The interior defect law (b+1)*prod(a_i-7) overshoots both ceilings at depth >= 13
(1080 at depth 13, 1715 at depth 14), so the measured defect is expected to plateau at one of the two candidate ceilings.
Compute depths 12, 13, 14 and read the plateau:
    plateau 480 -> finite-field computation matches the Thom-Porteous candidate in the b=4 test; the general b statement is proved in Theorem 5.1 for the acyclic Toeplitz range
    plateau 960 -> finite-field computation matches the factorial candidate (TP framing does not apply here;
                   investigate expected dimension breaks for this phi-family)
    neither      -> neither candidate matched; compute b=5 and re-model.

Sub-ceiling and ceiling ladder at equal charges:
    b=4 d=11: defect 296  -- sub-ceiling (two primes)
    b=4 d=12: defect 480  -- ceiling onset tau(4)=12 (two primes)
    b=4 d=13: defect 480  -- plateau, three primes

The archived b=4 ranks at depths 11-13 are reproduced by scripts/run_decisive_flint.py
using python-flint. The --validate flag here checks the b=4 matrix construction
against the proven b=2 closed-form formula; no python-flint is required for that.

Usage:
    python3 b4_engine.py --validate          # check b=4 construction vs b=2 formula
    python3 b4_engine.py --decisive          # run the b=4 depth sweep (needs python-flint)
    python3 b4_engine.py --depth 13          # single depth, two primes
"""
import argparse, itertools, sys, time
import numpy as np


def build_sparse(D, p, seed):
    """Koszul section map for D=(-a1,-a2,-a3,b), generic s in H^0(2,2,2,2), mod p."""
    from scipy.sparse import csr_matrix
    a = [-D[0], -D[1], -D[2]]; b = D[3]
    src_sizes = [a[0]+1, a[1]+1, a[2]+1, b-1]   # H^1 dims a_i+1; H^0(b-2) dim b-1
    tgt_sizes = [a[0]-1, a[1]-1, a[2]-1, b+1]   # H^1 dims a_i-1; H^0(b)  dim b+1
    src = np.array(list(itertools.product(*[np.arange(s) for s in src_sizes])),
                   dtype=np.int64)
    ntgt = int(np.prod(tgt_sizes))
    st = np.array([int(np.prod(tgt_sizes[k+1:])) for k in range(4)], dtype=np.int64)
    rng = np.random.default_rng(seed)
    Nmons = list(itertools.product(range(3), repeat=4))
    coeff = {d: int(rng.integers(1, p)) for d in Nmons}
    ci = np.arange(len(src))
    rows, cols, vals = [], [], []
    for d in Nmons:
        e0 = src[:, 0]-d[0]; e1 = src[:, 1]-d[1]; e2 = src[:, 2]-d[2]; f = src[:, 3]+d[3]
        ok = ((e0 >= 0) & (e0 < tgt_sizes[0]) & (e1 >= 0) & (e1 < tgt_sizes[1]) &
              (e2 >= 0) & (e2 < tgt_sizes[2]) & (f >= 0) & (f < tgt_sizes[3]))
        if not ok.any():
            continue
        ti = e0[ok]*st[0] + e1[ok]*st[1] + e2[ok]*st[2] + f[ok]*st[3]
        rows.append(ti); cols.append(ci[ok])
        vals.append(np.full(int(ok.sum()), coeff[d], dtype=np.int64))
    M = csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                   shape=(ntgt, len(src)))
    M.data %= p
    return M, len(src), ntgt


def rank_sparse(M_csr, p):
    """Pure-Python sparse Gaussian elimination over GF(p). Pivot by sparsity."""
    A = M_csr.tolil()
    rows = {i: dict(zip(A.rows[i], [int(x) % p for x in A.data[i]]))
            for i in range(A.shape[0]) if A.rows[i]}
    pivots, rank = {}, 0
    for i in sorted(rows.keys(), key=lambda k: len(rows[k])):
        row = rows[i]
        while row:
            c = min(row)
            if c in pivots:
                prow = rows[pivots[c]]
                fac = row[c]*pow(prow[c], p-2, p) % p
                for cc, v in prow.items():
                    nv = (row.get(cc, 0)-fac*v) % p
                    if nv: row[cc] = nv
                    elif cc in row: del row[cc]
            else:
                break
        if row:
            c = min(row); inv = pow(row[c], p-2, p)
            for cc in list(row): row[cc] = row[cc]*inv % p
            pivots[c] = i; rows[i] = row; rank += 1
    return rank


def rank_gfp(M_csr, p, backend="auto", depth_hint=None):
    """Exact rank over GF(p) using sparse elimination (sufficient for --validate).
    For archived decisive ranks at large depths, use run_decisive_flint.py.
    """
    return rank_sparse(M_csr, p)


def predict_b2(a):
    def h0(k): return 0 if any(t < 0 for t in k) else (k[0]+1)*(k[1]+1)*(k[2]+1)
    def h1(k):
        tot = 0
        for j in range(3):
            if k[j] <= -2 and all(k[i] >= 0 for i in range(3) if i != j):
                d = -k[j]-1
                for i in range(3):
                    if i != j: d *= (k[i]+1)
                tot += d
        return tot
    am4 = (a[0]-4, a[1]-4, a[2]-4); am6 = (a[0]-6, a[1]-6, a[2]-6)
    return 3*(a[0]-1)*(a[1]-1)*(a[2]-1) - (3*h0(am4) - h0(am6) + h1(am6))


def interior_b4(a):
    v = 1
    for x in a: v *= (x-7)
    return 5*v


def validate():
    print("Validation: sparse build + GF(p) rank vs proven b=2 formula")
    ok = True
    for a in [(4,4,4),(4,5,6),(5,5,5),(6,6,6),(7,7,7)]:
        D = tuple(-x for x in a)+(2,)
        M, s, t = build_sparse(D, 100003, 1)
        r = rank_sparse(M, 100003); f = predict_b2(a); ok &= (r == f)
        print(f"  a={a}: {t}x{s} nnz={M.nnz} rank={r} formula={f} "
              f"{'OK' if r == f else 'FAIL'}")
    print("VALIDATION", "PASS" if ok else "FAIL")
    return ok


def one_depth(depth, backend="auto"):
    a = (depth, depth, depth); D = tuple(-x for x in a)+(4,)
    t0 = time.time()
    M, s, t = build_sparse(D, 100003, 1)
    r1 = rank_gfp(M, 100003, backend, depth)
    M2, _, _ = build_sparse(D, 100019, 2)
    r2 = rank_gfp(M2, 100019, backend, depth)
    d1, d2 = min(s, t)-r1, min(s, t)-r2
    dt = time.time()-t0
    tag = "consistent" if d1 == d2 else "PRIME-DISAGREE"
    print(f"  a={a}: {t}x{s} nnz={M.nnz} rank={r1} defect={d1} "
          f"interior={interior_b4(a)} [{tag}] {dt:.1f}s")
    return d1 if d1 == d2 else None


def decisive():
    print("=== DECISIVE b=4 SWEEP (plateau decides 480 vs 960) ===")
    res = {}
    for depth in (12, 13, 14):
        res[depth] = one_depth(depth)
    print("\nVERDICT")
    plateau = res.get(14) or res.get(13)
    if plateau == 480:
        print("  defect plateau = 480 -> finite-field computation matches the Thom-Porteous candidate in the b=4 test; the general b statement is proved in Theorem 5.1 for the acyclic Toeplitz range.")
    elif plateau == 960:
        print("  defect plateau = 960 -> finite-field computation matches the factorial candidate (TP framing does not apply here).")
    elif plateau is not None:
        print(f"  defect plateau = {plateau} -> NEITHER law. Compute b=5, re-model.")
    else:
        print("  prime disagreement -> investigate before trusting.")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--validate", action="store_true")
    ap.add_argument("--decisive", action="store_true")
    ap.add_argument("--depth", type=int, default=None)
    ap.add_argument("--backend", default="auto", choices=["auto", "sparse"])
    args = ap.parse_args()
    if args.validate: validate()
    elif args.depth is not None: one_depth(args.depth, args.backend)
    elif args.decisive: decisive()
    else:
        validate()
