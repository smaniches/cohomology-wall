#!/usr/bin/env python3
"""Independent verification of the rank/defect formula (Theorem 3.2 of the note).

Engine: builds the actual Koszul section-map matrix in a sorted (colex-style)
basis with a fixed fresh seed, computes its rank by exact Gaussian elimination
over F_p (p = 2^31 - 1), and compares against the closed form, which has NO
fitted parameters.

Run:    python3 scripts/verify_rank_formula.py        (~1-2 min)
Output: comparison table on stdout; results/verification_table.csv;
        exit code 0 iff all 22 points match.
"""
import csv, os, random, sys, itertools
import numpy as np

random.seed(20260612)
P = 2**31 - 1

def atom_H0(k): return [(i, k - i) for i in range(k + 1)] if k >= 0 else []
def atom_H1(k):
    m = -k - 2
    return [(i, m - i) for i in range(m + 1)] if m >= 0 else []

def ambient(D):
    pieces = []
    for k in D:
        pieces.append([('H0', m) for m in atom_H0(k)] + [('H1', m) for m in atom_H1(k)])
    out = {}
    for combo in itertools.product(*pieces):
        q = sum(1 for t, _ in combo if t == 'H1')
        out.setdefault(q, []).append(combo)
    return out

DEG2 = [(2, 0), (1, 1), (0, 2)]
DEG2_4 = list(itertools.product(DEG2, repeat=4))

def mult_factor(typ, mon, dm, Dtgt):
    a, b = mon; da, db = dm
    if typ == 'H0':
        na, nb = a + da, b + db
        if Dtgt >= 0 and na + nb == Dtgt: return ('H0', (na, nb))
        return None
    na, nb = a - da, b - db
    if Dtgt <= -2 and na >= 0 and nb >= 0 and na + nb == (-Dtgt - 2):
        return ('H1', (na, nb))
    return None

def _rank_flint(M, p=P):
    """Exact rank over F_p via python-flint (compiled C; exact and BLAS-independent, fast on tested hardware)."""
    import flint
    R, C = M.shape
    A = flint.nmod_mat(R, C, p)
    Mp = (M % p).astype(object)
    for i in range(R):
        row = Mp[i]
        for j in range(C):
            v = int(row[j])
            if v:
                A[i, j] = v
    return A.rank()

def _rank_numpy(M, p=P):
    """Exact rank over F_p, pure numpy. Fallback when python-flint is absent.

    Reduces only the rows below the pivot (a contiguous block), in one vectorised
    op, instead of rebuilding the whole matrix each step. No per-step full-matrix
    allocation, so runtime does not blow up on slower numpy/BLAS builds.
    """
    M = (M % p).astype(np.int64)
    R, C = M.shape
    r = 0
    for c in range(C):
        nz = np.nonzero(M[r:R, c] % p)[0]
        if nz.size == 0:
            continue
        piv = r + int(nz[0])
        if piv != r:
            M[[r, piv]] = M[[piv, r]]
        inv = pow(int(M[r, c]), p - 2, p)
        M[r] = (M[r] * inv) % p
        below = M[r + 1:R]
        if below.shape[0]:
            factors = below[:, c].copy()
            nzb = np.nonzero(factors)[0]
            if nzb.size:
                below[nzb] = (below[nzb] - np.outer(factors[nzb], M[r])) % p
        r += 1
        if r == R:
            break
    return r

def rank_modp(M, p=P):
    """Exact rank over F_p: use flint if available (fast, exact, environment-agnostic),
    otherwise fall back to the pure-numpy elimination. Both are exact and agree."""
    try:
        import flint  # noqa: F401
        return _rank_flint(M, p)
    except Exception:
        return _rank_numpy(M, p)

def section_rank(D, trials=2):
    """Rank of the q=3 Koszul connecting map for D=(-a1,-a2,-a3,2), generic s."""
    Dm = tuple(d - 2 for d in D)
    q = sum(1 for d in D if d <= -2)
    srcA, tgtA = ambient(Dm), ambient(D)
    if q not in srcA or q not in tgtA: return None
    src = sorted(srcA[q]); tgt = sorted(tgtA[q])
    tidx = {t: i for i, t in enumerate(tgt)}
    best = 0
    for _ in range(trials):
        coeff = {tm: random.randint(1, P - 1) for tm in DEG2_4}
        M = np.zeros((len(tgt), len(src)), dtype=np.int64)
        for j, scol in enumerate(src):
            for tm in DEG2_4:
                new = []
                for f in range(4):
                    r = mult_factor(scol[f][0], scol[f][1], tm[f], D[f])
                    if r is None: new = None; break
                    new.append(r)
                if new is not None:
                    key = tuple(new)
                    if key in tidx:
                        M[tidx[key], j] = (M[tidx[key], j] + coeff[tm]) % P
        best = max(best, rank_modp(M))
    return len(src), len(tgt), best

# ---- Closed form (Theorem 3.2); all cohomology on (P^1)^3 via Kunneth ----
def h0_3(k): return 0 if any(t < 0 for t in k) else (k[0]+1)*(k[1]+1)*(k[2]+1)
def h1_3(k):
    tot = 0
    for j in range(3):
        if k[j] <= -2 and all(k[i] >= 0 for i in range(3) if i != j):
            d = -k[j] - 1
            for i in range(3):
                if i != j: d *= (k[i] + 1)
            tot += d
    return tot

def predict(a):
    src = (a[0]+1)*(a[1]+1)*(a[2]+1)            # h^3(A, O(D-N))
    tgt = 3*(a[0]-1)*(a[1]-1)*(a[2]-1)          # h^3(A, O(D))
    am4 = tuple(t-4 for t in a); am6 = tuple(t-6 for t in a)
    rank = tgt - (3*h0_3(am4) - h0_3(am6) + h1_3(am6))
    return src, tgt, rank, min(src, tgt) - rank

POINTS = [(4,4,4),(4,4,5),(4,4,6),(4,5,5),(4,5,6),(4,6,6),(4,6,7),(4,6,8),(4,7,7),
          (4,7,8),(5,5,5),(5,5,6),(5,6,6),(5,6,7),(5,6,8),(5,7,7),(5,7,8),(6,6,6),
          (6,6,7),(6,6,8),(6,7,7),(7,7,7)]

def main():
    here = os.path.dirname(os.path.abspath(__file__))
    out_csv = os.path.join(here, '..', 'results', 'verification_table.csv')
    rows, allok = [], True
    print(f"{'a':12s} {'src':>5s} {'tgt':>5s} {'meas.rank':>9s} {'pred.rank':>9s} "
          f"{'meas.def':>8s} {'pred.def':>8s}  match")
    for a in POINTS:
        D = tuple(-t for t in a) + (2,)
        ms, mt, mr = section_rank(D)
        ps, pt, pr, pd = predict(a)
        md = min(ms, mt) - mr
        ok = (mr == pr and md == pd and ms == ps and mt == pt)
        allok &= ok
        print(f"{str(a):12s} {ms:5d} {mt:5d} {mr:9d} {pr:9d} {md:8d} {pd:8d}  "
              f"{'YES' if ok else 'NO'}")
        rows.append([a[0], a[1], a[2], ms, mt, mr, pr, md, pd, ok])
    os.makedirs(os.path.dirname(out_csv), exist_ok=True)
    with open(out_csv, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['a1','a2','a3','src','tgt','meas_rank','pred_rank',
                    'meas_defect','pred_defect','match'])
        w.writerows(rows)
    print(f"\n{'ALL 22/22 MATCH' if allok else 'MISMATCH FOUND'}; "
          f"table written to results/verification_table.csv")
    sys.exit(0 if allok else 1)

if __name__ == '__main__':
    main()
