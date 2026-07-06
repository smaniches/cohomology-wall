#!/usr/bin/env python3
"""
toeplitz_minor_ideal_check.py -- machine verification of the two local claims
of Theorem 5.1 at small b.

Theorem 5.1 (general Toeplitz determinantal saturation) rests on two local
computations in S = k[x,y,z], with n = b-1:

  (i)  the ideal of maximal minors of the n x (n+2) banded Toeplitz matrix
       T_n(x,y,z) equals (x,y,z)^n;
  (ii) the Buchsbaum-Rim cokernel C_n = coker(S(-1)^(n+2) -> S^n) has finite
       length C(n+2,3) = C(b+1,3).

This script verifies both claims exactly for n = 1..8 (b = 2..9), several b
beyond any global finite-field computation in the archive.

  (i)  Every maximal minor is homogeneous of degree n, so I_n(T_n) is
       contained in (x,y,z)^n automatically, and the two ideals are equal iff
       the C(n+2,2) minors span the degree-n graded piece S_n, whose dimension
       is also C(n+2,2) -- that is, iff the minors are linearly independent.
       The minors are expanded by exact integer cofactor expansion and their
       coefficient matrix is checked for full rank.
  (ii) The graded piece of C_n in degree d has dimension n*C(d+2,2) minus the
       rank of the degree-d multiplication matrix of the Toeplitz map; the
       Hilbert series in Theorem 5.1 confines C_n to degrees 0..n-1, and the
       summed dimensions must equal C(n+2,3).

Both checks are deterministic (the matrix entries are literally x, y, z; no
random coefficients) and exact: ranks are taken over F_p with p = 2^31 - 1 via
python-flint, and full rank mod p certifies full rank over Q. Runtime: seconds.

This is a regression companion to the proof, for concrete n; it is not a
substitute for the general-n argument in the paper.

Usage:
    python3 scripts/toeplitz_minor_ideal_check.py     # n = 1..8; exit 0 iff all pass
"""
import sys
from itertools import combinations
from math import comb

try:
    import flint
except ImportError:
    raise SystemExit("python-flint is required for this script. "
                     "Install with: pip install python-flint")

P = 2**31 - 1


def toeplitz_columns(n):
    """Column j of T_n as a dict {row: variable}, variables coded 0=x, 1=y, 2=z."""
    cols = []
    for j in range(n + 2):
        col = {}
        for i in range(n):
            if j - i in (0, 1, 2):
                col[i] = j - i
        cols.append(col)
    return cols


def poly_mul_var(p, v):
    """Multiply a dict-polynomial {(ex,ey,ez): c} by the variable coded v."""
    out = {}
    for m, c in p.items():
        mm = list(m)
        mm[v] += 1
        out[tuple(mm)] = out.get(tuple(mm), 0) + c
    return out


def det_poly(cols, rows):
    """Exact determinant of the square submatrix of T_n given by (cols, rows),
    by cofactor expansion along the first remaining column. Entries are single
    variables or zero, so every intermediate coefficient is an exact integer."""
    if not rows:
        return {(0, 0, 0): 1}
    j = cols[0]
    total = {}
    for k, i in enumerate(rows):
        v = j.get(i)
        if v is None:
            continue
        sub = det_poly(cols[1:], rows[:k] + rows[k + 1:])
        term = poly_mul_var(sub, v)
        sign = 1 if k % 2 == 0 else -1
        for m, c in term.items():
            total[m] = total.get(m, 0) + sign * c
    return {m: c for m, c in total.items() if c}


def monomials(d):
    """All monomials (ex,ey,ez) of total degree d, in a fixed order."""
    return [(ex, ey, d - ex - ey) for ex in range(d + 1) for ey in range(d - ex + 1)]


def rank_fp(rows_of_entries, ncols):
    """Exact rank over F_P of a dense integer matrix given as a list of rows."""
    nr = len(rows_of_entries)
    M = flint.nmod_mat(nr, ncols, P)
    for i, row in enumerate(rows_of_entries):
        for j, c in row.items():
            M[i, j] = c % P
    return M.rank()


def check_minor_ideal(n):
    """Claim (i): the C(n+2,2) maximal minors of T_n are linearly independent
    in S_n, hence I_n(T_n) = (x,y,z)^n."""
    cols = toeplitz_columns(n)
    basis = {m: k for k, m in enumerate(monomials(n))}
    all_rows = list(range(n))
    rows = []
    for p, q in combinations(range(n + 2), 2):
        keep = [cols[j] for j in range(n + 2) if j not in (p, q)]
        d = det_poly(keep, all_rows)
        rows.append({basis[m]: c for m, c in d.items()})
    want = comb(n + 2, 2)
    got = rank_fp(rows, len(basis))
    return got, want


def check_coker_length(n):
    """Claim (ii): sum over d of dim (C_n)_d equals C(n+2,3). The degree-d
    presentation matrix sends the (j, m)-th source generator (column block j,
    monomial m of degree d-1) to x*m, y*m, z*m in the row blocks i with
    j - i = 0, 1, 2."""
    total = 0
    for d in range(n):
        tgt = {m: k for k, m in enumerate(monomials(d))}
        nrows_per_block = len(tgt)
        src_mons = monomials(d - 1) if d >= 1 else []
        cols = []
        for j in range(n + 2):
            for m in src_mons:
                col = {}
                for v in (0, 1, 2):
                    i = j - v
                    if 0 <= i < n:
                        mm = list(m)
                        mm[v] += 1
                        col[i * nrows_per_block + tgt[tuple(mm)]] = 1
                cols.append(col)
        rk = rank_fp(cols, n * nrows_per_block) if cols else 0
        total += n * nrows_per_block - rk
    return total, comb(n + 2, 3)


def main():
    ok = True
    print("Local Toeplitz model of Theorem 5.1: minor ideal and cokernel length")
    print(f"{'n':>2s} {'b':>2s}  {'minors indep.':>14s}  {'I=(x,y,z)^n':>11s}  "
          f"{'coker length':>12s}  {'C(b+1,3)':>8s}  verdict")
    for n in range(1, 9):
        gi, wi = check_minor_ideal(n)
        gl, wl = check_coker_length(n)
        good = (gi == wi) and (gl == wl)
        ok &= good
        print(f"{n:2d} {n+1:2d}  {f'{gi}/{wi}':>14s}  {'yes' if gi == wi else 'NO':>11s}  "
              f"{gl:12d}  {wl:8d}  {'PASS' if good else 'FAIL'}")
    print(f"\nTOEPLITZ LOCAL MODEL CHECK: {'PASS' if ok else 'FAIL'}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
