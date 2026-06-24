#!/usr/bin/env python3
"""Decisive b=4 run with the compiled flint backend (reproduces defect = 480).

Requires:  pip install python-flint
Memory-frugal: the rank routine fills a pre-zeroed flint.nmod_mat directly from
the sparse COO triples, with NO dense materialization and NO multi-million-entry
Python list, so the largest matrix (depth 13, 8640 x 8232) runs in well under
1 GB and a couple of minutes per prime on a single core.

The backend is validated against the proven b=2 closed form before the decisive
sweep. Outcome of record: defect = 480 at depths 12 and 13, three-prime
consistent (see results/b4_decisive_result.json). 480 = 48*C(5,3) is the
Thom-Porteous degree; the factorial candidate (b+1)!*2^3 = 960 disagrees with the computed finite-field value.

Usage:
    python3 scripts/run_decisive_flint.py                 # full: depths 11,12,13 x 3 primes
    python3 scripts/run_decisive_flint.py --depths 12,13  # subset of depths
    python3 scripts/run_decisive_flint.py --primes 100003 # subset of primes
"""
import argparse, sys, time
import numpy as np
try:
    import flint
except ImportError:
    flint = None
sys.path.insert(0, "scripts")
from b4_engine import build_sparse, predict_b2, interior_b4

PRIMES_DEFAULT = [100003, 100019, 100043]
DEPTHS_DEFAULT = [11, 12, 13]


def rank_flint(M_csr, p):
    """Exact GF(p) rank via flint.nmod_mat, built from COO without densifying."""
    if flint is None:
        raise SystemExit("python-flint is required for this script. Install with: pip install python-flint")
    nr, nc = M_csr.shape
    M = flint.nmod_mat(nr, nc, p)                 # pre-zeroed C buffer
    coo = M_csr.tocoo()
    r = coo.row; c = coo.col
    v = (coo.data % p).astype(np.int64)
    for i in range(v.shape[0]):
        M[int(r[i]), int(c[i])] = int(v[i])
    return M.rank()


def validate():
    print("Validate flint backend vs proven b=2 formula:")
    ok = True
    for a in [(4, 4, 4), (4, 5, 6), (5, 5, 5), (6, 6, 6), (7, 7, 7), (8, 9, 10)]:
        D = tuple(-x for x in a) + (2,)
        M, s, t = build_sparse(D, PRIMES_DEFAULT[0], 1)
        rk = rank_flint(M, PRIMES_DEFAULT[0]); f = predict_b2(a); ok &= (rk == f)
        print(f"  a={a}: {t}x{s} rank={rk} formula={f} {'OK' if rk == f else 'FAIL'}")
    print("VALIDATION", "PASS" if ok else "FAIL")
    assert ok, "flint backend failed validation -- abort"
    return ok


def sweep(depths, primes):
    print("\n=== DECISIVE b=4 SWEEP (flint) ===")
    print("Interior law 5*(a-7)^3 overshoots both ceilings at depth >= 13 (1080),")
    print("so the measured defect is expected to plateau at one of the two candidate ceilings: 480 or 960.\n")
    results = {}
    for depth in depths:
        a = (depth, depth, depth); D = tuple(-x for x in a) + (4,)
        defs = []
        t0 = time.time()
        for k, p in enumerate(primes, 1):
            M, s, t = build_sparse(D, p, k)
            defs.append(min(s, t) - rank_flint(M, p))
        dt = time.time() - t0
        consistent = len(set(defs)) == 1
        results[depth] = defs[0] if consistent else None
        tag = "consistent" if consistent else f"PRIME-DISAGREE {defs}"
        print(f"  a={a}: {t}x{s} defect={defs[0]} interior={interior_b4(a)} "
              f"[{tag}, {len(primes)} primes] {dt:.1f}s", flush=True)
    return results


def verdict(results):
    print("\n--- VERDICT ---")
    plateau = results.get(13) or results.get(max(results))
    if plateau == 480:
        print("  PLATEAU = 480 -> finite-field computation matches the Thom-Porteous candidate and disagrees with the factorial value 960.")
    elif plateau == 960:
        print("  PLATEAU = 960 -> finite-field computation matches the factorial candidate; Thom-Porteous framing does not apply here.")
    elif plateau is not None:
        print(f"  PLATEAU = {plateau} -> neither law; compute b=5 and re-model.")
    else:
        print("  prime disagreement at decisive depth -- investigate.")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--depths", default=",".join(map(str, DEPTHS_DEFAULT)))
    ap.add_argument("--primes", default=",".join(map(str, PRIMES_DEFAULT)))
    args = ap.parse_args()
    depths = [int(x) for x in args.depths.split(",")]
    primes = [int(x) for x in args.primes.split(",")]
    validate()
    verdict(sweep(depths, primes))
