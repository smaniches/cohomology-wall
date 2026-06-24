#!/usr/bin/env python3
"""
run_b7_decisive.py -- the b=7 corrected-onset witness (tau(7)=20 != rho(7)=21).

This is the decisive computation distinguishing the OLD cubic-crossing threshold
tau(b) from the corrected source-target onset rho(b). At b=7:

    d=20:  E = 8*19^3 = 54872,  F = 6*21^3 = 55566,  E < F  -> defect = C_7 + E - F = 1994
    d=21:  E = 8*20^3 = 64000,  F = 6*22^3 = 63888,  E >= F -> defect = C_7        = 2688

The old tau(7)=20 would predict the ceiling already at d=20; the measured defect
1994 (not 2688) refutes that and confirms rho(7)=21.

MEMORY. The dense F_p rank is taken on the full section-map matrix, of shape
E x F (about 55000 x 55000 at d=20). As a flint.nmod_mat (8 bytes/entry) that is
roughly 24 GB, plus LU workspace. Run the full computation on a host with
>= ~32 GB RAM. On a small machine, use --check, which builds the sparse matrix,
confirms its shape equals (E, F) and reports nnz, and does NOT attempt the rank.

Usage:
    python3 scripts/run_b7_decisive.py --check                 # dims only; runs anywhere
    python3 scripts/run_b7_decisive.py                         # full: d=20,21 x primes 100003,100019  (needs ~32 GB)
    python3 scripts/run_b7_decisive.py --depths 20 --primes 100003
    python3 scripts/run_b7_decisive.py --force                 # run the rank even if the RAM estimate says no
"""
import argparse, os, sys, time, platform
from datetime import datetime, timezone

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from b4_engine import build_sparse
from toeplitz_defect import E as Edim, F as Fdim, C_b, min_EF, defect_acyclic, in_acyclic_range

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
B = 7
PRIMES_DEFAULT = [100003, 100019]
DEPTHS_DEFAULT = [20, 21]


def _avail_bytes():
    try:
        return os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_AVPHYS_PAGES")
    except (ValueError, OSError, AttributeError):
        return None


def _dense_bytes(b, d):
    a = (d, d, d)
    return Edim(b, a) * Fdim(b, a) * 8  # flint.nmod_mat: 8 bytes/entry


def rank_flint(M_csr, p):
    import flint
    coo = M_csr.tocoo()
    nr, nc = coo.shape
    Mf = flint.nmod_mat(nr, nc, p)
    r, c = coo.row, coo.col
    v = (coo.data % p).astype(np.int64)
    for i in range(v.shape[0]):
        Mf[int(r[i]), int(c[i])] = int(v[i])
    return Mf.rank()


def check(depths):
    print("=== b=7 dimension / structure check (no dense rank) ===")
    ok = True
    for d in depths:
        a = (d, d, d)
        D = tuple(-d for _ in range(3)) + (B,)
        M, s, t = build_sparse(D, PRIMES_DEFAULT[0], 1)  # s=nsrc(cols), t=ntgt(rows)
        e, f, c = Edim(B, a), Fdim(B, a), C_b(B)
        pred = defect_acyclic(B, a)
        dims_ok = (t == e) and (s == f)
        ok &= dims_ok
        path = os.path.join(ROOT, "logs", f"b7_d{d}_check.txt")
        os.makedirs(os.path.join(ROOT, "logs"), exist_ok=True)
        with open(path, "w") as fh:
            fh.write(f"command        : python3 scripts/run_b7_decisive.py --check\n")
            fh.write(f"date_utc       : {datetime.now(timezone.utc).isoformat()}\n")
            fh.write(f"b              : {B}\n")
            fh.write(f"depth_d        : {d}\n")
            fh.write(f"built_rows     : {t}\n")
            fh.write(f"built_cols     : {s}\n")
            fh.write(f"nnz            : {M.nnz}\n")
            fh.write(f"E_source_dual  : {e}\n")
            fh.write(f"F_target_dual  : {f}\n")
            fh.write(f"min_src_tgt    : {min(s, t)}\n")
            fh.write(f"C_b            : {c}\n")
            fh.write(f"predicted_def  : {pred}\n")
            fh.write(f"predicted_rank : {min(e, f) - pred}\n")
            fh.write(f"branch         : {'C_b' if e >= f else 'C_b+E-F'}\n")
            fh.write(f"dense_bytes_est: {_dense_bytes(B, d)}\n")
            fh.write(f"rank           : PENDING (dense F_p rank not run; see header)\n")
            fh.write(f"acyclic_range  : {in_acyclic_range(B, a)} (need a_i >= {2*B+1})\n")
        print(f"  d={d}: built {t} x {s} (rows x cols), nnz={M.nnz}; "
              f"E={e} F={f} -> rows==E:{t==e} cols==F:{s==f}; "
              f"predicted defect={pred} [{'C_b' if e>=f else 'C_b+E-F'}]; wrote {os.path.relpath(path, ROOT)}")
        del M
    print("CHECK:", "PASS (dimensions match E x F)" if ok else "FAIL")
    return ok


def one(d, p, seed, force):
    a = (d, d, d)
    e, f, c = Edim(B, a), Fdim(B, a), C_b(B)
    pred = defect_acyclic(B, a)
    need = _dense_bytes(B, d)
    avail = _avail_bytes()
    if avail is not None and need > 0.85 * avail and not force:
        raise SystemExit(
            f"REFUSING d={d}: dense rank needs ~{need/2**30:.1f} GB but only "
            f"~{avail/2**30:.1f} GB available. Run on a larger host, or pass --force.")
    cmd = f"python3 scripts/run_b7_decisive.py --depths {d} --primes {p}"
    t0 = time.time()
    M, s, t = build_sparse(tuple(-d for _ in range(3)) + (B,), p, seed)
    rk = rank_flint(M, p)
    dt = time.time() - t0
    mn = min(s, t)
    defect = mn - rk
    import flint
    path = os.path.join(ROOT, "logs", f"b7_d{d}_p{p}.txt")
    os.makedirs(os.path.join(ROOT, "logs"), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(f"command       : {cmd}\n")
        fh.write(f"date_utc      : {datetime.now(timezone.utc).isoformat()}\n")
        fh.write(f"python        : {sys.version.split()[0]}\n")
        fh.write(f"platform      : {platform.platform()}\n")
        fh.write(f"numpy         : {np.__version__}\n")
        fh.write(f"python_flint  : {getattr(flint, '__version__', 'unknown')}\n")
        fh.write(f"charge_D      : (-{d},-{d},-{d},{B})\n")
        fh.write(f"b             : {B}\n")
        fh.write(f"depth_d       : {d}\n")
        fh.write(f"prime_p       : {p}\n")
        fh.write(f"seed          : {seed}\n")
        fh.write(f"matrix_rxc    : {t} x {s}\n")
        fh.write(f"E_source_dual : {e}\n")
        fh.write(f"F_target_dual : {f}\n")
        fh.write(f"min_src_tgt   : {mn}\n")
        fh.write(f"rank          : {rk}\n")
        fh.write(f"defect        : {defect}\n")
        fh.write(f"ceiling_48C   : {c}\n")
        fh.write(f"predicted_def : {pred}\n")
        fh.write(f"identity_ok   : {rk + defect == mn} (rank+defect==min(E,F)={mn})\n")
        fh.write(f"matches_pred  : {defect == pred}\n")
        fh.write(f"runtime_sec   : {dt:.2f}\n")
    flag = "MATCHES Theorem 5.1" if defect == pred else "DOES NOT MATCH prediction"
    print(f"  d={d} p={p}: matrix {t}x{s} rank={rk} defect={defect} "
          f"(predicted {pred}; {flag}); rank+defect=min(E,F): {rk+defect==mn}; {dt:.0f}s; wrote {os.path.relpath(path, ROOT)}")
    del M
    return defect


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true",
                    help="build and confirm dimensions only (runs anywhere)")
    ap.add_argument("--depths", nargs="+", type=int, default=DEPTHS_DEFAULT)
    ap.add_argument("--primes", nargs="+", type=int, default=PRIMES_DEFAULT)
    ap.add_argument("--force", action="store_true",
                    help="run the dense rank even if the RAM estimate says it will not fit")
    args = ap.parse_args()
    if args.check:
        sys.exit(0 if check(args.depths) else 1)
    print("=== b=7 DECISIVE RANK (full) ===")
    for d in args.depths:
        for k, p in enumerate(args.primes, 1):
            one(d, p, k, args.force)


if __name__ == "__main__":
    main()
