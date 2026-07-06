#!/usr/bin/env python3
"""Generate raw per-run logs for the decisive b=4 and b=5 exact-F_p rank computations.

Each log records the command, environment, prime, seed, matrix dimensions, rank,
defect, ceiling, runtime, and UTC date, so the heavy computations summarized in
results/b4_decisive_result.json and results/b5_decisive_result.json become fully
auditable execution artifacts rather than stored summaries.

Requires python-flint. The decisive ranks are heavy: b=5 depth 15 is about a
2.2 GB dense matrix over F_p and depth 16 about 3.2 GB, at minutes per prime on
a single core, so run on a machine with adequate memory. reproduce.py does NOT
run these by default; this script is the explicit, documented way to produce the
raw logs for an archival (e.g. Zenodo) release.

Usage:
    python3 scripts/make_decisive_logs.py            # the 7 shipped logs (b4 d13 x3, b5 d14 x2, b5 d15 x2)
    python3 scripts/make_decisive_logs.py --only 5 16 100003   # optional deeper plateau check (needs >4 GB RAM)
    python3 scripts/make_decisive_logs.py --only 5 15 100003
    python3 scripts/make_decisive_logs.py --test     # one cheap b=3 d=9 point, to check the log format
"""
import argparse, sys, time, os, platform
from datetime import datetime, timezone
from math import comb
import numpy as np
try:
    import flint
except ImportError:
    flint = None
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from b4_engine import build_sparse

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def ceiling(b):
    return 48 * comb(b + 1, 3)


def rank_flint(M_csr, p):
    if flint is None:
        raise SystemExit("python-flint is required for this script. Install with: pip install python-flint")
    nr, nc = M_csr.shape
    M = flint.nmod_mat(nr, nc, p)
    coo = M_csr.tocoo()
    r = coo.row; c = coo.col
    v = (coo.data % p).astype(np.int64)
    for i in range(v.shape[0]):
        M[int(r[i]), int(c[i])] = int(v[i])
    return M.rank()


def one(b, d, p, seed=1, outdir=None):
    if outdir is None:
        outdir = os.path.join(ROOT, "logs")
    os.makedirs(outdir, exist_ok=True)
    D = tuple(-d for _ in range(3)) + (b,)
    cmd = f"python3 scripts/make_decisive_logs.py --only {b} {d} {p}"
    t0 = time.time()
    M, s, t = build_sparse(D, p, seed)
    rk = rank_flint(M, p)
    dt = time.time() - t0
    mn = min(s, t); defect = mn - rk
    fv = getattr(flint, "__version__", "unknown")
    path = os.path.join(outdir, f"b{b}_d{d}_p{p}.txt")
    with open(path, "w") as f:
        f.write(f"command       : {cmd}\n")
        f.write(f"date_utc      : {datetime.now(timezone.utc).isoformat()}\n")
        f.write(f"python        : {sys.version.split()[0]}\n")
        f.write(f"platform      : {platform.platform()}\n")
        f.write(f"numpy         : {np.__version__}\n")
        f.write(f"python_flint  : {fv}\n")
        f.write(f"charge_D      : (-{d},-{d},-{d},{b})\n")
        f.write(f"b             : {b}\n")
        f.write(f"depth_d       : {d}\n")
        f.write(f"prime_p       : {p}\n")
        f.write(f"seed          : {seed}\n")
        f.write(f"matrix_rxc    : {t} x {s}\n")
        f.write(f"min_src_tgt   : {mn}\n")
        f.write(f"rank          : {rk}\n")
        f.write(f"defect        : {defect}\n")
        f.write(f"ceiling_48C   : {ceiling(b)}\n")
        f.write(f"runtime_sec   : {dt:.2f}\n")
    print(f"wrote {os.path.relpath(path, ROOT)}: defect={defect} rank={rk} matrix={t}x{s} runtime={dt:.1f}s")
    return defect


# Default set matches the raw logs shipped in logs/. Depth 16 is an optional
# deeper plateau check (run it explicitly with --only 5 16 PRIME on a >4 GB host).
DECISIVE = ([(4, 13, p) for p in (100003, 100019, 100043)]
            + [(5, 14, p) for p in (100003, 100019)]
            + [(5, 15, p) for p in (100003, 100019)])


def main():
    ap = argparse.ArgumentParser(description="Generate raw logs for the decisive b=4/b=5 exact-F_p ranks.")
    ap.add_argument("--only", nargs=3, type=int, metavar=("B", "D", "PRIME"), help="a single (b,d,prime)")
    ap.add_argument("--test", action="store_true", help="one cheap b=3 d=9 point, to check the log format")
    args = ap.parse_args()
    if args.test:
        one(3, 9, 100003); return
    if args.only:
        b, d, p = args.only; one(b, d, p); return
    for b, d, p in DECISIVE:
        one(b, d, p)


if __name__ == "__main__":
    main()
