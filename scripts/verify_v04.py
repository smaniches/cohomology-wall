#!/usr/bin/env python3
"""
verify_v04.py -- v0.4 acceptance verifier for the corrected source-target law.

Two parts, both run with no external dependency (pure arithmetic via
scripts/toeplitz_defect.py):

  PART A  Closed-form acceptance checks (Theorem 5.1 / Corollary 5.2):
            b=7, d=20  ->  defect 1994   (E < F branch: C_b + E - F)
            b=7, d=21  ->  defect 2688   (E >= F branch: C_b)
            b=4, d=13  ->  defect 480
            b=5, d=14  ->  defect 642
            b=5, d=15  ->  defect 960
          plus rho(b)==tau(b) for b=2..6, and rho(7)=21 != tau(7)=20.

  PART B  Raw-log cross-check: every logs/b{b}_d{d}_p{p}.txt is parsed, its
          internal identity rank + defect == min(E,F) is checked, and its
          defect is checked against the Theorem 5.1 formula for that (b,d).
          b=7 logs are OPTIONAL here: they require a dense rank of a roughly
          55000 x 55000 matrix over F_p (about 24 GB), so they are produced on
          an adequate-memory host by scripts/run_b7_decisive.py. If the b=7 raw
          logs are absent, PART B reports them PENDING (not a failure); if they
          are present, they are verified like the others.

Exit 0 iff PART A passes and every PRESENT log is consistent. Missing b=7 logs
do not cause a non-zero exit; they are reported as PENDING.

Usage:
    python3 scripts/verify_v04.py
    python3 scripts/verify_v04.py --require-b7    # also FAIL if b=7 logs absent
"""
import argparse, os, re, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from toeplitz_defect import E, F, C_b, min_EF, defect_acyclic, rho, tau

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Acceptance points: (b, d, expected_defect). The expected values are NOT
# hardcoded magic; they are re-derived from toeplitz_defect.defect_acyclic and
# only listed here as an independent transcription guard.
ACCEPT = [
    (7, 20, 1994),
    (7, 21, 2688),
    (4, 13, 480),
    (5, 14, 642),
    (5, 15, 960),
]


def part_a():
    print("=== PART A: closed-form acceptance (Theorem 5.1 / Corollary 5.2) ===")
    ok = True
    for b, d, exp in ACCEPT:
        a = (d, d, d)
        got = defect_acyclic(b, a)
        e, f = E(b, a), F(b, a)
        branch = "C_b" if e >= f else "C_b+E-F"
        good = (got == exp)
        ok &= good
        print(f"  b={b} d={d}: E={e} F={f} C_b={C_b(b)} min={min_EF(b,a)} "
              f"-> defect={got} [{branch}] expect={exp} {'PASS' if good else 'FAIL'}")
    print("  -- onset (rho) vs old cubic-crossing (tau) --")
    for b in range(2, 7):
        good = (rho(b) == tau(b))
        ok &= good
        print(f"  b={b}: rho={rho(b)} tau={tau(b)} agree={'PASS' if good else 'FAIL'}")
    good7 = (rho(7) == 21 and tau(7) == 20)
    ok &= good7
    print(f"  b=7: rho={rho(7)} tau={tau(7)}  (rho!=tau, tau superseded) "
          f"{'PASS' if good7 else 'FAIL'}")
    print(f"PART A: {'PASS' if ok else 'FAIL'}")
    return ok


def parse_log(path):
    d = {}
    with open(path) as f:
        for line in f:
            if ":" in line:
                k, _, v = line.partition(":")
                d[k.strip()] = v.strip()
    return d


def part_b(require_b7=False):
    print("\n=== PART B: raw-log cross-check vs Theorem 5.1 formula ===")
    logs_dir = os.path.join(ROOT, "logs")
    pat = re.compile(r"b(\d+)_d(\d+)_p(\d+)\.txt$")
    seen_b = {}
    ok = True
    if os.path.isdir(logs_dir):
        for fname in sorted(os.listdir(logs_dir)):
            m = pat.match(fname)
            if not m:
                continue
            b, d, p = int(m.group(1)), int(m.group(2)), int(m.group(3))
            seen_b.setdefault(b, 0)
            seen_b[b] += 1
            log = parse_log(os.path.join(logs_dir, fname))
            try:
                rk = int(log["rank"]); df = int(log["defect"]); mn = int(log["min_src_tgt"])
            except (KeyError, ValueError):
                print(f"  {fname}: FAIL (unparseable rank/defect/min_src_tgt)")
                ok = False
                continue
            a = (d, d, d)
            pred = defect_acyclic(b, a)
            id_ok = (rk + df == mn) and (mn == min_EF(b, a))
            df_ok = (df == pred)
            good = id_ok and df_ok
            ok &= good
            print(f"  {fname}: rank={rk} defect={df} min={mn} "
                  f"(formula def={pred}, min={min_EF(b,a)}) "
                  f"{'PASS' if good else 'FAIL'}")
    # b=7 status
    if seen_b.get(7, 0) == 0:
        msg = ("PENDING: no logs/b7_d*_p*.txt present. The b=7 dense rank is a "
               "~55000 x 55000 F_p matrix (~24 GB); run scripts/run_b7_decisive.py "
               "on an adequate-memory host to produce them.")
        if require_b7:
            print(f"  b=7 raw logs: FAIL ({msg})")
            ok = False
        else:
            print(f"  b=7 raw logs: {msg}")
    else:
        print(f"  b=7 raw logs: {seen_b[7]} present and checked above.")
    print(f"PART B: {'PASS' if ok else 'FAIL'}")
    return ok


def main():
    ap = argparse.ArgumentParser(description="v0.4 acceptance verifier.")
    ap.add_argument("--require-b7", action="store_true",
                    help="also FAIL if the b=7 raw logs are absent")
    args = ap.parse_args()
    a_ok = part_a()
    b_ok = part_b(require_b7=args.require_b7)
    print()
    if a_ok and b_ok:
        print("V0.4 VERIFIER: PASS")
        sys.exit(0)
    print("V0.4 VERIFIER: FAIL")
    sys.exit(1)


if __name__ == "__main__":
    main()
