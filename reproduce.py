#!/usr/bin/env python3
"""Single entry point: reproduce and verify the quantitative claims in the note.

Runs, in order:
  1. scripts/thom_porteous.py            (TP prediction; internal consistency)
  2. scripts/identity_checks.py          (R1-R6 regression identities)
  3. scripts/dims_audit.py               (source/target dimensions, defects)
  4. scripts/formula_sweep.py            (1,728-point sweep, ceiling stress)
  5. scripts/b4_engine.py --validate     (sparse F_p engine vs proven b=2 formula)
  6. scripts/verify_v04.py               (Theorem 5.1 acceptance: b=7 -> 1994/2688,
                                          rho vs tau, and b4/b5 raw-log cross-check)
  7. scripts/saturation_threshold.py --gate (onset rho(b); == tau(b) for b=2,3)
  8. scripts/verify_rank_formula.py      (exact from-scratch recomputation of all 22
                                          table points vs the closed form)
then verifies number provenance: the shipped results/verification_table.csv must
equal the paper's Table (Section 7) row for row.

  python3 reproduce.py
      Default. Runs all checks. FLINT-dependent stages use python-flint when present
      (fast, exact); they are reported SKIP (not FAIL) when absent, and the run
      still exits 0 on the remaining checks.

  python3 reproduce.py --smoke
      Lightweight check only. Skips FLINT-dependent stages with a WARNING.
      Suitable for CI or environments without python-flint.

  python3 reproduce.py --full
      Forces stage 7 to run even without python-flint (pure-numpy fallback).
      Exact but slower; progress streams row by row.

  python3 reproduce.py --archival
      Release-grade mode. Fails immediately if python-flint is absent. Runs all
      FLINT-dependent stages (does NOT rerun the decisive b=4 and b=5 heavy logs;
      verify those separately with scripts/verify_decisive_logs.py). Saves a
      transcript to logs/. Use before minting a Zenodo DOI or tagging a release.
      Full archival verification:
        python3 reproduce.py --archival
        python3 scripts/verify_decisive_logs.py

Stages run as child processes with single-threaded BLAS, unbuffered output, and a
hard per-stage timeout (env REPRODUCE_STAGE_TIMEOUT, default 600s; the no-flint
--full fallback uses REPRODUCE_FULL_TIMEOUT, default 5400s).
"""
import argparse, csv, importlib.util, os, subprocess, sys, tempfile

ROOT = os.path.dirname(os.path.abspath(__file__))
PY = sys.executable

# Paper Table (Section 7): (a1,a2,a3) -> (rank, defect)
PAPER_TABLE = {
    (4,4,4): (78, 3),   (4,4,5): (102, 6),  (4,4,6): (126, 9),
    (4,5,5): (132, 12), (4,5,6): (162, 18), (4,6,6): (197, 28),
    (4,6,7): (232, 38), (4,6,8): (267, 48), (4,7,7): (272, 48),
    (4,7,8): (312, 48), (5,5,5): (168, 24), (5,5,6): (204, 36),
    (5,6,6): (246, 48), (5,6,7): (288, 48), (5,6,8): (330, 48),
    (5,7,7): (336, 48), (5,7,8): (384, 48), (6,6,6): (295, 48),
    (6,6,7): (344, 48), (6,6,8): (393, 48), (6,7,7): (400, 48),
    (7,7,7): (464, 48),
}

# (name, rel, markers, exit_matters, needs_flint, full_only)
STAGES = [
    ("Thom-Porteous prediction", "scripts/thom_porteous.py",
     ["Internal agreement (symbolic == bracket): PASS"], True, False, False),
    ("Identity regression R1-R6", "scripts/identity_checks.py",
     ["R1 finite-difference identity, 28^3 lattice incl. negatives: PASS",
      "R2 min=src for a_i>=6 (34^3 sweep): PASS",
      "R3 H1(a-4)=0 on domain: PASS",
      "R4 TP det == 48[(b-1)^2+C(b-1,3)] for b=2..11: PASS",
      "R5 retrodiction: PASS",
      "PASS"], False, False, False),
    ("Dimension audit", "scripts/dims_audit.py",
     ["matches original measurements"], False, False, False),
    ("Formula sweep (1,728 pts)", "scripts/formula_sweep.py",
     ["formula violations (neg rank, neg defect, rank>min): 0",
      "max defect over [4,15]^3 = 48"], False, False, False),
    ("b=4 engine validation (vs proven b=2 formula)", "scripts/b4_engine.py --validate",
     ["VALIDATION PASS"], True, False, False),
    ("v0.4 acceptance (b=7 -> 1994/2688; b4/b5 log cross-check)", "scripts/verify_v04.py",
     ["V0.4 VERIFIER: PASS"], True, False, False),
    ("Saturation onset (rho source-target onset; b=2,3)", "scripts/saturation_threshold.py --gate",
     ["ONSET CHECK (b=2,3): PASS"], True, True, False),
    # Exact from-scratch recomputation of all 22 points. Uses flint when present
    # (~2s on tested hardware; exact and BLAS-independent) and skips cleanly when
    # flint is absent. Streams
    # progress row by row; hard timeout guards a pathological numpy fallback.
    ("Independent engine vs closed form (exact 22-pt recomputation)", "scripts/verify_rank_formula.py",
     ["ALL 22/22 MATCH"], True, True, False),
]

STAGE_TIMEOUT = int(os.environ.get("REPRODUCE_STAGE_TIMEOUT", "600"))
FULL_TIMEOUT  = int(os.environ.get("REPRODUCE_FULL_TIMEOUT", "5400"))

def _child_env():
    env = dict(os.environ)
    env.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
               MKL_NUM_THREADS="1", NUMEXPR_NUM_THREADS="1",
               PYTHONUNBUFFERED="1",
               # Children must emit UTF-8 on every platform. Without this, a default
               # Windows console (cp1252) encodes characters like the em-dash as byte
               # 0x97, which then breaks the captured-output decode below.
               PYTHONUTF8="1", PYTHONIOENCODING="utf-8")
    return env

def run_stage(name, rel, markers, exit_matters, live=False, timeout=STAGE_TIMEOUT):
    print(f"\n=== {name} ({rel}) ===", flush=True)
    parts = rel.split()
    cmd = [PY, "-u", os.path.join(ROOT, parts[0])] + parts[1:]
    if live:
        # Inherit stdout/stderr: live row-by-row progress, no pipe to deadlock.
        print(f"    running (heavy; streams below; timeout {timeout}s) ...", flush=True)
        try:
            rc = subprocess.run(cmd, cwd=ROOT, env=_child_env(), timeout=timeout).returncode
        except subprocess.TimeoutExpired:
            print(f"    TIMEOUT after {timeout}s", flush=True)
            print("--- stage verdict: FAIL", flush=True)
            return False
        ok = (rc == 0)
        print(f"    done ({'PASS' if ok else 'FAIL'})", flush=True)
        print(f"--- stage verdict: {'PASS' if ok else 'FAIL'}", flush=True)
        return ok
    print("    running ...", flush=True)
    with tempfile.TemporaryFile() as tf:  # binary capture; decoded explicitly below
        try:
            rc = subprocess.run(cmd, cwd=ROOT, stdout=tf, stderr=subprocess.STDOUT,
                                 timeout=timeout, env=_child_env()).returncode
        except subprocess.TimeoutExpired:
            print(f"    TIMEOUT after {timeout}s", flush=True)
            print("--- stage verdict: FAIL", flush=True)
            return False
        tf.seek(0)
        # Decode robustly: a stray non-UTF-8 byte must never abort the harness. The
        # stage's pass/fail is decided by its markers and exit code below, not by
        # decoding, so errors="replace" is safe and keeps Windows runs from crashing.
        out = tf.read().decode("utf-8", errors="replace")
    print(out.rstrip(), flush=True)
    ok = all(m in out for m in markers) and "FAIL" not in out and "MISMATCH FOUND" not in out
    if exit_matters:
        ok = ok and rc == 0
    print(f"    done ({'PASS' if ok else 'FAIL'})", flush=True)
    print(f"--- stage verdict: {'PASS' if ok else 'FAIL'}", flush=True)
    return ok

def check_provenance():
    print("\n=== Number provenance: CSV vs paper Table (Section 7) ===")
    path = os.path.join(ROOT, "results", "verification_table.csv")
    if not os.path.exists(path):
        print("FAIL: results/verification_table.csv not found"); return False
    ok, n = True, 0
    with open(path) as f:
        for row in csv.DictReader(f):
            a = (int(row["a1"]), int(row["a2"]), int(row["a3"]))
            exp = PAPER_TABLE.get(a)
            got = (int(row["meas_rank"]), int(row["meas_defect"]))
            pred = (int(row["pred_rank"]), int(row["pred_defect"]))
            if exp is None or got != exp or pred != exp:
                print(f"FAIL at {a}: paper={exp} measured={got} predicted={pred}")
                ok = False
            n += 1
    if n != len(PAPER_TABLE):
        print(f"FAIL: {n} rows in CSV, {len(PAPER_TABLE)} in paper table"); ok = False
    if ok:
        print(f"All {n} paper-table numbers trace to results/verification_table.csv: PASS")
    return ok

def main():
    ap = argparse.ArgumentParser(description="Reproduce and verify the note's quantitative claims.")
    ap.add_argument("--smoke", action="store_true",
                    help="lightweight: skip FLINT stages with WARNING.")
    ap.add_argument("--full", action="store_true",
                    help="force the 22-point recomputation even without python-flint.")
    ap.add_argument("--archival", action="store_true",
                    help="release-grade: fails if python-flint absent; saves transcript. "
                         "Pair with verify_decisive_logs.py for full archival check.")
    args = ap.parse_args()
    if args.archival and args.smoke:
        sys.exit("error: --archival and --smoke are mutually exclusive")

    flint_ok = importlib.util.find_spec("flint") is not None
    if args.archival and not flint_ok:
        sys.exit("ARCHIVAL MODE REQUIRES python-flint. "
                 "Install with: pip install python-flint")
    statuses, names = [], []
    for name, rel, markers, exit_matters, needs_flint, full_only in STAGES:
        names.append(name)
        is_verify = "verify_rank_formula" in rel
        if args.smoke and needs_flint:
            print("\n=== %s (%s) ===" % (name, rel))
            print("WARNING: skipped in --smoke mode (requires python-flint).")
            print("--- stage verdict: SKIP (smoke mode)")
            statuses.append("SKIP")
            continue
        skip = needs_flint and not flint_ok and not (is_verify and args.full) and not args.archival
        if skip:
            print("\n=== %s (%s) ===" % (name, rel))
            print("SKIP: python-flint is not installed. This stage runs an exact")
            print("finite-field (F_p) rank; install python-flint to enable it and")
            print("obtain the full verification (see requirements.txt).")
            if is_verify:
                print("(run 'python3 reproduce.py --full' to force the pure-numpy fallback.)")
            print("--- stage verdict: SKIP (optional dependency missing)")
            statuses.append("SKIP")
        else:
            slow_fallback = is_verify and not flint_ok
            ok = run_stage(name, rel, markers, exit_matters,
                           live=slow_fallback,
                           timeout=(FULL_TIMEOUT if slow_fallback else STAGE_TIMEOUT))
            statuses.append("PASS" if ok else "FAIL")
    statuses.append("PASS" if check_provenance() else "FAIL")
    names.append("Number provenance (CSV == paper table)")

    import datetime
    mode = ("ARCHIVAL" if args.archival else
            "SMOKE" if args.smoke else
            "FULL" if args.full else "DEFAULT")
    summary_lines = (
        ["", "=" * 56, "REPRODUCTION SUMMARY  [" + mode + "]"] +
        ["  %-54s %s" % (nm, st) for nm, st in zip(names, statuses)] +
        ["=" * 56]
    )
    for line in summary_lines: print(line)
    if args.archival:
        os.makedirs(os.path.join(ROOT, "logs"), exist_ok=True)
        date_str = datetime.date.today().isoformat()
        transcript = os.path.join(ROOT, "logs", "reproduce_archival_" + date_str + ".txt")
        with open(transcript, "w") as tf:
            tf.write("\n".join(summary_lines) + "\n")
        print("Archival transcript saved: " + transcript)
    failed = any(st == "FAIL" for st in statuses)
    if failed:
        print("FINAL: FAILURES PRESENT")
    elif "SKIP" in statuses:
        print("FINAL: CORE CHECKS PASS; OPTIONAL FLINT CHECKS SKIPPED")
        print("Note: install python-flint (see requirements.txt) for the exact finite-field")
        print("stages, or run 'python3 reproduce.py --full' to force the pure-numpy")
        print("recomputation of the 22-point table.")
    else:
        print("FINAL: ALL CHECKS PASS")
    sys.exit(1 if failed else 0)

if __name__ == "__main__":
    main()
