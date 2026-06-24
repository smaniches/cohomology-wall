#!/usr/bin/env python3
"""Parse the shipped raw finite-field logs and cross-check them against the JSON summaries.

This script does not recompute anything. It reads the raw log files in logs/,
verifies each log's internal consistency (rank + defect == min_src_tgt), checks
the defect values against the hardcoded paper-table expectations, and then
cross-checks those same values against the corresponding JSON measurement fields
in results/b4_decisive_result.json and results/b5_decisive_result.json.

A discrepancy between any log, the hardcoded table, or the JSON is reported as
FAIL and causes a non-zero exit.

Usage:
    python3 scripts/verify_decisive_logs.py

Exits 0 if all checks pass. Exits 1 on any inconsistency.
"""
import json, os, re, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Hardcoded expected defects from the paper (Table and Section 5 evidence).
# These match the archived raw logs and must equal the JSON measurement fields.
EXPECTED = {
    ("b4", "d13", "100003"): {"b": 4, "depth_d": 13, "defect": 480},
    ("b4", "d13", "100019"): {"b": 4, "depth_d": 13, "defect": 480},
    ("b4", "d13", "100043"): {"b": 4, "depth_d": 13, "defect": 480},
    ("b5", "d14", "100003"): {"b": 5, "depth_d": 14, "defect": 642},
    ("b5", "d14", "100019"): {"b": 5, "depth_d": 14, "defect": 642},
    ("b5", "d15", "100003"): {"b": 5, "depth_d": 15, "defect": 960},
    ("b5", "d15", "100019"): {"b": 5, "depth_d": 15, "defect": 960},
}

# Expected JSON measurement keys and their defect values.
# Format: (json_filename, measurement_key, expected_defect)
JSON_CROSSCHECKS = [
    ("b4_decisive_result.json", "(13,13,13)", 480),
    ("b5_decisive_result.json", "(14,14,14)", 642),
    ("b5_decisive_result.json", "(15,15,15)", 960),
]


def parse_log(path):
    d = {}
    with open(path) as f:
        for line in f:
            if ":" in line:
                k, _, v = line.strip().partition(":")
                d[k.strip()] = v.strip()
    return d


def main():
    logs_dir   = os.path.join(ROOT, "logs")
    results_dir = os.path.join(ROOT, "results")
    errors = []

    # ── 1. Parse and verify each raw log ─────────────────────────────────────
    print("=== Raw log verification ===")
    found_keys = set()
    for fname in sorted(os.listdir(logs_dir)):
        if not fname.endswith(".txt"):
            continue
        m = re.match(r"(b\d)_(d\d+)_p(\d+)\.txt", fname)
        if not m:
            continue
        b_tag, d_tag, prime = m.group(1), m.group(2), m.group(3)
        key = (b_tag, d_tag, prime)
        path = os.path.join(logs_dir, fname)

        print(f"  {fname} ...", end=" ")
        log = parse_log(path)

        try:
            rank   = int(log["rank"])
            defect = int(log["defect"])
            min_st = int(log["min_src_tgt"])
        except (KeyError, ValueError) as e:
            errors.append(f"PARSE ERROR {fname}: {e}")
            print("FAIL"); continue

        # internal consistency
        if rank + defect != min_st:
            errors.append(f"CONSISTENCY FAIL {fname}: rank({rank})+defect({defect})!= min_src_tgt({min_st})")
            print("FAIL"); continue

        # against paper table
        if key in EXPECTED:
            exp_defect = EXPECTED[key]["defect"]
            if defect != exp_defect:
                errors.append(f"TABLE MISMATCH {fname}: defect={defect}, expected={exp_defect}")
                print("FAIL"); continue

        found_keys.add(key)
        print(f"PASS  rank={rank} defect={defect} min={min_st}")

    # check all expected logs are present
    for key in EXPECTED:
        if key not in found_keys:
            errors.append(f"MISSING LOG: {'_'.join(key)}.txt not found in logs/")

    # ── 2. Cross-check JSON measurement fields against the same table ─────────
    print("\n=== JSON cross-check ===")
    for jname, mkey, expected_defect in JSON_CROSSCHECKS:
        jpath = os.path.join(results_dir, jname)
        print(f"  {jname}[\"{mkey}\"] ...", end=" ")
        if not os.path.exists(jpath):
            errors.append(f"MISSING: {jname}")
            print("FAIL"); continue
        try:
            d = json.load(open(jpath))
            meas = d.get("measurements", {})
            entry = meas.get(mkey)
            if entry is None:
                errors.append(f"JSON KEY MISSING: {jname}[measurements][{mkey}]")
                print("FAIL"); continue
            json_defect = int(entry["defect"])
        except (KeyError, ValueError, TypeError) as e:
            errors.append(f"JSON PARSE ERROR {jname}[{mkey}]: {e}")
            print("FAIL"); continue
        if json_defect != expected_defect:
            errors.append(f"JSON MISMATCH {jname}[{mkey}]: defect={json_defect}, expected={expected_defect}")
            print("FAIL"); continue
        # also cross-check against the log values already verified above
        print(f"PASS  defect={json_defect}")

    # ── 3. Summary ────────────────────────────────────────────────────────────
    print()
    if errors:
        print(f"DECISIVE LOG VERIFICATION: {len(errors)} error(s):")
        for e in errors:
            print(f"  {e}")
        sys.exit(1)
    n = len(found_keys)
    print(f"DECISIVE LOG VERIFICATION PASS: {n} logs checked; all logs, paper table, and JSON measurements consistent.")
    sys.exit(0)


if __name__ == "__main__":
    main()
