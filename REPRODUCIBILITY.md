# Reproducibility Guide

This document describes how to verify the quantitative claims in
"Anatomy of a Cohomology Wall: The Koszul section-map rank defect on the
tetraquadric Calabi-Yau threefold" (Santiago Maniches, 2026).

---

## Quick reference

| Command | Mode | python-flint required | Approximate runtime |
|---|---|---|---|
| `python3 reproduce.py --smoke` | smoke | No | < 10 s |
| `python3 reproduce.py` | default | Optional (skipped if absent) | < 30 s with flint; < 2 min without |
| `python3 reproduce.py --full` | full | No (NumPy fallback) | Minutes (machine-dependent) |
| `python3 reproduce.py --archival` | archival (stage checks) | **Yes** | < 2 min with flint |
| `python3 scripts/verify_decisive_logs.py` | log verification | No | < 5 s |
| `python3 scripts/verify_v04.py` | v0.4 acceptance (b=7 + log cross-check) | No | < 5 s |

---

## Modes

### Smoke (`--smoke`)
Runs only the lightweight algebraic and combinatorial checks.
Skips all FLINT-dependent stages with a printed WARNING.
Exits 0 even when stages are skipped.

```bash
python3 reproduce.py --smoke
```

Checks run: Thom-Porteous prediction, identity regression R1-R6,
dimension audit, formula sweep (1,728 points), number provenance (CSV == paper table).

Checks skipped: saturation onset gate, 22-point exact recomputation.

### Default (`python3 reproduce.py`)
Runs all stages. FLINT-dependent stages (saturation onset gate, 22-point
exact recomputation) use python-flint when present and are reported SKIP
(not FAIL) when absent. The run exits 0 either way.

### Full (`--full`)
Forces the 22-point exact recomputation to run even without python-flint,
using the exact modular NumPy fallback. This fallback is mathematically
exact but its runtime depends on the machine and NumPy/BLAS build.
Progress streams row by row to the console.

```bash
python3 reproduce.py --full
```

### Archival (`--archival`)
Release-grade mode. **Fails immediately if python-flint is absent.**
Runs all stages including all FLINT-dependent checks. Saves a dated
transcript to `logs/reproduce_archival_<date>.txt`.

Use this mode before minting a Zenodo DOI or tagging a release.
Note: --archival does NOT rerun the decisive b=4 and b=5 heavy logs.
For complete archival verification, also run the log verifier:

```bash
python3 reproduce.py --archival
python3 scripts/verify_decisive_logs.py
```

---

## Expected outputs

### Default (with python-flint)
```
REPRODUCTION SUMMARY  [DEFAULT]
  Thom-Porteous prediction                       PASS
  Identity regression R1-R6                      PASS
  Dimension audit                                PASS
  Formula sweep (1,728 pts)                      PASS
  b=4 engine validation (vs proven b=2 formula)  PASS
  v0.4 acceptance (b=7 -> 1994/2688; b4/b5 log cross-check) PASS
  Saturation onset (rho source-target onset; b=2,3)  PASS
  Independent engine vs closed form (exact ...)  PASS
  Number provenance (CSV == paper table)         PASS
FINAL: ALL CHECKS PASS
```

### Default (without python-flint)
```
FINAL: CORE CHECKS PASS; OPTIONAL FLINT CHECKS SKIPPED
```

---

## python-flint requirement

python-flint is required for:
- The saturation onset gate (stages run via `saturation_threshold.py --gate`)
- The exact 22-point recomputation (`verify_rank_formula.py`)
- The b=4 decisive sweep: `scripts/run_decisive_flint.py`
- The b=5 raw-log generation: `scripts/make_decisive_logs.py`

Install: `pip install python-flint`

python-flint is **not** required for the smoke mode or the b=4 construction
validation (`b4_engine.py --validate`).

---

## What each check verifies

### Verified by closed form (b=2 theorem)
- Rank formula (Theorem 3.2): the closed-form rank of the Koszul section map
  for `D = (-a1,-a2,-a3,2)`, all `a_i >= 4`.
- Saturation ceiling (Theorem 4.1): defect saturates at 48 for all `a_i >= 6`,
  by finite-difference identity and base-locus length.
- Identity regression R1-R6: six algebraic identities derived from the rank formula.
- 1,728-point formula sweep: stress-test of the closed form on a dense charge grid.

### Verified by exact finite-field computation
- 22-point verification table (Section 7): exact GF(p) ranks at 22 sampled
  charge vectors, matching the closed-form predictions to machine precision.
- Saturation onset gate (b=2,3): the onset equals the corrected source-target
  threshold rho(b) in the proven (b=2) and computed (b=3) cases. For b=2..6,
  rho(b) equals the old cubic-crossing tau(b); they first differ at b=7
  (tau(7)=20, rho(7)=21), which is the corrected statement of Theorem 5.1.
- b=4 construction validation: the b=4 section-map matrix passes the b=2
  closed-form consistency check.

### The b=4 and b=5 decisive finite-field computations
These are not run by `reproduce.py`; they require >1 GB RAM and python-flint.

b=4, depth 13, three primes (100003, 100019, 100043):
  defect = 480, matching the Thom-Porteous candidate.
  Raw logs: `logs/b4_d13_p100003.txt`, `logs/b4_d13_p100019.txt`, `logs/b4_d13_p100043.txt`
  Reproduce: `python3 scripts/run_decisive_flint.py` or `scripts/make_decisive_logs.py`

b=5, depth 14, two primes (100003, 100019):
  defect = 642, sub-ceiling (last depth below the onset rho(5)=15).
  Raw logs: `logs/b5_d14_p100003.txt`, `logs/b5_d14_p100019.txt`

b=5, depth 15, two primes (100003, 100019):
  defect = 960, first ceiling depth = onset rho(5)=15 (tau(5)=15 coincides here).
  Raw logs: `logs/b5_d15_p100003.txt`, `logs/b5_d15_p100019.txt`

### The b=7 corrected-onset witness (tau(7)=20 != rho(7)=21)
This is the decisive distinction between the old tau(b) and the corrected rho(b).
In the acyclic range, Theorem 5.1 predicts:

  b=7, d=20:  E=54872, F=55566, E<F  -> defect = C_7 + E - F = 1994  (not the ceiling 2688)
  b=7, d=21:  E=64000, F=63888, E>=F -> defect = C_7            = 2688

The matrix dimensions are confirmed by direct construction (rows x cols = E x F). The b=7 values 1994 and 2688 are theorem-derived validation targets. Dense finite-field b=7 rank logs are pending and are not claimed as part of the raw-log archive.
The dense F_p rank is a ~55000 x 55000 matrix over F_p (about 24 GB as a
flint.nmod_mat) and is NOT run by default. Produce it on a >= ~32 GB host:

  python3 scripts/run_b7_decisive.py --check    # dims only; runs anywhere
  python3 scripts/run_b7_decisive.py            # full rank; writes logs/b7_d{20,21}_p*.txt
  python3 scripts/verify_v04.py                 # checks predictions and any present logs
  python3 scripts/verify_v04.py --require-b7     # also fails if b=7 logs absent

These predictions follow from Theorem 5.1, whose two branches were reproduced over
F_p at b=4 (480) and b=5 (48,158,352,642,960).

---

## Stored prior evidence (not raw-log archived)

b=5, depth 16: defect 960 (stored plateau check; two primes).
This result is NOT part of the raw-log archive. It is recorded in
`results/b5_decisive_result.json` as stored prior evidence.
The onset claim requires only d=14 (sub-ceiling) and d=15 (ceiling),
both of which are raw-logged.

To regenerate: `python3 scripts/make_decisive_logs.py --only 5 16 100003`
(requires >4 GB RAM and python-flint).

---

## Note on finite-field transfer

The characteristic-zero proof for general b in the acyclic range a_i >= 2b+1 is
Theorem 5.1 (Bertini, the Toeplitz minor-ideal computation, and Buchsbaum-Rim
acyclicity). The finite-field computations certify the stated ranks for the
sampled sections and serve as validation and reproducibility checks of that
theorem; they are not themselves the proof.

---

## Runtimes (tested hardware)

All runtimes measured on a modern x86-64 machine with python-flint installed.

| Stage | Runtime |
|---|---|
| Smoke mode (full) | < 10 s |
| Default mode (full, with flint) | < 30 s |
| Archival mode | < 2 min |
| 22-point NumPy fallback (--full, no flint) | 5-30 min (machine-dependent) |
| b=4 d=13 per prime (run_decisive_flint.py) | ~2-5 min |
| b=5 d=15 per prime | ~8-15 min |
