# Anatomy of a Cohomology Wall

**The Koszul section-map rank defect on the tetraquadric Calabi–Yau threefold**

[![verify](https://github.com/smaniches/cohomology-wall/actions/workflows/verify.yml/badge.svg?branch=main)](https://github.com/smaniches/cohomology-wall/actions/workflows/verify.yml)

Santiago Maniches (ORCID [0009-0005-6480-1987](https://orcid.org/0009-0005-6480-1987)) · Independent research · July 2026

This repository contains the note (`paper/`), all verification code (`scripts/`),
and the machine-readable verification table (`results/`). `reproduce.py` (default) verifies
the b=2 results (closed form and measured F_p ranks), the identity regressions, the
1,728-point sweep, the table, and validates the b=4 engine against the proven b=2
formula, and runs the v0.4 acceptance verifier (`scripts/verify_v04.py`). The b=5 ladder is recorded in `results/b5_decisive_result.json` and
`results/saturation_onset.json`. Raw logs are shipped for d=14 and d=15; depths 11,
12, and 13 are stored/reproducible values rather than raw-log archived runs.
Regenerate individual b=5 logs with `python3 scripts/make_decisive_logs.py --only 5 D PRIME`. The b=5 ceiling is fully logged at depth 15 over two primes. A deeper depth-16 plateau check is recorded as prior evidence (not part of the raw-log archive) in
`results/b5_decisive_result.json`. Recompute it with
`python3 scripts/saturation_threshold.py --bd 5 15` and `--bd 5 16` (add a second prime,
e.g. `--bd 5 15 100019`, to repeat the two-prime check); this is heavier than the default
chain, since the b=5 matrices are about 16k by 16k at depth 15 and larger at depth 16.

---

## Results, with status labels

The note uses explicit claim labels throughout: **[proven]**, **[computed]**, and **[open]**.

| # | Statement | Status |
|---|-----------|--------|
| 1 | Slices of a generic section of `O(2,2,2,2)` form a regular sequence on `(P^1)^3` with finite base scheme of length **48** (slice map is a linear isomorphism `V⊗W ≅ V^⊕3`; Bertini; Cohen–Macaulay) | **proven** (Lemma 3.1) |
| 2 | Closed-form rank/defect law for the family `D=(−a1,−a2,−a3,2)`, all `a_i ≥ 4`, no fitted parameters (decisive vanishing: `H^1(O(a−4)) = 0` identically on the domain) | **proven** (Theorem 3.2) |
| 3 | Saturation ceiling `defect = 48 = 3!·2^3` for all `a_i ≥ 6` — two independent proofs: a finite-difference identity (third difference of a cubic, step 2), and `coker μ ≅ H^0(O_Z)` | **proven** (Theorem 4.1) |
| 4 | Boundary `H^1` corrections arise where individual ambient `H^1` groups switch on (`a_i = 4` entries); the source-target minimum flip is a separate saturation effect | **proven** (Remark 3.4) |
| 5 | Rising sub-ceiling branch below the ceiling, with exact acyclic-range formula `defect = C_b + E − F` when `E < F`; low-b computations match the formula | **proven** in the acyclic Toeplitz range (Theorem 5.1), **computed** as finite-field checks |
| 6 | General Toeplitz determinantal saturation: maximal-minor ideal equals the `(b−1)`st power of the 48-point base-locus ideal; local cokernel length is `48·C(b+1,3)` | **proven** (Theorem 5.1) |
| 7 | Corrected ceiling condition in the acyclic range: `defect = 48·C(b+1,3)` exactly when `(b+1)·Π(a_i−1) ≥ (b−1)·Π(a_i+1)`; otherwise `defect = 48·C(b+1,3)+E−F` | **proven** (Theorem 5.1) |
| 8 | Bicubic reproduction: pencil of cubics saturates at 9 = Bézout; nets give defect 0 | **computed** |

## The finite-field distinction at b=4 and b=5

Two ceiling laws fit the hand-computable data (b=2: 48; b=3: 192):

- **Thom–Porteous:** `48·[(b−1)² + C(b−1,3)]` → 48, 192, **480**, **960**, 1680, …
- **Factorial pattern:** `(b+1)!·2³` → 48, 192, **960**, **5760**, …

They first diverge at **b = 4: 480 vs 960**, again at **b = 5: 960 vs 5760**. Exact
`F_p` rank computation distinguishes the two candidate ceilings and matches the Thom–Porteous values in the computed b=4 and b=5 cases. The interior law
overshoots the ceiling past the onset, so the measured defect plateaus at the Thom--Porteous ceiling in the computed cases:

| b | `a` | matrix | measured defect | interior law |
|---|---|---|---|---|
| 4 | (11,11,11) | 5000×5184 | 296 | 320 |
| 4 | (12,12,12) | 6655×6591 | **480** | 625 |
| 4 | (13,13,13) | 8640×8232 | **480** | 1080 |
| 5 | (13,13,13) | 10368×10976 | 352 | 384 |
| 5 | (15,15,15) | 16464×16384 | **960** | 1296 |
| 5 | (16,16,16) | 20250×19652 | **960** | 2058 |

b=4 plateaus at **480** (three primes, depth 13); b=5 plateaus at **960** (two primes).
**Exact finite-field rank computations at b=4 and b=5 match the proven Toeplitz source-target formula and disagree with the factorial values 960 and 5760 at the first two divergent cases.**
See `results/b4_decisive_result.json` and `results/b5_decisive_result.json`. Raw per-run logs
(command, environment, matrix size, rank, defect, prime, seed, runtime, date) for these finite-field
ranks are produced by `python3 scripts/make_decisive_logs.py` and written to `logs/`.
Logs included in this archive:
- `logs/b4_d13_p100003.txt`, `logs/b4_d13_p100019.txt`, `logs/b4_d13_p100043.txt`
  (b=4, depth 13, three primes: finite-field saturation test, defect 480)
- `logs/b5_d14_p100003.txt`, `logs/b5_d14_p100019.txt`
  (b=5, depth 14, two primes: last sub-ceiling depth below the onset, defect 642)
- `logs/b5_d15_p100003.txt`, `logs/b5_d15_p100019.txt`
  (b=5, depth 15, two primes: first ceiling depth = onset rho(5)=15, defect 960)

The depth-16 plateau check (b=5, defect 960) is recorded as stored prior evidence in
`results/b5_decisive_result.json` and is not part of the raw-log archive. The onset
claim requires only d=14 (sub-ceiling, 642) and d=15 (ceiling, 960), both raw-logged.
Recompute d=16 with: `python3 scripts/make_decisive_logs.py --only 5 16 100003`

`reproduce.py` does not rerun the decisive heavy logs by default.

### Corrected equal-charge onset

The corrected equal-charge threshold is rho(b), defined by the source-target condition E >= F in the acyclic range. For b=2,3,4,5, the old cubic-crossing value tau(b) happens to agree with rho(b). From b=7 onward it can differ; tau(7)=20 while rho(7)=21.

In the acyclic range `a_i >= 2b+1`, Theorem 5.1 proves the defect equals the ceiling
`C_b = 48*C(b+1,3)` exactly when `E >= F`, with `E = (b+1)*prod(a_i-1)` and
`F = (b-1)*prod(a_i+1)`; for `E < F` it is `C_b + E - F`. The corrected equal-charge
onset is therefore

> **rho(b) = min{ d >= 2b+1 : (b+1)(d-1)^3 >= (b-1)(d+1)^3 },**

the first depth at which `E >= F`. The old cubic-crossing reference
`tau(b) = (2b-1) + ceil(2*(b(b-1))^(1/3))` is **not** the saturation threshold in
general. It agrees with `rho(b)` through `b=6` but fails from `b=7`: `tau(7)=20`
while `rho(7)=21`. At `b=7, d=20` the source-target law gives defect 1994
(`E=54872`, `F=55566`, `E<F`), strictly below the ceiling 2688, so `tau(7)=20` is
refuted; saturation begins at `d=21=rho(7)`.

The equal-charge ladders below are exact `F_p` computations
(`scripts/saturation_threshold.py`); they are validation of the proven formula, not
its proof. For `b=2..6` the onset `rho(b)` coincides with `tau(b)`; the `b=7` row is
the proven-formula prediction (Theorem 5.1).

| b | ladder (depth: defect) | ceiling C_b | 2b | 2b+1 | onset rho(b) | tau(b) |
|---|---|---|---|---|---|---|
| 2 | 4:3, 5:24, 6:48, 7:48 | 48 | 4 | 5 | **6** | 6 |
| 3 | 6:4, 7:32, 8:106, 9:192, 10:192 | 192 | 6 | 7 | **9** | 9 |
| 4 | 11:296, 12:480, 13:480 | 480 | 8 | 9 | **12** | 12 |
| 5 | 11:48, 12:158, 13:352, 14:642, 15:960 | 960 | 10 | 11 | **15** | 15 |
| 7 | 20:1994, 21:2688 (confirmed by randomized sparse finite-field computation) | 2688 | 14 | 15 | **21** | 20 |

At `b=5`, `a_i = 2b+1 = 11` gives defect 48, twenty times below the ceiling 960: the
`2b+1` bound is where the acyclic-range formula begins to apply, not where the
ceiling is reached. The general-`b` Toeplitz determinantal theorem is **proven** in
the acyclic range (Theorem 5.1): the maximal-minor ideal of the section-induced
Toeplitz map is `(x,y,z)^(b-1)` locally, the local Buchsbaum-Rim cokernel length is
`48*C(b+1,3)`, and `a_i >= 2b+1` is exactly the term-by-term acyclicity that makes
the global-section complex exact. The `b=7` dense deterministic run (about 24-32 GB,
`scripts/run_b7_decisive.py`) remains pending on an adequate-memory host. In the
interim, the randomized sparse Wiedemann/Kaltofen-Saunders computation
(`scripts/wiedemann_rank.py`; exact arithmetic mod p with randomized rank recovery),
validated against the archived dense b=3-5 ranks and repeated at independent primes,
confirms the predicted defects: 1994 at d=20 and 2688 at d=21, each at two primes
(`logs/b7_d{20,21}_sparse_p{100003,100019}.txt`). See `paper/`
(Theorem 5.1, Corollary 5.2), `results/saturation_onset.json`, and
`results/b7_decisive_prediction.json`.

## Layout

```
cohomology-wall/
├── README.md                  ← this file
├── REPRODUCIBILITY.md         ← modes, runtimes, and what each check verifies
├── CHANGELOG.md               ← release history
├── reproduce.py               ← entry point: all checks; flint makes the exact stages fast; --full forces the numpy fallback
├── requirements.txt           ← deps: numpy, scipy, python-flint (F_p onset gate)
├── LICENSE                    ← MIT (code; see Licensing below)
├── LICENSE-DOCS               ← CC BY 4.0 (paper, text, figures)
├── NOTICE                     ← license split + math-not-copyrightable note
├── CITATION.cff               ← citation metadata
├── MANIFEST.md                ← file inventory with SHA-256 checksums
├── checksums.sha256           ← machine-verifiable (sha256sum -c)
├── .zenodo.json               ← Zenodo deposit metadata
├── .github/workflows/verify.yml ← CI: checksums + reproduce.py + log verifier
├── logs/                      ← raw F_p run logs (b=4, b=5), b=7 dimension checks, b=7 sparse confirmation, archival transcripts
├── paper/
│   ├── tetraquadric_cohomology_wall.tex
│   └── tetraquadric_cohomology_wall.pdf
├── scripts/
│   ├── verify_rank_formula.py ← exact F_p rank (flint if present, numpy fallback) vs closed form; writes CSV
│   ├── identity_checks.py     ← regression: finite-difference identity, vanishings, TP determinant
│   ├── formula_sweep.py       ← 1,728-point sweep, ceiling stress, intersection numbers
│   ├── dims_audit.py          ← source/target dimensions, min-flip table
│   ├── thom_porteous.py       ← TP prediction, symbolic + closed form, b = 2…11
│   ├── b4_engine.py           ← vectorized sparse build + exact F_p rank (decisive b=4)
│   ├── saturation_threshold.py ← exact F_p defect ladders; onset rho(b) vs tau(b), 2b, 2b+1 (--gate)
│   ├── toeplitz_defect.py      ← canonical evaluator (E, F, C_b, rho, tau, defect_acyclic)
│   ├── toeplitz_minor_ideal_check.py ← local model of Theorem 5.1: minor ideal + cokernel length, b=2..9
│   ├── verify_v04.py           ← v0.4 acceptance: b=7 -> 1994/2688; cross-checks b4/b5 logs
│   ├── run_b7_decisive.py      ← b=7 corrected-onset runner (--check anywhere; full dense rank needs ~32-40 GB)
│   └── wiedemann_rank.py       ← randomized sparse F_p rank (Wiedemann/Kaltofen-Saunders); runs on commodity hardware; --validate revalidates against archived dense ranks first
└── results/
    ├── verification_table.csv ← generated by verify_rank_formula.py
    ├── b4_decisive_result.json ← the b=4 ceiling result (480, three primes)
    ├── b5_decisive_result.json ← the b=5 ceiling result (960, two primes)
    ├── saturation_onset.json   ← defect ladders, rho(b) and tau(b) (b=2..7), with provenance
    └── b7_decisive_prediction.json ← b=7 (1994/2688): proven by Theorem 5.1, confirmed dims, confirmed by randomized sparse computation at two primes per depth; dense deterministic rank pending
```

## Quickstart

```bash
pip install -r requirements.txt
python3 reproduce.py          # core checks; FLINT-dependent exact stages run only when python-flint is installed.
python3 reproduce.py --full   # forces the 22-point exact recomputation through the exact NumPy fallback when python-flint is absent.
python3 scripts/verify_v04.py # v0.4 acceptance: b=7 -> 1994/2688 (proven-formula) + b4/b5 raw-log cross-check.
```

Individual scripts:

| Script | Runtime | Expected outcome | Exit code |
|---|---|---|---|
| `scripts/thom_porteous.py` | < 1 s | symbolic == bracket for b=2…11; table incl. 480-vs-960 line | 0 on pass |
| `scripts/identity_checks.py` | < 5 s | six `PASS` lines (R1–R6), no `FAIL` | n/a (visual) |
| `scripts/dims_audit.py` | < 1 s | all 10 defects match; min-flip table | n/a (visual) |
| `scripts/formula_sweep.py` | < 10 s | 0 violations on 1,728 points; max defect = 48 | n/a (visual) |
| `scripts/verify_rank_formula.py` | ~2 s with python-flint; slower pure-numpy fallback otherwise | `ALL 22/22 MATCH`; writes `results/verification_table.csv` | 0 on pass |
| `scripts/saturation_threshold.py` | ~20 s | onset table for rho(b) and tau(b); `--gate` checks the low-b onset calculation and the corrected rho/tau comparison | 0 on pass (`--gate`) |
| `scripts/toeplitz_minor_ideal_check.py` | < 10 s | the two local claims of Theorem 5.1 (minor ideal `(x,y,z)^(b-1)`; per-point cokernel length `C(b+1,3)`) verified exactly for b=2..9 | 0 on pass |
| `scripts/wiedemann_rank.py --validate` | ~5 min | randomized sparse rank reproduces all five archived exact dense ranks (b=3,4,5) | 0 on pass |
| `scripts/wiedemann_rank.py --b7` | ~40 min | randomized sparse confirmation of the b=7 defects at d=20,21, two primes each; runs on commodity hardware (also wired as a manual `.github/workflows/b7-sparse.yml` job on free GitHub runners) | 0 on pass |

## Provenance and independence

All rank computations are exact finite-field linear algebra; no floating point enters
any reported rank. The 22-point verification table uses p = 2³¹−1. The saturation and
decisive b=4/b=5 checks use the smaller independent primes recorded in the result files
and scripts, including 100003, 100019, and 100043 where applicable. The table in the paper was computed by the shipped implementation
(`scripts/verify_rank_formula.py`) using fixed finite-field arithmetic and a fixed
seed (`20260612`). The closed form contains
**no fitted parameters**: it is the output of Theorem 3.2, derived before comparison.

## Honest scope

The cohomology *values* on this manifold are implicit in the region formulae of
Constantin–Lukas (arXiv:1808.09992) and computable by cohomCalg (arXiv:1003.5217);
proofs of cohomology formulae by different methods exist for surfaces and
elliptically fibered threefolds (arXiv:1906.08769, arXiv:2009.01275). The
contribution of this note is a self-contained, elementary, Koszul-level **proof
of the mechanism** (defect = base-locus data; ceiling = intersection number =
third finite difference; boundary `H^1` corrections arise where individual ambient `H^1` groups switch on, the source-target minimum flip being a separate saturation effect) and the proof of the general-`b` Toeplitz determinantal saturation theorem in the acyclic range (Theorem 5.1), with the finite-field computations serving as validation and reproducibility checks rather than as proof.

## Licensing

This repository is dual-licensed by content type:

- **Code** — everything under `scripts/`, plus `reproduce.py`: **MIT** (see `LICENSE`).
- **Paper, text, and figures** — the manuscript under `paper/`, this README, `CHANGELOG.md`, and the figures: **CC BY 4.0** (see `LICENSE-DOCS`).

The mathematical content itself — theorem statements, formulas, the numerical
invariants, and the proof ideas — is not copyrightable and is not restricted by
either license; see `NOTICE`.

Because the paper is CC BY 4.0, it can be included in a citable archive. `CITATION.cff` records the software repository as MIT and provides a `preferred-citation` for the article marked CC-BY-4.0. The repository itself remains split-licensed: code under MIT and paper/docs/figures under CC BY 4.0. See `LICENSE`, `LICENSE-DOCS`, and `NOTICE`.

## Citation

See `CITATION.cff`. Suggested: *S. Maniches, “Anatomy of a Cohomology Wall: the
Koszul section-map rank defect on the tetraquadric Calabi–Yau threefold,” 2026.*
