# Abstraction-discovery benchmark (experimental v1)

This is an executable benchmark for **discovering algebraic structure**, not an
AI theorem prover. It tests predictions on group-held-out parameters, coordinate
changes, and deliberate counterexamples. The local measurements use an exact,
bounded finite-field rank computation. The global answers come from a published
theorem. Neither a high prediction score nor code execution establishes that an
agent rediscovered the theorem or supplied a convincing mathematical explanation.

This directory is an **additional experiment**, not a modification of the
historical proofs, archived rank logs, or originally published research record.
The standard repository manifest/checksum gate covers these new files as well.

## Two explicitly different sources of truth

**Local module (COMPUTED).** For `1 <= n <= 8`, the matrix `T_n(L0,L1,L2)` has
`n` rows and `n+2` columns, with entry `L_(j-i)` for `j-i` in `{0,1,2}` and zero
otherwise. The three linear forms in `x,y,z` have coefficients in the prime
field `F_101`. The target is the degree-`d` dimension of

    C_n = coker(T_n : F_101[x,y,z](-1)^(n+2) -> F_101[x,y,z]^n).

The oracle enumerates degree-`d` monomials, constructs the induced multiplication
matrix, and takes its exact modular rank with sparse column elimination. A
separately written dense *row* elimination is used in the tests. `forms` is a
3-by-3 array, one row of coefficients per linear form. The identity forms are
`[x,y,z]`. Invertible changes of the variable basis preserve the graded
Hilbert function; singular triples need not. Computations over `F_101` do not by
themselves establish a characteristic-zero generalization.

**Global defect (THEOREM-DERIVED).** For `D=(-a1,-a2,-a3,b)` with `b >= 2`,
`a_i >= 2b+1`, and `b <= 12`, `a_i <= 64`, the reference returns the
source-target defect from the existing *Theorem 5.1* in `cohomology-wall`.
The returned integer is **not** a newly measured Koszul rank. It must never be
counted as an independent reproduction of the theorem. The verified boundary
values `b=7, a=(20,20,20) -> 1994` and `b=7, a=(21,21,21) -> 2688` are
included as oracle regression tests.

Both oracles are deterministic and use only the Python standard library.

## Run

From the repository root:

```bash
python3 discovery/benchmark.py oracle --family local --n 3 --degree 2
python3 discovery/benchmark.py oracle --family local --n 1 --degree 3 --forms singular
python3 discovery/benchmark.py oracle --family global --b 7 --a 20,20,20

python3 discovery/benchmark.py generate --seed 20261122 --out /tmp/discovery-demo
python3 discovery/benchmark.py baseline \
  --train /tmp/discovery-demo/train.json \
  --challenge /tmp/discovery-demo/challenge.json \
  --out /tmp/discovery-demo/baseline.json
python3 discovery/benchmark.py score \
  --manifest /tmp/discovery-demo/manifest.json \
  --challenge /tmp/discovery-demo/challenge.json \
  --predictions /tmp/discovery-demo/baseline.json \
  --out /tmp/discovery-demo/report.json

python3 -m unittest discover -s tests -p 'test_*.py'
```

`generate` creates `train.json` (with answers), `challenge.json` (without
answers), and `manifest.json` (seed, counts, SHA-256 of both canonicalized sets).
For an honest prospective test, choose the *challenge seed after freezing*
the candidate submission and store the manifest digest externally. Changing the
seed changes case order, not the groups or difficulty. A later version needs
new parameter generators for genuine multiple independent test sets.

Submission `predictions.json` is a JSON list containing **exactly one** record
for each challenge case:

```json
[
  {"id": "<copy a case id>", "prediction": 0}
]
```

The illustrative list above is **not** a complete submission. Predictions must
be nonnegative integers, IDs must match the supplied challenge, and every case
must be covered exactly once. The evaluator rejects duplicate JSON keys,
additional case fields, contaminated challenge answers, unexpected IDs,
missing predictions, invalid value types, and mismatched challenge manifests.
It **does not execute submitted model code**. It reports correctness by regime
and deliberately leaves `mechanism_assessment` as `NOT_EVALUATED`.

## The split

| Regime | Train | Challenge | What it probes |
| --- | --- | --- | --- |
| Local module | n=1,2,3; all degrees through n+1, standard x/y/z | n=4,5,6; all degrees through n+1 | Extrapolation in presentation size |
| Local change of coordinates | No changed-basis training | Invertible changes for n=4,5,6 | Representation invariance |
| Local singular counterexample | No degenerate training | Rank-deficient linear-form triples at n=1,2,3 | Detecting a missing hypothesis |
| Global dimension/defect | b=2,3,4, equal and unequal charges | b=5,6,7,8, around the onset and unequal charges | Nonlinear ceiling and source-target transition |

Default seed 20261122 creates **26 training cases and 72 challenge cases**.
The deliberately weak reference baseline memorizes local values and fits a
quadratic ceiling in `b` to the global training labels while retaining the
known source-target dimension correction. It gets **12/72** challenge answers
right on the version checked here. This is a harness sanity check, **not** a
state-of-the-art mathematical discovery baseline. A model can obtain full
numerical accuracy by computing the same ranks or reading the existing theorem.

## What numerical scoring cannot establish

Mathematical explanation requires a separate review. A credible mechanism
submission should state:

1. The graded presentation and the exact assumptions on the linear forms.
2. A general structural claim that predicts the Hilbert function, its support,
   and its total length, with a proof or independently checkable certificate.
3. Why invertible variable changes preserve that claim and an explicit
   counterexample when the forms become singular.
4. For the global result, the relationship between source-target dimensions,
   the local cokernel, and onset. A formula fitted to `b=2,3,4` is insufficient.
5. At least one **new, independently tested consequence** not used to construct
   the hypothesis, and a comparison with direct elimination and interpolation.

These are review criteria, not implemented automatic proof checks. If a model
simply reads Theorem 5.1 or the evaluator source, then predicts its values,
that is legitimate **software reproduction** but not evidence of independent
abstraction discovery. A participant able to inspect this repository can also
query the oracle on challenge parameters; consequently these challenge labels
are only *protocol-blinded*, not technically secret. Genuine blinded testing
requires an independently operated evaluator, withheld source/answers and
post-freeze challenge construction. Public mathematical results may already
appear in model training data, a separate contamination question.

The first scientific outcome of this benchmark is the reproducible **evaluation
infrastructure**, not a positive claim about any AI model's discovery ability.
