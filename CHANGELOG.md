# Changelog

## v0.4.1 - 2026-06-23 - Release-clean patch

Release-engineering patch over v0.4.0. No change to the theorem or to any numerical result; the general Toeplitz determinantal theorem (Theorem 5.1) and the corrected source-target onset are unchanged.

- README: removed file-tree fragments that had been spliced into the opening paragraph; restored the prose and the `scripts/` and `results/` trees. The equal-charge onset is stated as `rho(b)` (the source-target condition `E >= F`); `tau(b)` is described only as the old cubic-crossing value that coincides with `rho(b)` for `b=2,3,4,5` and can differ from `b=7` (`tau(7)=20`, `rho(7)=21`).
- Paper, Section 9: the `[computed]` status line no longer lists `tau(b)=6,9,12,15` as a saturation onset; it records the low-b coincidence `rho(b)=tau(b)` for `b=2,3,4,5`, with `rho(b)` the corrected source-target threshold.
- Theorem 5.1 proof hardening: the term order is specified as a weight-refined order (`y`-degree, then lex in `x,z`) and the displayed monomial is shown to be the unique, non-cancelling leading term; faithfully-flat descent from the completed local ring to the local ring is made explicit for the Fitting-ideal equality; Buchsbaum-Rim exactness is cited to the Buchsbaum-Eisenbud acyclicity criterion (Eisenbud, Commutative Algebra, Thm. 20.9) via the expected codimension 3; and the graded Hilbert-series length is reconciled with the completed local length at each reduced base point.
- Claim tags: Corollary 5.2 tagged `[proven]`; Proposition 5.3 (the b=4 finite-field distinction) tagged `[computed]` and described explicitly as a finite-field computation rather than a theorem. Cross-references to it now read Proposition, not Theorem.
- First use of `tau` and `rho` is annotated to state that both are equal-charge thresholds (defined for `a_1=a_2=a_3`), so `tau` is not applied to asymmetric charges.
- Cosmetics: a single release date 2026-06-23 across the PDF, CITATION.cff, CHANGELOG, and .zenodo.json; contact email lowercased to `santiago.maniches@proton.me`.
- `.zenodo.json`: kept `upload_type: software` with top-level `license: MIT` (internally consistent), and made the split license explicit in the description (code MIT under LICENSE; paper and docs CC BY 4.0 under LICENSE-DOCS; mathematical results not copyrightable).
- Paper recompiled with `lmodern` and T1 fonts: zero Type 3 fonts, clean text extraction, 15 pages.

## v0.4.0 - 2026-06-23 - General Toeplitz theorem and corrected source-target onset

- Replaced the former universal `tau(b)` ceiling conjecture with Theorem 5.1, a proven all-`b >= 2` Toeplitz determinantal saturation theorem in the acyclic range.
- Proved that the section-induced Toeplitz maximal-minor ideal is exactly `I_Z^(b-1)`, and computed the local Buchsbaum-Rim cokernel length as `48*C(b+1,3)`.
- Corrected the global defect statement: for `a_i >= 2b+1`, the defect is `C_b` when `E >= F`, and `C_b + E - F` when `E < F`.
- Replaced the old equal-charge onset by `rho(b) = min{d >= 2b+1 : (b+1)(d-1)^3 >= (b-1)(d+1)^3}`. The old `tau(b)` agrees through `b=6` but fails at `b=7`.
- Updated the paper, README, Zenodo metadata, and machine-readable onset data to reflect the corrected theorem.
- Audit hardening (release-engineering pass): added `scripts/toeplitz_defect.py` (canonical closed-form evaluator for E, F, C_b, rho, tau, and the acyclic-range defect), `scripts/verify_v04.py` (acceptance verifier: asserts the Theorem 5.1 defect is 1994 at b=7,d=20 and 2688 at b=7,d=21, asserts rho==tau for b=2..6 and rho(7)=21 != tau(7)=20, and cross-checks every shipped b=4/b=5 raw log against the formula), and `scripts/run_b7_decisive.py` (the b=7 corrected-onset runner; `--check` confirms the matrix dimensions equal E x F on any machine, and the full dense F_p rank, about 24 GB, runs on an adequate-memory host).
- The b=7 source/target dimensions were confirmed by direct construction (54872 x 55566 at d=20; 64000 x 63888 at d=21); the dense b=7 F_p rank is recorded as a pending big-RAM run in `results/b7_decisive_prediction.json`, with the predicted 1994/2688 following from Theorem 5.1. The Theorem 5.1 piecewise formula was independently reproduced over F_p at b=4 (defect 480) and b=5 (48, 158, 352, 642, 960), covering both the C_b and the C_b+E-F branches.
- README and REPRODUCIBILITY coherence pass: the equal-charge onset is now stated as rho(b) (the source-target condition E>=F) with tau(b) retained only as the old cubic-crossing reference that fails from b=7; the stale 'general-b proof reduced to two open inputs' wording was removed (Theorem 5.1 proves both inputs in the acyclic range). `reproduce.py` now runs `verify_v04.py` as a core, flint-free stage.

## v0.3.0 — 2026-06-22 — Saturation onset correction; b=5 computation; figures; CC BY relicense

Substantive correction to the deep-region threshold, plus a second divergent-case
confirmation. No previously reported *value* is retracted: the b=2 closed form, the
22-point table, the 1,728-point sweep, and the b=4 ceiling (480, three primes) stand.
What changed is the **onset** at which the ceiling is reached, which earlier drafts
stated inconsistently as `a_i ≥ 2b`.

### Correction (the published threshold was wrong)
- The deep-region ceiling `48·C(b+1,3)` is reached at the **cubic crossing**
  `τ(b) = (2b−1) + ⌈(8·b·(b−1))^(1/3)⌉ = 6, 9, 12, 15` for b = 2,3,4,5, **not** at
  `a_i ≥ 2b` (Conjecture 5.2 had `2b`; internally inconsistent, since the proven b=2
  ceiling needs `a_i ≥ 6`, not `2b = 4`) and **not** at `a_i ≥ 2b+1` (skyscraper
  acyclicity). Between `2b` and `τ(b)` the defect lies in a rising sub-ceiling regime close to
  the cubic interior reference, with small boundary/H^1 corrections. Exact `F_p` ladders, recorded or recomputed this release: b=2 `3,24,48,…`;
  b=3 `4,32,106,192,192`; b=4 `296` at depth 11; b=5 depths 11-13 values `48,158,352` are stored results; depth 14 has raw two-prime logs, defect `642` (sub-ceiling); depth 15 has raw two-prime logs, defect `960` (ceiling); depth 16 is a stored prior plateau check, not part of the raw-log archive.
- Remark 5.4 rewritten: the Eagon–Northcott/Buchsbaum–Rim passage identifies the expected
  local-length mechanism, but for `b≥3` the equality `length(coker φ) = 48·C(b+1,3)` (constant
  in a) remains conditional on the expected-codimension/genericity input; the **defect** equals `h^0`
  of that cokernel only after a second vanishing — `H^1` of the rank-2 kernel and
  rank-(b−1) image sheaves under the `O(a−2)` twist — which is the remaining vanishing input at `τ(b)`, strictly
  above the `2b+1` skyscraper bound. The conjecture now has **two** open inputs
  (expected codimension; the kernel/image `H^1` vanishing at τ(b)), not one.

### b=5 (second divergent case)
- New Remark 5.7: defect = **960** at (15,15,15) and (16,16,16) (two primes), matching
  `48·C(6,3)`, disagreeing with the factorial value `6!·2³ = 5760`. Onset at depth 15 = τ(5).
- Abstract, Section 5 intro, finite-field-test remark, and the open-problems list updated
  to carry b=5 and τ(b).

### Code and reproducibility
- New `scripts/saturation_threshold.py`: exact-`F_p` defect ladders and onset table;
  `--gate` asserts onset == τ(b) for b=2,3 and that 2b, 2b+1 are sub-ceiling.
- `reproduce.py`: added the onset-gate stage; `run_stage` accepts script arguments;
  the onset gate is reported SKIP (not failure) when `python-flint` is absent, so
  `pip install -r requirements.txt && python3 reproduce.py` exits 0 on any machine.
- `requirements.txt`: relaxed to lower bounds (`numpy>=1.21`, `scipy>=1.7`) with
  `python-flint>=0.6.0` listed for the exact-`F_p` onset gate and the b=5 runs.
- New `results/saturation_onset.json`: ladders and τ(b) with per-value provenance.
- `scripts/saturation_threshold.py` gained a `--bd B D` flag: recomputes the defect at
  any single charge `(b,d)` (closed form at b=2, exact `F_p` otherwise) and prints its
  provenance and status against the ceiling, interior law, and τ(b).

### Figures and significance
- Three figures added to the paper: a mechanism schematic (Koszul map → Serre-dual
  contraction φ → codimension-3 degeneracy locus → local Artinian length `48·C(b+1,3)`);
  a two-panel ladder figure (b=3 and b=5) showing the defect tracking the interior cubic
  and locking onto the ceiling exactly at τ(b), with measured `F_p` markers; and an onset
  figure plotting τ(b) against the false candidates 2b and 2b+1 across b=2..7.
- New "Significance and transferability" section: frames the result as a mechanistic
  derivation of a wall (proven at b=2; for b≥3 the local length is conditional on expected codimension and
  the global identification is a program), its use as a benchmark for fitted or machine-learned cohomology
  formulae (the Thom–Porteous vs factorial divergence at b=4 and b=5), and the
  string-model-building and degeneracy-locus context. Adds the
  Brodie–Constantin–Deen–Lukas machine-learning reference (arXiv:1906.08730).
- PDF metadata (`pdftitle`/`pdfauthor`/`pdfsubject`/`pdfkeywords`) and a
  "Code and data availability" section added.
- PDF recompiled: now 13 pages.

### Licensing
- Documentation relicensed from "all rights reserved" to **CC BY 4.0**: new `LICENSE-DOCS`
  (full legal code) and `NOTICE` (the MIT-code / CC-BY-docs split, the
  math-not-copyrightable disclaimer, and the preferred citation). `CITATION.cff` records the software repository as MIT and provides a preferred article
  citation marked CC-BY-4.0. Repository code (`scripts/`, `reproduce.py`) remains MIT;
  paper, documentation, README, CHANGELOG, and figures remain CC BY 4.0. README licensing
  section and file tree updated.

## v0.2.0 — 2026-06-21 — Public-release pass: reproducibility fix, proof-gap closures, metadata

Public first release. No quantitative result changed: the b=2 closed form, the
1,728-point sweep, the table, and the b=4 finite-field ceiling check (defect = 480, three
primes) are all unchanged and independently re-verified. This pass closes two
small proof-writeup gaps, fixes a reproducibility defect in the b=4 computation script,
and makes the release metadata internally consistent.

### Reproducibility
- `scripts/run_decisive_flint.py` rewritten to be memory-frugal. The previous
  version densified the matrix and built a multi-million-entry Python list
  (`A.todense().flatten().tolist()`), which exhausts memory at the large-depth
  depth (depth 13, 8640 x 8232). The rank routine now fills a pre-zeroed
  `flint.nmod_mat` directly from the sparse COO triples; depth 13 runs in well
  under 1 GB. A `--depths`/`--primes` CLI is added; defaults run depths 11,12,13
  across the three primes. Validated against the proven b=2 formula and confirmed
  to reproduce defect = 480 at the large b=4 depths.
- `README.md` no longer states that `reproduce.py` reproduces every number. It
  now states precisely what the default chain checks (b=2 closed form and
  measured ranks, identities, sweep, table, b=4-engine validation) and that the
  b=4 ceiling itself is reproduced by `scripts/run_decisive_flint.py`.

### Paper (paper/tetraquadric_cohomology_wall.tex; PDF recompiled, 9 pages)
- Proposition (Dichotomy), 2+2 item: the step "empty base locus ⇒ defect 0" is
  completed. Emptiness gives surjectivity of the sheaf map; surjectivity on H^0
  (hence defect 0) is now derived from the vanishing of H^1 of the kernel sheaf,
  which holds for all a_i ≥ 4 on the surface (P^1)^2 by Künneth non-negativity,
  the same mechanism as Proof 2 of the ceiling theorem with one fewer Koszul step.
- Lemma (Slices of a generic section): an explicit transversality/reducedness
  line is added, recording that for generic s the 48 base points are reduced lci
  points with regular local rings. This is the input the thickening-length count
  in the saturation remark already relied upon.
- Saturation remark: "single missing input" softened to "main missing input",
  since granting expected codimension still leaves the Eagon–Northcott /
  Buchsbaum–Rim cohomology passage to be written.
- Removed the uncited Griffiths–Harris bibliography entry.

### Metadata and packaging
- `CITATION.cff`: marked `type: software` and annotated that the declared MIT
  license is the code license; the paper PDF is not part of a software deposit.
- `README.md` Licensing: a Zenodo software DOI from this repo covers code only;
  to include the paper in a citable deposit, relicense it (e.g. CC-BY-4.0) first.
- Added `.gitignore` (`__pycache__/`, `*.pyc`); removed a stray
  `scripts/__pycache__/` that had been packaged into the v1.4 tarball.
- `checksums.sha256` and `MANIFEST.md` regenerated against the v0.2.0 file set.
- Version label moved from v1.4 (internal) to v0.2.0 for public release.

## v1.4 — 2026-06-12 — Correction: Saturation Principle restored to open status

This revision corrects a single epistemological defect in the v1.3 paper. The
repository metadata (`README.md`) was already
honest — it labelled the general-`b` ceiling identification as **open**
(Conjecture 5.3) — but the TeX source `paper/tetraquadric_cohomology_wall.tex`
had **upgraded that statement from open to proven without a complete argument**.
The paper and the repository now agree.

### What was wrong (v1.3 paper)

A `[proven]` *Theorem* "Saturation Principle" asserted, for all `a_i ≥ 2b`,
`def(D) = 48·C(b+1,3) = deg D_{b-2}(φ_s)`. Its short proof silently assumed the
conclusion at three points:

1. **Genericity within the realised family.** The Thom–Porteous degree formula
   (Prop. "Thom–Porteous degree") holds *only if* `D_{b-2}(φ_s)` has the
   expected codimension 3. The maps `φ_s` range only over the Toeplitz
   subvariety realised by a single section `s`, not over all bundle maps, so
   this hypothesis is not automatic for `b ≥ 3`. It was assumed, not proven.
2. **Wrong resolution.** For `b ≥ 3` the contraction is genuinely partial; the
   cokernel is resolved by Eagon–Northcott / Buchsbaum–Rim, not by the Koszul
   complex used in the `b = 2` "Proof 2". The squeeze identifying the cokernel
   with `H⁰` of a structure sheaf was not re-derived through the correct terms.
3. **A nonexistent lemma.** The proof of "Adjugate annihilation" cited a
   "Lemma B" that does not appear in the paper and sketched, without proof, the
   claim that the maximal minors of `φ_s` *generate* `I_Z^{b-1}` (a
   determinantal-to-symbolic-power statement equivalent in strength to the
   conjecture itself). Only the easy containment direction is actually true.

The numerology was never in question: `48·C(b+1,3) = 48·[(b−1)²+C(b−1,3)]`
identically in `b`, and the value is confirmed computationally at `b = 4`
(480, three primes). The defect was purely in the *proof status* of the
general-`b` identification.

### What changed (v1.4 paper)

- The abstract no longer claims a general-`b` proof and drops a dangling
  "(Theorem 7.1)" reference; it scopes the proven contribution to `b = 2` plus
  a falsifiable conjecture for `b ≥ 3`.
- "Corollary (Ceiling identification, formerly Conjecture 5.3)" `[proven]` →
  **Conjecture (Ceiling identification)** `[open]`.
- "Theorem (Saturation Principle)" `[proven]` and its proof are **removed**;
  the statement now lives only as the conjecture above.
- "Proposition (Adjugate annihilation)" `[proven]` → **Lemma (Fitting
  annihilation; one-directional)** `[proven]`, stating and proving only the two
  facts that are genuinely true: (1) the maximal minors annihilate the cokernel
  (Fitting), and (2) `Fitt(φ_s) ⊆ I_Z^{b-1}`, i.e. `D_{b-2}(φ_s) ⊇ V(I_Z^{b-1})`.
  The nonexistent "Lemma B" and the unproven generation claim are gone.
- A new **Remark (Saturation Principle: status, obstruction, strategy)** records
  exactly what is proven unconditionally, isolates the single missing input
  (expected codimension of `D_{b-2}(φ_s)` within the section-induced family),
  and shows that granting it the conjecture follows via Thom–Porteous +
  Eagon–Northcott + the equal-length containment.
- The "Thom–Porteous degree" proposition keeps its `[proven]` status as a class
  computation but now states that its expected-dimension hypothesis is itself
  the open point for `b ≥ 3`.
- The "Status of claims" ledger is corrected to match the body.

No code, no results, and no numbers changed. `reproduce.py` and all scripts are
untouched; the verification table and the `b = 4` finite-field result are unchanged.

### Verification of this revision

- LaTeX recompiles cleanly (9 pages, 0 undefined references, 0 errors).
- No `\ref` is left dangling; every reference resolves to a label.
- `sha256sum -c checksums.sha256` passes against the files in this archive.
