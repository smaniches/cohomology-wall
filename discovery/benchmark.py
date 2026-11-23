#!/usr/bin/env python3
"""Bounded, stdlib-only abstraction-discovery benchmark for cohomology-wall.

The local oracle constructs and eliminates the graded Toeplitz presentation over
F_101. The separate global oracle evaluates the PUBLISHED acyclic-range theorem;
it is not an independent matrix-rank measurement or proof oracle.

No model inference, untrusted-code execution, network access, or dependency on
an archived results file is required.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from itertools import product
from math import comb, prod
from pathlib import Path
import random

P = 101
IDENTITY = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
CHANGE = ((1, 2, 0), (0, 1, 1), (1, 0, 1))
CHANGE2 = ((2, 0, 0), (0, 3, 0), (0, 0, 5))
SINGULAR = ((1, 0, 0), (0, 1, 0), (1, 1, 0))
RANK_ONE = ((1, 0, 0), (1, 0, 0), (1, 0, 0))
MAX_N, MAX_B, MAX_A = 8, 12, 64


def _integer(value: object, name: str, low: int, high: int) -> int:
    if type(value) is not int or not low <= value <= high:
        raise ValueError(f"{name} must be an integer in [{low},{high}]")
    return value


def _form_matrix(value: object) -> tuple[tuple[int, ...], ...]:
    if not isinstance(value, (tuple, list)) or len(value) != 3:
        raise ValueError("forms must be a 3-by-3 integer matrix")
    out = []
    for row in value:
        if not isinstance(row, (tuple, list)) or len(row) != 3:
            raise ValueError("forms must be a 3-by-3 integer matrix")
        out.append(tuple(_integer(x, "form coefficient", -100, 100) % P for x in row))
    return tuple(out)


def monomials(d: int) -> list[tuple[int, int, int]]:
    if d < 0:
        return []
    return [(x, y, d - x - y) for x in range(d + 1) for y in range(d - x + 1)]


def modular_column_rank(columns: list[dict[int, int]]) -> int:
    """Exact column rank over F_101, with normalized sparse pivot columns."""
    pivots: dict[int, dict[int, int]] = {}
    for original in columns:
        vec = {r: v % P for r, v in original.items() if v % P}
        while vec:
            k = min(vec)
            if k not in pivots:
                inv = pow(vec[k], -1, P)
                pivots[k] = {r: (v * inv) % P for r, v in vec.items()}
                break
            multiple = vec[k]
            for r, v in pivots[k].items():
                result = (vec.get(r, 0) - multiple * v) % P
                if result:
                    vec[r] = result
                else:
                    vec.pop(r, None)
    return len(pivots)


def local_hilbert(n: int, degree: int, forms: object = IDENTITY) -> int:
    """Compute dim coker[T_n(L0,L1,L2): S(-1)^(n+2) -> S^n]_degree.

    T_ij = L_(j-i) for 0 <= j-i <= 2, otherwise zero, and each L_k is
    the linear form in x,y,z given by forms[k]. Computes degree components
    through finite-field Gaussian elimination, NOT a Hilbert-series formula.
    """
    n = _integer(n, "n", 1, MAX_N)
    degree = _integer(degree, "degree", 0, MAX_N + 2)
    forms = _form_matrix(forms)
    source_monomials = monomials(degree - 1)
    target_monomials = monomials(degree)
    index = {m: k for k, m in enumerate(target_monomials)}
    columns: list[dict[int, int]] = []
    for j in range(n + 2):
        for mon in source_monomials:
            column: dict[int, int] = {}
            for i in range(n):
                t = j - i
                if t < 0 or t > 2:
                    continue
                for variable, coefficient in enumerate(forms[t]):
                    if not coefficient:
                        continue
                    changed = list(mon)
                    changed[variable] += 1
                    row = i * len(target_monomials) + index[tuple(changed)]
                    column[row] = (column.get(row, 0) + coefficient) % P
            columns.append(column)
    target_dim = n * len(target_monomials)
    return target_dim - modular_column_rank(columns)


def global_defect(b: int, a: object) -> int:
    """Theorem 5.1 reference, acyclic range ONLY; not a computed matrix rank."""
    b = _integer(b, "b", 2, MAX_B)
    if not isinstance(a, (list, tuple)) or len(a) != 3:
        raise ValueError("a must contain exactly three integer charges")
    charges = tuple(_integer(x, "charge", 2 * b + 1, MAX_A) for x in a)
    e = (b + 1) * prod(x - 1 for x in charges)
    f = (b - 1) * prod(x + 1 for x in charges)
    return 48 * comb(b + 1, 3) + min(0, e - f)


def onset(b: int) -> int:
    b = _integer(b, "b", 2, MAX_B)
    for d in range(2 * b + 1, MAX_A + 1):
        if (b + 1) * (d - 1) ** 3 >= (b - 1) * (d + 1) ** 3:
            return d
    raise ValueError("the onset is outside the supported charge range")


def answer(case: dict) -> int:
    family = case["family"]
    args = case["inputs"]
    if family == "local":
        return local_hilbert(args["n"], args["degree"], args["forms"])
    if family == "global":
        return global_defect(args["b"], args["a"])
    raise ValueError(f"unknown family: {family}")


def canonical(obj: object) -> bytes:
    return (json.dumps(obj, sort_keys=True, separators=(",", ":"),
                       ensure_ascii=True, allow_nan=False) + "\n").encode("utf-8")


def digest(obj: object) -> str:
    return hashlib.sha256(canonical(obj)).hexdigest()


def _case(family: str, regime: str, inputs: dict) -> dict:
    key = {"family": family, "regime": regime, "inputs": inputs}
    return {"id": hashlib.sha256(canonical(key)).hexdigest()[:16], **key}


def make_cases(seed: int) -> tuple[list[dict], list[dict]]:
    """Group-disjoint parameter holdout; no challenge answer is serialized."""
    rng = random.Random(seed)
    train: list[dict] = []
    challenge: list[dict] = []
    for n in (1, 2, 3):
        for d in range(n + 2):
            train.append(_case("local", "local-training", {
                "n": n, "degree": d, "forms": IDENTITY}))
    # Changed coordinate systems with unobserved n, plus a rank-deficient foil.
    for n in (4, 5, 6):
        for d in range(n + 2):
            challenge.append(_case("local", "local-new-n", {
                "n": n, "degree": d, "forms": IDENTITY}))
            challenge.append(_case("local", "local-coordinate-change", {
                "n": n, "degree": d, "forms": CHANGE if n % 2 else CHANGE2}))
    for n in (1, 2, 3):
        for d in (n, n + 1):
            challenge.append(_case("local", "local-degenerate-foil", {
                "n": n, "degree": d, "forms": SINGULAR if d == n else RANK_ONE}))
    for b in (2, 3, 4):
        low = 2 * b + 1
        for d in sorted(set((low, onset(b) - 1, onset(b), onset(b) + 2))):
            if d >= low:
                train.append(_case("global", "global-training", {
                    "b": b, "a": (d, d, d)}))
        train.append(_case("global", "global-training", {
            "b": b, "a": (low, low + 1, low + 2)}))
    for b in (5, 6, 7, 8):
        low = 2 * b + 1
        for d in sorted(set((low, onset(b) - 1, onset(b), onset(b) + 1))):
            if d >= low:
                challenge.append(_case("global", "global-new-b", {
                    "b": b, "a": (d, d, d)}))
        challenge.append(_case("global", "global-asymmetric", {
            "b": b, "a": (low, low + 2, onset(b) + 3)}))
        challenge.append(_case("global", "global-asymmetric", {
            "b": b, "a": (onset(b), onset(b) + 2, onset(b) + 4)}))
    # Challenge cases may be shuffled. The freeze should precede seed selection.
    rng.shuffle(train)
    rng.shuffle(challenge)
    assert len({x["id"] for x in train + challenge}) == len(train) + len(challenge)
    return [{**x, "answer": answer(x)} for x in train], challenge


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(json.dumps(value, indent=2, sort_keys=True,
                                ensure_ascii=True, allow_nan=False).encode("utf-8") + b"\n")


def _unique_pairs(pairs):
    result = {}
    for k, v in pairs:
        if k in result:
            raise ValueError(f"duplicate JSON key: {k}")
        result[k] = v
    return result


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"),
                      object_pairs_hook=_unique_pairs,
                      parse_constant=lambda v: (_ for _ in ()).throw(ValueError(f"nonfinite JSON: {v}")))


def score(challenge: list[dict], predictions: list[dict]) -> dict:
    if not isinstance(challenge, list) or not isinstance(predictions, list):
        raise ValueError("challenge and predictions must be JSON lists")
    cases = {}
    for case in challenge:
        if not isinstance(case, dict) or set(case) != {"id", "family", "regime", "inputs"}:
            raise ValueError("challenge must contain unlabeled well-formed cases")
        if case["id"] != _case(case["family"], case["regime"], case["inputs"])["id"]:
            raise ValueError("challenge case ID does not match its content")
        if case["id"] in cases:
            raise ValueError("duplicate challenge case")
        cases[case["id"]] = case
    seen = {}
    for item in predictions:
        if not isinstance(item, dict) or set(item) != {"id", "prediction"}:
            raise ValueError("each prediction must contain exactly id and prediction")
        cid = item["id"]
        if cid not in cases or cid in seen:
            raise ValueError("unknown or duplicate prediction id")
        value = item["prediction"]
        if type(value) is not int or value < 0:
            raise ValueError("predictions must be nonnegative integers")
        seen[cid] = value
    if set(seen) != set(cases):
        raise ValueError(f"missing {len(set(cases) - set(seen))} predictions")
    groups: dict[str, dict] = {}
    incorrect = []
    for cid, case in cases.items():
        expected = answer(case)
        group = case["regime"]
        counts = groups.setdefault(group, {"correct": 0, "total": 0})
        counts["total"] += 1
        if seen[cid] == expected:
            counts["correct"] += 1
        else:
            incorrect.append({"id": cid, "regime": group,
                              "predicted": seen[cid]})
    return {"numerical_status": "PASS" if not incorrect else "FAIL",
            "mechanism_assessment": "NOT_EVALUATED",
            "independent_discovery": "NOT_ESTABLISHED",
            "correct": sum(g["correct"] for g in groups.values()),
            "total": len(cases), "by_regime": groups,
            "incorrect": incorrect, "challenge_sha256": digest(challenge)}


def validate_manifest(challenge: list[dict], manifest: dict) -> None:
    """Rebuild declared deterministic challenge before comparing frozen IDs."""
    if not isinstance(manifest, dict) or manifest.get("benchmark") != "cohomology-abstraction-discovery-v1":
        raise ValueError("unsupported benchmark manifest")
    if type(manifest.get("seed")) is not int:
        raise ValueError("manifest seed must be an integer")
    expected_train, expected_challenge = make_cases(manifest["seed"])
    if (manifest.get("challenge_sha256") != digest(challenge)
            or manifest.get("challenge_sha256") != digest(expected_challenge)
            or manifest.get("train_sha256") != digest(expected_train)
            or manifest.get("challenge_cases") != len(challenge)
            or manifest.get("train_cases") != len(expected_train)
            or canonical(challenge) != canonical(expected_challenge)):
        raise ValueError("challenge disagrees with its deterministic manifest")


def _lagrange(x: int, points: list[tuple[int, int]]) -> Fraction:
    return sum((Fraction(y) * prod((Fraction(x - a, b - a)
                                    for b, _ in points if b != a), start=Fraction(1))
                for a, y in points), Fraction(0))


def simple_baseline(train: list[dict], challenge: list[dict]) -> list[dict]:
    """Memorized local table + degree-2 global ceiling extrapolation.

    The global baseline is intentionally competitive with dimension-aware
    interpolation: it uses exact E-F from the published problem definition,
    but learns its ceiling on b=2,3,4 solely from training examples.
    """
    local = {(x["inputs"]["n"], x["inputs"]["degree"]): x["answer"]
             for x in train if x["family"] == "local"}
    caps: dict[int, int] = {}
    for x in train:
        if x["family"] == "global":
            b = x["inputs"]["b"]
            caps[b] = max(caps.get(b, 0), x["answer"])
    if sorted(caps) != [2, 3, 4]:
        raise ValueError("baseline expects b=2,3,4 in training")
    fit = sorted(caps.items())
    result = []
    for x in challenge:
        q = x["inputs"]
        if x["family"] == "local":
            value = local.get((q["n"], q["degree"]), 0)
        else:
            b, a = q["b"], q["a"]
            e = (b + 1) * prod(v - 1 for v in a)
            f = (b - 1) * prod(v + 1 for v in a)
            fitted = _lagrange(b, fit)
            value = max(0, int(fitted) + min(0, e - f))
        result.append({"id": x["id"], "prediction": value})
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    oracle = sub.add_parser("oracle", help="Query the bounded local arithmetic or published global reference")
    oracle.add_argument("--family", choices=("local", "global"), required=True)
    oracle.add_argument("--n", type=int)
    oracle.add_argument("--degree", type=int)
    oracle.add_argument("--forms", default="identity", choices=("identity", "change", "singular", "rank-one"))
    oracle.add_argument("--b", type=int)
    oracle.add_argument("--a", help="three comma-separated charges")
    gen = sub.add_parser("generate", help="Generate train labels and unlabeled challenge")
    gen.add_argument("--seed", type=int, required=True)
    gen.add_argument("--out", type=Path, required=True)
    baseline = sub.add_parser("baseline", help="Generate deliberately naive predictions from the training set")
    baseline.add_argument("--train", type=Path, required=True)
    baseline.add_argument("--challenge", type=Path, required=True)
    baseline.add_argument("--out", type=Path, required=True)
    evaluate = sub.add_parser("score", help="Score complete prediction JSON without executing user code")
    evaluate.add_argument("--challenge", type=Path, required=True)
    evaluate.add_argument("--predictions", type=Path, required=True)
    evaluate.add_argument("--manifest", type=Path, required=True)
    evaluate.add_argument("--out", type=Path)
    args = parser.parse_args()
    try:
        if args.command == "oracle":
            if args.family == "local":
                if args.n is None or args.degree is None:
                    raise ValueError("local oracle needs --n and --degree")
                forms = {"identity": IDENTITY, "change": CHANGE,
                         "singular": SINGULAR, "rank-one": RANK_ONE}[args.forms]
                val = local_hilbert(args.n, args.degree, forms)
                response = {"family": "local", "value": val,
                            "source": "computed exact F_101 matrix rank"}
            else:
                if args.b is None or args.a is None:
                    raise ValueError("global oracle needs --b and --a")
                val = global_defect(args.b, tuple(int(v) for v in args.a.split(",")))
                response = {"family": "global", "value": val,
                            "source": "derived from published Theorem 5.1 (NOT measured rank)"}
            print(json.dumps(response, sort_keys=True))
        elif args.command == "generate":
            train, challenge = make_cases(args.seed)
            write_json(args.out / "train.json", train)
            write_json(args.out / "challenge.json", challenge)
            write_json(args.out / "manifest.json", {
                "benchmark": "cohomology-abstraction-discovery-v1",
                "seed": args.seed, "train_cases": len(train),
                "challenge_cases": len(challenge),
                "train_sha256": digest(train),
                "challenge_sha256": digest(challenge),
                "caveat": "self-service protocol blinding only; published theorem and oracle source are accessible"})
            print(f"train={len(train)} challenge={len(challenge)} out={args.out}")
        elif args.command == "baseline":
            pred = simple_baseline(read_json(args.train), read_json(args.challenge))
            write_json(args.out, pred)
            print(f"baseline predictions={len(pred)} out={args.out}")
        elif args.command == "score":
            challenge = read_json(args.challenge)
            validate_manifest(challenge, read_json(args.manifest))
            result = score(challenge, read_json(args.predictions))
            if args.out:
                write_json(args.out, result)
            print(json.dumps(result, sort_keys=True, indent=2))
    except (ValueError, KeyError, TypeError, OSError, json.JSONDecodeError) as exc:
        parser.exit(2, f"benchmark error: {exc}\n")


if __name__ == "__main__":
    main()
