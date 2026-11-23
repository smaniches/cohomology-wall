"""Independent numerical checks and adversarial protocol tests for discovery v1."""
from __future__ import annotations

from collections import defaultdict
from itertools import product
from math import comb
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from discovery import benchmark as bm  # noqa: E402


def independent_dense_cokernel(n, degree, forms):
    """Separate ROW-reduction implementation; uses no benchmark rank routine."""
    targets = [(i, e) for i in range(n) for e in bm.monomials(degree)]
    sources = [(j, e) for j in range(n + 2) for e in bm.monomials(degree - 1)]
    mat = []
    for i, target in targets:
        row = []
        for j, source in sources:
            diff = tuple(target[v] - source[v] for v in range(3))
            k = j - i
            if k not in (0, 1, 2) or diff not in ((1, 0, 0), (0, 1, 0), (0, 0, 1)):
                row.append(0)
            else:
                row.append(forms[k][diff.index(1)] % bm.P)
        mat.append(row)
    r, c = len(mat), len(sources)
    rank = 0
    for j in range(c):
        pivot = next((i for i in range(rank, r) if mat[i][j] % bm.P), None)
        if pivot is None:
            continue
        mat[rank], mat[pivot] = mat[pivot], mat[rank]
        inv = pow(mat[rank][j], -1, bm.P)
        mat[rank] = [(v * inv) % bm.P for v in mat[rank]]
        for i in range(rank + 1, r):
            fac = mat[i][j] % bm.P
            if fac:
                mat[i] = [(u - fac * v) % bm.P
                          for u, v in zip(mat[i], mat[rank])]
        rank += 1
        if rank == r:
            break
    return r - rank


class ExactOracleTests(unittest.TestCase):
    def test_known_hilbert_lengths_and_finite_support(self):
        for n in range(1, 7):
            with self.subTest(n=n):
                vals = [bm.local_hilbert(n, d) for d in range(n + 2)]
                self.assertEqual(sum(vals), comb(n + 2, 3))
                self.assertEqual(vals[0], n)
                self.assertEqual(vals[n:], [0, 0])

    def test_linear_change_of_coordinates(self):
        for n in range(1, 6):
            for degree in range(n + 2):
                with self.subTest(n=n, degree=degree):
                    target = bm.local_hilbert(n, degree)
                    self.assertEqual(target, bm.local_hilbert(n, degree, bm.CHANGE))
                    self.assertEqual(target, bm.local_hilbert(n, degree, bm.CHANGE2))

    def test_singular_counterexamples(self):
        for forms in (bm.SINGULAR, bm.RANK_ONE):
            self.assertGreater(bm.local_hilbert(1, 3, forms), 0)
            self.assertEqual(bm.local_hilbert(1, 3, bm.IDENTITY), 0)

    def test_independent_dense_rank_matches_sparse_column_rank(self):
        for n, degree, forms in product((1, 2, 3), (0, 1, 2, 3),
                                        (bm.IDENTITY, bm.CHANGE, bm.SINGULAR)):
            with self.subTest(n=n, d=degree, forms=forms):
                self.assertEqual(bm.local_hilbert(n, degree, forms),
                                 independent_dense_cokernel(n, degree, forms))

    def test_global_reference_boundary_and_validation(self):
        self.assertEqual(bm.global_defect(7, (20, 20, 20)), 1994)
        self.assertEqual(bm.global_defect(7, (21, 21, 21)), 2688)
        self.assertEqual(bm.onset(7), 21)
        self.assertEqual(bm.onset(2), 6)
        for b in range(2, 11):
            d = bm.onset(b)
            self.assertEqual(bm.global_defect(b, (d, d, d)),
                             48 * comb(b + 1, 3))
            self.assertLess(bm.global_defect(b, (d - 1,) * 3),
                            48 * comb(b + 1, 3))
        with self.assertRaises(ValueError):
            bm.global_defect(7, (13, 13, 13))
        with self.assertRaises(ValueError):
            bm.local_hilbert(True, 1)
        with self.assertRaises(ValueError):
            bm.local_hilbert(1, -1)
        with self.assertRaises(ValueError):
            bm.local_hilbert(1, 1, [[1, 0], [0, 1], [1, 1]])
        with self.assertRaises(ValueError):
            bm.global_defect(13, (27, 27, 27))


class ProtocolTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.train, cls.challenge = bm.make_cases(20261122)

    def test_deterministic_unlabeled_split(self):
        train2, challenge2 = bm.make_cases(20261122)
        self.assertEqual(bm.digest(self.train), bm.digest(train2))
        self.assertEqual(bm.digest(self.challenge), bm.digest(challenge2))
        self.assertEqual(len(self.train), 26)
        self.assertEqual(len(self.challenge), 72)
        self.assertTrue(all("answer" in x for x in self.train))
        self.assertTrue(all("answer" not in x for x in self.challenge))
        ids = [x["id"] for x in self.train + self.challenge]
        self.assertEqual(len(ids), len(set(ids)))
        train_b = {x["inputs"]["b"] for x in self.train if x["family"] == "global"}
        test_b = {x["inputs"]["b"] for x in self.challenge if x["family"] == "global"}
        self.assertTrue(train_b.isdisjoint(test_b))
        self.assertEqual(train_b, {2, 3, 4})
        self.assertEqual(test_b, {5, 6, 7, 8})

    def test_complete_reference_passes_by_regime(self):
        reference = [{"id": x["id"], "prediction": bm.answer(x)}
                     for x in self.challenge]
        result = bm.score(self.challenge, reference)
        self.assertEqual(result["mechanism_assessment"], "NOT_EVALUATED")
        self.assertEqual(result["numerical_status"], "PASS")
        self.assertEqual(result["correct"], 72)
        self.assertEqual(result["incorrect"], [])
        self.assertEqual(set(result["by_regime"]), {
            "local-new-n", "local-coordinate-change", "local-degenerate-foil",
            "global-new-b", "global-asymmetric"})

    def test_baseline_extrapolation_does_not_pass(self):
        prediction = bm.simple_baseline(self.train, self.challenge)
        result = bm.score(self.challenge, prediction)
        self.assertEqual(result["numerical_status"], "FAIL")
        self.assertLess(result["correct"], 72)
        self.assertEqual(result["correct"], 12)

    def test_score_fail_closed_on_corruption(self):
        good = [{"id": x["id"], "prediction": bm.answer(x)}
                for x in self.challenge]
        corrupted = [dict(x) for x in good]
        corrupted[0]["prediction"] += 1
        result = bm.score(self.challenge, corrupted)
        self.assertEqual(result["numerical_status"], "FAIL")
        self.assertEqual(result["correct"], 71)
        with self.assertRaises(ValueError):
            bm.score(self.challenge, good[:-1])
        with self.assertRaises(ValueError):
            bm.score(self.challenge, good + [good[0]])
        with self.assertRaises(ValueError):
            bm.score(self.challenge, good + [{"id": "fake", "prediction": 0}])
        bad_boolean = [dict(x) for x in good]
        bad_boolean[0]["prediction"] = True
        with self.assertRaises(ValueError):
            bm.score(self.challenge, bad_boolean)
        contaminated = [dict(x) for x in self.challenge]
        contaminated[0]["answer"] = bm.answer(contaminated[0])
        with self.assertRaises(ValueError):
            bm.score(contaminated, good)
        changed_case = [dict(x) for x in self.challenge]
        changed_case[0]["inputs"] = dict(changed_case[0]["inputs"])
        changed_case[0]["inputs"]["something_new"] = 7
        with self.assertRaises(ValueError):
            bm.score(changed_case, good)

    def test_manifest_rejects_mutation(self):
        manifest = {"benchmark": "cohomology-abstraction-discovery-v1",
                    "seed": 20261122, "train_cases": len(self.train),
                    "challenge_cases": len(self.challenge),
                    "train_sha256": bm.digest(self.train),
                    "challenge_sha256": bm.digest(self.challenge)}
        bm.validate_manifest(self.challenge, manifest)
        tampered = dict(manifest)
        tampered["challenge_sha256"] = "0" * 64
        with self.assertRaises(ValueError):
            bm.validate_manifest(self.challenge, tampered)
        mutated = [dict(x) for x in self.challenge]
        mutated[0]["regime"] = "wrong"
        with self.assertRaises(ValueError):
            bm.validate_manifest(mutated, manifest)

    def test_duplicate_json_keys_refused(self):
        with tempfile.TemporaryDirectory() as directory:
            p = Path(directory) / "duplicate.json"
            p.write_text('{"a":1,"a":2}', encoding="utf-8")
            with self.assertRaises(ValueError):
                bm.read_json(p)

    def test_cli_round_trip_without_external_packages(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            script = str(ROOT / "discovery" / "benchmark.py")
            def call(*argv):
                return subprocess.run([sys.executable, script, *map(str, argv)],
                                      text=True, capture_output=True)
            self.assertEqual(call("oracle", "--family", "global", "--b", 7,
                                  "--a", "20,20,20").returncode, 0)
            self.assertEqual(call("generate", "--seed", 20261122,
                                  "--out", root).returncode, 0)
            self.assertEqual(call("baseline", "--train", root / "train.json",
                                  "--challenge", root / "challenge.json",
                                  "--out", root / "predictions.json").returncode, 0)
            self.assertEqual(call("score", "--challenge", root / "challenge.json",
                                  "--predictions", root / "predictions.json",
                                  "--manifest", root / "manifest.json",
                                  "--out", root / "score.json").returncode, 0)
            self.assertEqual(bm.read_json(root / "score.json")["correct"], 12)
            self.assertEqual(call("oracle", "--family", "global", "--b", 7,
                                  "--a", "13,13,13").returncode, 2)


if __name__ == "__main__":
    unittest.main()
