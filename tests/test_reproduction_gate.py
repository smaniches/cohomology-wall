"""Regression tests for the reproduction orchestrator's fail-closed stage verdicts.

These tests execute temporary child processes. They do not recompute or alter
the historical mathematical experiment or any archived result.
"""

from __future__ import annotations

import contextlib
import io
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import reproduce


class ReproductionStageGateTests(unittest.TestCase):
    def run_fixture(
        self,
        script: str,
        *,
        markers: tuple[str, ...] = ("CHECK PASSED",),
        live: bool = False,
        timeout: float = 5.0,
    ) -> bool:
        with tempfile.TemporaryDirectory() as directory:
            Path(directory, "stage.py").write_text(script, encoding="utf-8")
            with patch.object(reproduce, "ROOT", directory):
                with contextlib.redirect_stdout(io.StringIO()):
                    return reproduce.run_stage(
                        "synthetic gate fixture",
                        "stage.py",
                        list(markers),
                        live=live,
                        timeout=timeout,
                    )

    def test_zero_exit_with_expected_markers_passes(self) -> None:
        self.assertTrue(
            self.run_fixture("print('CHECK PASSED'); print('SECOND MARKER')\n",
                             markers=("CHECK PASSED", "SECOND MARKER"))
        )

    def test_nonzero_exit_after_expected_markers_fails(self) -> None:
        # Previously gave a false PASS when exit_matters=False.
        self.assertFalse(
            self.run_fixture("print('CHECK PASSED'); raise SystemExit(7)\n")
        )

    def test_missing_marker_fails_even_with_zero_exit(self) -> None:
        self.assertFalse(self.run_fixture("print('UNRELATED OUTPUT')\n"))

    def test_explicit_failure_marker_fails_even_with_zero_exit(self) -> None:
        self.assertFalse(self.run_fixture("print('CHECK PASSED'); print('FAIL')\n"))

    def test_timeout_fails(self) -> None:
        self.assertFalse(
            self.run_fixture(
                "import time\nprint('CHECK PASSED', flush=True)\ntime.sleep(5)\n",
                timeout=0.1,
            )
        )

    def test_live_mode_requires_zero_exit(self) -> None:
        self.assertFalse(
            self.run_fixture("print('CHECK PASSED'); raise SystemExit(8)\n", live=True)
        )


if __name__ == "__main__":
    unittest.main()
