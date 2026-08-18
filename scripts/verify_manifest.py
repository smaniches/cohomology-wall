#!/usr/bin/env python3
"""Verify the checksum registry and human-readable MANIFEST as one invariant."""

from __future__ import annotations

import hashlib
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CHECKSUMS = ROOT / "checksums.sha256"
MANIFEST = ROOT / "MANIFEST.md"
EXCLUDED = {"MANIFEST.md", "checksums.sha256"}
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
MANIFEST_ROW_RE = re.compile(
    r"^\| `(?P<path>[^`]+)` \| `(?P<prefix>[0-9a-f]{16})` \| (?P<size>\d+) \|$"
)


def fail(message: str) -> None:
    raise SystemExit(f"MANIFEST CHECK FAILED: {message}")


def tracked_files() -> set[str]:
    proc = subprocess.run(
        ["git", "ls-files", "-z"],
        cwd=ROOT,
        check=True,
        stdout=subprocess.PIPE,
    )
    return {
        path
        for path in proc.stdout.decode("utf-8").split("\0")
        if path and path not in EXCLUDED
    }


def parse_checksums() -> dict[str, str]:
    entries: dict[str, str] = {}
    for line_no, raw in enumerate(CHECKSUMS.read_text(encoding="utf-8").splitlines(), 1):
        if not raw.strip():
            continue
        try:
            digest, path = raw.split("  ./", 1)
        except ValueError:
            fail(f"checksums.sha256:{line_no}: malformed row")
        if not SHA256_RE.fullmatch(digest):
            fail(f"checksums.sha256:{line_no}: invalid SHA-256")
        if path in entries:
            fail(f"checksums.sha256:{line_no}: duplicate path {path}")
        entries[path] = digest
    return entries


def parse_manifest() -> dict[str, tuple[str, int]]:
    entries: dict[str, tuple[str, int]] = {}
    for line_no, raw in enumerate(MANIFEST.read_text(encoding="utf-8").splitlines(), 1):
        match = MANIFEST_ROW_RE.fullmatch(raw)
        if not match:
            continue
        path = match.group("path")
        if path in entries:
            fail(f"MANIFEST.md:{line_no}: duplicate path {path}")
        entries[path] = (match.group("prefix"), int(match.group("size")))
    return entries


def main() -> None:
    checksums = parse_checksums()
    manifest = parse_manifest()
    tracked = tracked_files()

    checksum_paths = set(checksums)
    manifest_paths = set(manifest)
    if checksum_paths != tracked:
        fail(
            "checksum registry does not match tracked files: "
            f"missing={sorted(tracked - checksum_paths)}, "
            f"extra={sorted(checksum_paths - tracked)}"
        )
    if manifest_paths != checksum_paths:
        fail(
            "MANIFEST table does not match checksum registry: "
            f"missing={sorted(checksum_paths - manifest_paths)}, "
            f"extra={sorted(manifest_paths - checksum_paths)}"
        )

    for path in sorted(checksums):
        data = (ROOT / path).read_bytes()
        digest = hashlib.sha256(data).hexdigest()
        expected_digest = checksums[path]
        prefix, expected_size = manifest[path]
        if digest != expected_digest:
            fail(f"{path}: SHA-256 mismatch")
        if prefix != digest[:16]:
            fail(f"{path}: MANIFEST SHA prefix mismatch")
        if expected_size != len(data):
            fail(f"{path}: MANIFEST byte count mismatch")

    print(f"MANIFEST CHECK PASS: {len(checksums)} tracked files verified")


if __name__ == "__main__":
    main()
