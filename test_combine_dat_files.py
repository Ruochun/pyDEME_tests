#!/usr/bin/env python3
"""
Tests for combine_dat_files.py
"""

import os
import sys
import tempfile
import traceback
from pathlib import Path

try:
    from combine_dat_files import combine_dat_files, read_dat_lines
except ImportError as e:
    print(f"Error: Could not import combine_dat_files module: {e}")
    print("Make sure combine_dat_files.py is in the same directory.")
    sys.exit(1)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def write_file(path: Path, content: str) -> None:
    path.write_text(content, encoding="utf-8")


def run_test(name: str, fn) -> bool:
    try:
        fn()
        print(f"  PASS  {name}")
        return True
    except Exception:
        print(f"  FAIL  {name}")
        traceback.print_exc()
        return False


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_read_dat_lines_skips_comments():
    with tempfile.TemporaryDirectory() as tmp:
        p = Path(tmp) / "data.dat"
        write_file(p, "# comment\n1.0 2.0 3.0\n# another comment\n4.0 5.0 6.0\n")
        lines = read_dat_lines(p)
        assert lines == ["1.0 2.0 3.0", "4.0 5.0 6.0"], f"Got: {lines}"


def test_read_dat_lines_skips_empty_lines():
    with tempfile.TemporaryDirectory() as tmp:
        p = Path(tmp) / "data.dat"
        write_file(p, "\n1.0 2.0 3.0\n\n4.0 5.0 6.0\n\n")
        lines = read_dat_lines(p)
        assert lines == ["1.0 2.0 3.0", "4.0 5.0 6.0"], f"Got: {lines}"


def test_basic_combine():
    with tempfile.TemporaryDirectory() as tmp:
        coords = Path(tmp) / "coords.dat"
        radii = Path(tmp) / "radii.dat"
        output = Path(tmp) / "output.dat"

        write_file(coords, "1.0 2.0 3.0\n4.0 5.0 6.0\n")
        write_file(radii, "0.5\n1.5\n")

        combine_dat_files(coords, radii, output)

        lines = output.read_text().splitlines()
        assert lines[0] == "0 1.0 2.0 3.0 0.5", f"Got: {lines[0]!r}"
        assert lines[1] == "1 4.0 5.0 6.0 1.5", f"Got: {lines[1]!r}"


def test_combine_with_comment_lines():
    with tempfile.TemporaryDirectory() as tmp:
        coords = Path(tmp) / "coords.dat"
        radii = Path(tmp) / "radii.dat"
        output = Path(tmp) / "output.dat"

        write_file(coords, "# xyz file\n1.0 2.0 3.0\n4.0 5.0 6.0\n")
        write_file(radii, "# radius file\n0.5\n1.5\n")

        combine_dat_files(coords, radii, output)

        lines = output.read_text().splitlines()
        assert len(lines) == 2
        assert lines[0] == "0 1.0 2.0 3.0 0.5"
        assert lines[1] == "1 4.0 5.0 6.0 1.5"


def test_index_is_zero_based():
    with tempfile.TemporaryDirectory() as tmp:
        coords = Path(tmp) / "coords.dat"
        radii = Path(tmp) / "radii.dat"
        output = Path(tmp) / "output.dat"

        write_file(coords, "0.1 0.2 0.3\n0.4 0.5 0.6\n0.7 0.8 0.9\n")
        write_file(radii, "1.0\n2.0\n3.0\n")

        combine_dat_files(coords, radii, output)

        lines = output.read_text().splitlines()
        assert lines[0].startswith("0 ")
        assert lines[1].startswith("1 ")
        assert lines[2].startswith("2 ")


def test_mismatched_counts_raises():
    with tempfile.TemporaryDirectory() as tmp:
        coords = Path(tmp) / "coords.dat"
        radii = Path(tmp) / "radii.dat"
        output = Path(tmp) / "output.dat"

        write_file(coords, "1.0 2.0 3.0\n4.0 5.0 6.0\n")
        write_file(radii, "0.5\n")  # only one radius for two coords

        try:
            combine_dat_files(coords, radii, output)
            assert False, "Expected ValueError was not raised"
        except ValueError:
            pass  # expected


def test_single_point():
    with tempfile.TemporaryDirectory() as tmp:
        coords = Path(tmp) / "coords.dat"
        radii = Path(tmp) / "radii.dat"
        output = Path(tmp) / "output.dat"

        write_file(coords, "10.0 20.0 30.0\n")
        write_file(radii, "5.0\n")

        combine_dat_files(coords, radii, output)

        lines = output.read_text().splitlines()
        assert len(lines) == 1
        assert lines[0] == "0 10.0 20.0 30.0 5.0"


# ---------------------------------------------------------------------------
# Runner
# ---------------------------------------------------------------------------

def main():
    tests = [
        ("read_dat_lines skips #-comment lines", test_read_dat_lines_skips_comments),
        ("read_dat_lines skips empty lines", test_read_dat_lines_skips_empty_lines),
        ("basic combine produces correct output", test_basic_combine),
        ("combine with comment lines in both files", test_combine_with_comment_lines),
        ("output index is 0-based", test_index_is_zero_based),
        ("mismatched line counts raises ValueError", test_mismatched_counts_raises),
        ("single-point files", test_single_point),
    ]

    passed = 0
    failed = 0
    for name, fn in tests:
        if run_test(name, fn):
            passed += 1
        else:
            failed += 1

    print(f"\n{passed} passed, {failed} failed out of {len(tests)} tests.")
    if failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
