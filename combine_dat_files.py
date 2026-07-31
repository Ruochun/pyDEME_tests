#!/usr/bin/env python3
"""
Combine two dat files into one.

One input dat file has xyz coordinates of points (space separated, one point
per line). The other input dat file has the radius of each point (one number
per line). Both files may contain #-marked comment lines which are skipped.

Output format (space separated, one point per line):
  <0-based index> <x> <y> <z> <radius>

Example:
  python combine_dat_files.py coords.dat radii.dat output.dat
"""

from __future__ import annotations

import argparse
from pathlib import Path


def read_dat_lines(path: Path) -> list[str]:
    """Return non-comment, non-empty lines from a dat file."""
    lines = []
    with path.open(encoding="utf-8") as f:
        for line in f:
            stripped = line.strip()
            if stripped and not stripped.startswith("#"):
                lines.append(stripped)
    return lines


def combine_dat_files(coords_file: Path, radii_file: Path, output_file: Path) -> None:
    """Combine xyz coords and radius dat files into a single output dat file."""
    coord_lines = read_dat_lines(coords_file)
    radius_lines = read_dat_lines(radii_file)

    if len(coord_lines) != len(radius_lines):
        raise ValueError(
            f"Number of data points does not match: "
            f"{len(coord_lines)} coords vs {len(radius_lines)} radii"
        )

    with output_file.open("w", encoding="utf-8") as out:
        for i, (coord_line, radius_line) in enumerate(zip(coord_lines, radius_lines)):
            parts = coord_line.split()
            if len(parts) < 3:
                raise ValueError(
                    f"Line {i} of coords file has fewer than 3 values: {coord_line!r}"
                )
            x, y, z = parts[0], parts[1], parts[2]
            r = float(radius_line.strip()) # corrected, the input is now indeed radius
            out.write(f"{i} {x} {y} {z} {r}\n")


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Combine xyz-coords dat file and radius dat file into one."
    )
    ap.add_argument("coords_dat", type=Path, help="Input dat file with x y z (one point per line)")
    ap.add_argument("radii_dat", type=Path, help="Input dat file with radius (one value per line)")
    ap.add_argument("output_dat", type=Path, help="Output dat file (index x y z radius)")
    args = ap.parse_args()

    combine_dat_files(args.coords_dat, args.radii_dat, args.output_dat)
    print(f"Done. Output written to {args.output_dat}")


if __name__ == "__main__":
    main()
