#!/usr/bin/env python3
"""
Convert CSV with header columns X, Y, Z, r into XYZ(+radius) format.

Output format:
  line 1: number of particles
  line 2: empty line
  line 3+: x y z r

Example:
  python csv_to_xyz_with_radius.py particles.csv particles.xyz
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def csv_to_xyz(csv_file: Path, xyz_file: Path) -> None:
    xs, ys, zs, rs = [], [], [], []

    with csv_file.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)

        required = {"X", "Y", "Z", "r"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(
                f"CSV must contain columns {required}, found {reader.fieldnames}"
            )

        for row in reader:
            xs.append(float(row["X"]))
            ys.append(float(row["Y"]))
            zs.append(float(row["Z"]))
            rs.append(float(row["r"]))

    n = len(xs)

    with xyz_file.open("w") as f:
        f.write(f"{n}\n")
        f.write("\n")  # XYZ comment line
        for x, y, z, r in zip(xs, ys, zs, rs):
            f.write(f"{x} {y} {z} {r}\n")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("input_csv", type=Path, help="Input CSV file with X,Y,Z,r columns")
    ap.add_argument("output_xyz", type=Path, help="Output XYZ file (x y z r)")
    args = ap.parse_args()

    csv_to_xyz(args.input_csv, args.output_xyz)


if __name__ == "__main__":
    main()
