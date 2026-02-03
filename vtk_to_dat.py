#!/usr/bin/env python3
"""
Convert legacy ASCII VTK (DATASET ... + POINTS + POINT_DATA/SCALARS radius)
to Chrono-style particle .dat lines:  n x y z r

Example:
  python vtk_to_dat.py Lbox_va_A225200.vtk out.dat
  python vtk_to_dat.py Lbox_va_A225200.vtk out.dat --radius-name radius
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import List, Tuple


def _parse_points(lines: List[str]) -> Tuple[List[Tuple[float, float, float]], int]:
    npoints = None
    points: List[Tuple[float, float, float]] = []

    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.upper().startswith("POINTS"):
            parts = line.split()
            if len(parts) < 3:
                raise ValueError(f"Malformed POINTS line: {line}")
            npoints = int(parts[1])

            needed = 3 * npoints
            tokens: List[str] = []
            i += 1
            while i < len(lines) and len(tokens) < needed:
                tokens.extend(lines[i].split())
                i += 1

            if len(tokens) < needed:
                raise ValueError(f"POINTS section ended early: got {len(tokens)} floats, need {needed}")

            vals = list(map(float, tokens[:needed]))
            points = [(vals[j], vals[j + 1], vals[j + 2]) for j in range(0, needed, 3)]
            return points, npoints
        i += 1

    raise ValueError("POINTS section not found.")


def _parse_scalar_under_point_data(
    lines: List[str], npoints: int, scalar_name: str
) -> List[float]:
    # Find POINT_DATA, then scan SCALARS blocks
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.upper().startswith("POINT_DATA"):
            parts = line.split()
            if len(parts) < 2:
                raise ValueError(f"Malformed POINT_DATA line: {line}")
            pd_n = int(parts[1])
            # Usually pd_n == npoints, but we won't hard-fail if it isn't.
            count = pd_n

            i += 1
            while i < len(lines):
                l = lines[i].strip()
                if l.upper().startswith("SCALARS"):
                    sp = l.split()
                    if len(sp) < 2:
                        raise ValueError(f"Malformed SCALARS line: {l}")
                    name = sp[1]

                    # Next line often "LOOKUP_TABLE default"
                    i += 1
                    if i < len(lines) and lines[i].strip().upper().startswith("LOOKUP_TABLE"):
                        i += 1

                    # Read scalar values (count of them)
                    tokens: List[str] = []
                    while i < len(lines) and len(tokens) < count:
                        nxt = lines[i].strip()
                        if not nxt:
                            i += 1
                            continue

                        # If a new section keyword appears (rare mid-way), stop
                        # (We keep it conservative: keyword-like + not starting with number)
                        if re.match(r"^[A-Z_]", nxt) and not re.match(r"^[-+0-9.]", nxt):
                            break

                        tokens.extend(lines[i].split())
                        i += 1

                    if len(tokens) < count:
                        raise ValueError(
                            f"SCALARS {name} ended early: got {len(tokens)} values, need {count}"
                        )

                    if name == scalar_name:
                        return list(map(float, tokens[:count]))

                    # else: continue scanning for other SCALARS blocks
                    continue

                i += 1

            break
        i += 1

    raise ValueError(f"SCALARS '{scalar_name}' not found under POINT_DATA.")


def vtk_to_dat(in_vtk: Path, out_dat: Path, radius_name: str = "radius") -> None:
    lines = in_vtk.read_text(errors="ignore").splitlines(True)

    points, npoints = _parse_points(lines)
    radius = _parse_scalar_under_point_data(lines, npoints, radius_name)

    if len(points) != len(radius):
        raise ValueError(f"Points ({len(points)}) and radius ({len(radius)}) length mismatch.")

    with out_dat.open("w") as f:
        for idx, ((x, y, z), r) in enumerate(zip(points, radius)):
            # Match your format: n x y z r
            f.write(f"{idx} {x} {y} {z} {r}\n")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("input_vtk", type=Path, help="Input legacy ASCII .vtk file")
    ap.add_argument("output_dat", type=Path, help="Output .dat file (n x y z r per line)")
    ap.add_argument("--radius-name", default="radius", help="SCALARS array name for radius (default: radius)")
    args = ap.parse_args()

    vtk_to_dat(args.input_vtk, args.output_dat, radius_name=args.radius_name)


if __name__ == "__main__":
    main()
