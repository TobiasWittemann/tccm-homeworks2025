#!/usr/bin/env python3

"""
Parse traj.xyz produced by the Fortran MD code and extract E(tot) vs step.
Outputs:
  - energy.csv  (step, time_ps, T_kJmol, V_kJmol, E_kJmol)
  - energy.png  
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import matplotlib.pyplot as plt


LINE_RE = re.compile(
    r"Step=\s*(\d+).*?"
    r"T\(kJ/mol\)=\s*([+-]?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)\s*.*?"
    r"V\(kJ/mol\)=\s*([+-]?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)\s*.*?"
    r"E\(kJ/mol\)=\s*([+-]?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)"
)

def parse_traj_xyz(path: Path):
    steps, Tvals, Vvals, Evals = [], [], [], []

    with path.open("r", encoding="utf-8", errors="replace") as f:
        while True:
            first = f.readline()
            if not first:
                break

            first = first.strip()
            if not first:
                continue

            try:
                natoms = int(first.split()[0])
            except ValueError:
                # Not a frame header; skip and keep scanning
                continue

            comment = f.readline()
            if not comment:
                break

            m = LINE_RE.search(comment)
            if m:
                step = int(m.group(1))
                T = float(m.group(2))
                V = float(m.group(3))
                E = float(m.group(4))
                steps.append(step)
                Tvals.append(T)
                Vvals.append(V)
                Evals.append(E)

            # Skip atom lines of this frame
            for _ in range(natoms):
                if not f.readline():
                    break

    return steps, Tvals, Vvals, Evals


def write_csv(outpath: Path, steps, Tvals, Vvals, Evals, dt: float):
    with outpath.open("w", newline="", encoding="utf-8") as g:
        w = csv.writer(g)
        w.writerow(["step", "time_ps", "T_kJmol", "V_kJmol", "E_kJmol"])
        for s, T, V, E in zip(steps, Tvals, Vvals, Evals):
            time_ps = s * dt
            w.writerow([s, f"{time_ps:.10g}", f"{T:.10g}", f"{V:.10g}", f"{E:.10g}"])


def plot_energy(outpng: Path, steps, Evals, dt: float):
    times = [s * dt for s in steps]

    plt.figure()
    plt.plot(times, Evals, marker="o", linestyle="-")
    plt.xlabel("Time (ps)")
    plt.ylabel("Total energy E (kJ/mol)")
    plt.title("MD energy drift: E(t)")
    plt.tight_layout()
    plt.savefig(outpng, dpi=200)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("traj", type=Path, help="Path to traj.xyz written by the MD code")
    ap.add_argument("--dt", type=float, default=0.02, help="Time step in ps (default: 0.02)")
    ap.add_argument("--csv", type=Path, default=Path("energy.csv"), help="Output CSV (default: energy.csv)")
    ap.add_argument("--png", type=Path, default=Path("energy.png"), help="Output plot PNG (default: energy.png)")
    ap.add_argument("--no-plot", action="store_true", help="Do not generate plot")
    args = ap.parse_args()

    steps, Tvals, Vvals, Evals = parse_traj_xyz(args.traj)
    if not steps:
        raise SystemExit(
            "No frames parsed. Check that traj.xyz contains lines like "
            "'Step= ... T(kJ/mol)= ... V(kJ/mol)= ... E(kJ/mol)= ...'."
        )

    write_csv(args.csv, steps, Tvals, Vvals, Evals, args.dt)

    if not args.no_plot:
        plot_energy(args.png, steps, Evals, args.dt)

    print(f"Parsed {len(steps)} frames.")
    print(f"Wrote: {args.csv}")
    if not args.no_plot:
        print(f"Wrote: {args.png}")


if __name__ == "__main__":
    main()
