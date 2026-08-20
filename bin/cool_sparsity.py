#!/usr/bin/env python3
"""Report contact-matrix density and sparsity without jq/bc dependencies."""

from __future__ import annotations

import argparse
from pathlib import Path

import cooler


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("cool_uri")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    matrix = cooler.Cooler(args.cool_uri)
    nbins = int(matrix.info["nbins"])
    nnz = int(matrix.info["nnz"])
    possible = nbins * (nbins + 1) // 2
    density = nnz / possible if possible else 0.0

    with args.output.open("w", encoding="utf-8") as handle:
        handle.write("cool_uri\tnbins\tnnz\tpossible_upper_triangle\tdensity\tsparsity\n")
        handle.write(
            f"{args.cool_uri}\t{nbins}\t{nnz}\t{possible}\t{density:.8g}\t{1.0 - density:.8g}\n"
        )


if __name__ == "__main__":
    main()

