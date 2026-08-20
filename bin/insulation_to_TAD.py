#!/usr/bin/env python3
"""Convert cooltools insulation boundaries into intervals between boundaries."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def insulation_to_tad(input_file: Path, output_file: Path, window_size: int) -> None:
    table = pd.read_csv(input_file, sep="\t")
    boundary_col = f"is_boundary_{window_size}"
    required = {"chrom", "start", boundary_col}
    missing = sorted(required - set(table.columns))
    if missing:
        raise ValueError(f"Missing insulation columns: {', '.join(missing)}")

    boundaries = table.loc[table[boundary_col].fillna(False)].copy()
    tads: list[list[object]] = []
    for chrom, group in boundaries.groupby("chrom", sort=False):
        starts = sorted(group["start"].astype(int).unique())
        for left, right in zip(starts, starts[1:]):
            if right > left:
                tads.append([chrom, left, right, f"{chrom}:{left}-{right}", 0, "."])

    output = pd.DataFrame(
        tads,
        columns=["chrom", "start", "end", "name", "score", "strand"],
    )
    output.to_csv(output_file, sep="\t", index=False, header=False)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Convert cooltools insulation boundaries to a BED6 TAD file"
    )
    parser.add_argument("input", type=Path)
    parser.add_argument("-o", "--output", type=Path, default=Path("tads.bed"))
    parser.add_argument("-w", "--window", type=int, required=True)
    args = parser.parse_args()
    insulation_to_tad(args.input, args.output, args.window)


if __name__ == "__main__":
    main()
