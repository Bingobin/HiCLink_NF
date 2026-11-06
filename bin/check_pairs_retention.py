#!/usr/bin/env python3
"""Check that selection preserves all UU/UR/RU records from parse or dedup stats."""

import argparse
from pathlib import Path


PAIR_TYPES = ("UU", "UR", "RU")


def read_counts(path):
    counts = {}
    for line in Path(path).read_text().splitlines():
        key, value = line.split("\t", 1)
        if key == "total" or key.startswith("pair_types/"):
            counts[key] = int(value)
    if "total" not in counts:
        raise ValueError(f"Missing total count in {path}")
    return counts


def validate_retention(before_path, selected_path):
    before = read_counts(before_path)
    selected = read_counts(selected_path)
    rows = [
        (pair_type, before.get(f"pair_types/{pair_type}", 0),
         selected.get(f"pair_types/{pair_type}", 0))
        for pair_type in PAIR_TYPES
    ]
    rows.append(("TOTAL", sum(row[1] for row in rows), selected["total"]))
    mismatches = [f"{kind}: expected {before}, retained {after}"
                  for kind, before, after in rows if before != after]
    unexpected = [key.removeprefix("pair_types/") for key, count in selected.items()
                  if key.startswith("pair_types/") and count
                  and key.removeprefix("pair_types/") not in PAIR_TYPES]
    if mismatches or unexpected:
        raise ValueError(
            "UU/UR/RU retention check failed: "
            + "; ".join(mismatches + [f"unexpected pair type: {kind}" for kind in unexpected])
        )
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--before", "--dedup", dest="before", required=True, type=Path,
                        help="Stats immediately before selection: parse if dedup is skipped, otherwise dedup")
    parser.add_argument("--selected", required=True, type=Path)
    args = parser.parse_args()
    try:
        rows = validate_retention(args.before, args.selected)
    except ValueError as exc:
        parser.error(str(exc))
    print("pair_type\texpected\tretained")
    for row in rows:
        print(*row, sep="\t")


if __name__ == "__main__":
    main()
