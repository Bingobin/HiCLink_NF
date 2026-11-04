#!/usr/bin/env python3
"""Validate and normalize a HiCLink_NF samplesheet."""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path


FIELDS = [
    "ID",
    "R1",
    "R2",
    "assay",
    "condition",
    "replicate",
    "control",
    "peak_file",
]
REQUIRED = {"ID", "assay"}
SAMPLE_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
URI_RE = re.compile(r"^[A-Za-z][A-Za-z0-9+.-]*://")


def normalize_assay(value: str) -> str:
    key = re.sub(r"[-_\s]", "", value.strip().lower())
    aliases = {
        "hic": "hic",
        "hichip": "hichip",
    }
    if key not in aliases:
        raise ValueError("assay must be one of: hic, hichip")
    return aliases[key]


def normalize_path(value: str, base_dir: Path) -> str:
    value = value.strip()
    if not value or URI_RE.match(value):
        return value
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = base_dir / path
    return str(path.resolve(strict=False))


def validate(
    input_path: Path, output_path: Path, base_dir: Path, start_from: str = "fastq"
) -> None:
    errors: list[str] = []
    normalized: list[dict[str, str]] = []
    seen_samples: set[str] = set()

    with input_path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        headers = set(reader.fieldnames or [])
        required = REQUIRED | ({"R1", "R2"} if start_from == "fastq" else set())
        missing = sorted(required - headers)
        if missing:
            raise ValueError(f"missing required columns: {', '.join(missing)}")

        for line_number, raw in enumerate(reader, start=2):
            row = {field: (raw.get(field) or "").strip() for field in FIELDS}
            sample = row["ID"]
            if not SAMPLE_RE.fullmatch(sample):
                errors.append(
                    f"line {line_number}: invalid sample '{sample}'; use letters, digits, '.', '_' or '-'"
                )
            elif sample in seen_samples:
                errors.append(f"line {line_number}: duplicate sample '{sample}'")
            seen_samples.add(sample)

            try:
                row["assay"] = normalize_assay(row["assay"])
            except ValueError as exc:
                errors.append(f"line {line_number}: {exc}; received '{row['assay']}'")

            for field in ("R1", "R2"):
                if start_from == "fastq" and not row[field]:
                    errors.append(f"line {line_number}: {field} is empty")
                row[field] = normalize_path(row[field], base_dir)

            row["peak_file"] = normalize_path(row["peak_file"], base_dir)
            row["condition"] = row["condition"] or sample
            row["replicate"] = row["replicate"] or "1"
            normalized.append(row)

    if not normalized:
        errors.append("samplesheet has no data rows")
    if errors:
        raise ValueError("\n".join(errors))

    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(normalized)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--base-dir", required=True, type=Path)
    parser.add_argument("--start-from", choices=("fastq", "features"), default="fastq")
    args = parser.parse_args()
    try:
        validate(args.input, args.output, args.base_dir, args.start_from)
    except (OSError, ValueError) as exc:
        print(f"Samplesheet validation failed:\n{exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
