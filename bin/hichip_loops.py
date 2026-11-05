#!/usr/bin/env python3
"""Input discovery, validation, and format handling for the HiChIP loop subworkflow."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import importlib.metadata
import json
import re
import shlex
import subprocess
import sys
from pathlib import Path


SAMPLE_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
SAMPLE_FIELDS = ["ID", "assay", "condition", "replicate", "control", "pairs", "matrix"]
COMPARISON_FIELDS = ["id", "treatment", "control"]


def require(test, message):
    if not test:
        raise ValueError(message)


def validate_settings(cfg):
    require(bool(cfg["callers"]) and set(cfg["callers"]) <= {"cloops2", "mustache", "fithichip"},
            "Loop callers must be selected from cloops2,mustache,fithichip")
    if 'fithichip' in cfg['callers']:
        fit = cfg['fithichip']
        require(fit['interaction_type'] in (1, 2, 3, 4), 'FitHiChIP supports one interaction type per run (1-4)')
        require(fit['fold_change'] >= 2 and fit['bcv'] > 0, 'FitHiChIP fold_change must be >=2 and BCV >0')
    for key in ("resolution", "min_distance", "max_distance", "balance_max_iters"):
        require(isinstance(cfg[key], int) and cfg[key] > 0, f"{key} must be a positive integer")
    require(cfg["min_distance"] < cfg["max_distance"], "min_distance must be smaller than max_distance")
    if "mustache" in cfg["callers"]:
        res, dist = cfg["resolution"], cfg["max_distance"]
        require(200 * res <= dist <= min(2000 * res, 2_000_000),
                "Mustache requires 200*resolution <= max_distance <= min(2000*resolution, 2Mb)")
    for key in ("eps", "minpts"):
        require(re.fullmatch(r"[1-9][0-9]*(,[1-9][0-9]*)*", str(cfg[key])),
                f"{key} must contain comma-separated positive integers")
    for key in ("diff_pcut", "diff_ma_fdr", "mustache_pt", "mustache_diff_pt", "sparsity"):
        require(0 < float(cfg[key]) < 1, f"{key} must be strictly between 0 and 1")
    chroms = cfg["chroms"]
    require(bool(chroms) and len(set(chroms)) == len(chroms)
            and all(re.fullmatch(r"chr([1-9]|1[0-9]|2[0-2]|X)", c) for c in chroms),
            "chroms must contain unique entries from chr1-22/chrX")


def find_pairs(root, sample):
    # Prefer published valid pairs over less explicitly named pairs files.
    for suffix in (".valid.pairs.gz", ".valid.pairs", ".pairs.gz", ".pairs"):
        path = root / "03_pairs" / (sample + suffix)
        if path.is_file():
            return path.resolve()
    raise ValueError(f"Missing pairs for {sample} under {root / '03_pairs'}")


def find_matrix(root, sample, resolution):
    multi = root / "04_matrices" / f"{sample}.mcool"
    if multi.is_file():
        return multi.resolve()
    exact = root / "04_matrices" / f"{sample}.{resolution}.cool"
    if exact.is_file():
        return exact.resolve()
    plain = root / "04_matrices" / f"{sample}.cool"
    if plain.is_file():
        return plain.resolve()
    candidates = []
    pattern = re.compile(re.escape(sample) + r"\.(\d+)\.cool$")
    for path in (root / "04_matrices").glob("*.cool"):
        match = pattern.fullmatch(path.name)
        if match:
            binsize = int(match[1])
            if binsize > 0 and binsize <= resolution and resolution % binsize == 0:
                candidates.append((binsize, path))
    if candidates:
        return max(candidates, key=lambda item: item[0])[1].resolve()
    raise ValueError(f"Missing compatible cool/mcool for {sample} under {root / '04_matrices'}")


def find_peaks(root, sample, supplied=''):
    if supplied:
        return Path(supplied).resolve(strict=True)
    path = root / '06_hichip/peaks' / f'{sample}.peaks.bed'
    if path.is_file():
        return path.resolve()
    raise ValueError(f'Missing peaks for {sample}: provide peak_file or 06_hichip/peaks/{sample}.peaks.bed')


def discover(sheet, root, cfg, output_dir):
    validate_settings(cfg)
    with sheet.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        require({"ID", "assay"} <= set(reader.fieldnames or []), "Samplesheet requires ID and assay columns")
        records = []
        seen = set()
        for raw in reader:
            row = {key: (value or "").strip() for key, value in raw.items() if key}
            sid = row["ID"]
            require(SAMPLE_RE.fullmatch(sid) and sid not in seen, f"Invalid or duplicated sample ID: {sid}")
            seen.add(sid)
            assay = re.sub(r"[-_\s]", "", row["assay"].lower())
            require(assay in ("hic", "hichip"), f"Unknown assay for {sid}: {row['assay']}")
            if assay != "hichip":
                continue
            records.append({
                "ID": sid, "assay": assay, "condition": row.get("condition") or sid,
                "replicate": row.get("replicate") or "1", "control": row.get("control", ""),
                "pairs": str(find_pairs(root, sid)),
                "matrix": str(find_matrix(root, sid, cfg["resolution"])),
            })
            if 'fithichip' in cfg['callers']:
                records[-1]['peak_file'] = str(find_peaks(root, sid, row.get('peak_file', '')))
    require(bool(records), "No HiChIP samples found in the samplesheet")
    by_id = {row["ID"]: row for row in records}
    comparisons = []
    for row in records:
        sid, control = row["ID"], row["control"]
        if not control or not cfg["run_diff"]:
            continue
        require(control in by_id and control != sid, f"Invalid HiChIP control for {sid}: {control}")
        # The current HMGA2 IDs encode cell line and IP target in their first two fields.
        # Generic IDs are allowed; explicit matching is the user's responsibility.
        tparts, cparts = sid.split("."), control.split(".")
        if len(tparts) >= 3 and len(cparts) >= 3:
            require(tparts[:2] == cparts[:2], f"Cross-cell-line or cross-IP comparison: {sid} vs {control}")
        comparisons.append({"id": f"{sid}_vs_{control}", "treatment": sid, "control": control})
    signature = {
        "settings": cfg, "results_dir": str(root.resolve()),
        "samplesheet_sha256": hashlib.sha256(sheet.read_bytes()).hexdigest(),
        "helper_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "input_stats": [
            {"path": row[field], "size": Path(row[field]).stat().st_size,
             "mtime_ns": Path(row[field]).stat().st_mtime_ns}
            for row in records for field in (("pairs", "matrix", "peak_file") if 'fithichip' in cfg['callers'] else ("pairs", "matrix"))
        ],
    }
    previous_signature = Path(cfg["outdir"]) / "pipeline_info/loop_inputs_signature.json"
    if previous_signature.exists():
        require(json.loads(previous_signature.read_text()) == signature,
                "Loop settings, samplesheet, helper, or input files changed; use a new hichip_loop_outdir")
    output_dir.mkdir(parents=True, exist_ok=True)
    sample_fields = SAMPLE_FIELDS + (['peak_file'] if 'fithichip' in cfg['callers'] else [])
    for name, fields, rows in (("loop_samples.tsv", sample_fields, records),
                               ("loop_comparisons.tsv", COMPARISON_FIELDS, comparisons)):
        with (output_dir / name).open("x", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
    (output_dir / "loop_settings.json").write_text(json.dumps(cfg, indent=2) + "\n")
    (output_dir / "loop_inputs_signature.json").write_text(json.dumps(signature, indent=2) + "\n")
    print(f"Discovered {len(records)} HiChIP samples and {len(comparisons)} directed comparisons")
    print('FASTQ, BAM and reference FASTA were not accessed; peak paths are required only for FitHiChIP')


def matrix_source(source, resolution):
    import cooler
    import h5py

    with h5py.File(source, "r") as handle:
        if "resolutions" in handle:
            available = [int(key) for key in handle["resolutions"] if key.isdigit()]
            candidates = [r for r in available if r <= resolution and resolution % r == 0]
            require(candidates, f"No resolution can generate {resolution} bp in {source}")
            selected = max(candidates)
            uri = f"{source}::/resolutions/{selected}"
        else:
            uri = str(source)
    clr = cooler.Cooler(uri)
    require(clr.binsize and clr.binsize <= resolution and resolution % clr.binsize == 0,
            f"Cannot derive {resolution} bp from {source}; no upsampling or fractional coarsening is allowed")
    return uri, clr


def check_sample(meta, pairs, matrix, cfg, output):
    validate_settings(cfg)
    require(SAMPLE_RE.fullmatch(meta["id"]), "Invalid sample ID")
    require(meta.get("assay") == "hichip", "This subworkflow accepts HiChIP samples only")
    report = {"sample": meta, "pairs": str(pairs), "matrix": str(matrix), "settings": cfg}
    if "cloops2" in cfg["callers"]:
        opener = gzip.open if pairs.suffix == ".gz" else open
        columns = None
        first = None
        with opener(pairs, "rt") as handle:
            for line in handle:
                if line.startswith("#columns:"):
                    columns = line.split(":", 1)[1].split()
                if not line.startswith("#"):
                    first = line.rstrip().split("\t")
                    break
        expected = ["readID", "chrom1", "pos1", "chrom2", "pos2", "strand1", "strand2", "pair_type"]
        require(columns and columns[:8] == expected, f"Unsupported pairs column order: {pairs}")
        require(first and len(first) >= 8 and first[7] in ("UU", "UR", "RU", "MR", "MU"),
                f"Expected nonempty upstream-filtered valid pairs: {pairs}")
        report["pairs_check"] = "header and first record only; not a complete file scan"
    if set(cfg['callers']) & {'mustache', 'fithichip'}:
        uri, clr = matrix_source(matrix, cfg["resolution"])
        require(set(cfg["chroms"]) <= set(clr.chromnames), f"Configured chromosomes missing: {matrix}")
        report["source_uri"] = uri
        report["source_resolution"] = int(clr.binsize)
        report["chromsizes"] = {str(c): int(v) for c, v in clr.chromsizes.items()}
    packages = set()
    if 'cloops2' in cfg['callers']:
        packages.add('cLoops2')
    if set(cfg['callers']) & {'mustache', 'fithichip'}:
        packages.add('cooler')
    if 'mustache' in cfg['callers']:
        packages.add('mustache-hic')
    report["software"] = {pkg: importlib.metadata.version(pkg) for pkg in sorted(packages)}
    output.write_text(json.dumps(report, indent=2) + "\n")


def check_cool(path, resolution, chroms):
    import cooler
    import h5py
    import numpy as np

    clr = cooler.Cooler(str(path))
    require(clr.binsize == resolution, f"Cool resolution mismatch: {path}")
    with h5py.File(path, "r") as handle:
        require("weight" in handle["bins"], f"Missing weights: {path}")
        weights = handle["bins/weight"]
        require(bool(np.all(weights.attrs.get("converged", False))), f"Weights did not converge: {path}")
        require(not bool(weights.attrs.get("divisive_weights", False)), "Multiplicative ICE weights are required")
    for chrom in chroms:
        values = clr.bins().fetch(chrom)["weight"].to_numpy()
        require(np.count_nonzero(np.isfinite(values) & (values > 0)) >= 10,
                f"Too few usable bins on {chrom}: {path}")


def run_tool(args):
    print(shlex.join([str(item) for item in args]), flush=True)
    subprocess.run([str(item) for item in args], check=True)


def prepare_cool(source, target, cfg, cpus):
    import h5py
    import numpy as np

    require(not target.exists() and source.resolve() != target.resolve(), "Working cool must be a new file, not the source")
    uri, clr = matrix_source(source, cfg["resolution"])
    if clr.binsize == cfg["resolution"]:
        run_tool(["cooler", "cp", uri, target])
    else:
        run_tool(["cooler", "coarsen", uri, "--factor", cfg["resolution"] // clr.binsize,
                  "--nproc", cpus, "--chunksize", 1_000_000, "--out", target])
    with h5py.File(target, "r") as handle:
        converged = ("weight" in handle["bins"] and bool(np.all(handle["bins/weight"].attrs.get("converged", False))))
    if not converged:
        # Only the newly created working matrix is modified; preserve source weights.
        run_tool(["cooler", "balance", target, "--force", "--nproc", cpus, "--chunksize", 1_000_000,
                  "--ignore-diags", 2, "--mad-max", 5, "--min-nnz", 10, "--min-count", 0,
                  "--tol", "1e-5", "--max-iters", cfg["balance_max_iters"], "--convergence-policy", "error"])
    check_cool(target, cfg["resolution"], cfg["chroms"])


def check_pre(path):
    meta = json.loads(path.read_text())
    require(meta.get("Unique PETs", 0) > 0 and meta.get("data", {}).get("cis"), f"No usable cis PETs: {path}")
    for record in meta["data"]["cis"].values():
        require(Path(record["ixy"]).is_file(), f"Missing ixy file: {record['ixy']}")


def require_loops(paths):
    for path in paths:
        with path.open() as handle:
            require(next(handle, "").split("\t")[0] == "loopId", f"Invalid cLoops2 header: {path}")
            require(any(line.strip() for line in handle), f"No candidate loops: {path}; inspect sequencing depth and parameters")


def filter_mustache(source, target, minimum, maximum):
    total = kept = 0
    with source.open() as src, target.open("x", newline="") as dst:
        reader = csv.reader(src, delimiter="\t")
        header = next(reader, [])
        require(len(header) == 8 and header[:3] == ["BIN1_CHR", "BIN1_START", "BIN1_END"],
                f"Missing or invalid Mustache output: {source}")
        writer = csv.writer(dst, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        for row in reader:
            if not row:
                continue
            require(len(row) == 8, f"Invalid Mustache output columns: {source}")
            total += 1
            if row[0] == row[3] and minimum <= abs(int(row[4]) - int(row[1])) <= maximum:
                writer.writerow(row)
                kept += 1
    print(f"{target.name}: native={total}, distance_filtered={kept}")
    if not kept:
        print("WARNING: 0 loops; this is a valid empty result, not evidence of differential interaction", file=sys.stderr)


def check_pair_reports(first, second):
    a, b = json.loads(first.read_text()), json.loads(second.read_text())
    require(a.get("chromsizes") == b.get("chromsizes"), "Comparison matrices have different chromosome sizes")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("discover")
    p.add_argument("--samplesheet", type=Path, required=True)
    p.add_argument("--results-dir", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    p = sub.add_parser("check-sample")
    p.add_argument("--meta", required=True)
    p.add_argument("--pairs", type=Path, required=True)
    p.add_argument("--matrix", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p = sub.add_parser("prepare-cool")
    p.add_argument("--source", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--cpus", type=int, required=True)
    p = sub.add_parser("check-pre")
    p.add_argument("path", type=Path)
    p = sub.add_parser("require-loops")
    p.add_argument("paths", nargs="+", type=Path)
    p = sub.add_parser("filter-mustache")
    p.add_argument("source", type=Path)
    p.add_argument("target", type=Path)
    p.add_argument("minimum", type=int)
    p.add_argument("maximum", type=int)
    p = sub.add_parser("check-pair-reports")
    p.add_argument("first", type=Path)
    p.add_argument("second", type=Path)
    for name in ("discover", "check-sample", "prepare-cool"):
        sub.choices[name].add_argument("--settings", required=True, help="JSON-encoded settings")
    args = parser.parse_args()
    if args.command == "discover":
        discover(args.samplesheet, args.results_dir.resolve(), json.loads(args.settings), args.output_dir)
    elif args.command == "check-sample":
        check_sample(json.loads(args.meta), args.pairs, args.matrix, json.loads(args.settings), args.output)
    elif args.command == "prepare-cool":
        prepare_cool(args.source, args.output, json.loads(args.settings), args.cpus)
    elif args.command == "check-pre":
        check_pre(args.path)
    elif args.command == "require-loops":
        require_loops(args.paths)
    elif args.command == "filter-mustache":
        filter_mustache(args.source, args.target, args.minimum, args.maximum)
    elif args.command == "check-pair-reports":
        check_pair_reports(args.first, args.second)


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, subprocess.CalledProcessError) as exc:
        print(f"HiChIP loop workflow error: {exc}", file=sys.stderr)
        raise SystemExit(2)
