"""Real gzip merge and single-run symlink tests; no biological data are used."""

import csv
import gzip
import hashlib
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
NEXTFLOW = shutil.which("nextflow")
if not NEXTFLOW and (Path.home() / "bin/nextflow").is_file():
    NEXTFLOW = str(Path.home() / "bin/nextflow")


@unittest.skipUnless(NEXTFLOW, "Nextflow is required for real merge checks")
class MergeFastqTests(unittest.TestCase):
    def test_gzip_members_and_mate_order_are_preserved_without_changing_sources(self):
        with tempfile.TemporaryDirectory(prefix="merge-fastq-") as directory:
            root = Path(directory)
            rows, expected, hashes = [], {"MULTI": ["A", "B"], "SINGLE": ["C"]}, {}
            for library, runs in expected.items():
                for run in runs:
                    paths = []
                    for mate in (1, 2):
                        path = root / f"{run}.R{mate}.fastq.gz"
                        with gzip.open(path, "wt") as handle:
                            handle.write(f"@{run}/{mate}\nACGT\n+\nIIII\n")
                        hashes[path] = hashlib.sha256(path.read_bytes()).hexdigest()
                        paths.append(str(path))
                    rows.append([library, run, *paths])
            sheet = root / "inputs.csv"
            with sheet.open("w", newline="") as handle:
                writer = csv.writer(handle)
                writer.writerow(["ID", "run_id", "R1", "R2"])
                writer.writerows(rows)
            result = subprocess.run(
                [NEXTFLOW, "-c", str(ROOT / "nextflow.config"), "run", str(ROOT / "tests/merge_fastq.nf"),
                 "-profile", "test", "-ansi-log", "false", "-w", str(root / "work"),
                 "--input", str(sheet), "--outdir", str(root / "results")],
                cwd=root, capture_output=True, text=True,
                env={**os.environ, "NXF_OFFLINE": "true"}, timeout=60,
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            outputs = [line.split("\t") for line in result.stdout.splitlines() if line.startswith("MERGED\t")]
            self.assertEqual({row[1] for row in outputs}, set(expected))
            for _, library, r1, r2 in outputs:
                for mate, filename in enumerate((r1, r2), start=1):
                    with gzip.open(filename, "rt") as handle:
                        names = handle.read().splitlines()[::4]
                    self.assertEqual(names, [f"@{run}/{mate}" for run in expected[library]])
                    self.assertEqual(Path(filename).is_symlink(), library == "SINGLE")
                manifest = root / f"results/pipeline_info/sequencing_runs/{library}.sequencing_runs.tsv"
                with manifest.open() as handle:
                    runs = list(csv.DictReader(handle, delimiter="\t"))
                self.assertEqual([row["run_id"] for row in runs], expected[library])
            for path, digest in hashes.items():
                self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), digest)


if __name__ == "__main__":
    unittest.main()
