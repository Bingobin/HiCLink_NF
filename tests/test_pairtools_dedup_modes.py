"""Exercise both dedup modes using stub tasks and tiny fixtures only."""

import csv
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


@unittest.skipUnless(NEXTFLOW, "Nextflow is required for stub workflow checks")
class PairtoolsDedupModeTests(unittest.TestCase):
    def run_stub(self, skip, repeated=False, run_loops=True):
        with tempfile.TemporaryDirectory(prefix="pairtools-dedup-") as directory:
            root = Path(directory)
            sheet = root / "samples.csv"
            with (ROOT / "tests/data/samplesheet.csv").open() as handle:
                reader = csv.DictReader(handle)
                fields = reader.fieldnames
                rows = list(reader)
            for row in rows:
                for key in ("R1", "R2", "peak_file"):
                    if row.get(key):
                        row[key] = str((ROOT / "tests/data" / row[key]).resolve())
            if repeated:
                second = dict(rows[-1])
                extra = root / "run2"
                extra.mkdir()
                for key in ("R1", "R2"):
                    target = extra / Path(second[key]).name
                    shutil.copyfile(second[key], target)
                    second[key] = str(target)
                rows.append(second)
            with sheet.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields)
                writer.writeheader()
                writer.writerows(rows)
            output = root / "results"
            result = subprocess.run(
                [NEXTFLOW, "run", str(ROOT / "main.nf"), "-profile", "test", "-stub-run",
                 "-ansi-log", "false", "-w", str(root / "work"), "--input", str(sheet),
                 "--outdir", str(output), "--skip_pairtools_dedup", str(skip).lower(),
                 "--run_loops", str(run_loops).lower()],
                cwd=root, capture_output=True, text=True,
                env={**os.environ, "NXF_OFFLINE": "true"}, timeout=60,
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            qc = output / "01_qc/pairtools"
            for row in rows:
                sample = row["ID"]
                self.assertTrue((output / f"03_pairs/{sample}.valid.pairs.gz").exists())
                self.assertTrue((qc / f"{sample}.parse.stats").exists())
                self.assertTrue((qc / f"{sample}.pairtools.stats").exists())
                self.assertEqual((qc / f"{sample}.dedup.stats").exists(), not skip)
                self.assertEqual((qc / f"{sample}.dedup_policy.tsv").read_text(),
                                 f"pairtools_dedup\t{'skipped' if skip else 'enabled'}\n")
            self.assertTrue((output / "pipeline_info/multiqc/multiqc_report.html").exists())
            with (output / "pipeline_info/samplesheet.validated.csv").open() as handle:
                validated = list(csv.DictReader(handle))
            for row in validated:
                self.assertTrue((output / f"01_qc/fastp/{row['run_id']}.fastp.json").exists())
                self.assertTrue((output / f"01_qc/fastqc_raw/{row['run_id']}.R1_fastqc.zip").exists())
            libraries = {row["ID"] for row in validated}
            for sample in libraries:
                manifest = output / f"pipeline_info/sequencing_runs/{sample}.sequencing_runs.tsv"
                self.assertEqual(manifest.read_text().splitlines()[1:],
                                 sorted(row["run_id"] for row in validated if row["ID"] == sample))
            with (output / "pipeline_info/execution_trace.txt").open() as handle:
                tasks = list(csv.DictReader(handle, delimiter="\t"))
            for process in ("FASTQC_RAW", "FASTP"):
                self.assertEqual(sum(task["name"].startswith(process + " (") for task in tasks), len(rows))
            for process in ("MERGE_FASTQ", "BWA_ALIGN", "PAIRTOOLS", "COOLER_MATRIX"):
                self.assertEqual(sum(task["name"].startswith(process + " (") for task in tasks), len(libraries))
            self.assertEqual(any(task["name"].startswith("COOLTOOLS_DOTS (") for task in tasks), run_loops)
            self.assertFalse(any(task["name"].startswith("FITHICHIP (") for task in tasks))
            self.assertTrue(any(task["name"].startswith("HICHIP_PEAKS (") for task in tasks))

    def test_default_dedup_mode(self):
        self.run_stub(False)

    def test_skip_dedup_mode(self):
        self.run_stub(True)

    def test_repeated_library_runs_with_dedup(self):
        self.run_stub(False, repeated=True)

    def test_repeated_library_runs_without_dedup(self):
        self.run_stub(True, repeated=True)

    def test_repeated_runs_can_skip_all_loop_calling(self):
        self.run_stub(True, repeated=True, run_loops=False)

    def test_standalone_loop_entry_can_be_disabled_without_inputs(self):
        with tempfile.TemporaryDirectory(prefix="loops-disabled-") as directory:
            result = subprocess.run(
                [NEXTFLOW, "-c", str(ROOT / "conf/hichip_loops.config"), "run", str(ROOT / "hichip_loops.nf"),
                 "-ansi-log", "false", "--run_loops", "false", "--input", "/missing/samplesheet.csv"],
                cwd=directory, capture_output=True, text=True,
                env={**os.environ, "NXF_OFFLINE": "true"}, timeout=60,
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertIn("Loop calling disabled", result.stdout)
            self.assertNotIn("Submitted process", result.stdout)


if __name__ == "__main__":
    unittest.main()
