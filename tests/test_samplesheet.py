"""Samplesheet contracts for raw-read runs and feature restarts."""

import csv
import importlib.util
from pathlib import Path
import tempfile
import unittest


SPEC = importlib.util.spec_from_file_location(
    "validate_samplesheet", Path(__file__).resolve().parents[1] / "bin/validate_samplesheet.py"
)
SAMPLESHEET = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(SAMPLESHEET)


class SamplesheetTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name).resolve()
        self.source = self.root / "samples.csv"
        self.output = self.root / "validated.csv"

    def validate(self, contents, mode="fastq"):
        self.source.write_text(contents, encoding="utf-8")
        SAMPLESHEET.validate(self.source, self.output, self.root, mode)
        with self.output.open(newline="", encoding="utf-8") as handle:
            return list(csv.DictReader(handle))

    def test_feature_restart_accepts_metadata_without_reads(self):
        rows = self.validate("ID,assay\nS1,HiChIP\n", "features")
        self.assertEqual(rows[0]["assay"], "hichip")
        self.assertEqual(rows[0]["condition"], "S1")
        self.assertEqual(rows[0]["replicate"], "1")
        self.assertEqual(rows[0]["R1"], "")
        self.assertEqual(rows[0]["R2"], "")

    def test_raw_read_mode_still_requires_read_columns_and_values(self):
        with self.assertRaisesRegex(ValueError, "missing required columns: R1, R2"):
            self.validate("ID,assay\nS1,hic\n")
        with self.assertRaisesRegex(ValueError, "R1 is empty"):
            self.validate("ID,R1,R2,assay\nS1,,reads_R2.fq.gz,hic\n")

    def test_feature_restart_preserves_metadata_and_resolves_external_peak(self):
        rows = self.validate(
            "ID,assay,condition,replicate,control,peak_file\n"
            "S1,hichip,knockdown,2,S0,chip.narrowPeak.gz\n",
            "features",
        )
        self.assertEqual(rows[0]["condition"], "knockdown")
        self.assertEqual(rows[0]["replicate"], "2")
        self.assertEqual(rows[0]["control"], "S0")
        self.assertEqual(rows[0]["peak_file"], str(self.root / "chip.narrowPeak.gz"))

    def test_feature_restart_rejects_duplicate_ids_and_invalid_assays(self):
        with self.assertRaisesRegex(ValueError, "duplicate sample"):
            self.validate("ID,assay\nS1,hic\nS1,hichip\n", "features")
        with self.assertRaisesRegex(ValueError, "assay must be"):
            self.validate("ID,assay\nS1,rnaseq\n", "features")


if __name__ == "__main__":
    unittest.main()
