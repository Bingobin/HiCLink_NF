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

    def test_feature_restart_rejects_conflicting_ids_and_invalid_assays(self):
        with self.assertRaisesRegex(ValueError, "conflicting library metadata"):
            self.validate("ID,assay\nS1,hic\nS1,hichip\n", "features")
        with self.assertRaisesRegex(ValueError, "assay must be"):
            self.validate("ID,assay\nS1,rnaseq\n", "features")

    def test_same_library_runs_receive_distinct_run_ids(self):
        rows = self.validate(
            "ID,R1,R2,assay,condition,replicate\n"
            "S1,runA/R1.fq.gz,runA/R2.fq.gz,hichip,SCR,1\n"
            "S1,runB/R1.fq.gz,runB/R2.fq.gz,hichip,SCR,1\n"
            "S2,R1.fq.gz,R2.fq.gz,hic,SH,2\n"
        )
        self.assertEqual([row["ID"] for row in rows], ["S1", "S1", "S2"])
        self.assertEqual([row["run_id"] for row in rows], ["S1.run001", "S1.run002", "S2"])
        self.assertEqual(rows[0]["R1"], str(self.root / "runA/R1.fq.gz"))

    def test_repeated_id_requires_consistent_library_metadata(self):
        for field in SAMPLESHEET.LIBRARY_FIELDS:
            base = {"ID": "S1", "R1": "a1.fq", "R2": "a2.fq", "assay": "hichip",
                    "condition": "SCR", "replicate": "1", "control": "CTRL", "peak_file": "a.bed"}
            changed = {**base, "R1": "b1.fq", "R2": "b2.fq", field: "hic" if field == "assay" else "different"}
            fields = list(base)
            contents = ",".join(fields) + "\n" + ",".join(base.values()) + "\n" + ",".join(changed.values()) + "\n"
            with self.subTest(field=field), self.assertRaisesRegex(ValueError, "conflicting library metadata"):
                self.validate(contents)

    def test_same_fastq_cannot_be_counted_twice_within_library(self):
        with self.assertRaisesRegex(ValueError, "reused FASTQ"):
            self.validate("ID,R1,R2,assay\nS1,a.fq,b.fq,hic\nS1,a.fq,b.fq,hic\n")
        with self.assertRaisesRegex(ValueError, "reused FASTQ"):
            self.validate("ID,R1,R2,assay\nS1,a.fq,a.fq,hic\n")

    def test_generated_qc_ids_cannot_collide_with_another_library(self):
        with self.assertRaisesRegex(ValueError, "run ID collision"):
            self.validate(
                "ID,R1,R2,assay\nS1,a1.fq,a2.fq,hic\n"
                "S1,b1.fq,b2.fq,hic\nS1.run001,c1.fq,c2.fq,hic\n"
            )

    def test_feature_restart_collapses_sequencing_runs_without_accessing_fastqs(self):
        rows = self.validate(
            "ID,R1,R2,assay\nS1,/missing/a1,/missing/a2,hichip\n"
            "S1,/missing/b1,/missing/b2,hichip\n", "features"
        )
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["ID"], "S1")


if __name__ == "__main__":
    unittest.main()
