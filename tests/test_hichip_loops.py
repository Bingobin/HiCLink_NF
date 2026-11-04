"""Metadata, discovery, and format tests; no loop caller is executed."""
import csv
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest


SPEC = importlib.util.spec_from_file_location(
    "hichip_loops", Path(__file__).resolve().parents[1] / "bin/hichip_loops.py"
)
LOOPS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(LOOPS)


def settings():
    return {"callers": ["cloops2", "mustache"], "run_diff": True,
            "resolution": 10000, "min_distance": 20000, "max_distance": 2000000,
            "chroms": ["chr1"], "eps": "5000,10000", "minpts": "10",
            "diff_pcut": 0.01, "diff_ma_fdr": 0.1, "mustache_pt": 0.05,
            "mustache_diff_pt": 0.05, "sparsity": 0.88, "balance_max_iters": 1000,
            "outdir": "test-output"}


class LoopInputTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        (self.root / "03_pairs").mkdir()
        (self.root / "04_matrices").mkdir()
        self.sheet = self.root / "samples.csv"

    def files(self, sample, matrix="mcool"):
        (self.root / "03_pairs" / f"{sample}.valid.pairs.gz").touch()
        (self.root / "04_matrices" / f"{sample}.{matrix}").touch()

    def discover(self, contents, cfg=None):
        self.sheet.write_text(contents)
        out = self.root / "discovered"
        cfg = cfg or settings()
        cfg["outdir"] = str(self.root / "published")
        LOOPS.discover(self.sheet, self.root, cfg, out)
        with (out / "loop_samples.tsv").open() as handle:
            samples = list(csv.DictReader(handle, delimiter="\t"))
        with (out / "loop_comparisons.tsv").open() as handle:
            comparisons = list(csv.DictReader(handle, delimiter="\t"))
        return samples, comparisons

    def test_discovers_exact_ids_and_ignores_missing_server_fastqs(self):
        for sample in ("A549.RAD21.SCR", "A549.RAD21.SH2", "A549.HMGA2.OE_HA"):
            self.files(sample)
        samples, comparisons = self.discover(
            "ID,R1,R2,assay,condition,replicate,control\n"
            "A549.RAD21.SCR,/missing/server_R1,/missing/server_R2,hichip,SCR,1,\n"
            "A549.RAD21.SH2,,,hichip,SH2,1,A549.RAD21.SCR\n"
            "A549.HMGA2.OE_HA,,,hichip,OE_HA,1,\n")
        self.assertEqual(len(samples), 3)
        self.assertEqual(len(comparisons), 1)
        self.assertEqual(comparisons[0]["treatment"], "A549.RAD21.SH2")
        self.assertTrue(samples[0]["matrix"].endswith("A549.RAD21.SCR.mcool"))

    def test_reuses_one_control_for_multiple_treatments(self):
        for sample in ("CTRL", "SH1", "SH2"):
            self.files(sample)
        _, comparisons = self.discover("ID,assay,control\nCTRL,hichip,\nSH1,hichip,CTRL\nSH2,hichip,CTRL\n")
        self.assertEqual([row["control"] for row in comparisons], ["CTRL", "CTRL"])

    def test_rejects_unknown_self_and_cross_ip_controls(self):
        for sample in ("A549.HMGA2.OE", "A549.RAD21.SCR"):
            self.files(sample)
        for control in ("missing", "A549.HMGA2.OE", "A549.RAD21.SCR"):
            with self.subTest(control=control), self.assertRaisesRegex(ValueError, "control|Cross-cell"):
                self.discover("ID,assay,control\nA549.RAD21.SCR,hichip,\n"
                              f"A549.HMGA2.OE,hichip,{control}\n")

    def test_rejects_missing_pairs_without_selecting_similar_id(self):
        self.files("S10")
        with self.assertRaisesRegex(ValueError, "Missing pairs for S1 "):
            self.discover("ID,assay\nS1,hichip\n")

    def test_matrix_priority_and_compatible_finer_cool_fallback(self):
        self.files("S1", "5000.cool")
        selected = LOOPS.find_matrix(self.root, "S1", 10000)
        self.assertEqual(selected.name, "S1.5000.cool")
        (self.root / "04_matrices/S1.10000.cool").touch()
        self.assertEqual(LOOPS.find_matrix(self.root, "S1", 10000).name, "S1.10000.cool")
        (self.root / "04_matrices/S1.mcool").touch()
        self.assertEqual(LOOPS.find_matrix(self.root, "S1", 10000).name, "S1.mcool")

    def test_rejects_coarser_named_matrix(self):
        self.files("S1", "25000.cool")
        with self.assertRaisesRegex(ValueError, "Missing compatible"):
            LOOPS.find_matrix(self.root, "S1", 10000)

    def test_hic_rows_are_not_sent_to_the_hichip_subworkflow(self):
        self.files("HC")
        samples, comparisons = self.discover("ID,assay\nHI,hic\nHC,hichip\n")
        self.assertEqual([row["ID"] for row in samples], ["HC"])
        self.assertEqual(comparisons, [])

    def test_diff_disabled_keeps_individual_calls_only(self):
        self.files("S1")
        cfg = settings()
        cfg["run_diff"] = False
        _, comparisons = self.discover("ID,assay,control\nS1,hichip,not_present\n", cfg)
        self.assertEqual(comparisons, [])

    def test_rejects_duplicate_ids_and_silent_mustache_distance_adjustment(self):
        self.files("S1")
        with self.assertRaisesRegex(ValueError, "duplicated"):
            self.discover("ID,assay\nS1,hichip\nS1,hichip\n")
        cfg = settings()
        cfg["max_distance"] = 1000000
        with self.assertRaisesRegex(ValueError, "Mustache requires"):
            LOOPS.validate_settings(cfg)

    def test_distance_filter_preserves_native_fields_and_valid_empty_results(self):
        source = self.root / "native.tsv"
        header = "BIN1_CHR\tBIN1_START\tBIN1_END\tBIN2_CHROMOSOME\tBIN2_START\tBIN2_END\tFDR\tDETECTION_SCALE\n"
        source.write_text(header + "chr1\t0\t10000\tchr1\t20000\t30000\t0.001\t1.6\n"
                          + "chr1\t0\t10000\tchr1\t10000\t20000\t0.002\t1.6\n")
        target = self.root / "filtered.tsv"
        LOOPS.filter_mustache(source, target, 20000, 2000000)
        self.assertEqual(len(target.read_text().splitlines()), 2)
        self.assertIn("0.001\t1.6", target.read_text())
        empty = self.root / "empty.tsv"
        LOOPS.filter_mustache(source, empty, 100000, 2000000)
        self.assertEqual(empty.read_text(), header)

    def test_comparison_reports_must_have_matching_chromosome_lengths(self):
        first, second = self.root / "a.json", self.root / "b.json"
        first.write_text(json.dumps({"chromsizes": {"chr1": 100}}))
        second.write_text(json.dumps({"chromsizes": {"chr1": 101}}))
        with self.assertRaisesRegex(ValueError, "different chromosome sizes"):
            LOOPS.check_pair_reports(first, second)

    def test_signature_allows_identical_retry_and_rejects_changed_inputs(self):
        self.files("S1")
        self.sheet.write_text("ID,assay\nS1,hichip\n")
        cfg = settings()
        cfg["outdir"] = str(self.root / "published")
        first = self.root / "first"
        LOOPS.discover(self.sheet, self.root, cfg, first)
        info = Path(cfg["outdir"]) / "pipeline_info"
        info.mkdir(parents=True)
        (info / "loop_inputs_signature.json").write_bytes((first / "loop_inputs_signature.json").read_bytes())
        LOOPS.discover(self.sheet, self.root, cfg, self.root / "retry")
        (self.root / "04_matrices/S1.mcool").write_bytes(b"changed")
        with self.assertRaisesRegex(ValueError, "use a new hichip_loop_outdir"):
            LOOPS.discover(self.sheet, self.root, cfg, self.root / "changed")


if __name__ == "__main__":
    unittest.main()
