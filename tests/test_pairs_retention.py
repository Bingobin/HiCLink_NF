"""Regression tests for valid-pair selection and retention auditing."""

import importlib.util
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "bin/check_pairs_retention.py"
SPEC = importlib.util.spec_from_file_location("check_pairs_retention", SCRIPT)
RETENTION = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(RETENTION)


class PairsRetentionTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.dedup = self.root / "sample.dedup.stats"
        self.selected = self.root / "sample.pairtools.stats"

    def write_stats(self, path, counts):
        path.write_text(
            f"total\t{sum(counts.values())}\n"
            + "".join(f"pair_types/{kind}\t{count}\n" for kind, count in counts.items())
            + "summary/frac_cis\t0.5\n",
            encoding="utf-8",
        )

    def test_all_unique_and_rescued_contacts_are_preserved(self):
        self.write_stats(self.dedup, {"UU": 10, "UR": 3, "RU": 2, "DD": 50, "MM": 4})
        self.write_stats(self.selected, {"UU": 10, "UR": 3, "RU": 2})
        self.assertEqual(RETENTION.validate_retention(self.dedup, self.selected),
                         [("UU", 10, 10), ("UR", 3, 3), ("RU", 2, 2), ("TOTAL", 15, 15)])

    def test_uu_only_selection_fails_when_rescued_contacts_exist(self):
        self.write_stats(self.dedup, {"UU": 10, "UR": 3, "RU": 2, "DD": 50})
        self.write_stats(self.selected, {"UU": 10})
        with self.assertRaisesRegex(ValueError, "UR: expected 3, retained 0"):
            RETENTION.validate_retention(self.dedup, self.selected)

    def test_genuinely_uu_only_sample_passes(self):
        self.write_stats(self.dedup, {"UU": 10, "DD": 50})
        self.write_stats(self.selected, {"UU": 10})
        self.assertEqual(RETENTION.validate_retention(self.dedup, self.selected)[-1],
                         ("TOTAL", 10, 10))

    def test_each_type_is_checked_even_when_total_matches(self):
        self.write_stats(self.dedup, {"UU": 10, "UR": 3, "RU": 2})
        self.write_stats(self.selected, {"UU": 10, "UR": 2, "RU": 3})
        with self.assertRaisesRegex(ValueError, "UR: expected 3, retained 2"):
            RETENTION.validate_retention(self.dedup, self.selected)

    def test_unwanted_contacts_fail(self):
        self.write_stats(self.dedup, {"UU": 10, "UR": 3, "RU": 2, "DD": 50})
        self.write_stats(self.selected, {"UU": 10, "UR": 3, "RU": 2, "DD": 1})
        with self.assertRaisesRegex(ValueError, "unexpected pair type: DD"):
            RETENTION.validate_retention(self.dedup, self.selected)

    def test_missing_total_is_not_treated_as_empty_sample(self):
        self.dedup.touch()
        self.write_stats(self.selected, {"UU": 10})
        with self.assertRaisesRegex(ValueError, "Missing total count"):
            RETENTION.validate_retention(self.dedup, self.selected)

    def test_nextflow_selector_keeps_exactly_uu_ur_ru(self):
        source = (ROOT / "modules/local/preprocessing.nf").read_text()
        expression = re.search(r"pairtools select[^\n]*\n\s*'([^']+)'", source).group(1)
        types = ("UU", "UR", "RU", "DD", "MU", "MR", "MM", "NU", "NR", "NN", "WW", "XX")
        kept = [kind for kind in types
                if eval(expression, {"__builtins__": {}}, {"pair_type": kind})]
        self.assertEqual(kept, ["UU", "UR", "RU"])
        self.assertIn("check_pairs_retention.py", source)
        self.assertIn('path("${meta.id}.retention.tsv"), emit: retention', source)

    def test_cli_reports_counts_and_exits_nonzero_on_loss(self):
        self.write_stats(self.dedup, {"UU": 10, "UR": 3, "RU": 2})
        self.write_stats(self.selected, {"UU": 10, "UR": 3, "RU": 2})
        command = [sys.executable, str(SCRIPT), "--dedup", str(self.dedup),
                   "--selected", str(self.selected)]
        result = subprocess.run(command, capture_output=True, text=True, check=True)
        self.assertIn("pair_type\texpected\tretained\n", result.stdout)
        self.assertIn("TOTAL\t15\t15\n", result.stdout)
        self.write_stats(self.selected, {"UU": 10})
        result = subprocess.run(command, capture_output=True, text=True)
        self.assertEqual(result.returncode, 2)
        self.assertIn("UU/UR/RU retention check failed", result.stderr)
        self.assertEqual(result.stdout, "")

    def test_skipped_dedup_checks_all_parse_counts_without_requiring_dedup_stats(self):
        parsed = self.root / "sample.parse.stats"
        self.write_stats(parsed, {"UU": 100, "UR": 30, "RU": 20, "MM": 4})
        self.write_stats(self.selected, {"UU": 100, "UR": 30, "RU": 20})
        self.assertFalse(self.dedup.exists())
        result = subprocess.run(
            [sys.executable, str(SCRIPT), "--before", str(parsed),
             "--selected", str(self.selected)],
            capture_output=True, text=True, check=True,
        )
        self.assertIn("TOTAL\t150\t150\n", result.stdout)
        self.write_stats(self.selected, {"UU": 10, "UR": 3, "RU": 2})
        with self.assertRaisesRegex(ValueError, "UU: expected 100, retained 10"):
            RETENTION.validate_retention(parsed, self.selected)

    def test_cooler_preserves_repeated_coordinate_counts(self):
        from click.testing import CliRunner
        import cooler
        from cooler.cli import cli

        chroms = self.root / "chrom.sizes"
        chroms.write_text("chr1\t100000\n")
        pairs = self.root / "repeated.pairs"
        pairs.write_text(
            "#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type\n"
            "r1\tchr1\t100\tchr1\t50100\t+\t-\tUU\n"
            "r2\tchr1\t100\tchr1\t50100\t+\t-\tUU\n"
            "r3\tchr1\t200\tchr1\t50200\t+\t-\tUR\n"
        )
        matrix = self.root / "contacts.cool"
        result = CliRunner().invoke(
            cli, ["cload", "pairs", "-c1", "2", "-p1", "3", "-c2", "4", "-p2", "5",
                  f"{chroms}:5000", str(pairs), str(matrix)],
        )
        self.assertEqual(result.exit_code, 0, result.output)
        self.assertEqual(cooler.Cooler(str(matrix)).info["sum"], 3)


if __name__ == "__main__":
    unittest.main()
