"""Validate external-output boundary without running RepeatMasker."""

import importlib.util
import tempfile
import unittest
from pathlib import Path

MODULE = (
    Path(__file__).resolve().parents[1] / "workflow/scripts/repeatmasker_outputs.py"
)


class RepeatMaskerOutputsTests(unittest.TestCase):
    def setUp(self):
        spec = importlib.util.spec_from_file_location("repeatmasker_outputs", MODULE)
        self.helper = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(self.helper)
        tmp = tempfile.TemporaryDirectory(prefix="te-rm-boundary-")
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)
        self.raw = self.root / "raw.out"
        self.gff = self.root / "raw.gff"
        self.out = self.root / "normalized.out"
        self.out_gff = self.root / "normalized.gff"

    def normalize(self, **kwargs):
        self.helper.normalize_outputs(
            self.raw, self.gff, self.out, self.out_gff, **kwargs
        )

    def test_explicit_no_hits_and_empty_input(self):
        self.raw.write_text("There were no repetitive sequences detected.\n")
        self.normalize()
        self.assertIn("SW", self.out.read_text())
        self.assertIn("gff-version", self.out_gff.read_text())
        self.raw.unlink()
        self.normalize(empty=True)
        self.assertIn("SW", self.out.read_text())

    def test_missing_output_is_failure(self):
        with self.assertRaises((ValueError, FileNotFoundError, RuntimeError)):
            self.normalize()

    def test_hits_require_gff(self):
        self.raw.write_text(
            "   SW  perc perc perc\n score div. del. ins. query ID\n\n 100 0 0 0 tx1 1 3 (3) + Alu SINE/Alu 1 3 (0) 1\n"
        )
        with self.assertRaises((ValueError, FileNotFoundError, RuntimeError)):
            self.normalize()

    def test_malformed_or_empty_external_output_is_failure(self):
        for content in ["", "unexpected failure message\n", " 100 0 0 bad\n"]:
            with self.subTest(content=content):
                self.raw.write_text(content)
                with self.assertRaises(ValueError):
                    self.normalize()

    def test_headerless_first_hit_is_preserved(self):
        self.raw.write_text(" 100 0 0 0 tx1 1 3 (3) + Alu SINE/Alu 1 3 (0) 1\n")
        self.gff.write_text(
            '##gff-version 2\ntx1\tRepeatMasker\tsimilarity\t1\t3\t100\t+\t.\tTarget "Alu"\n'
        )
        self.normalize()
        self.assertEqual(self.out.read_text().count("tx1"), 1)

    def test_merge_preserves_first_hit_and_writes_one_header(self):
        header = "   SW  perc perc perc\n score div. del. ins. query ID\n\n"
        first = self.root / "first.out"
        second = self.root / "second.out"
        first.write_text(header + " 100 0 0 0 tx1 1 3 (3) + Alu SINE/Alu 1 3 (0) 1\n")
        second.write_text(header)
        gff1 = self.root / "first.gff"
        gff2 = self.root / "second.gff"
        gff1.write_text(
            '##gff-version 2\ntx1\tRepeatMasker\tsimilarity\t1\t3\t100\t+\t.\tTarget "Alu"\n'
        )
        gff2.write_text("##gff-version 2\n")
        self.helper.merge_outputs([first, second], [gff1, gff2], self.out, self.out_gff)
        self.assertEqual(self.out.read_text().count("tx1"), 1)
        self.assertEqual(self.out.read_text().count("SW"), 1)
        self.assertEqual(self.out_gff.read_text().count("gff-version"), 1)
        # A no-hit GFF2 placeholder must not downgrade a real GFF3 hit file.
        gff1.write_text(gff1.read_text().replace("gff-version 2", "gff-version 3"))
        self.helper.merge_outputs([second, first], [gff2, gff1], self.out, self.out_gff)
        self.assertTrue(self.out_gff.read_text().startswith("##gff-version 3\n"))
        self.assertEqual(self.out_gff.read_text().count("gff-version"), 1)


if __name__ == "__main__":
    unittest.main()
