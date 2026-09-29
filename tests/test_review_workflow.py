"""Isolated smoke tests; RepeatMasker is a deterministic stub, never a biology test.

Run with Python unittest in the analysis environment. SNAKEMAKE may name an
external executable; all outputs and synthetic reference data live in /tmp.
"""

import csv
import json
import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SNAKEMAKE = os.environ.get("SNAKEMAKE", shutil.which("snakemake") or "")
HEADER = "   SW  perc perc perc  query position in query matching repeat position in repeat\n score div. del. ins. sequence begin end (left) repeat class/family begin end (left) ID\n\n"


@unittest.skipUnless(SNAKEMAKE, "Snakemake executable required (set SNAKEMAKE)")
class WorkflowSmokeTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="te-review-")
        self.addCleanup(self.tmp.cleanup)
        self.work = Path(self.tmp.name)
        shutil.copy(ROOT / "Snakefile", self.work / "Snakefile")
        shutil.copytree(ROOT / "workflow", self.work / "workflow")
        (self.work / "pc.txt").write_text("tx1\n")
        (self.work / "lnc.txt").write_text("tx2\n")
        (self.work / "transcripts.fa").write_text(">tx1\nACGTAC\n>tx2\nACGTACGT\n")
        (self.work / "genome.fa").write_text(">chr1\n" + "ACGT" * 12 + "\n")
        records = []
        for tid, start, end, strand, biotype, exons in [
            ("tx1", 1, 8, "+", "protein_coding", [(1, 3), (6, 8)]),
            ("tx2", 13, 24, "-", "lncRNA", [(13, 16), (21, 24)]),
        ]:
            attrs = f'gene_id "g{tid}"; transcript_id "{tid}"; transcript_type "{biotype}"; gene_type "{biotype}";'
            for feature, lo, hi in [("transcript", start, end)] + [
                ("exon", a, b) for a, b in exons
            ]:
                records.append(
                    f"chr1\tfixture\t{feature}\t{lo}\t{hi}\t.\t{strand}\t.\t{attrs}\n"
                )
        (self.work / "annotation.gtf").write_text("".join(records))
        self.config = dict(
            gencode_gtf="annotation.gtf",
            gencode_fasta="transcripts.fa",
            genome_fasta="genome.fa",
            classification_mode="id_file",
            pc_transcript_ids_file="pc.txt",
            lncrna_transcript_ids_file="lnc.txt",
            sequence_mode="spliced",
            n_chunks=3,
            datasets=["toy"],
            threads=1,
        )
        bindir = self.work / "bin"
        bindir.mkdir()
        stub = bindir / "RepeatMasker"
        stub.write_text(
            "#!/usr/bin/env python\n"
            + """import os, pathlib, sys
args = sys.argv[1:]
source = pathlib.Path(args[-1])
base = pathlib.Path(args[args.index('-dir') + 1]) / source.name
ids = [line[1:].split()[0] for line in source.read_text().splitlines() if line.startswith('>')]
hits = [] if os.environ.get('TE_TEST_NO_HITS') else [tid for tid in ids if tid.startswith('tx1')]
header = """
            + repr(HEADER)
            + """
if hits:
    pathlib.Path(str(base) + '.out').write_text(header + ''.join(f' 100 0.0 0.0 0.0 {tid} 1 3 (3) + Alu SINE/Alu 1 3 (0) 1\\n' for tid in hits))
    pathlib.Path(str(base) + '.out.gff').write_text('##gff-version 2\\n' + ''.join(f'{tid}\\tRepeatMasker\\tsimilarity\\t1\\t3\\t100\\t+\\t.\\tTarget "Alu"\\n' for tid in hits))
else:
    pathlib.Path(str(base) + '.out').write_text('There were no repetitive sequences detected.\\n')
"""
        )
        stub.chmod(0o755)
        self.env = dict(
            os.environ,
            PATH=str(bindir) + os.pathsep + os.environ["PATH"],
            XDG_CACHE_HOME=str(self.work / "cache"),
        )

    def run_workflow(self, *extra, success=True):
        (self.work / "config.json").write_text(json.dumps(self.config))
        command = [
            SNAKEMAKE,
            "--snakefile",
            "Snakefile",
            "--configfile",
            "config.json",
            "--cores",
            "1",
            "--workflow-profile",
            "none",
            "--profile",
            "none",
            "pipeline_all_features",
            *extra,
        ]
        result = subprocess.run(
            command,
            cwd=self.work,
            env=self.env,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            timeout=120,
        )
        if success:
            self.assertEqual(result.returncode, 0, result.stdout)
        else:
            self.assertNotEqual(result.returncode, 0, result.stdout)
        return result.stdout

    def test_mode_dependencies_and_feature_target(self):
        for mode, fasta in [("spliced", "transcripts.fa"), ("unspliced", "genome.fa")]:
            with self.subTest(mode=mode):
                self.config["sequence_mode"] = mode
                out = self.run_workflow("--dry-run")
                self.assertIn("extract_all_features", out)
                self.assertIn(fasta, out)
                self.assertIn("pc.txt", out)
                self.assertIn("lnc.txt", out)
                self.assertNotIn("rule statistical_analysis:", out)

    def test_invalid_chunk_counts_rejected(self):
        for count in [0, -1, 1.5, True, "2"]:
            with self.subTest(n_chunks=count):
                self.config["n_chunks"] = count
                self.assertIn("n_chunks", self.run_workflow("--dry-run", success=False))

    def test_source_gtf_subset_with_fasta_classification(self):
        source = self.work / "annotation.gtf"
        source.write_text(
            source.read_text()
            + 'chr1\tfixture\ttranscript\t30\t40\t.\t+\t.\tgene_id "g3"; transcript_id "excluded"; transcript_type "lncRNA";\n'
        )
        original_source = source.read_text()
        (self.work / "pc.fa").write_text(">tx1|gene1 description\nACGTAC\n")
        (self.work / "lnc.fa").write_text(">tx2 description\nACGTACGT\n")
        self.config.update(
            source_gtf="annotation.gtf",
            gencode_gtf="subset.gtf",
            classification_mode="fasta",
            pc_transcripts_fasta="pc.fa",
            lncrna_transcripts_fasta="lnc.fa",
        )
        self.run_workflow()
        self.assertEqual(source.read_text(), original_source)
        self.assertNotIn("excluded", (self.work / "subset.gtf").read_text())
        ids = list(
            (self.work / "resources/annotation").glob("subset_transcript_ids.*.txt")
        )
        self.assertEqual(len(ids), 1)
        self.assertEqual(set(ids[0].read_text().splitlines()), {"tx1", "tx2"})
        with (
            self.work / "results/toy/features/all_transcripts_te_features.csv"
        ).open() as handle:
            self.assertEqual(
                {row["transcript_id"] for row in csv.DictReader(handle)}, {"tx1", "tx2"}
            )

    def test_all_no_hit_features_retain_transcripts(self):
        self.env["TE_TEST_NO_HITS"] = "1"
        self.run_workflow()
        with (
            self.work / "results/toy/features/all_transcripts_te_features.csv"
        ).open() as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual({row["transcript_id"] for row in rows}, {"tx1", "tx2"})
        self.assertTrue(all(float(row["global_rm_count"]) == 0 for row in rows))

    @unittest.skipUnless(
        shutil.which("bedtools"), "bedtools required for integration fixture"
    )
    def test_feature_extraction_both_modes_with_empty_and_no_hit_chunks(self):
        for mode, lengths in [
            ("spliced", {"tx1": 6, "tx2": 8}),
            ("unspliced", {"tx1": 8, "tx2": 12}),
        ]:
            with self.subTest(mode=mode):
                self.config["sequence_mode"] = mode
                self.run_workflow()
                dataset = "toy" if mode == "spliced" else "toy_unspliced"
                with (
                    self.work
                    / "results"
                    / dataset
                    / "features/all_transcripts_te_features.csv"
                ).open() as handle:
                    rows = {row["transcript_id"]: row for row in csv.DictReader(handle)}
                self.assertEqual(set(rows), set(lengths))
                for tid, length in lengths.items():
                    self.assertEqual(float(rows[tid]["transcript_length"]), length)
                self.assertEqual(float(rows["tx1"]["te_count"]), 1)
                self.assertEqual(float(rows["tx2"]["te_count"]), 0)
                self.assertFalse(
                    (self.work / "results" / dataset / "analysis").exists()
                )
                chunks = self.work / "results" / dataset / "repeatmasker/chunks"
                self.assertEqual((chunks / "chunk_2.fa").read_text(), "")
                self.assertTrue((chunks / "chunk_1.out").exists())
                self.assertTrue((chunks / "chunk_2.out.gff").exists())
                # A same-path input content change must invalidate completed jobs.
                fasta = self.work / (
                    "transcripts.fa" if mode == "spliced" else "genome.fa"
                )
                fasta.write_text(fasta.read_text().replace("ACGT", "TCGT", 1))
                rerun = self.run_workflow("--dry-run")
                self.assertIn("rule prepare_transcript_fasta:", rerun)
                self.assertIn("rule extract_all_features:", rerun)
                pc_ids = self.work / "pc.txt"
                pc_ids.write_text(pc_ids.read_text() + "additional_fixture_id\n")
                rerun = self.run_workflow("--dry-run")
                self.assertIn("rule extract_ids_from_fasta:", rerun)
                self.assertIn("rule extract_all_features:", rerun)


if __name__ == "__main__":
    unittest.main()
