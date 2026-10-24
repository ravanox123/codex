import tempfile
import unittest
from pathlib import Path

from dna_lab.core import Sequence, analyze, read_fasta, reverse_complement, translate


class AnalysisTests(unittest.TestCase):
    def test_translation(self):
        self.assertEqual(translate("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"), "MAIVMGR*KGAR*")
        self.assertEqual(translate("ATGNNNTAAA"), "MX*")
        self.assertEqual(translate("AATGGCC", 1), "MA")

    def test_ambiguity_complements(self):
        self.assertEqual(reverse_complement("ACGTRYSWKMBDHVN"), "NBDHVKMWSRYACGT")

    def test_gc_and_overlapping_motif(self):
        report = analyze(Sequence("test", "AAAAGCNN"), "AA")
        self.assertEqual(report["motif_positions"], [1, 2, 3])
        self.assertEqual(report["gc_percent"], 33.33)
        self.assertEqual(report["ambiguous_bases"], 2)
        self.assertIsNone(analyze(Sequence("unknown", "NNN"))["gc_percent"])

    def test_fasta_parsing_and_errors(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "input.fasta"
            path.write_text(">first\natg\nccn\n\n>second\nTAA\n")
            self.assertEqual(read_fasta(path), [Sequence("first", "ATGCCN"), Sequence("second", "TAA")])
            for invalid in ("ATG", ">empty\n", ">bad\nATZ", "", ">\nATG"):
                path.write_text(invalid)
                with self.assertRaises(ValueError):
                    read_fasta(path)


if __name__ == "__main__":
    unittest.main()
