#! /usr/bin/env python

import unittest
from mutyper.ancestor import Ancestor


class TestAncestor(unittest.TestCase):
    def setUp(self):
        self.anc = Ancestor("tests/test_data/ancestor.fa")
        self.anc_stranded = Ancestor(
            "tests/test_data/ancestor.fa",
            strand_file="tests/test_data/strandedness.bed",
        )

    def test_seq(self):
        anc_seq = self.anc["foo"][:].seq
        anc_seq_true = "AAACCCgggTTT"
        self.assertEqual(anc_seq, anc_seq_true)

    def test_context(self):
        self.assertEqual(
            list(self.anc.region_contexts("foo")),
            [
                None,
                "AAA",
                "AAC",
                "ACC",
                "CCC",
                None,
                None,
                None,
                None,
                None,
                "AAA",
                None,
            ],
        )

        # same as above but using strand polarized ancestor
        self.assertEqual(
            list(self.anc_stranded.region_contexts("foo")),
            [
                None,
                "AAA",
                "AAC",
                "ACC",
                "GGG",
                None,
                None,
                None,
                None,
                None,
                "AAA",
                None,
            ],
        )

        self.assertEqual(list(self.anc.region_contexts("foo", 1, 3)), ["AAA", "AAC"])

        self.assertEqual(list(self.anc.region_contexts("foo", 9, 11)), [None, "AAA"])

    def test_mutation_type(self):
        self.assertEqual(self.anc.mutation_type("foo", 1, "A", "T"), ("AAA", "ATA"))
        # infinite sites violation
        self.assertEqual(self.anc.mutation_type("foo", 1, "C", "T"), (None, None))
        # low confidence (lower case) ancestral state
        self.assertEqual(self.anc.mutation_type("foo", 7, "G", "T"), (None, None))
        # reverse complement
        self.assertEqual(self.anc.mutation_type("foo", 10, "G", "T"), ("AAA", "ACA"))
        # same as above but using strand polarized ancestor
        self.assertEqual(
            self.anc_stranded.mutation_type("foo", 10, "G", "T"), ("AAA", "ACA")
        )


if __name__ == "__main__":
    unittest.main()
