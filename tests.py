import unittest
from bioinformatics_stronghold import valid_dna, valid_rna


class TestValidDNASequences(unittest.TestCase):

    def test_valid_dna(self):
        dna_sequence = "ATTGCCATG"
        msg = "Strings containing only A, T, C, and G are valid"
        self.assertTrue(valid_dna(dna_sequence), msg)

    def test_invalid_dna_bases(self):
        dna_sequence = "JTIENGKDI"
        msg = "Strings with bases other than A, T, C, and G are invalid"
        self.assertFalse(valid_dna(dna_sequence), msg)

    def test_invalid_dna_type(self):
        dna_sequence = 2829
        msg = "Sequences that are not strings are invalid"
        self.assertFalse(valid_dna(dna_sequence), msg)

    def test_valid_rna(self):
        rna_sequence = "AUUGCCAUG"
        msg = "Strings containing only A, U, C, and G are valid"
        self.assertTrue(valid_rna(rna_sequence), msg)

    def test_invalid_rna_bases(self):
        rna_sequence = "JTIENGKDI"
        msg = "Strings with bases other than A, U, C, and G are invalid"
        self.assertFalse(valid_rna(rna_sequence), msg)

    def test_invalid_rna_type(self):
        rna_sequence = 2829
        msg = "Sequences that are not strings are invalid"
        self.assertFalse(valid_rna(rna_sequence), msg)


if __name__ == '__main__':
    unittest.main()
