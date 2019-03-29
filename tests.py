import unittest
from utilities import UniprotClient
from exceptions import PopulationError, ProteinNotFoundError
from bioinformatics_stronghold import (valid_dna, valid_rna,
                                       calculate_expected_offspring)


class TestUniprotClient(unittest.TestCase):
    """
    Tests the UniprotClient class
    """
    def setUp(self):
        self.valid_id = "B5ZC00"
        self.invalid_id = "BA0123456789AB"
        self.u = UniprotClient()

    def test_invalid_accession_id(self):
        """
        Invalid protein IDs should raise ProteinNotFoundError
        """
        with self.assertRaises(ProteinNotFoundError):
            self.u.get_protein(self.invalid_id)

    def test_valid_return(self):
        """
        Valid protein IDs should return a dict of non-empty strings
        """
        protein_data = self.u.get_protein(self.valid_id)
        self.assertIs(type(protein_data), dict)

        for key, val in protein_data.items():
            self.assertIs(type(val), str)
            self.assertTrue(val)


class TestValidPopulation(unittest.TestCase):

    def test_calculate_expected_offspring_input_handling(self):
        pop_list_short = [0, 2]
        pop_list_long = [24, 6, 1, 1, 4, 33, 75, 192, 12]

        with self.assertRaises(PopulationError):
            calculate_expected_offspring(pop_list_short)
            self.fail("population cannot be fewer than six values")

        with self.assertRaises(PopulationError):
            calculate_expected_offspring(pop_list_short)
            self.fail("population cannot be more than six values")


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
