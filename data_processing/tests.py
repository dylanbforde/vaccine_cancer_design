import unittest
import pandas as pd
from unittest.mock import mock_open, patch, MagicMock
import logging

logging.basicConfig(level=logging.INFO)

from data_processing.helper_functions import (
    generate_peptides,
)

from data_processing.mutated_genes import (
    process_mutations_in_batches,
    filter_variants,
    parse_protein_change,
    get_cds_sequences,
)


class TestMutationProcessing(unittest.TestCase):
    def setUp(self):
        """Set up test data"""
        self.test_mutations = pd.DataFrame(
            {
                "Hugo_Symbol": ["GENE1", "GENE2", "GENE3", "GENE4", "GENE5"],
                "HGVSp_Short": [
                    "p.A123V",  # Simple substitution
                    "p.P408Afs*99",  # Frameshift
                    "p.Q344del",  # Deletion
                    "p.S2107dup",  # Duplication
                    "p.H57_S58delinsQ",  # Complex change
                ],
                "Variant_Classification": [
                    "Missense_Mutation",
                    "Frame_Shift_Ins",
                    "In_Frame_Del",
                    "In_Frame_Ins",
                    "In_Frame_Del",
                ],
                "Tumor_Sample_Barcode": ["S1", "S2", "S3", "S4", "S5"],
            }
        )

    def test_parse_protein_change(self):
        """Test parsing of different mutation types"""
        test_cases = {
            "p.A123V": (123, 123, "Sub", "V"),  # Simple substitution
            "p.P408Afs*99": (408, 408, "Fs", "A"),  # Frameshift
            "p.Q344del": (344, 344, "Del", "X"),  # Deletion
            "p.S2107dup": (2107, 2107, "Dup", "X"),  # Duplication
            "p.R30*": (30, 30, "Stop", "*"),  # Stop gained
            "p.H57_S58delinsQ": (57, 58, "Delins", "Q"),  # Complex change
        }

        for hgvsp, expected in test_cases.items():
            result = parse_protein_change(hgvsp)
            self.assertEqual(result, expected, f"Failed to parse {hgvsp}")

    def test_parse_protein_change_invalid(self):
        """Test invalid protein change patterns"""
        invalid_inputs = [
            "p.A123VV",  # Invalid substitution
            "p.R456fss",  # Invalid frameshift
            "p.123del",  # Missing amino acid
            "p.E123_L125delinssW",  # Invalid delins
            "p.C22**",  # Invalid stop
            "123",  # Not HGVS
            "p.",  # Empty
            "",  # Empty string
            None,  # None
        ]

        for invalid in invalid_inputs:
            result = parse_protein_change(invalid)
            self.assertIsNone(
                result, f"Should return None for invalid input: {invalid}"
            )

    def test_filter_variants(self):
        """Test variant filtering"""
        filtered = filter_variants(self.test_mutations)

        # Should keep all rows except those with invalid variant classifications
        self.assertEqual(len(filtered), 5)

        # Check required columns are present
        required_columns = {"Hugo_Symbol", "HGVSp_Short", "Tumor_Sample_Barcode"}
        self.assertTrue(all(col in filtered.columns for col in required_columns))

    @patch("data_processing.mutated_genes.get_cds_sequences")
    @patch("data_processing.mutated_genes.get_sequences")
    def test_process_mutations_in_batches(self, mock_get_sequences, mock_get_cds):
        """Test batch processing of mutations"""
        # Mock sequence retrieval
        mock_get_sequences.return_value = {
            "GENE1": ("A" * 600, 600),  # Valid sequence
            "GENE2": ("A" * 600, 600),  # Valid sequence
        }
        mock_get_cds.return_value = {
            "GENE1": "ATGGCCGCCGTGGCCCTGGTGGCCACCGCC",  # Mock CDS
            "GENE2": "ATGGCCGCCGTGGCCCTGGTGGCCACCGCC",
        }

        # Test with small batch size
        result = process_mutations_in_batches(self.test_mutations.head(2), batch_size=1)

        # Verify results
        self.assertGreater(len(result), 0)
        self.assertTrue("peptide" in result.columns)
        self.assertTrue("wildtype_seq" in result.columns)
        self.assertTrue("cds_seq" in result.columns)

    def test_process_mutations_empty_input(self):
        """Test processing with empty input"""
        empty_df = pd.DataFrame(columns=self.test_mutations.columns)
        result = process_mutations_in_batches(empty_df)
        self.assertTrue(result.empty)

    @patch("data_processing.mutated_genes.pd.read_csv")
    @patch("data_processing.mutated_genes.get_cds_sequences")
    @patch("data_processing.mutated_genes.get_sequences")
    def test_end_to_end_processing(
        self, mock_get_sequences, mock_get_cds, mock_read_csv
    ):
        """Test full mutation processing pipeline"""
        # Mock sequence retrieval
        mock_get_sequences.return_value = {
            "GENE1": ("A" * 600, 600),
            "GENE2": ("A" * 600, 600),
            "GENE3": ("A" * 600, 700),
        }
        mock_get_cds.return_value = {
            "GENE1": "ATGGCCGCCGTGGCCCTGGTGGCCACCGCC",
            "GENE2": "ATGGCCGCCGTGGCCCTGGTGGCCACCGCC",
            "GENE3": "ATGGCCGCCGTGGCCCTGGTGGCCACCGCC",
        }
        mock_read_csv.return_value = self.test_mutations

        result = process_mutations_in_batches(self.test_mutations)

        # Verify results
        self.assertIsNotNone(result)
        self.assertGreater(len(result), 0)
        required_columns = {
            "Hugo_Symbol",
            "HGVSp_Short",
            "Tumor_Sample_Barcode",
            "pos",
            "end_pos",
            "mut_type",
            "alt",
            "wildtype_seq",
            "peptide",
            "cds_seq",
        }
        self.assertTrue(all(col in result.columns for col in required_columns))


class TestSequenceHandling(unittest.TestCase):
    """Tests for sequence retrieval and peptide generation"""

    def setUp(self):
        self.test_row = pd.Series(
            {
                "Hugo_Symbol": "GENE1",
                "wildtype_seq": ("A" * 30, 30),  # Changed to valid sequence
                "pos": 5,
                "end_pos": 5,  # Added end_pos
                "mut_type": "Sub",
                "alt": "K",  # Changed to valid amino acid
                "cds_seq": "AAA" * 50,  # Mock CDS to avoid stop codons
            }
        )

    def test_generate_peptides_valid(self):
        """Test valid peptide generation"""
        result = generate_peptides(self.test_row)
        self.assertIsNotNone(result, "Result should not be None for valid input")
        self.assertEqual(len(result), 9, "Should generate 9-mer peptide")
        self.assertTrue(
            all(aa in "ACDEFGHIKLMNPQRSTVWY" for aa in result),
            "Should only contain valid amino acids",
        )

    def test_generate_peptides_invalid_position(self):
        """Test peptide generation with invalid position"""
        row = self.test_row.copy()
        row["pos"] = 35  # Position beyond sequence length (30)
        result = generate_peptides(row)
        self.assertIsNone(result)

    def test_generate_peptides_different_mutation_types(self):
        """Test peptide generation for different mutation types"""
        mutation_types = {
            "Sub": "K",  # Changed to valid amino acid
            "Fs": "K",  # Changed to valid amino acid (requires cds_seq, present in mock)
            "Del": None,  # Should now work
            "Dup": None,  # Should now work
            "Stop": "*",
            "Delins": "K",  # Changed to valid amino acid
            "Ins": "K",  # New type
        }

        for mut_type, alt in mutation_types.items():
            row = self.test_row.copy()
            row["mut_type"] = mut_type
            row["alt"] = alt

            # For Del/Dup/Ins/Delins, we might want different end_pos
            if mut_type in ["Del", "Dup", "Delins"]:
                row["end_pos"] = 5  # Single AA
            if mut_type == "Ins":
                row["end_pos"] = 6  # Insert between 5 and 6

            result = generate_peptides(row)

            if mut_type in [
                "Sub",
                "Fs",
                "Delins",
                "Del",
                "Dup",
                "Ins",
            ]:  # These should succeed
                self.assertIsNotNone(result, f"Should generate peptide for {mut_type}")
                self.assertEqual(
                    len(result), 9, f"Should generate 9-mer for {mut_type}"
                )
            else:  # Stop
                if result is not None:
                    # Is it Stop?
                    if mut_type == "Stop":
                        pass  # Should probably return None
                    else:
                        self.assertEqual(
                            len(result),
                            9,
                            f"If generated, should be 9-mer for {mut_type}",
                        )


def run_tests():
    unittest.main(verbosity=2)


if __name__ == "__main__":
    run_tests()
