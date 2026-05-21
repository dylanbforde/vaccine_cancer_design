import os
import tempfile
import unittest

import torch
from torch_geometric.data import Data

from model.training import process_iedb_data
from vaccine_design.pipeline import (
    PeptideEncoder,
    PeptideMHCPredictor,
    SinkhornPeptideAlignment,
)


class TestIEDBProcessing(unittest.TestCase):
    def test_process_iedb_data_reads_two_header_export(self):
        rows = [
            [
                "Epitope",
                "MHC Restriction",
                "Assay",
                "Assay",
                "Other",
            ],
            [
                "Name",
                "Name",
                "Qualitative Measurement",
                "Method",
                "Ignored",
            ],
            [
                "SLYNTVATL",
                "HLA-A*02:01",
                "Positive",
                "MHC binding assay",
                "x",
            ],
            [
                "TOO-LONG-PEPTIDE",
                "HLA-A*02:01",
                "Positive",
                "MHC binding assay",
                "x",
            ],
            [
                "GILGFVFTL",
                "H2-Kb",
                "Negative",
                "MHC binding assay",
                "x",
            ],
        ]

        with tempfile.NamedTemporaryFile("w", delete=False, newline="") as handle:
            path = handle.name
            for row in rows:
                handle.write(",".join(row) + "\n")

        try:
            data = process_iedb_data(path)
        finally:
            os.unlink(path)

        self.assertEqual(len(data), 1)
        self.assertEqual(data.iloc[0]["Epitope"], "SLYNTVATL")
        self.assertEqual(data.iloc[0]["binding"], 1)


class TestSinkhornAlignment(unittest.TestCase):
    def test_transport_plan_is_doubly_stochastic_for_nine_mer(self):
        torch.manual_seed(0)
        layer = SinkhornPeptideAlignment(
            hidden_dim=8,
            num_motifs=2,
            motif_length=9,
            n_iters=20,
            epsilon=0.7,
        )
        node_features = torch.randn(9, 8)
        motifs = layer._motifs_for_length(9)
        plan, _scores = layer.transport_plan(node_features, motifs)

        expected = torch.ones(2, 9)
        self.assertTrue(torch.allclose(plan.sum(dim=-1), expected, atol=1e-4))
        self.assertTrue(torch.allclose(plan.sum(dim=-2), expected, atol=1e-4))

    def test_predictor_forward_uses_alignment_features(self):
        encoder = PeptideEncoder()
        x = encoder.encode_peptide("SLYNTVATL")
        edge_index = encoder.create_edge_index(x.size(0))
        data = Data(
            x=x,
            edge_index=edge_index,
            batch=torch.zeros(x.size(0), dtype=torch.long),
        )
        model = PeptideMHCPredictor(hidden_dim=8, num_alignment_motifs=2)
        model.eval()

        with torch.no_grad():
            out = model(data)

        self.assertEqual(tuple(out.shape), (1, 1))
        self.assertGreaterEqual(float(out.item()), 0.0)
        self.assertLessEqual(float(out.item()), 1.0)


if __name__ == "__main__":
    unittest.main()
