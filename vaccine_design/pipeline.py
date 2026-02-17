import torch
import torch.nn as nn
from torch_geometric.data import Data
from torch_geometric.nn import GCNConv, global_mean_pool
import pandas as pd
import logging


class PeptideEncoder:
    """Encodes peptide sequences into numerical features for GNN input"""

    def __init__(self):
        # Amino acid properties: hydrophobicity, size, charge, etc.
        self.aa_features = {
            "A": [0.62, 88.6, 0.0],  # Alanine
            "R": [-2.53, 173.4, 1.0],  # Arginine
            "N": [-0.78, 114.1, 0.0],  # Asparagine
            "D": [-0.90, 111.1, -1.0],  # Aspartic Acid
            "C": [0.29, 108.5, 0.0],  # Cysteine
            "Q": [-0.85, 143.8, 0.0],  # Glutamine
            "E": [-0.74, 138.4, -1.0],  # Glutamic Acid
            "G": [0.48, 60.1, 0.0],  # Glycine
            "H": [-0.40, 153.2, 0.1],  # Histidine
            "I": [1.38, 166.7, 0.0],  # Isoleucine
            "L": [1.06, 166.7, 0.0],  # Leucine
            "K": [-1.50, 168.6, 1.0],  # Lysine
            "M": [0.64, 162.9, 0.0],  # Methionine
            "F": [1.19, 189.9, 0.0],  # Phenylalanine
            "P": [0.12, 112.7, 0.0],  # Proline
            "S": [-0.18, 89.0, 0.0],  # Serine
            "T": [-0.05, 116.1, 0.0],  # Threonine
            "W": [0.81, 227.8, 0.0],  # Tryptophan
            "Y": [0.26, 193.6, 0.0],  # Tyrosine
            "V": [1.08, 140.0, 0.0],  # Valine
            "X": [0.0, 0.0, 0.0],  # Padding/Unknown
        }

        # Normalize features (Z-score) to ensure model stability
        # Extract all feature vectors (excluding padding 'X' for calculation stats if desired,
        # but including it effectively assumes it's "average-ish" or 0. Let's exclude X for stats).
        valid_aas = [k for k in self.aa_features if k != "X"]
        raw_values = torch.tensor([self.aa_features[k] for k in valid_aas])

        mean = raw_values.mean(dim=0)
        std = raw_values.std(dim=0)

        # Avoid division by zero
        std[std == 0] = 1.0

        for k, v in self.aa_features.items():
            if k == "X":
                self.aa_features[k] = [0.0, 0.0, 0.0]  # Keep padding as 0
            else:
                t = torch.tensor(v)
                normalized = (t - mean) / std
                self.aa_features[k] = normalized.tolist()

    def encode_peptide(self, peptide):
        """Convert peptide sequence to node features"""
        features = []
        for aa in peptide:
            if aa in self.aa_features:
                features.append(self.aa_features[aa])
            else:
                features.append(self.aa_features["X"])
        return torch.tensor(features, dtype=torch.float)

    def create_edge_index(self, peptide_length):
        """Create edge connections between amino acids"""
        # Create edges between adjacent residues
        edges = []
        for i in range(peptide_length):
            for j in range(i + 1, peptide_length):
                # Connect residues within 3 positions of each other
                if j - i <= 3:
                    edges.append([i, j])
                    edges.append([j, i])  # Bidirectional edges
        return torch.tensor(edges, dtype=torch.long).t()


class PeptideMHCPredictor(nn.Module):
    """GNN model for predicting peptide-MHC binding"""

    def __init__(self, input_dim=3, hidden_dim=64, num_layers=2, dropout_rate=0.2):
        super().__init__()

        self.convs = nn.ModuleList()
        self.convs.append(GCNConv(input_dim, hidden_dim))
        for _ in range(num_layers - 1):
            self.convs.append(GCNConv(hidden_dim, hidden_dim))

        self.batch_norms = nn.ModuleList()
        for _ in range(num_layers):
            self.batch_norms.append(nn.BatchNorm1d(hidden_dim))

        self.fc = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout_rate),
            nn.Linear(hidden_dim, 1),
            nn.Sigmoid(),
        )

    def forward(self, data):
        x, edge_index = data.x, data.edge_index

        for conv, batch_norm in zip(self.convs, self.batch_norms):
            x = torch.relu(batch_norm(conv(x, edge_index)))

        x = global_mean_pool(x, data.batch)
        return self.fc(x)


class VaccineDesignPipeline:
    """End-to-end pipeline for neoantigen prediction and vaccine design"""

    def __init__(self, model_path):
        self.peptide_encoder = PeptideEncoder()
        self.model = PeptideMHCPredictor()  # Use default or configurable params
        self.model.load_state_dict(torch.load(model_path))
        self.model.eval()

    def process_mutations(self, mutations_df):
        """Process mutation data and generate peptide candidates"""
        processed = []
        for _, row in mutations_df.iterrows():
            peptide = row["peptide"]
            if pd.isna(peptide):
                continue

            try:
                x = self.peptide_encoder.encode_peptide(peptide)
                edge_index = self.peptide_encoder.create_edge_index(len(peptide))
                graph_data = Data(
                    x=x,
                    edge_index=edge_index,
                    peptide=peptide,
                    mutation_info={
                        "gene": row["Hugo_Symbol"],
                        "sample": row["Tumor_Sample_Barcode"],
                        "position": row["pos"],
                        "mutation": row["alt"],
                    },
                )
                processed.append(graph_data)
            except Exception as e:
                logging.error(f"Error processing peptide {peptide}: {e}")

        return processed

    def predict_binding(self, graph_data_list):
        """Predict MHC binding probability for peptide candidates"""
        predictions = []

        with torch.no_grad():
            for data in graph_data_list:
                # Add batch dimension
                data.batch = torch.zeros(data.x.size(0), dtype=torch.long)

                # Get prediction
                pred = self.model(data)

                predictions.append(
                    {
                        "peptide": data.peptide,
                        "binding_score": pred.item(),
                        **data.mutation_info,
                    }
                )

        return pd.DataFrame(predictions)

    def design_vaccine(self, mutations_df, binding_threshold=0.5):
        """Design personalized vaccine from mutation data"""
        # Process mutations into graph format
        graph_data = self.process_mutations(mutations_df)

        # Predict binding scores
        predictions = self.predict_binding(graph_data)

        # Filter and rank candidates
        candidates = predictions[predictions["binding_score"] >= binding_threshold]
        candidates = candidates.sort_values("binding_score", ascending=False)

        return candidates


def main(
    mutations_file="data/mutations_processed.csv",
    model_path="model/trained_model.pt",
    candidates_file="data/vaccine_candidates.csv",
    report_file="docs/vaccine_report.txt",
    binding_threshold=0.5,
):
    # Set up logging
    logging.basicConfig(level=logging.INFO)

    # Load processed mutations
    mutations = pd.read_csv(mutations_file)

    # Initialize pipeline
    pipeline = VaccineDesignPipeline(model_path=model_path)

    # Design vaccine candidates
    candidates = pipeline.design_vaccine(mutations, binding_threshold=binding_threshold)

    # Save results
    candidates.to_csv(candidates_file, index=False)
    logging.info(f"Generated {len(candidates)} vaccine candidates")

    # Generate report
    top_candidates = candidates.head(10)
    report = f"""
    Vaccine Design Report
    ---------------------
    Mutations File: {mutations_file}
    Model Path: {model_path}
    Binding Threshold: {binding_threshold:.3f}
    Total candidates analyzed: {len(mutations)}
    Promising candidates: {len(candidates)}

    Top 10 Candidates:
    {top_candidates.to_string()}

    Statistics:
    - Mean binding score: {candidates["binding_score"].mean():.3f}
    - Median binding score: {candidates["binding_score"].median():.3f}
    - Unique genes targeted: {len(candidates["gene"].unique())}
    - Binding Score Distribution:
        {candidates["binding_score"].describe().to_string()}
    """

    with open(report_file, "w") as f:
        f.write(report)


if __name__ == "__main__":
    main()
