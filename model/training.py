import pandas as pd
import torch
from torch_geometric.loader import DataLoader
from torch.optim import Adam
import numpy as np
from sklearn.model_selection import train_test_split
import logging
from tqdm import tqdm
from torch_geometric.data import Data

from vaccine_design.pipeline import PeptideEncoder, PeptideMHCPredictor


def process_iedb_data(file_path="data/tcell_full_v3.csv"):
    """Process IEDB T cell assay data for training"""
    df = pd.read_csv(file_path)
    data = df[["Epitope", "Qualitative Measurement", "MHC Restriction", "Assay"]].copy()

    # More efficient filtering
    data = data[
        (data["Epitope"].str.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]{9}"))
        & (data["MHC Restriction"].str.startswith("HLA", na=False))
    ]

    data["binding"] = data["Qualitative Measurement"].map(
        {"Positive": 1, "Negative": 0}
    )

    return data.dropna(subset=["binding", "Epitope"])


def create_training_data(iedb_data, peptide_encoder):
    """Convert peptides to graph data format"""
    graph_data = []

    # Optimization: use itertuples for faster row iteration
    for row in iedb_data.itertuples(index=False):
        peptide = getattr(row, "Epitope")
        x = peptide_encoder.encode_peptide(peptide)
        edge_index = peptide_encoder.create_edge_index(len(peptide))
        data = Data(
            x=x,
            edge_index=edge_index,
            y=torch.tensor([getattr(row, "binding")], dtype=torch.float),
        )
        graph_data.append(data)

    return graph_data


def train_model(train_loader, val_loader, model, epochs=50):
    """Train the GNN model"""
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = model.to(device)
    optimizer = Adam(model.parameters(), lr=0.001)
    criterion = torch.nn.BCELoss()

    best_val_loss = float("inf")
    best_model = None

    for epoch in range(epochs):
        # Training
        model.train()
        train_loss = 0
        for batch in tqdm(train_loader, desc=f"Epoch {epoch + 1}/{epochs}"):
            batch = batch.to(device)
            optimizer.zero_grad()

            out = model(batch)
            loss = criterion(out, batch.y.view(-1, 1))

            loss.backward()
            optimizer.step()
            train_loss += loss.item()

        # Validation
        model.eval()
        val_loss = 0
        predictions = []
        true_labels = []

        with torch.no_grad():
            for batch in val_loader:
                batch = batch.to(device)
                out = model(batch)
                val_loss += criterion(out, batch.y.view(-1, 1)).item()

                predictions.extend(out.cpu().numpy())
                true_labels.extend(batch.y.cpu().numpy())

        # Calculate metrics
        train_loss /= len(train_loader)
        val_loss /= len(val_loader)

        predictions = np.array(predictions)
        true_labels = np.array(true_labels)

        # Calculate AUC-ROC
        from sklearn.metrics import roc_auc_score

        auc_score = roc_auc_score(true_labels, predictions)

        logging.info(
            f"Epoch {epoch + 1}: "
            f"Train Loss: {train_loss:.4f}, "
            f"Val Loss: {val_loss:.4f}, "
            f"AUC-ROC: {auc_score:.4f}"
        )

        # Save best model
        if val_loss < best_val_loss:
            best_val_loss = val_loss
            best_model = model.state_dict()

    return best_model


def main():
    # Set up logging
    logging.basicConfig(level=logging.INFO)

    # Load and process IEDB data
    logging.info("Loading IEDB data...")
    iedb_data = process_iedb_data("data/tcell_full_v3.csv")

    # Create graph data
    logging.info("Creating graph representations...")
    peptide_encoder = PeptideEncoder()
    graph_data = create_training_data(iedb_data, peptide_encoder)

    # Split data
    train_data, val_data = train_test_split(graph_data, test_size=0.2, random_state=42)

    # Create data loaders
    train_loader = DataLoader(train_data, batch_size=32, shuffle=True)
    val_loader = DataLoader(val_data, batch_size=32)

    # Initialize and train model
    logging.info("Training model...")
    model = PeptideMHCPredictor()
    best_model = train_model(train_loader, val_loader, model)

    # Save trained model
    torch.save(best_model, "model/trained_model.pt")
    logging.info("Model saved successfully")


if __name__ == "__main__":
    main()
