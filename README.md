# Personalized Cancer Vaccine Design with Graph AI

This project implements an end-to-end pipeline for designing personalized cancer vaccines using Graph Neural Networks (GNNs). It takes patient-specific mutation data, identifies valid neoantigens, and predicts their binding affinity to MHC molecules to recommend optimal vaccine candidates.

## The Pipeline

The workflow consists of three main stages:

### 1. Data Processing & Neoantigen Generation
Raw mutation data (from TCGA) is processed to generate candidate neoantigen peptides.
- **Variant Filtering**: Selects relevant somatic mutations (Missense, Frameshift, In-Frame Del/Ins).
- **Sequence Retrieval**: Fetches Wildtype protein sequences from **UniProt**.
- **Biologically Accurate Mutation Simulation**:
    - **Frameshifts**: Fetches **Coding DNA Sequences (CDS)** from **Ensembl** to simulate frameshifts at the genetic level, ensuring accurate downstream protein translation.
    - **Range Mutations**: Correctly parses and handles complex mutations like large deletions, insertions, and duplications.
- **Peptide Generation**: Extracts precise 9-mer peptides centered around the mutation site, ensuring no artificial padding is introduced.

### 2. GNN-based MHC Binding Prediction
A Graph Neural Network (GNN) model predicts the binding affinity between peptide sequences and MHC molecules.
- **Graph Construction**: Peptides are represented as graphs where nodes are amino acids (with biophysical features like hydrophobicity, charge) and edges represent spatial/sequential relationships.
- **Training**: The model is trained on validated epitope-MHC binding data from the **IEDB** database.

### 3. Vaccine Candidate Selection
The trained model evaluates the generated neoantigens.
- **Encoding**: Patient neoantigens are converted into graph structures.
- **Prediction**: The model predicts the probability of MHC binding.
- **Ranking**: Candidates are filtered and ranked by binding probability to select the most promising targets for vaccine formulation.

## Prerequisites

### Data
This pipeline requires external datasets which must be placed in the `data/` directory:
1.  **Mutation Data**: `data/brca_tcga_pan_can_atlas_2018/data_mutations.txt` (TCGA PanCan Atlas)
2.  **IEDB Data**: `data/tcell_full_v3.csv` (Immune Epitope Database)

A helper script is provided to automatically download and extract these files:
```bash
uv run download_data.py
```

### Environment
The project is managed using `uv` for dependency management and execution.

## Usage

To run the full pipeline, execute the following steps in order:

1.  **Download Data**:
    ```bash
    uv run download_data.py
    ```

2.  **Process Data** (Parses mutations and generates peptides):
    ```bash
    uv run process_data.py
    ```
    *Output: `data/mutations_processed.csv`*

3.  **Run Pipeline** (Train Model & Predict Candidates):
    ```bash
    uv run main.py
    ```
    *Output: `model/trained_model.pt`, `data/vaccine_candidates.csv`, `docs/vaccine_report.txt`*

## Project Structure
- `download_data.py`: Helper script to fetch TCGA and IEDB datasets.
- `process_data.py`: ETL script for mutation parsing and peptide generation.
- `data_processing/`: Logic for parsing HGVS mutations, fetching CDS/Protein sequences.
- `model/`: GNN model definition (`training.py`) and training logic.
- `vaccine_design/`: Inference pipeline and candidate ranking.
- `main.py`: Entry point for model training and vaccine design.
