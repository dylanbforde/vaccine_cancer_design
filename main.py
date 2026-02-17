import os
import logging
import sys
from model import training
from vaccine_design import pipeline

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)


def check_data_availability():
    """Check if required data files exist"""
    required_files = {
        "Mutation Data": "data/brca_tcga_pan_can_atlas_2018/data_mutations.txt",
        "IEDB Data": "data/tcell_full_v3.csv",
    }

    missing = []
    for name, path in required_files.items():
        if not os.path.exists(path):
            missing.append((name, path))

    return missing


def main():
    logging.info("Starting Vaccine Design Pipeline...")

    # 1. Data Check
    missing_data = check_data_availability()
    if missing_data:
        logging.error("Missing required data files:")
        for name, path in missing_data:
            logging.error(f"  - {name}: {path}")
        logging.info("Please download the data and place it in the 'data' directory.")
        logging.info("Refer to README.md for data sources.")
        return

    # 2. Data Processing Check
    processed_mutations_file = "data/mutations_processed.csv"
    if not os.path.exists(processed_mutations_file):
        logging.error(
            "Processed mutations file not found: data/mutations_processed.csv"
        )
        logging.info("Please run the data processing script first:")
        logging.info("  uv run process_data.py")
        return
    else:
        logging.info(f"Using processed mutations from: {processed_mutations_file}")

    # 3. Model Training (GNN for MHC Binding)
    model_path = "model/trained_model.pt"
    if not os.path.exists(model_path):
        logging.info("Step 2/3: Training MHC Binding Model...")
        try:
            training.main()
        except Exception as e:
            logging.error(f"Failed during model training: {e}")
            return
    else:
        logging.info(
            f"Step 2/3: Model already trained ({model_path} exists). Skipping."
        )

    # 4. Vaccine Design (Prediction & Selection)
    logging.info("Step 3/3: Designing Vaccine Candidates...")
    try:
        pipeline.main()
    except Exception as e:
        logging.error(f"Failed during vaccine design: {e}")
        return

    logging.info("Pipeline completed successfully!")
    logging.info(
        "Check 'data/vaccine_candidates.csv' and 'docs/vaccine_report.txt' for results."
    )


if __name__ == "__main__":
    main()
