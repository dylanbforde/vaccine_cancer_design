import logging
import sys
import os
from data_processing import mutated_genes

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)


def main():
    """
    Entry point for data processing pipeline.
    Parses mutations, fetches sequences, and generates peptides.
    """
    logging.info("Starting Data Processing...")

    # Check if input data exists
    input_file = "data/brca_tcga_pan_can_atlas_2018/data_mutations.txt"
    if not os.path.exists(input_file):
        logging.error(f"Input file not found: {input_file}")
        logging.error("Please run 'uv run download_data.py' first.")
        return

    try:
        mutated_genes.main()
        logging.info("Data processing completed successfully.")
    except Exception as e:
        logging.error(f"Data processing failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
