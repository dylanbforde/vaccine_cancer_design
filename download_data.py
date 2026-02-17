import os
import requests
import zipfile
import logging
import sys
from tqdm import tqdm

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)

DATA_DIR = "data"
os.makedirs(DATA_DIR, exist_ok=True)

TRACK_URLS = {
    "TCGA_BRCA": {
        "url": "https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/brca_tcga_pan_can_atlas_2018/data_mutations.txt",
        "filename": "brca_tcga_pan_can_atlas_2018/data_mutations.txt",  # Nested path
        "extract": False,
        "type": "txt",
    },
    "IEDB_TCELL": {
        "url": "https://www.iedb.org/downloader.php?file_name=doc/tcell_full_v3.zip",
        "filename": "tcell_full_v3.zip",
        "extract": True,
        "type": "zip",
    },
}


def download_file(url, filepath):
    """Download a file with progress bar and User-Agent"""
    session = requests.Session()
    session.headers.update(
        {
            "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36"
        }
    )

    response = session.get(url, stream=True)
    response.raise_for_status()  # Check for errors

    total_size = int(response.headers.get("content-length", 0))
    block_size = 1024  # 1 Kibibyte

    # Ensure directory exists
    os.makedirs(os.path.dirname(filepath), exist_ok=True)

    logging.info(f"Downloading {os.path.basename(filepath)}...")
    with (
        open(filepath, "wb") as file,
        tqdm(
            desc=os.path.basename(filepath),
            total=total_size,
            unit="iB",
            unit_scale=True,
            unit_divisor=1024,
        ) as bar,
    ):
        for data in response.iter_content(block_size):
            size = file.write(data)
            bar.update(size)


def extract_file(filepath, extract_to, archive_type):
    """Extract archive file"""
    logging.info(f"Extracting {os.path.basename(filepath)}...")
    try:
        if archive_type == "zip":
            with zipfile.ZipFile(filepath, "r") as zip_ref:
                zip_ref.extractall(path=extract_to)
        logging.info("Extraction complete.")
    except Exception as e:
        logging.error(f"Failed to extract {filepath}: {e}")


def main():
    for name, info in TRACK_URLS.items():
        filepath = os.path.join(DATA_DIR, info["filename"])

        # Download
        if not os.path.exists(filepath):
            try:
                download_file(info["url"], filepath)
            except Exception as e:
                logging.error(f"Failed to download {name}: {e}")
                continue
        else:
            logging.info(f"{info['filename']} already exists. Skipping download.")

        # Extract
        if info["extract"]:
            if name == "IEDB_TCELL":
                # Check for target file to avoid re-extraction
                expected_file = os.path.join(DATA_DIR, "tcell_full_v3.csv")
                if not os.path.exists(expected_file):
                    extract_file(filepath, DATA_DIR, info["type"])
                else:
                    logging.info(f"IEDB data already extracted to {expected_file}")

    logging.info("\nData download complete! You can now run 'uv run main.py'")


if __name__ == "__main__":
    main()
