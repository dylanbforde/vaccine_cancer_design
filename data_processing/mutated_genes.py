import pandas as pd
import re
import logging
from typing import Optional, Tuple, Dict
from data_processing.helper_functions import (
    generate_peptides,
    get_sequences,
    get_cds_sequences,
)

logging.basicConfig(level=logging.INFO)


def parse_protein_change(hgvsp_short: str) -> Optional[Tuple[int, int, str, str]]:
    """
    Parse HGVSp notation into (start_pos, end_pos, mutation_type, alt_aa)
    Returns None if the format is invalid
    """
    if pd.isna(hgvsp_short) or not hgvsp_short.startswith("p."):
        return None

    # Handle stop codon (e.g., p.R30*)
    stop_match = re.match(r"^p\.([A-Z])(\d+)\*$", hgvsp_short)
    if stop_match:
        pos = int(stop_match.group(2))
        return (pos, pos, "Stop", "*")

    patterns = [
        # Simple substitutions (e.g., p.E678K)
        (r"^p\.([A-Z])(\d+)([A-Z])$", "Sub"),
        # Frameshifts (e.g., p.P408Afs*99)
        (r"^p\.([A-Z])(\d+)([A-Z])fs\*\d+$", "Fs"),
        # Deletions (e.g., p.Q344del, p.K2_L3del)
        (r"^p\.([A-Z])(\d+)(_([A-Z])(\d+))?del$", "Del"),
        # Duplications (e.g., p.S2107dup, p.K2_L3dup)
        (r"^p\.([A-Z])(\d+)(_([A-Z])(\d+))?dup$", "Dup"),
        # Insertions (e.g., p.K2_L3insQ)
        (r"^p\.([A-Z])(\d+)_([A-Z])(\d+)ins([A-Z]+)$", "Ins"),
        # Complex changes (e.g., p.H57_S58delinsQ)
        (r"^p\.([A-Z])(\d+)(_([A-Z])(\d+))?delins([A-Z]+)$", "Delins"),
    ]

    for pattern, mut_type in patterns:
        match = re.match(pattern, hgvsp_short)
        if match:
            start_pos = int(match.group(2))

            # Extract end position
            if mut_type in ["Del", "Dup", "Delins"]:
                # Group 4/5 are typically the end of range if present
                # Standard regex groups: 1=StartAA, 2=StartPos, 3=_EndAA EndPos part, 4=EndAA, 5=EndPos
                # (For Delins, group 6 is Alt)
                if match.group(5):
                    end_pos = int(match.group(5))
                else:
                    end_pos = start_pos
            elif mut_type == "Ins":
                # Group: 1=StartAA, 2=StartPos, 3=EndAA, 4=EndPos, 5=Alt
                end_pos = int(match.group(4))
            else:
                end_pos = start_pos

            # Extract Alt
            if mut_type == "Delins":
                alt = match.group(6)
            elif mut_type == "Ins":
                alt = match.group(5)
            elif mut_type == "Fs":
                # Frameshift alt is the first new AA
                alt = match.group(3)
            elif mut_type == "Sub":
                alt = match.group(3)
            else:
                alt = "X"  # Del/Dup usually don't have 'alt' sequence in HGVSp (implied from WT)

            return (start_pos, end_pos, mut_type, alt)

    return None


def filter_variants(df: pd.DataFrame) -> pd.DataFrame:
    """Filter for processable mutation types"""
    valid_variants = [
        "Missense_Mutation",
        "Frame_Shift_Ins",
        "Frame_Shift_Del",
        "In_Frame_Del",
        "In_Frame_Ins",
    ]

    filtered = df[
        df["Variant_Classification"].isin(valid_variants) & df["HGVSp_Short"].notna()
    ][["Hugo_Symbol", "HGVSp_Short", "Tumor_Sample_Barcode"]]

    logging.info(f"Filtered {len(filtered)} variants from {len(df)} total mutations")
    return filtered


def process_mutations_in_batches(
    mutations_df: pd.DataFrame, batch_size: int = 1000
) -> pd.DataFrame:
    """Process mutations in batches to handle large datasets"""
    all_results = []
    total_batches = len(mutations_df) // batch_size + (
        1 if len(mutations_df) % batch_size > 0 else 0
    )

    for batch_num in range(total_batches):
        start_idx = batch_num * batch_size
        end_idx = min((batch_num + 1) * batch_size, len(mutations_df))
        batch_df = mutations_df.iloc[start_idx:end_idx]

        logging.info(
            f"Processing batch {batch_num + 1}/{total_batches} ({len(batch_df)} mutations)"
        )

        # Parse protein changes
        # Vectorized parsing for performance
        parsed_series = batch_df["HGVSp_Short"].apply(parse_protein_change)
        valid_mask = parsed_series.notna()

        if valid_mask.any():
            # Create parsed_df from valid results
            # The parsed_series contains tuples: (start_pos, end_pos, mut_type, alt)
            parsed_data = pd.DataFrame(
                parsed_series[valid_mask].tolist(),
                columns=["pos", "end_pos", "mut_type", "alt"],
                index=batch_df.loc[valid_mask].index,
            )

            # Combine with original data
            parsed_df = pd.concat(
                [
                    batch_df.loc[
                        valid_mask, ["Hugo_Symbol", "HGVSp_Short", "Tumor_Sample_Barcode"]
                    ],
                    parsed_data,
                ],
                axis=1,
            )

            # Identify genes needing CDS
            genes_needing_cds = set(
                parsed_df[parsed_df["mut_type"] == "Fs"]["Hugo_Symbol"]
            )

            # Fetch protein sequences for batch
            sequences = get_sequences(parsed_df["Hugo_Symbol"].unique())
            parsed_df["wildtype_seq"] = parsed_df["Hugo_Symbol"].map(sequences)

            # Fetch CDS sequences for frameshifts
            cds_sequences = {}
            if genes_needing_cds:
                logging.info(
                    f"Fetching CDS for {len(genes_needing_cds)} genes with frameshifts"
                )
                cds_sequences = get_cds_sequences(list(genes_needing_cds))

            parsed_df["cds_seq"] = parsed_df["Hugo_Symbol"].map(cds_sequences)

            # Generate peptides
            parsed_df["peptide"] = parsed_df.apply(generate_peptides, axis=1)

            # Keep only valid results
            valid_results = parsed_df.dropna(subset=["peptide"])
            all_results.append(valid_results)

            logging.info(
                f"Batch {batch_num + 1}: Generated {len(valid_results)} valid peptides"
            )

    if all_results:
        final_df = pd.concat(all_results, ignore_index=True)
        logging.info("\nMutation type distribution:")
        print(final_df["mut_type"].value_counts())
        return final_df
    return pd.DataFrame()


def main():
    """Main pipeline for processing mutations into vaccine candidates"""
    logging.info("Loading mutation data...")
    mutations = pd.read_csv(
        "data/brca_tcga_pan_can_atlas_2018/data_mutations.txt",
        sep="\t",
        comment="#",
        low_memory=False,
    )

    logging.info(f"Total mutations loaded: {len(mutations)}")

    # Filter variants
    filtered = filter_variants(mutations)
    logging.info(f"Filtered to {len(filtered)} relevant mutations")

    # Process in batches
    results = process_mutations_in_batches(filtered)

    if len(results) > 0:
        results.to_csv("data/mutations_processed.csv", index=False)
        logging.info(f"Successfully processed {len(results)} mutations")
        logging.info("Results saved to mutations_processed.csv")

        return results
    else:
        logging.warning("No valid mutations processed")
        return pd.DataFrame()


if __name__ == "__main__":
    main()
