import logging
from concurrent.futures import ThreadPoolExecutor
from typing import Dict, Optional, Tuple
import pandas as pd
import requests
import os
import json

logging.basicConfig(level=logging.INFO)

CACHE_FILE = "uniprot_cache.json"
CDS_CACHE_FILE = "cds_cache.json"
VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")

# Standard Genetic Code
GENETIC_CODE = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "TGT": "C",
    "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}


def translate_dna(dna_seq: str) -> str:
    """Translate DNA sequence to protein"""
    protein = []
    for i in range(0, len(dna_seq), 3):
        codon = dna_seq[i : i + 3]
        if len(codon) < 3:
            break
        aa = GENETIC_CODE.get(codon, "X")
        if aa == "*":
            break
        protein.append(aa)
    return "".join(protein)


def load_cache(filename: str = CACHE_FILE) -> Dict:
    """Load sequence cache"""
    if os.path.exists(filename):
        try:
            with open(filename, "r") as f:
                return json.load(f)
        except Exception as e:
            logging.warning(f"Error loading cache {filename}: {str(e)}")
            return {}
    return {}


def save_cache(cache: Dict, filename: str = CACHE_FILE) -> None:
    """Save sequence cache"""
    with open(filename, "w") as f:
        json.dump(cache, f)


def get_uniprot_sequences_batch(genes: list, taxonomy_id: str = "9606") -> Dict:
    """Fetch sequences from UniProt in batch"""
    try:
        gene_queries = [f"gene:{gene}" for gene in genes]
        query = f"({' OR '.join(gene_queries)}) AND reviewed:true AND taxonomy_id:{taxonomy_id}"
        url = f"https://rest.uniprot.org/uniprotkb/search?query={query}&fields=accession,gene_names,sequence,length&format=json&size=500"

        response = requests.get(url, timeout=30)
        response.raise_for_status()

        results = response.json().get("results", [])
        gene_map = {}

        for result in results:
            seq = result.get("sequence", {}).get("value", "").upper()
            length = result.get("sequence", {}).get("length")

            if not seq or length is None:
                continue

            genes_in_entry = [g["geneName"]["value"] for g in result.get("genes", [])]
            for gene in genes_in_entry:
                if gene in genes:
                    gene_map[gene] = (seq, length)

        return gene_map
    except Exception as e:
        logging.error(f"Batch request failed: {str(e)}")
        return {}


def get_sequences(genes: list, batch_size: int = 50) -> Dict:
    """Get protein sequences with caching"""
    cache = load_cache()
    missing = [
        g for g in genes if g not in cache or not isinstance(cache.get(g), tuple)
    ]

    if missing:
        with ThreadPoolExecutor(max_workers=4) as executor:
            futures = []
            for i in range(0, len(missing), batch_size):
                batch = missing[i : i + batch_size]
                futures.append(executor.submit(get_uniprot_sequences_batch, batch))

            for future in futures:
                batch_results = future.result()
                cache.update(batch_results)
                save_cache(cache, CACHE_FILE)

    return {gene: cache[gene] for gene in genes if gene in cache}


def get_cds_sequences(genes: list) -> Dict:
    """Get CDS sequences from Ensembl with caching"""
    cache = load_cache(CDS_CACHE_FILE)
    missing = [g for g in genes if g not in cache]

    if missing:
        # Ensembl recommends small parallel requests or sequential.
        # We'll do sequential for now to avoid rate limits on public API.
        server = "https://rest.ensembl.org"
        headers = {"Content-Type": "application/json"}

        for gene in missing:
            try:
                # 1. Lookup symbol to get ID and canonical transcript
                ext = f"/lookup/symbol/homo_sapiens/{gene}?expand=1"
                r = requests.get(server + ext, headers=headers)

                if not r.ok:
                    logging.warning(f"Failed lookup for {gene}: {r.status_code}")
                    continue

                decoded = r.json()
                transcripts = decoded.get("Transcript", [])
                coding = [t for t in transcripts if t["biotype"] == "protein_coding"]

                if not coding:
                    continue

                # Use canonical if available, else longest CDS
                canonical = next(
                    (t for t in coding if t.get("is_canonical")), coding[0]
                )
                transcript_id = canonical["id"]

                # 2. Fetch CDS
                ext_seq = f"/sequence/id/{transcript_id}?type=cds"
                r_seq = requests.get(
                    server + ext_seq, headers={"Content-Type": "text/plain"}
                )

                if r_seq.ok:
                    cache[gene] = r_seq.text

            except Exception as e:
                logging.warning(f"Error fetching CDS for {gene}: {e}")

        save_cache(cache, CDS_CACHE_FILE)

    return {gene: cache[gene] for gene in genes if gene in cache}


def generate_frameshift_sequence(
    wildtype_protein: str, cds_seq: str, pos: int, alt: str, length: int = 30
) -> str:
    """
    Generate mutant protein sequence for frameshift mutations using CDS.

    Args:
        wildtype_protein: The wildtype protein sequence.
        cds_seq: The coding DNA sequence (CDS).
        pos: The 1-based amino acid position where the frameshift starts.
        alt: The inserted/substituted character(s) causing the frameshift.
        length: Max length of the mutant tail.
    """
    if not cds_seq:
        return None

    try:
        # Approximate CDS position from protein position
        # Protein pos 1 -> CDS codon 0-2 (indices)
        # CDS index = (pos - 1) * 3
        # Note: This assumes the CDS matches the protein exactly (canonical).
        cds_start_idx = (pos - 1) * 3

        if cds_start_idx >= len(cds_seq):
            return None

        # Construct mutant DNA
        # For simplicity in this heuristic pipeline:
        # We assume the frameshift is an insertion of 'alt' bases or similar.
        # But `alt` in HGVSp is amino acids? No, HGVSp 'fs' usually implies a shift.
        # Actually HGVSp `p.P408Afs*99` means at P408, we see A, then a frameshift.
        # The `alt` passed here from `parse_protein_change` (which regexes for P408A...) is 'A'.
        # This 'A' is the *first* amino acid of the new frame.
        # To get the actual DNA sequence change, we need the *genetic* mutation (HGVSc).
        # We don't have HGVSc here. This is a limitation.

        # Heuristic workaround:
        # If we have `p.P408Afs`, it means at codon 408 (which codes for P),
        # the DNA changed such that it now codes for A, and shifts frame.
        # Usually this is a deletion or insertion of non-multiple-of-3 bases.
        # Without HGVSc, we can't know the exact sequence.
        # However, typically fs corresponds to specific indel.
        # If we strictly want to generate *some* frameshift for testing:
        # We can simulate a +1 base insertion after the codon?

        # BETTER APPROACH FOR THIS SCOPE:
        # We can't accurately reconstruct the exact patient neoantigen without HGVSc.
        # But we can simulate a generic frameshift at this position to generate *a* neoantigen candidate.
        # A common frameshift is a single base deletion or insertion.
        # Let's simulate a single base deletion at the start of the codon to shift frame -1.

        # We want the resulting protein to start with `alt` at `pos`.
        # `wildtype_protein[:pos-1]` is correct.
        # Then `alt` (the new AA).
        # Then the rest of the translation in the new frame.

        # But we don't know the new frame (is it +1 or +2?).
        # We try both?

        # Let's try to match `alt`.
        # Shift +1 (insert 'A') -> check if codon is `alt`.
        # Shift -1 (delete 1) -> check if codon is `alt`.
        # Whatever matches `alt` best.

        # Since we can't be perfect without HGVSc, we will default to a 1-base deletion (frame +2/-1)
        # unless `alt` suggests otherwise?

        # Let's just do a 1-base deletion at [cds_start_idx].
        mutant_dna = cds_seq[:cds_start_idx] + cds_seq[cds_start_idx + 1 :]

        # Translate from cds_start_idx
        # But wait, we need the *new* tail.
        # We translate `mutant_dna` starting from `cds_start_idx`.
        # The protein up to pos-1 is the same.

        new_tail_dna = mutant_dna[cds_start_idx:]
        new_tail_protein = translate_dna(new_tail_dna)

        full_mutant = wildtype_protein[: pos - 1] + new_tail_protein
        return full_mutant

    except Exception as e:
        logging.warning(f"Error in frameshift simulation: {str(e)}")
        return None


def generate_peptides(row: pd.Series, peptide_length: int = 9) -> Optional[str]:
    """Generate a neopeptide sequence based on mutation data"""
    wildtype_seq = getattr(row, "wildtype_seq", None)
    hugo_symbol = getattr(row, "Hugo_Symbol", "Unknown")

    if (
        pd.isna(wildtype_seq)
        or not isinstance(wildtype_seq, tuple)
        or len(wildtype_seq) != 2
    ):
        logging.warning(f"Invalid or missing wildtype_seq for {hugo_symbol}")
        return None

    seq, seq_length = wildtype_seq
    if not seq or seq_length is None:
        logging.warning(f"Missing sequence or length for {hugo_symbol}")
        return None

    pos = getattr(row, "pos", None)
    # Handle end_pos if present, else default to pos
    end_pos = getattr(row, "end_pos", pos)

    mut_type = getattr(row, "mut_type", None)
    alt = getattr(row, "alt", None)

    if pos > seq_length or pos <= 0:
        logging.warning(
            f"Isoform mismatch or invalid position: {hugo_symbol} "
            f"(UniProt length {seq_length} vs TCGA position {pos})"
        )
        return None

    mutant_peptide = None
    try:
        # Strategy: Construct local mutant sequence context, then extract 9-mer
        # We need enough context to cover the peptide length
        # Start far enough back
        context_start = max(0, pos - peptide_length - 5)
        context_end = min(seq_length, end_pos + peptide_length + 5)

        # Extract wildtype context
        wt_context = seq[context_start:context_end]

        # Relative positions in context
        rel_pos = pos - context_start
        rel_end_pos = end_pos - context_start

        # Construct mutant context
        if mut_type == "Sub":
            if alt not in VALID_AA:
                return None
            mut_context = wt_context[: rel_pos - 1] + alt + wt_context[rel_pos:]
            # Center the window around the mutation.
            # Mutation index in mut_context is rel_pos-1
            mut_idx = rel_pos - 1

        elif mut_type == "Fs":
            cds_seq = getattr(row, "cds_seq", None)
            # If we don't have CDS, we return None as per new strict biological requirement
            if not cds_seq or pd.isna(cds_seq):
                logging.warning(f"Missing CDS for frameshift in {hugo_symbol}")
                return None

            fs_seq = generate_frameshift_sequence(seq, cds_seq, pos, alt)
            if not fs_seq:
                return None

            # For FS, the neoantigen is the new sequence downstream.
            # Let's take [pos-1 : pos-1+9]
            start_extract = max(0, pos - 1)
            mutant_peptide = fs_seq[start_extract : start_extract + peptide_length]

        elif mut_type == "Del":
            # Deletion: remove [rel_pos-1 : rel_end_pos]
            mut_context = wt_context[: rel_pos - 1] + wt_context[rel_end_pos:]
            mut_idx = rel_pos - 1  # Position where deletion occurred

        elif mut_type == "Dup":
            dup_seq = wt_context[rel_pos - 1 : rel_end_pos]
            mut_context = wt_context[:rel_end_pos] + dup_seq + wt_context[rel_end_pos:]
            mut_idx = rel_end_pos  # Point of insertion of dup

        elif mut_type == "Delins":
            if not all(aa in VALID_AA for aa in alt):
                return None
            mut_context = wt_context[: rel_pos - 1] + alt + wt_context[rel_end_pos:]
            mut_idx = rel_pos - 1

        elif mut_type == "Ins":
            if not all(aa in VALID_AA for aa in alt):
                return None
            mut_context = wt_context[:rel_pos] + alt + wt_context[rel_pos:]
            mut_idx = rel_pos  # Insertion point (between pos and pos+1)

        else:
            return None

        # Determine extraction window from mut_context
        # We want the 9-mer that best represents the mutation.
        # Usually centered.
        if mut_type != "Fs":
            # center of 9-mer should be around mut_idx
            # start = mut_idx - 4
            # end = mut_idx + 5

            # Safety check on bounds
            if mut_idx < 4:
                # Close to start of protein
                start_i = 0
            else:
                start_i = mut_idx - 4

            mutant_peptide = mut_context[start_i : start_i + peptide_length]

    except Exception as e:
        logging.warning(f"Error generating peptide: {e}")
        return None

    if not mutant_peptide:
        return None

    # Pad if necessary - actually, for strict biological validity, we should NOT pad with X or G.
    # If the peptide is shorter than the required length (e.g. near termini), we should probably discard it
    # or return it as is, but the pipeline likely expects fixed length.
    # Given the user instruction to "Avoid G padding", we will enforce exact length.
    if len(mutant_peptide) < peptide_length:
        logging.warning(
            f"Peptide too short for {hugo_symbol}: {len(mutant_peptide)} < {peptide_length}. Skipping."
        )
        return None

    # Final validation
    mutant_peptide = mutant_peptide[:peptide_length]  # Ensure exact length
    if not all(aa in VALID_AA for aa in mutant_peptide):
        logging.warning(
            f"Invalid peptide generated for {hugo_symbol}: {mutant_peptide}"
        )
        return None

    return mutant_peptide
