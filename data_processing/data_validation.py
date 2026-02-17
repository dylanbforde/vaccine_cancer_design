import pandas as pd
import logging
import re
from typing import Tuple, Dict
from mutated_genes import parse_protein_change

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')

VALID_AA = set('ACDEFGHIKLMNPQRSTVWY')

def validate_mutations(df: pd.DataFrame) -> Tuple[bool, Dict]:
    """Validate mutation data using the same parser as the main pipeline"""
    errors = {}
    
    # Use the same parser as the main pipeline for consistency
    valid_variants = df['HGVSp_Short'].apply(lambda x: parse_protein_change(x) is not None)
    
    invalid_count = (~valid_variants).sum()
    valid_count = valid_variants.sum()
    
    if invalid_count > 0:
        errors['invalid_protein_changes'] = invalid_count
    errors['valid_protein_changes'] = valid_count
    
    logging.info(f"Valid mutations: {valid_count}")
    logging.info(f"Invalid mutations: {invalid_count}")
    
    # Consider validation successful if we have valid mutations
    return valid_count > 0, errors

def validate_peptides(peptides: pd.Series) -> Tuple[bool, Dict]:
    """
    Validates generated peptides:
    - Contains only valid amino acids
    - Has correct length (9-mer)
    """
    errors = {}
    
    if peptides.empty:
        errors['empty_peptides'] = 0
        return False, errors
    
    # Check for valid amino acids
    invalid_aa = peptides.apply(
        lambda x: False if pd.isna(x) else not all(aa in VALID_AA for aa in x)
    )
    invalid_aa_count = invalid_aa.sum()
    
    if invalid_aa_count > 0:
        errors['invalid_amino_acids'] = invalid_aa_count
    
    # Check length (should be 9-mer)
    invalid_length = peptides.apply(
        lambda x: False if pd.isna(x) else len(x) != 9
    )
    invalid_length_count = invalid_length.sum()
    
    if invalid_length_count > 0:
        errors['invalid_length'] = invalid_length_count
    
    return len(errors) == 0, errors

def cross_validate_mutations(raw_df: pd.DataFrame, processed_df: pd.DataFrame) -> Dict[str, int]:
    """Cross validate raw mutations against processed results"""
    stats = {
        'total_raw_mutations': len(raw_df),
        'total_processed': len(processed_df),
        'valid_mutations': 0,
        'valid_peptides': 0,
        'processing_efficiency': 0.0
    }
    
    # Count valid mutations (those we can parse)
    valid_mutations = raw_df['HGVSp_Short'].apply(lambda x: parse_protein_change(x) is not None)
    stats['valid_mutations'] = valid_mutations.sum()
    
    # Count valid peptides
    if 'peptide' in processed_df.columns:
        stats['valid_peptides'] = processed_df['peptide'].notna().sum()
    
    # Calculate efficiency
    if stats['valid_mutations'] > 0:
        stats['processing_efficiency'] = (stats['valid_peptides'] / stats['valid_mutations']) * 100
    
    logging.info("Cross-validation stats:")
    for key, value in stats.items():
        if key == 'processing_efficiency':
            logging.info(f"{key}: {value:.1f}%")
        else:
            logging.info(f"{key}: {value}")
    
    return stats

def generate_validation_report(validation_results: Dict[str, Tuple[bool, Dict]]) -> str:
    """Generate a comprehensive validation report"""
    report = ["Data Validation Report", "========================"]
    
    for stage, (is_valid, errors) in validation_results.items():
        report.append(f"\nStage: {stage.upper()}")
        report.append(f"Status: {'VALID' if is_valid else 'INVALID'}")
        
        if errors:
            report.append("Statistics:")
            for error, count in errors.items():
                if isinstance(count, float):
                    report.append(f"  - {error}: {count:.1f}%")
                else:
                    report.append(f"  - {error}: {count:,} instances")
        else:
            report.append("No errors detected")
    
    return "\n".join(report)

def main():
    """Run full validation pipeline"""
    validation_results = {}

    # Load data
    print("Loading raw mutations...")
    raw_mutations = pd.read_csv(
        'data/brca_tcga_pan_can_atlas_2018/data_mutations.txt', 
        sep='\t', 
        comment='#',
        low_memory=False
    )
    
    # Validate raw mutations
    is_valid, errors = validate_mutations(raw_mutations)
    validation_results['raw_mutations'] = (is_valid, errors)

    # Validate processed peptides
    try:
        print("Loading processed mutations...")
        processed = pd.read_csv('data/mutations_processed.csv')
        is_valid, errors = validate_peptides(processed['peptide'])
        validation_results['processed_peptides'] = (is_valid, errors)
        
        # Add cross-validation
        print("\nPerforming cross-validation...")
        cross_stats = cross_validate_mutations(raw_mutations, processed)
        validation_results['cross_validation'] = (True, cross_stats)
        
    except FileNotFoundError:
        validation_results['processed_peptides'] = (False, {'file_not_found': 1})

    # Generate and display report
    report = generate_validation_report(validation_results)
    print("\n" + report)

    # Save report if successful
    if all(result[0] for result in validation_results.values()):
        with open('docs/validation_report.txt', 'w') as f:
            f.write(report)
        print("\nValidation report saved successfully.")
    else:
        print("\nValidation errors detected. Report not saved.")

if __name__ == "__main__":
    main()