import pandas as pd
import re

def analyze_protein_changes(df: pd.DataFrame):
    """Analyzes the protein changes in the given DataFrame."""

    print("Analyzing protein changes...\n")

    # 1. Frequency of Variant Classifications
    print("Variant Classification Frequencies:")
    print(df['Variant_Classification'].value_counts())
    print("\n" + "-" * 40 + "\n")

    # 2. Unique HGVSp_Short Patterns
    print("Unique HGVSp_Short patterns (first 50):")
    unique_patterns = df['HGVSp_Short'].dropna().unique()
    for pattern in unique_patterns[:50]:
        print(pattern)
    print("\n" + "-" * 40 + "\n")

    # 3. Examples of HGVSp_Short for each Variant Classification
    print("Examples of HGVSp_Short for each Variant Classification:")
    for var_class in df['Variant_Classification'].unique():
        print(f"\n{var_class}:")
        examples = df[df['Variant_Classification'] == var_class]['HGVSp_Short'].dropna().head(10)
        for example in examples:
            print(example)

    # 4. Basic Statistics on the Length of HGVSp_Short
    df['HGVSp_Short_Length'] = df['HGVSp_Short'].str.len()
    print("\n" + "-" * 40 + "\n")
    print("Basic Statistics on the Length of HGVSp_Short:")
    print(df['HGVSp_Short_Length'].describe())

def main():
    """Main function to load data and perform analysis."""

    # Load data
    try:
        raw_mutations = pd.read_csv('./data/brca_tcga_pan_can_atlas_2018/data_mutations.txt', sep='\t', comment='#', low_memory=False)
    except FileNotFoundError:
        print("Error: data_mutations.txt not found. Please check the file path.")
        return

    # Analyze protein changes
    analyze_protein_changes(raw_mutations)

if __name__ == "__main__":
    main()