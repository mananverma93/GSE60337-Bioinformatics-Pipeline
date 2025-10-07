import pandas as pd
import os

# Check if running via Snakemake or standalone
try:
    input_de_results = snakemake.input.de_results
    input_annotation = snakemake.input.annotation
    output_annotated = snakemake.output[0]
except NameError:
    print("Running in standalone mode")
    input_de_results = "data/processed/DE_full_C57BL_6J_vs_CBA_CaJ.csv"
    input_annotation = "data/raw/GPL6246-18741.txt"
    output_annotated = "data/processed/DE_annotated_C57BL_6J_vs_CBA_CaJ.csv"

# Ensure output directory exists
os.makedirs(os.path.dirname(output_annotated), exist_ok=True)

print("Loading DE results...")
de_results = pd.read_csv(input_de_results)
print(f"DE results shape: {de_results.shape}")

print("\nLoading annotation file...")
# Read the file, skipping comment lines that start with #
annotation = pd.read_csv(input_annotation, sep='\t', comment='#', low_memory=False)
print(f"Annotation shape: {annotation.shape}")
print(f"Annotation columns: {annotation.columns.tolist()}")

# Keep only necessary columns from annotation
# ID = probe ID, gene_assignment = gene info
if 'ID' in annotation.columns and 'gene_assignment' in annotation.columns:
    annotation_subset = annotation[['ID', 'gene_assignment']].copy()
    
    # Extract gene symbol from gene_assignment
    # Format: "NM_008866 // Lypla1 // lysophospholipase 1 // ..."
    def extract_gene_symbol(gene_assignment):
        if pd.isna(gene_assignment) or gene_assignment == '---':
            return 'Unknown'
        parts = str(gene_assignment).split('//')
        if len(parts) >= 2:
            return parts[1].strip()
        return 'Unknown'
    
    annotation_subset['gene_symbol'] = annotation_subset['gene_assignment'].apply(extract_gene_symbol)
    annotation_subset = annotation_subset[['ID', 'gene_symbol']]
    annotation_subset.columns = ['gene_id', 'gene_symbol']
    
    print(f"\nAnnotation subset shape: {annotation_subset.shape}")
    print(f"Sample annotations:")
    print(annotation_subset.head(10))
else:
    print("ERROR: Expected columns not found in annotation file")
    print(f"Available columns: {annotation.columns.tolist()}")
    exit(1)

# Merge DE results with annotations
print("\nMerging DE results with gene annotations...")
de_annotated = de_results.merge(annotation_subset, on='gene_id', how='left')

# Fill missing gene symbols
de_annotated['gene_symbol'] = de_annotated['gene_symbol'].fillna('Unknown')

# Reorder columns to put gene_symbol first
cols = ['gene_id', 'gene_symbol'] + [col for col in de_annotated.columns if col not in ['gene_id', 'gene_symbol']]
de_annotated = de_annotated[cols]

# Save annotated results
de_annotated.to_csv(output_annotated, index=False)

print(f"\n{'='*60}")
print("ANNOTATION COMPLETE")
print(f"{'='*60}")
print(f"Total genes: {len(de_annotated)}")
print(f"Annotated genes: {(de_annotated['gene_symbol'] != 'Unknown').sum()}")
print(f"Unknown genes: {(de_annotated['gene_symbol'] == 'Unknown').sum()}")

if 'significant' in de_annotated.columns:
    sig_annotated = de_annotated[de_annotated['significant']]
    print(f"\nSignificant genes: {len(sig_annotated)}")
    print(f"\nTop 10 significant genes with annotations:")
    print(sig_annotated[['gene_id', 'gene_symbol', 'log2FC', 'p_adjusted']].head(10))

print(f"\nAnnotated results saved to: {output_annotated}")