import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import os
from collections import Counter
from statsmodels.stats.multitest import multipletests

# Check if running via Snakemake or standalone
try:
    input_data = snakemake.input.data
    input_metadata = snakemake.input.metadata
    output_results = snakemake.output.results
    output_volcano = snakemake.output.volcano
    group1_name = snakemake.params.group1
    group2_name = snakemake.params.group2
except NameError:
    print("Running in standalone mode")
    input_data = "data/processed/GSE60337_full_normalized.csv"
    input_metadata = "data/raw/GSE60337_series_matrix.txt"
    output_results = "data/processed/DE_full_C57BL_6J_vs_CBA_CaJ.csv"
    output_volcano = "data/processed/volcano_plot_full_C57BL_6J_vs_CBA_CaJ.png"
    group1_name = "C57BL/6J"
    group2_name = "CBA/CaJ"

# Ensure output directory exists
os.makedirs(os.path.dirname(output_results), exist_ok=True)

# Parse metadata
sample_to_strain = {}
try:
    with open(input_metadata, 'r', encoding='utf-8') as f:
        lines = f.readlines()
except UnicodeDecodeError:
    with open(input_metadata, 'r', encoding='latin-1') as f:
        lines = f.readlines()

for line in lines:
    if line.startswith('!Sample_geo_accession'):
        samples = [s.strip('"') for s in line.strip().split('\t')[1:]]
    if line.startswith('!Sample_characteristics_ch1') and 'genetic background:' in line:
        strains = [s.strip('"').replace('genetic background: ', '') for s in line.strip().split('\t')[1:]]
        sample_to_strain = dict(zip(samples, strains))
        break

# Load normalized data
df_norm = pd.read_csv(input_data, index_col=0)

# Define comparison groups
group1_samples = [s for s, strain in sample_to_strain.items() if strain == group1_name and s in df_norm.columns]
group2_samples = [s for s, strain in sample_to_strain.items() if strain == group2_name and s in df_norm.columns]

print(f"{'='*60}")
print(f"DIFFERENTIAL EXPRESSION ANALYSIS")
print(f"{'='*60}")
print(f"Comparison: {group1_name} vs {group2_name}")
print(f"Group 1: {len(group1_samples)} samples")
print(f"Group 2: {len(group2_samples)} samples")

# Extract data
group1_data = df_norm[group1_samples]
group2_data = df_norm[group2_samples]

# Calculate statistics
gene_names = []
fold_changes = []
log2_fold_changes = []
p_values = []
group1_means = []
group2_means = []

print("\nCalculating statistics...")
for gene in df_norm.index:
    g1_values = group1_data.loc[gene].values
    g2_values = group2_data.loc[gene].values
    
    g1_mean = np.mean(g1_values)
    g2_mean = np.mean(g2_values)
    
    fold_change = g1_mean / g2_mean
    log2_fc = np.log2(fold_change)
    
    t_stat, p_val = stats.ttest_ind(g1_values, g2_values)
    
    gene_names.append(gene)
    group1_means.append(g1_mean)
    group2_means.append(g2_mean)
    fold_changes.append(fold_change)
    log2_fold_changes.append(log2_fc)
    p_values.append(p_val)

# Create results dataframe
de_results = pd.DataFrame({
    'gene_id': gene_names,
    f'mean_{group1_name}': group1_means,
    f'mean_{group2_name}': group2_means,
    'fold_change': fold_changes,
    'log2FC': log2_fold_changes,
    'p_value': p_values
})

de_results = de_results.sort_values('p_value')

# Clean and apply FDR correction
de_results_clean = de_results.dropna(subset=['p_value'])
de_results_clean = de_results_clean[np.isfinite(de_results_clean['p_value'])]

reject, p_adjusted, _, _ = multipletests(de_results_clean['p_value'], 
                                         alpha=0.05, method='fdr_bh')

de_results_clean['p_adjusted'] = p_adjusted
de_results_clean['significant'] = reject

sig_genes = de_results_clean[de_results_clean['significant']]

print(f"\n{'='*60}")
print(f"RESULTS")
print(f"{'='*60}")
print(f"Total genes tested: {len(de_results_clean)}")
print(f"Significant genes (FDR < 0.05): {len(sig_genes)}")

# Save results
de_results_clean.to_csv(output_results, index=False)

# Volcano plot
plt.figure(figsize=(10, 7))

plt.scatter(de_results_clean['log2FC'], 
           -np.log10(de_results_clean['p_value'].replace(0, 1e-300)),
           alpha=0.4, s=10, c='gray', label='Not significant')

if len(sig_genes) > 0:
    plt.scatter(sig_genes['log2FC'], 
               -np.log10(sig_genes['p_value'].replace(0, 1e-300)),
               alpha=0.8, s=30, c='red', label=f'Significant (n={len(sig_genes)})')

plt.axhline(y=-np.log10(0.05), color='blue', linestyle='--', alpha=0.5, label='p=0.05')
plt.axvline(x=-1, color='green', linestyle='--', alpha=0.5)
plt.axvline(x=1, color='green', linestyle='--', alpha=0.5, label='|Log2FC| = 1')

plt.xlabel(f'Log2 Fold Change ({group1_name} vs {group2_name})')
plt.ylabel('-Log10(p-value)')
plt.title(f'Volcano Plot: {len(group1_samples)} vs {len(group2_samples)} samples')
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(output_volcano, dpi=300)
plt.close()

print(f"Results saved to {output_results}")
print(f"Volcano plot saved to {output_volcano}")
print("Analysis complete!")