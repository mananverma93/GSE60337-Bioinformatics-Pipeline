import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from matplotlib.patches import Patch
import os

# Check if running via Snakemake or standalone
try:
    input_data = snakemake.input.data
    input_metadata = snakemake.input.metadata
    output_pca = snakemake.output.pca
    output_heatmap = snakemake.output.heatmap
except NameError:
    print("Running in standalone mode")
    input_data = "data/processed/GSE60337_full_normalized.csv"
    input_metadata = "data/raw/GSE60337_series_matrix.txt"
    output_pca = "data/processed/pca_by_strain_full.png"
    output_heatmap = "data/processed/top_variable_genes_heatmap_full.png"

# Ensure output directory exists
os.makedirs(os.path.dirname(output_pca), exist_ok=True)

# Load data
df_norm = pd.read_csv(input_data, index_col=0)
print(f"Data shape: {df_norm.shape}")

# Parse metadata for strain mapping
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

# PCA Analysis
df_transposed = df_norm.T
gene_variances = df_transposed.var(axis=0)
df_filtered = df_transposed.loc[:, gene_variances > 0]

scaler = StandardScaler()
df_scaled = scaler.fit_transform(df_filtered)

pca = PCA()
pca_result = pca.fit_transform(df_scaled)

print(f"PC1: {pca.explained_variance_ratio_[0]:.1%}")
print(f"PC2: {pca.explained_variance_ratio_[1]:.1%}")

# Generate colors for strains
unique_strains = sorted(set(sample_to_strain.values()))
colors = cm.tab20(np.linspace(0, 1, 20))
colors2 = cm.tab20b(np.linspace(0, 1, 20))
all_colors = np.vstack([colors, colors2])
strain_colors = {strain: all_colors[i % len(all_colors)] for i, strain in enumerate(unique_strains)}

# PCA plot
plt.figure(figsize=(14, 10))
for i, sample in enumerate(df_filtered.index):
    strain = sample_to_strain.get(sample, 'Unknown')
    color = strain_colors.get(strain, 'gray')
    plt.scatter(pca_result[i, 0], pca_result[i, 1], s=80, alpha=0.7, 
               color=color, edgecolors='black', linewidth=0.5)

legend_elements = [Patch(facecolor=strain_colors[strain], label=strain, edgecolor='black') 
                  for strain in unique_strains]
plt.legend(handles=legend_elements, title='Mouse Strains (39 total)', 
          bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=7, ncol=2)

plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%} variance)')
plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%} variance)')
plt.title('PCA: CD4+ T Cell Expression by Mouse Strain\n(82 samples across 39 inbred strains)')
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(output_pca, dpi=300, bbox_inches='tight')
plt.close()
print(f"Saved PCA plot to {output_pca}")

# Top variable genes heatmap
import seaborn as sns
top_variable_genes = gene_variances.sort_values(ascending=False).head(50).index
top_gene_data = df_norm.loc[top_variable_genes]
top_var_zscore = (top_gene_data.T - top_gene_data.mean(axis=1)) / top_gene_data.std(axis=1)

plt.figure(figsize=(16, 10))
sns.heatmap(top_var_zscore.T, cmap='RdYlBu_r', center=0, 
            xticklabels=False, yticklabels=False, cbar_kws={'label': 'Z-score'})
plt.title('Top 50 Most Variable Genes Across 82 Samples\n(Z-score normalized)')
plt.xlabel('Samples (82 total)')
plt.ylabel('Genes (50 most variable)')
plt.tight_layout()
plt.savefig(output_heatmap, dpi=300, bbox_inches='tight')
plt.close()
print(f"Saved heatmap to {output_heatmap}")