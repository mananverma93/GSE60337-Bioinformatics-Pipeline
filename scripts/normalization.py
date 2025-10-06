import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# Check if running via Snakemake or standalone
try:
    # Running via Snakemake
    input_file = snakemake.input[0]
    output_data = snakemake.output.data
    output_plot = snakemake.output.plot
except NameError:
    # Running standalone for testing
    print("Running in standalone mode (not via Snakemake)")
    input_file = "data/processed/GSE60337_full_qc.csv"
    output_data = "data/processed/GSE60337_full_normalized.csv"
    output_plot = "data/processed/normalization_comparison_full.png"

# Ensure output directory exists
os.makedirs(os.path.dirname(output_data), exist_ok=True)

# Load QC-filtered data
df = pd.read_csv(input_file, index_col=0)

print(f"Original data shape: {df.shape}")
print(f"Original data range: {df.min().min():.2f} to {df.max().max():.2f}")

# Check if data needs log transformation
if df.min().min() > 0 and df.max().max() > 100:
    print("Applying log2 transformation")
    df_log = np.log2(df + 1)
else:
    print("Data appears already log-transformed or normalized")
    df_log = df.copy()


# Quantile normalization function
def quantile_normalize(df):
    """
    quantile normalization implementation
    """
    sorted_df = pd.DataFrame(np.sort(df.values, axis=0), 
                           index=df.index, columns=df.columns)
    
    quantile_means = sorted_df.mean(axis=1)
    ranks = df.rank(method='average')
    normalized = pd.DataFrame(index=df.index, columns=df.columns)
    
    for col in df.columns:
        col_ranks = ranks[col]
        sorted_means = np.sort(quantile_means.values)
        normalized[col] = sorted_means[col_ranks.astype(int) - 1]
    
    return normalized

# Apply quantile normalization
df_norm = quantile_normalize(df_log)

# Save normalized data
df_norm.to_csv(output_data)

print("Normalization complete!")
print(f"Normalized shape: {df_norm.shape}")
print(f"Normalized data range: {df_norm.min().min():.2f} to {df_norm.max().max():.2f}")
print(f"Saved normalized data to {output_data}")

# Create figure with subplots
fig, axes = plt.subplots(2, 2, figsize=(15, 12))

# Sample columns for visualization
sample_indices = list(range(0, len(df_log.columns), 5))
sample_cols = df_log.columns[sample_indices]

# Boxplots
axes[0, 0].boxplot([df_log[col].dropna() for col in sample_cols], tick_labels=range(len(sample_cols)))
axes[0, 0].set_title("Before Normalization (Log2)")
axes[0, 0].set_xlabel("Samples (every 5th)")
axes[0, 0].set_ylabel("Log2 Expression")
axes[0, 0].tick_params(axis='x', rotation=45, labelsize=8)

axes[0, 1].boxplot([df_norm[col].dropna() for col in sample_cols], tick_labels=range(len(sample_cols)))
axes[0, 1].set_title("After Quantile Normalization")
axes[0, 1].set_xlabel("Samples (every 5th)")
axes[0, 1].set_ylabel("Normalized Expression")
axes[0, 1].tick_params(axis='x', rotation=45, labelsize=8)

# Quantile plots
for col in df_log.columns:
    axes[1, 0].plot(df_log[col].dropna().sort_values(), 
                   np.linspace(0, 1, len(df_log[col].dropna())), 
                   alpha=0.3, linewidth=0.5)
axes[1, 0].set_title("Before Normalization - Quantile Plot (all 82 samples)")
axes[1, 0].set_xlabel("Expression Value")
axes[1, 0].set_ylabel("Quantile")

for col in df_norm.columns:
    axes[1, 1].plot(df_norm[col].dropna().sort_values(), 
                   np.linspace(0, 1, len(df_norm[col].dropna())), 
                   alpha=0.3, linewidth=0.5)
axes[1, 1].set_title("After Normalization - Quantile Plot (all 82 samples)")
axes[1, 1].set_xlabel("Normalized Expression Value")
axes[1, 1].set_ylabel("Quantile")

plt.tight_layout()
plt.savefig(output_plot, dpi=300, bbox_inches='tight')
plt.close()

print("\n--- Normalization QC ---")
print("Before normalization:")
print(f"  Mean of means: {df_log.mean().mean():.3f}")
print(f"  Std of means: {df_log.mean().std():.3f}")
print(f"  Mean of stds: {df_log.std().mean():.3f}")

print("After normalization:")
print(f"  Mean of means: {df_norm.mean().mean():.3f}")
print(f"  Std of means: {df_norm.mean().std():.3f}")
print(f"  Mean of stds: {df_norm.std().mean():.3f}")

print(f"\nQuantile alignment check:")
print(f"All samples have identical distributions: {df_norm.apply(lambda x: x.sort_values().values).nunique(axis=1).max() == 1}")