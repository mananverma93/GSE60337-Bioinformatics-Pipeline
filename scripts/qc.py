import pandas as pd
import os

# Access Snakemake inputs/outputs/params
input_file = snakemake.input[0]
output_file = snakemake.output[0]
threshold = snakemake.params.threshold

# Ensure output directory exists
os.makedirs(os.path.dirname(output_file), exist_ok=True)

# Load expression matrix
df = pd.read_csv(input_file, index_col=0)
print("Original shape (probes x samples):", df.shape)
print(f"Number of samples: {df.shape[1]}")

# Quality control
df_qc = df[df.mean(axis=1) > threshold]
print("After QC shape (probes x samples):", df_qc.shape)
print(f"Probes removed: {df.shape[0] - df_qc.shape[0]}")

# Save the QCed expression matrix
df_qc.to_csv(output_file)
print("Saved QCed expression matrix to:", output_file)