# Analysis of Mouse Gene Expression Divergence (GSE60337)

This repository provides the reproducible bioinformatics workflow (built with Snakemake) used to analyze the **GSE60337** microarray dataset, which investigates gene expression differences between two common inbred mouse strains.

The primary goal is to identify **Differentially Expressed (DE) genes** that distinguish the **C57BL/6J** ("Black Six") and **CBA/CaJ** strains.

## Analysis Workflow: From Raw Data to Differential Expression

The pipeline processes the raw expression data through four distinct, reproducible stages:

### 1. Quality Control (QC)
The raw expression matrix is initially filtered to ensure data integrity.
Gene probes with low overall expression are removed to reduce noise and increase statistical power. This step ensures that only reliably detected genes proceed to normalization.

### 2. Normalization
Data is adjusted across all samples to ensure that differences in expression values are due to biological variation, not technical artifacts.
The filtered data is normalized (e.g., using a method like RMA or quantile normalization, depending on your script) to make sample distributions comparable.

### 3. Exploratory Data Analysis (EDA)
After normalization, the data is explored to verify batch-effect removal and identify major sources of variation.
A Principal Component Analysis (PCA) plot is generated (`pca_by_strain_full.png`) to visualize sample clustering and confirm separation based on the biological factor of interest (mouse strain).
A Heatmap of the top variable genes is generated to visualize global expression patterns.

### 4. Differential Expression (DE) Analysis
The core statistical test to identify genes whose expression levels are significantly different between the two strains.
Comparison: C57BL/6J vs. CBA/CaJ.
Method: A two-sample Student's t-test is applied for every gene to calculate raw p-values and log2fold changes.
Correction: P-values are adjusted using the FDR method to control the False Discovery Rate, resulting in final adjusted p-values for significance testing.
