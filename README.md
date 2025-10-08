# Analysis of Mouse Gene Expression Divergence (GSE60337) 🧬

This repository provides the **reproducible bioinformatics workflow** (built with **Snakemake**) used to analyze the GSE60337 microarray dataset, which investigates gene expression differences between two common inbred mouse strains.

The primary goal is to identify **Differentially Expressed (DE) genes** that distinguish the C57BL/6J ("Black Six") and CBA/CaJ strains.

---

## 📊 Project Overview

| Category | Detail |
| :--- | :--- |
| **Dataset** | GSE60337 - CD4+ T cell gene expression across 82 samples |
| **Platform** | Affymetrix Mouse Gene 1.0 ST Array (24,922 probes) |
| **Result** | **65 significantly DE genes** ($\text{FDR} < 0.05$) |
| **Key Finding** | Identified immune-related genes ($\text{H2}$-$\text{D1}$, $\text{H2}$-$\text{Ob}$, $\text{Apobec3}$) consistent with T-cell biology |

---

## 🛠 Technical Stack

* **Language:** Python 3.x
* **Libraries:** pandas, numpy, scipy, matplotlib, seaborn, scikit-learn, statsmodels
* **Workflow:** **Snakemake**
* **Analysis:** Differential expression, $\text{PCA}$, hierarchical clustering, $\text{FDR}$ correction

---

## Analysis Workflow: From Raw Data to Differential Expression

The automated $\text{Snakemake}$ pipeline processes the raw expression data through four distinct, reproducible stages:

### 1. Quality Control ($\text{QC}$)
The raw expression matrix is initially filtered to ensure data integrity and reduce noise.
* **Action:** Gene probes with low overall expression ($\text{mean}$ $\text{expression} > 20$) are removed.
* **Result:** Reduced from 24,922 to **23,236 probes**.

### 2. Normalization
Data is adjusted across all samples to ensure that differences in expression values are due to biological variation, not technical artifacts.
* **Action:** Filtered data is normalized using **log2 transformation** and **quantile normalization**.
* **Validation:** Confirmed via boxplots and correlation analysis ($>0.98$ correlation across samples).

### 3. Exploratory Data Analysis ($\text{EDA}$)
After normalization, the data is explored to verify batch-effect removal and identify major sources of variation.
* **Analysis:** **Principal Component Analysis ($\text{PCA}$)** is performed to visualize sample clustering and confirm separation based on the biological factor of interest (mouse strain).
    * $\text{PC1}$ explained $13.2\%$ variance, $\text{PC2}$ explained $7.4\%$.
    * **Strain-specific clustering observed.**
* **Analysis:** A **Heatmap** of the top 50 variable genes is generated to visualize global expression patterns.

### 4. Differential Expression ($\text{DE}$) Analysis
The core statistical test to identify genes whose expression levels are significantly different between the two strains.
* **Comparison:** C57BL/6J ($\text{n}=5$) vs. CBA/CaJ ($\text{n}=4$).
* **Method:** A two-sample **Student's $\text{t}$-test** is applied for every gene.
* **Correction:** $\text{P}$-values are adjusted using the **$\text{FDR}$ (False Discovery Rate) method (Benjamini-Hochberg)** to control the family-wise error rate.
* **Result:** **65 genes significant** at $\text{FDR} < 0.05$.
* **Annotation:** Mapped **$99.2\%$ of probes** to gene symbols using $\text{GPL}6246$ annotation.

---

## 📈 Key Statistical Insights

* **Multiple Testing Correction is Critical:** $\text{1,528}$ genes with $\text{p} < 0.05$ were reduced to just **65** with $\text{FDR} < 0.05$.
* **Sample Size Matters:** Comparing 5 vs. 4 samples yielded **65** significant genes, versus only **1** gene when comparing 2 vs. 2 samples, demonstrating the impact of statistical power.
* **High Correlation $\neq$ No Differences:** Despite high inter-sample correlation ($>0.98$), $\text{PCA}$ revealed distinct strain-specific expression patterns.

---

## 🔗 Data Source

* **GEO Accession:** [GSE60337](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60337)
* **Platform:** [GPL6246](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL6246) (Affymetrix Mouse Gene 1.0 ST Array)
