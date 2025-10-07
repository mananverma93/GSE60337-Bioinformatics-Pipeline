# Snakefile for GSE60337 Bioinformatics Pipeline

# Configuration
RAW_DATA = "data/raw/GSE60337_expression.csv"
METADATA = "data/raw/GSE60337_series_matrix.txt"
PROCESSED_DIR = "data/processed"

# Define all output files
rule all:
    input:
        qc = f"{PROCESSED_DIR}/GSE60337_full_qc.csv",
        normalized = f"{PROCESSED_DIR}/GSE60337_full_normalized.csv",
        norm_plots = f"{PROCESSED_DIR}/normalization_comparison_full.png",
        pca_plot = f"{PROCESSED_DIR}/pca_by_strain_full.png",
        heatmap = f"{PROCESSED_DIR}/top_variable_genes_heatmap_full.png", 
        de_results = f"{PROCESSED_DIR}/DE_full_C57BL_6J_vs_CBA_CaJ.csv",
        volcano = f"{PROCESSED_DIR}/volcano_plot_full_C57BL_6J_vs_CBA_CaJ.png",
        annotated = f"{PROCESSED_DIR}/DE_annotated_C57BL_6J_vs_CBA_CaJ.csv"  


# Rule 1: Quality Control
rule quality_control:
    input:
        RAW_DATA
    output:
        f"{PROCESSED_DIR}/GSE60337_full_qc.csv"
    params:
        threshold = 20.0
    script:
        "scripts/qc.py"

# Rule 2: Normalization
rule normalization:
    input:
        f"{PROCESSED_DIR}/GSE60337_full_qc.csv"
    output:
        data = f"{PROCESSED_DIR}/GSE60337_full_normalized.csv",
        plot = f"{PROCESSED_DIR}/normalization_comparison_full.png"
    script:
        "scripts/normalization.py"

# Rule 3: Exploratory Data Analysis
rule eda:
    input:
        data = f"{PROCESSED_DIR}/GSE60337_full_normalized.csv",
        metadata = METADATA
    output:
        pca = f"{PROCESSED_DIR}/pca_by_strain_full.png",
        heatmap = f"{PROCESSED_DIR}/top_variable_genes_heatmap_full.png"
    script:
        "scripts/eda.py"

# Rule 4: Differential Expression Analysis
rule differential_expression:
    input:
        data = f"{PROCESSED_DIR}/GSE60337_full_normalized.csv",
        metadata = METADATA
    output:
        results = f"{PROCESSED_DIR}/DE_full_C57BL_6J_vs_CBA_CaJ.csv",
        volcano = f"{PROCESSED_DIR}/volcano_plot_full_C57BL_6J_vs_CBA_CaJ.png"
    params:
        group1 = "C57BL/6J",
        group2 = "CBA/CaJ"
    script:
        "scripts/dea.py"

# Rule 5: Gene Annotation
rule gene_annotation:
    input:
        de_results = f"{PROCESSED_DIR}/DE_full_C57BL_6J_vs_CBA_CaJ.csv",
        annotation = "data/raw/GPL6246-18741.txt"  # Updated filename
    output:
        f"{PROCESSED_DIR}/DE_annotated_C57BL_6J_vs_CBA_CaJ.csv"
    script:
        "scripts/gene_annotation.py"