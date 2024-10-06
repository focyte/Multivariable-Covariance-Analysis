# Covariance Modeling
Covariance modeling is a method used to extract insights from multivariate biological data. It models and analyses the relationships between multiple predictors and a response variable, enabling a better understanding of how different factors contribute to an outcome.

For instance, in the paper "FOXA1 regulates alternative splicing in prostate cancer," we used a multivariable covariance approach to (1) fit cumulative SRG expression as a function of transcription factor expression using generalized linear regression, and (2) measure the relative contributions of these factors in the model.

Here, I apply covariance modeling to explore the relationship between key protein components of the HIF1A hypoxic pathway and their role in predicting gene expression levels related to hypoxia.


---

## Table of Contents
1. [Overview](#overview)
2. [Pipeline Stages](#pipeline-stages)
    - [Data Sources](#data-sources)
    - [Data Preparation](#data-preparation)
    - [Analysis Pipeline](#analysis-pipeline)
3. [Results](#results)

---

## Overview

This project draws on data from my research during my postdoctoral position in the group of Professor Tyson V Sharp, where we investigated the role of the tumor suppressor gene LIMD1 in cancer. LIMD1, along with PHD2 and VHL, plays a role in the proteasomal degradation of HIF1A protein. HIF1A is a transcription factor responsible for the cellular response to low oxygen levels (hypoxia). Loss of LIMD1 leads to stabilization of HIF1A protein and increased transcription of HIF1A-regulated genes. I contributed to a publication that revealed a HIF-1–LIMD1 negative regulatory feedback mechanism in lung cancer:

![HIF‐1–LIMD1-model](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6079541/bin/EMMM-10-e8304-g013.jpg)

A HIF–LIMD1 negative feedback mechanism mitigates the pro‐tumorigenic effects of hypoxia - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6079541/

In this project, we investigate how proteins in this pathway—HIF1A, PHD2, and LIMD1—contribute to hypoxic gene expression in cell lines from the Cancer Cell Line Encyclopedia (CCLE).

---

## Pipeline Stages

### Data Sources

Proteomics Data
Data on protein expression for HIF1A, PHD2, and LIMD1 was downloaded from DepMap. The data was derived from quantitative profiling by mass spectrometry across 375 cell lines from the Gygi lab. Notably, data for VHL protein was missing in several cell lines and was therefore omitted.

Transcriptomics Data
CCLE RNAseq gene expression data (RSEM, transcript) for 1019 cell lines was used: CCLE_RNAseq_rsem_transcripts_tpm_20180929.txt. To generate a hypoxia score, genes were selected from the Buffa et al. 2010 publication (link).

## Data Preparation
The following steps were taken to prepare the data:

1. Loading Data
- RNAseq data was loaded from the provided file, and relevant columns corresponding to selected cell lines were filtered.
- Proteomics data for HIF1A, PHD2, and LIMD1 was loaded from individual CSV files, with columns renamed for consistency and merged into a combined dataframe.
- The protein VHL was not included in the model because there were a large number of cell lines where data was missing for this protein.
  
2. Filtering and Merging RNAseq Data
- The RNAseq data was filtered to retain only relevant genes and cell lines, and merged with gene annotations to facilitate downstream analyses.
- RNAseq data was reshaped into a form suitable for modeling, and finally merged with proteomics data to create the full dataset.
  
3. Hypoxia Score Calculation
- A "hypoxia score" was generated using the Buffa gene list. Expression data for these genes was log-transformed (adding 1 to handle zeros), and the median expression was calculated to determine the hypoxia score for each cell line.
- Cell lines were classified as CASE (high hypoxia) or CNTR (control) based on their hypoxia score.
  
4. Output Files
- x.tsv: Protein expression values for modeling.
- y.tsv: Calculated hypoxia scores.
- merged_df.csv: Combined RNAseq and proteomics data.
- Buffa.csv: Processed hypoxia scores with classifications.
- Buffa-Score-notes.txt: Missing genes from the Buffa list.
- BuffaScore_violin_plot.pdf: Violin plot visualizing hypoxia score distribution.

## Analysis Pipeline

1. Model Fitting
- Generalized Linear Model (GLM)
- A GLM was fit to model the relationship between the predictors (protein expression levels) and the Buffa hypoxia score.
- Summary statistics and coefficients were saved for further analysis.
  
2. Coefficient Extraction and Plotting
- Visualization
- Model coefficients were extracted (excluding the intercept) and plotted.
- Significant coefficients (p-value < 0.05) were highlighted in the resulting plot (coefficients_plot.pdf).
  
3. Bootstrapping for Variable Importance
- Bootstrapping
- Bootstrapping with 1000 iterations was performed to estimate the relative importance of each predictor in explaining the hypoxia score.
- Relative contributions (as %R²) were calculated and saved, providing insight into the influence of each protein.
  
4. Combining and Analyzing Results
- Merging Results
- GLM coefficients were merged with bootstrapping results to create a comprehensive output file (results.csv).
- The R² values were adjusted based on the sign of the GLM coefficients to highlight positive and negative contributions.
  
5. Output Files
- GLM Summary: glm_model_summary_hypoxia.txt
- GLM Coefficients: glm_coefficients_hypoxia.csv
- Bootstrapping Results: relimp_summary_hypoxia.csv
- Combined Results: results.csv
- Plots: coefficients_plot.pdf, R2_plot.pdf

## Results

The results reflect the known biological roles of the proteins HIF1A, LIMD1, and PHD2 in the regulation of hypoxic gene expression.
HIF1A protein expression is positively associated with hypoxic gene expression (Buffa) score, whereas LIMD1 and PHD2 are negatively associated.
The coefficients are not large, but this could be explained by several factors including:
 - Not all genes in the Buffa Score are HIF1A dependent
 - Other transcription factors, such as HIF2A, regulate gene expression in response to hypoxia
 - The model didn't include VHL protein, an important part of this regulatory process

The covariance method shows that HIF1A has a very strong contribution to R2, explaining 74.6% of the model.

<img src="https://raw.githubusercontent.com/focyte/Multivariable-Covariance-Analysis/refs/heads/main/R2Graph.jpg" alt="R2_graph" width="500"/>




