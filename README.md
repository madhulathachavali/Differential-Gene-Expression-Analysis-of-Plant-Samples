# Differential Gene Expression Analysis of RNA-Seq Data

This repository presents an RNA-Seq differential gene expression (DEG) analysis of rice roots using bulk RNA sequencing data. It includes data processing, differential expression analysis, and visualization with the DESeq2 R package.

## Dataset Information
This analysis uses RNA sequencing data from the GEO dataset GSE283509, which includes bulk RNA-Seq data from intact and digested rice roots (Oryza sativa). The data focuses on single-cell RNA-seq protoplasting induced gene filtering.

Key Details:
Conditions:
Intact rice roots: Harvested from gel-grown conditions.
Digested rice roots: Protoplasts and undigested chunks.
Platform: Illumina NextSeq 500

RNA Extraction: Zymo MagBead RNA Isolation Kit, followed by Lexogen QuantSeq 3′ FWD RNA-Seq library preparation.

Alignment: Reads aligned to MSU Rice Genome v7, counted with HTSeq-Count.

Samples: 8 total samples, with replicates from each condition.

Reference: Zhu M, Hsu C, Lucas OP, Taylor IW, Mijar M, Nolan TM, Sadanandom A, Bennett MJ, Benfey PN, Pandey BK


## Project Overview

1. **Data Preprocessing**: Importing raw gene expression count data from text files and organizing it into a count matrix.
2. **Metadata Preparation**: Annotating the samples with experimental conditions (e.g., Condition1, Condition2) to facilitate comparison across samples.
3. **Differential Gene Expression (DEG) Analysis**: Using the **DESeq2** package in R to analyze gene expression data and identify differentially expressed genes (DEGs) between experimental conditions.
4. **Visualization**: Generating visualizations such as heatmaps, volcano plots, and results tables to interpret the results of DEG analysis.
5. **Significant Genes**: Extract a list of genes that show significant differential expression between experimental conditions.

---

## Prerequisites

- **R** (version 4.0 or higher)
- **DESeq2** package
- **ggplot2** package for visualization
- **pheatmap** package for heatmap generation

install the necessary R packages using the following commands:

```r
install.packages("DESeq2")
install.packages("ggplot2")
install.packages("pheatmap")

