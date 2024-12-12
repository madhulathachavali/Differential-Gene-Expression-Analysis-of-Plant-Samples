# Differential Gene Expression Analysis of RNA-Seq Data

## Project Overview

This repository provides tools and scripts for:

1. **Data Preprocessing**: Importing raw gene expression count data from text files and organizing it into a count matrix.
2. **Metadata Preparation**: Annotating the samples with experimental conditions (e.g., Condition1, Condition2) to facilitate comparison across samples.
3. **Differential Gene Expression (DEG) Analysis**: Using the **DESeq2** package in R to analyze gene expression data and identify differentially expressed genes (DEGs) between experimental conditions.
4. **Visualization**: Generating visualizations such as heatmaps, volcano plots, and results tables to interpret the results of DEG analysis.
5. **Significant Genes**: Extract a list of genes that show significant differential expression between experimental conditions.

---

## Prerequisites

To run this analysis, you will need:

- **R** (version 4.0 or higher)
- **DESeq2** package
- **ggplot2** package for visualization
- **pheatmap** package for heatmap generation

You can install the necessary R packages using the following commands:

```r
install.packages("DESeq2")
install.packages("ggplot2")
install.packages("pheatmap")

