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
```
## Differential Gene Expression Analysis 

Differential Gene Expression analysis is used to identify genes whose expression levels significantly differ between experimental conditions.

Step 1: Data Preparation
Before performing DEG analysis, RNA-Seq data is loaded into the count matrix and metadata. The metadata describes each sample's experimental conditions (e.g., Condition1, Condition2, etc.)

```r

Step 1: Data Preparation

col_data <- data.frame(
  condition = factor(c("Condition1", "Condition1", "Condition2", "Condition2", 
                       "Condition3", "Condition3", "Condition4", "Condition4")),
  row.names = colnames(count_matrix)  # Sample names from the count matrix
)
```
Step 2: Create DESeq2 Dataset
With the count matrix and metadata, the DESeqDataSet object is created, which is required by DESeq2 for analysis.
```r
# Create DESeq2 dataset from count matrix and metadata

dds <- DESeqDataSetFromMatrix(countData = count_matrix, 
                              colData = col_data, 
                              design = ~condition)
```
Step 3: Run DESeq2 Analysis
```r
# Run the DESeq2 differential expression analysis
dds <- DESeq(dds)
```
Step 4: Results
```r
# Extract results from DESeq2 analysis
res <- results(dds)

# View summary of results
summary(res)
```
## Visualization

Step 1: Volcano Plot

A volcano plot defines the relationship between the log2 fold change (x-axis) and the -log10 p-value (y-axis) and it shows significantly upregulated or downregulated genes.
```r
library(ggplot2)

# Set thresholds for p-value and log2 fold change
pval_threshold <- 0.05
lfc_threshold <- 1

# Add significance categories based on thresholds
res$significance <- "Not Significant"
res$significance[which(res$pvalue < pval_threshold & res$log2FoldChange > lfc_threshold)] <- "Upregulated"
res$significance[which(res$pvalue < pval_threshold & res$log2FoldChange < -lfc_threshold)] <- "Downregulated"

# Generate volcano plot

ggplot(res, aes(x = log2FoldChange, y = -log10(pvalue), color = significance)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("Not Significant" = "gray", "Upregulated" = "red", "Downregulated" = "blue")) +
  xlim(c(-5, 5)) + ylim(c(0, 10)) +  # Adjust axis limits
  geom_vline(xintercept = c(-lfc_threshold, lfc_threshold), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed", color = "black") +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 P-value") +
  theme_minimal() +
  theme(legend.title = element_blank())

Step 2: Heatmap
visualize gene expression patterns across samples
```r
library(pheatmap)

# Select the top significant genes (padj < 0.1)
top_genes <- rownames(res[which(res$padj < 0.1),])

# Generate the heatmap for the top differentially expressed genes
pheatmap(assay(dds)[top_genes,], 
         scale = "row", 
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean")

```
