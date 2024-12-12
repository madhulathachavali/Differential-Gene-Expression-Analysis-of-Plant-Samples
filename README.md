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

![Volcano plot-1](https://github.com/user-attachments/assets/96e637c2-4edb-4651-af0c-c7a615654fff)


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

## Results

### Summary of Differential Gene Expression Analysis

After performing the differential gene expression (DEG) analysis using the **DESeq2** package, differentially expressed between the various conditions (e.g., Intact vs. Digested rice roots) is as follows:

- **Total genes analyzed**: 55,991
- **Significant genes (adjusted p-value < 0.1)**: 4,193 genes
- **Upregulated genes**: 1,892 genes (7.5%)
- **Downregulated genes**: 1,301 genes (5.2%)
- **Low count genes (mean count < 10)**: 11,459 genes (45%)

### Table: Condition4 vs. Condition1

| Gene ID         | Base Mean | Log2 Fold Change | p-value  | Adjusted p-value (padj) |
|-----------------|-----------|------------------|----------|------------------------|
| LOC_Os01g01010  | 7.32      | -3.52            | 0.021    | 0.092                  |
| LOC_Os01g01030  | 3.31      | -2.15            | 0.199    | NA                     |
| LOC_Os01g01060  | 246.04    | 0.65             | 0.059    | 0.189                  |
| LOC_Os01g01295  | 22.15     | 1.71             | 0.001    | 0.010                  |
| LOC_Os01g01520  | 13.60     | -4.22            | 0.005    | 0.031                  |

> **Note**: The `log2FoldChange` represents the change in gene expression between conditions. A positive value indicates upregulation, while a negative value indicates downregulation.

### Significant Genes

Genes with an **adjusted p-value** < 0.1 were considered statistically significant. For example:

- **Upregulated genes**: LOC_Os01g01295, LOC_Os01g01520, etc.
- **Downregulated genes**: LOC_Os01g01010, LOC_Os01g01030, etc.

### Volcano Plot

A **volcano plot** was generated to visualize the relationship between the **log2 fold change** and **p-value** for each gene. The plot highlights genes that are significantly upregulated or downregulated.

![Volcano Plot](https://github.com/madhulathachavali/Differential-Gene-Expression-Analysis-of-Plant-Samples/blob/main/images/volcano_plot.png?raw=true)

> **Thresholds used**:  
> - **p-value threshold**: 0.05  
> - **Log2 fold change threshold**: ±1

### Heatmap of Significantly Differentially Expressed Genes

A heatmap was created to visualize the expression levels of significantly differentially expressed genes across the different conditions. The heatmap provides an overview of the expression patterns, with rows representing genes and columns representing samples.

![Heatmap](https://github.com/madhulathachavali/Differential-Gene-Expression-Analysis-of-Plant-Samples/blob/main/images/heatmap.png?raw=true)

---

- **Upregulated Genes**: Genes that are **upregulated** in one condition (e.g., Condition4) relative to the other (e.g., Condition1) indicate an increase in expression, which might be associated with biological processes specific to that condition.
  
- **Downregulated Genes**: Genes that are **downregulated** suggest a decrease in expression, which could reflect the suppression of certain pathways or regulatory mechanisms in the condition.

- **Key Findings**:
  - Certain genes, such as **LOC_Os01g01010** and **LOC_Os01g01520**, show significant changes in expression, potentially indicating their involvement in biological responses to the conditions analyzed (e.g., intact vs. digested rice roots).
  - Genes like **LOC_Os01g01295** show upregulation, which may be associated with stress responses or adaptation to the experimental conditions.

---

## Summary of DEG Analysis
In this analysis, we successfully identified genes with significant changes in expression between the conditions. The combination of **volcano plots** and **heatmaps** helped to visualize the differentially expressed genes, and the detailed results table provided further insight into the magnitude and significance of these changes.

