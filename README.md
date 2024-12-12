# Differential-Gene-Expression-Analysis-of-Plant-Samples

# load .txt files into R
# Filter out the .txt files from the list
txt_files <- txt_files[grepl("\\.txt$", txt_files)]

# Initialize an empty list to store the count data
count_data_list <- list()

# Loop through each .txt file and load its count data (column V2)
for (file in txt_files) {
  count_data <- read.table(file, header = FALSE, sep = "\t")
  count_data_list[[file]] <- count_data$V2  # Store the count values (V2)
}

# Check the structure of the list (first few items)
head(count_data_list)

# create a count matrix 
# Convert the list to a matrix
count_matrix <- do.call(cbind, count_data_list)

# Set row names to gene names (the first column in the .txt files)
rownames(count_matrix) <- read.table(txt_files[1], header = FALSE, sep = "\t")$V1

# Check the matrix dimensions and first few rows
dim(count_matrix)
head(count_matrix)

# Prepare meta data 
# Create the colData with sample information (you can customize this based on your data)
col_data <- data.frame(
    condition = factor(c("Condition1", "Condition1", "Condition2", "Condition2", "Condition3", "Condition3", "Condition4", "Condition4")),
    row.names = colnames(count_matrix)  # Sample names from the count matrix
)

# Check the colData structure
head(col_data)

# Check if the sample names in col_data match the column names in count_matrix
all(rownames(col_data) == colnames(count_matrix))  # Should return TRUE


# Run DESeq2 for differential expression analysis:

# Load DESeq2 library
library(DESeq2)

# Create the DESeqDataSet object
dds <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = col_data,
    design = ~ condition  # Design formula (adjust if needed)
)

# Run DESeq2 analysis
dds <- DESeq(dds)

# Get the results of the differential expression analysis
res <- results(dds)

# View summary of the results
summary(res)

#  MA plot
plotMA(res)

# Volcano plot (log2FoldChange vs -log10(p-value))
volcano_data <- data.frame(log2FC = res$log2FoldChange, pval = -log10(res$pvalue))
plot(volcano_data$log2FC, volcano_data$pval, 
     pch = 20, col = ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > 1, "red", "black"),
     xlab = "Log2 Fold Change", ylab = "-log10(p-value)",
     main = "Volcano Plot")


# Results 

# Filter for significant genes (adjust p-value and log2FoldChange thresholds)
significant_genes <- res[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1),]

# View the top significant genes
head(significant_genes)

