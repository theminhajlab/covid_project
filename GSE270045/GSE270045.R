# Data Download & Pre-processing

# Load required packages
library(tidyverse)
library(DESeq2)
library(conflicted)
library(readr)
library(dplyr)
library(GEOquery)
library(EnhancedVolcano)
library(ComplexHeatmap)

# Load count data, specifying that the first column contains gene names
count_data <- read.table("GSE270045_counts.tsv", header = TRUE, row.names = 1)

# Load metadata, specifying that the first column contains sample names
meta_data <- read.csv("GSE270045_metadata.csv", header = TRUE, row.names = 1)

# Verify the sample order is now the same
all(colnames(count_data) == rownames(meta_data))

# Convert your data frame into a matrix
counts <- as.matrix(count_data) 

# Filtering out low expression genes less than 10
keep <- rowSums(counts >= 10) >= 3
counts_filtered <- counts[keep, ]

# Explicitly ensure the matrix is of integer type.
counts_filtered <- as.matrix(round(counts_filtered))
storage.mode(counts_filtered) <- "integer"

# Ensure the 'Conditions' column is a factor
meta_data$Conditions <- factor(meta_data$Conditions)

# Reorder the levels to set 'normal' as the reference
meta_data$Conditions <- relevel(meta_data$Conditions, ref = "normal")

# 2. Differential Expression Analysis

# Create a DESeqDataSet object 
dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = meta_data,
  design = ~ Conditions
)

# Run DESeq2 analysis
dds <- DESeq(dds)

# Extract results 
results <- results(dds, contrast = c("Conditions", "covid", "normal"))

# View a summary of the results
base::summary(results)

# 02. Identify significantly DEGs
significant_degs <- subset(results, padj < 0.05 & abs(log2FoldChange) > 1)

# Sort the results by adjusted p-value
significant_degs_sorted <- significant_degs[order(significant_degs$padj), ]

# 03. Export the list of DEGs with gene names
significant_degs_df <- as.data.frame(significant_degs_sorted)
significant_degs_df <- cbind(Gene = rownames(significant_degs_df), significant_degs_df)
write.csv(significant_degs_df,
          file = "significant_degs.csv",
          row.names = FALSE)

# Generate the volcano plot
volcano <- EnhancedVolcano(
  results,
  lab = rownames(results),
  x = 'log2FoldChange',
  y = 'padj',  
  pointSize = 2.0,
  labSize = 5.0,
  col = c('#636363', '#9ecae1', '#a1d99b', '#e6550d'),
  colAlpha = 1,
  legendLabels = c('NS', 'Log2FC', 'p-value', 'p-Value & Log2FC'),
  legendPosition = 'top',
  legendLabSize = 16,
  legendIconSize = 5.0,
  title = "Volcano Plot: COVID vs Normal",
  titleLabSize = 16,
  subtitle = "",
  subtitleLabSize = 18,
  pCutoff = 0.001, 
  FCcutoff = 2,
  cutoffLineType = "dashed",
  border = "partial"
)

# Export the volcano plot
ggsave(
  filename = "figures/volcano_plot.png",
  plot = volcano,
  width = 8.5,
  height = 8,
  dpi = 300
)

# --- HEATMAP CODE USING COMPLEXHEATMAP ---
# Get the top 10 up-regulated genes based on log2FoldChange
results_df <- as.data.frame(results)
up_regulated_sorted <- results_df[order(results_df$log2FoldChange, decreasing = TRUE), ]
top_up_names <- head(rownames(up_regulated_sorted), 10)

# Get the top 10 down-regulated genes based on log2FoldChange
down_regulated_sorted <- results_df[order(results_df$log2FoldChange, decreasing = FALSE), ]
top_down_names <- head(rownames(down_regulated_sorted), 10)

# Combine the gene names to plot
genes_to_plot <- c(top_up_names, top_down_names)

# Ensure the gene names are present in the normalized counts matrix
normalized_counts <- counts(dds, normalized = TRUE)
genes_to_plot <- base::intersect(genes_to_plot, rownames(normalized_counts))

# If there are genes to plot, create the heatmap
if (length(genes_to_plot) > 0) {
  # Extract normalized counts for the genes to plot
  top_genes_normalized_counts <- normalized_counts[genes_to_plot, ]
  
  # Z-score scaling for row
  mat <- t(scale(t(top_genes_normalized_counts)))
  
  # Create a column annotation object
  ha <- HeatmapAnnotation(Condition = meta_data$Conditions)
  
  # Create the heatmap and save to file
  png("heatmap.png", width = 800, height = 1000, res = 100)
  Heatmap(
    mat,
    name = "Z-score",
    row_names_gp = gpar(fontsize = 8),
    top_annotation = ha,
    show_column_names = FALSE,
    column_title = "Heatmap of Top 20 DEGs"
  )
  dev.off()
  
} else {
  print("No genes were found in your normalized counts matrix to create a heatmap.")
}

# --- EXPORT TOP 15 GENES ---
# Sort the results by log2(fold change) in descending order to find the top up-regulated genes
up_regulated_sorted <- results[order(results$log2FoldChange, decreasing = TRUE), ]

# Subset to get the top 15 up-regulated genes
top15_up_genes <- as.data.frame(head(up_regulated_sorted, 15))

# Add a 'Gene' column and export
top15_up_genes <- cbind(Gene = rownames(top15_up_genes), top15_up_genes)
write.csv(top15_up_genes, file = "top15_up_regulated_genes.csv", row.names = FALSE)

# Get the top 15 down-regulated genes
down_regulated_sorted <- results[order(results$log2FoldChange, decreasing = FALSE), ]
top15_down_genes <- as.data.frame(head(down_regulated_sorted, 15))

# Add a 'Gene' column and export
top15_down_genes <- cbind(Gene = rownames(top15_down_genes), top15_down_genes)
write.csv(top15_down_genes, file = "top15_down_regulated_genes.csv", row.names = FALSE)

# Read the two CSV files into R
up_genes <- read.csv("top15_up_regulated_genes.csv")
down_genes <- read.csv("top15_down_regulated_genes.csv")

# Combine the two data frames into a single one
combined_genes <- rbind(up_genes, down_genes)

# View the first few rows of the combined data frame to confirm
head(combined_genes)

# Export the combined data frame to a new CSV file
write.csv(combined_genes, file = "top30_combined_genes.csv", row.names = FALSE)







