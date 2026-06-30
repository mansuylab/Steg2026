# Prep_Hermann2018.R
# This script processes single-cell RNA-seq data from Hermann et al., 2018, to generate cell type reference files for CIBERSORTx.
# Clustering information is loaded directly from provided loupe files rather than being generated within this script.

# Reference:
# Hermann BP, Cheng K, Singh A, et al. The Mammalian Spermatogenesis Single-Cell Transcriptome, from Spermatogonial Stem Cells to Spermatids.
# Cell Rep. 2018;25(6):1650-1667.e8. doi: 10.1016/j.celrep.2018.10.026. PMID: 30404016; PMCID: PMC6384825.

# Loupe files with cluster IDs were downloaded from:
# Hermann, Brian (2018), “Queryable single-cell RNA-seq (10x Genomics) datasets of Human and Mouse spermatogenic cells”, 
# Mendeley Data, V1, doi: 10.17632/kxd5f8vpt4.1 

# Load required libraries
library(Seurat)
library(dplyr)
library(Matrix)
library(gplots)
library(data.table)
library(ggplot2)

# Clear workspace
rm(list = ls())

# Define paths for saving plots and outputs
path_to_plots <- "Steg2026/deconvolution/plots"
path_to_outputs <- "/Steg2026/deconvolution/output_files"

# Load cluster IDs for Adult samples
cluster_IDs <- read.csv("/Steg2026/Hermann2018_Adu_ID4sorted/Hermann2018_Adu_ID4sorted.csv", row.names = 1)
# information taken from the .cloupe files in doi: 10.17632/kxd5f8vpt4.1

# Load expression data for Adult samples
barcodes_Adu_ID4sorted <- readLines("/Steg2026/Hermann2018_Adu_ID4sorted/GSE109033_Ad-Id4GFP-bright-dim-CD9bright_filteredmatrixbarcodes.tsv")
features_Adu_ID4sorted <- read.delim("/Steg2026/Hermann2018_Adu_ID4sorted/GSE109033_Ad-Id4GFP-bright-dim-CD9bright_filteredmatrixgenes.tsv", header = FALSE, row.names = 1)
matrix_Adu_ID4sorted <- readMM(file = "/Steg2026/Hermann2018_Adu_ID4sorted/GSE109033_Ad-Id4GFP-bright-dim-CD9bright_filteredmatrix.mtx")

# Set barcodes and feature names as column and row names of the matrix
colnames(matrix_Adu_ID4sorted) <- barcodes_Adu_ID4sorted
rownames(matrix_Adu_ID4sorted) <- features_Adu_ID4sorted$V2

# Define a function to assign cell types based on cluster information
categorize_cell <- function(cell_id, cluster_IDs) {
  if (cell_id %in% rownames(cluster_IDs)) {
    return(cluster_IDs$Cells[which(rownames(cluster_IDs) == cell_id)])
  } else {
    return("Unknown")
  }
}

# Update the column names of the matrix using the provided cluster IDs
new_colnames_merged <- vector("character", length = length(colnames(matrix_Adu_ID4sorted)))
for (i in seq_along(new_colnames_merged)) {
  new_colnames_merged[i] <- categorize_cell(colnames(matrix_Adu_ID4sorted)[i], cluster_IDs)
}

# Assign new column names based on cell type
colnames(matrix_Adu_ID4sorted) <- new_colnames_merged

# Create a data frame that includes GeneSymbol as the first column
matrix_df_Adu_ID4sorted <- data.frame(GeneSymbol = rownames(matrix_Adu_ID4sorted), matrix_Adu_ID4sorted, check.names = FALSE)

# Write the data frame to a text file for CIBERSORTx
write.table(matrix_df_Adu_ID4sorted, file = paste(path_to_outputs, "/Hermann2018_Adu_ID4sorted.txt", sep = ""), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
