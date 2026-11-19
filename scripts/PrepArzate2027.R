# PrepArzate2027.R
# This script downloads the Seurat object from Arzate-Mejia et al. and creates
# a single-cell reference file for CIBERSORTx.
#
# Reference:
# Arzate-Mejia RG, Schopp T, Uzel K, Lazar-Contes I, Mansuy IM.
# Early life stress affects the transcription and chromatin accessibility of
# spermatogonial cells. bioRxiv. 2026. doi: 10.64898/2026.08.20.745941.

# Load required libraries
library(Seurat)
library(Matrix)

# Clear workspace
rm(list = ls())

# Define paths for downloaded data and output files
path_to_data <- "Steg2026/deconvolution/Arzate2027"
path_to_outputs <- "Steg2026/deconvolution/output_files"

dir.create(path_to_data, recursive = TRUE, showWarnings = FALSE)
dir.create(path_to_outputs, recursive = TRUE, showWarnings = FALSE)

# GEO accession will be added when the dataset is public
geo_accession <- "GSEXXXXX"
seurat_file_name <- "seurat_integrated.rds"

geo_series_folder <- paste0(
  substr(geo_accession, 1, nchar(geo_accession) - 3),
  "nnn"
)
seurat_url <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/",
  geo_series_folder,
  "/",
  geo_accession,
  "/suppl/",
  seurat_file_name
)
seurat_file <- file.path(path_to_data, seurat_file_name)

if (geo_accession == "GSEXXXXX") {
  stop("Replace geo_accession with the final GEO accession number before running this script.")
}

# Download and load the Seurat object
if (!file.exists(seurat_file)) {
  download.file(seurat_url, destfile = seurat_file, mode = "wb")
}

seurat_object <- readRDS(seurat_file)

# Use raw RNA counts and the Control cells as the reference
DefaultAssay(seurat_object) <- "RNA"
seurat_object <- JoinLayers(seurat_object, assay = "RNA")
seurat_object_control <- subset(seurat_object, subset = condition == "Control")

# Get raw counts and replace cell barcodes with the refined cell type labels
count_matrix <- GetAssayData(seurat_object_control, assay = "RNA", layer = "counts")
colnames(count_matrix) <- seurat_object_control$CellType

# Create a data frame that includes GeneSymbol as the first column
count_df <- data.frame(GeneSymbol = rownames(count_matrix), count_matrix, check.names = FALSE)

# Write the reference file for CIBERSORTx
write.table(
  count_df,
  file = file.path(path_to_outputs, "Arzate2027_Control.txt"),
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)
