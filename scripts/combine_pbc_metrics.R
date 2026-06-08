# Activating the sink function to redirect output and messages to the log file
log_file <- snakemake@log[['log']]
sink(log_file, append=TRUE, split=TRUE)

# Reading in the PBC metric files
pbc_files <- snakemake@input[['pbc_metrics']]

# Log the input files
cat("Processing PBC metric files:\n")
cat(paste(pbc_files, collapse="\n"), "\n")

# Combine the PBC metrics
pbc_combined <- do.call(rbind, lapply(pbc_files, function(x) {
  data <- read.table(x, header=FALSE, sep="\t", stringsAsFactors=FALSE)
  sample_name <- gsub("\\.pbc.metrics$", "", basename(x))
  return(data.frame(Sample=sample_name, NRF=data$V5, PBC1=data$V6, PBC2=data$V7))
}))

# Writing the combined PBC metrics to file
write.table(pbc_combined, snakemake@output[['combined_metrics']], row.names=FALSE, sep="\t", quote=FALSE)
