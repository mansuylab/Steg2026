# Activating the sink function to redirect output and messages to the log file
log_file <- snakemake@log[['log']]
sink(log_file, append=TRUE, split=TRUE)

# Reading in the TSS score files
tss_files <- snakemake@input[['tss_scores']]
values_files <- snakemake@input[['tss_values']]

# Log the input files
cat("Processing TSS score files:\n")
cat(paste(tss_files, collapse="\n"), "\n")

cat("Processing TSS values files:\n")
cat(paste(values_files, collapse="\n"), "\n")

# Combine the TSS scores
tss_combined <- do.call(rbind, lapply(tss_files, function(x) {
  read.table(x, header=TRUE, sep="\t", stringsAsFactors=FALSE)
}))

# Writing the combined TSS scores to file
write.table(tss_combined, snakemake@output[['combined_scores']], row.names=FALSE, sep="\t", quote=FALSE)

# Reading and plotting the TSS values
cat("Generating TSS plots...\n")

pdf(snakemake@output[['combined_plot']])

# Set up the initial plot with the first dataset
first_values <- read.table(values_files[1], header=TRUE, sep="\t", stringsAsFactors=FALSE)
plot(100*(-9:10-.5), first_values$x, type="b", 
     main="Combined TSS Values",
     xlab="distance to TSS",
     ylab="aggregate TSS score",
     col=1, # Color for the first dataset
     ylim=c(min(sapply(values_files, function(f) min(read.table(f, header=TRUE, sep="\t", stringsAsFactors=FALSE)$x))),
            max(sapply(values_files, function(f) max(read.table(f, header=TRUE, sep="\t", stringsAsFactors=FALSE)$x)))))

# Add the rest of the datasets on top
for (i in 2:length(values_files)) {
  values <- read.table(values_files[i], header=TRUE, sep="\t", stringsAsFactors=FALSE)
  lines(100*(-9:10-.5), values$x, type="b", col=i) # Different color for each dataset
}
legend_labels <- gsub("_TSS_Values.txt", "", basename(values_files))

legend("topright", legend=legend_labels, fill=1:length(values_files), title="Samples")
dev.off()

# Stop the sink function to release the connection to the log file
sink()
