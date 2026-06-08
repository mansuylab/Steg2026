library(ATACseqQC)

# Read input arguments from Snakemake object
bamFile <- snakemake@input[['bam']]
output_file <- snakemake@output[['plot']]
bamFile.label <- gsub(".bam", "", basename(bamFile), fixed=TRUE)

# Generate fragment size distribution plot
pdf(output_file)
fragSizeDist(bamFile, bamFile.label)
dev.off()
