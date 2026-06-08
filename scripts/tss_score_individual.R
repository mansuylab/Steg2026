# Activating the sink function to redirect output and messages to the log file
log_file <- snakemake@log[['log']]
sink(log_file,append=TRUE, split=TRUE)

library(GenomicRanges)
library(ATACseqQC)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

bam.file <- snakemake@input[['bam']]
tss_score_output <- snakemake@output[['tss_score']]
tss_values_output <- snakemake@output[['tss_values']]
tss_plot_output <- snakemake@output[['tss_plot']]

# Extracting the sample name from the bam file
sample_name <- basename(bam.file)
sample_name <- sub(".bam", "", sample_name)

txs <- transcripts(TxDb.Mmusculus.UCSC.mm10.knownGene)

bam <- readBamFile(bamFile=bam.file, tag=character(0), which=GRanges("chr1", IRanges(1, 1e6)), asMates=FALSE, bigFile=TRUE)
tsse <- TSSEscore(bam, txs)

# Saving TSS score and summary
cat("Writing results to:", tss_score_output, "\n")
results <- data.frame(Sample = sample_name, TSS_Score = tsse$TSSEscore)
write.table(results, file = tss_score_output, row.names=FALSE, sep="\t", quote=FALSE)

# Save the values for plotting
cat("Writing results to:", tss_values_output, "\n")
write.table(tsse$values, file = tss_values_output, row.names=TRUE, sep="\t", quote=FALSE)

# Generating and saving the plot
pdf(tss_plot_output)
plot(100*(-9:10-.5), tsse$values, type="b", 
     xlab="distance to TSS",
     ylab="aggregate TSS score")
dev.off()

# Stop the sink function to release the connection to the log file
sink()