# Activating the sink function to redirect output and messages to the log file
log_file <- snakemake@log[['log']]
sink(log_file,append=TRUE, split=TRUE)

# Load necessary libraries for data processing and visualization
library(GenomicRanges)
library(ATACseqQC)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(Rsamtools)
library(BSgenome.Mmusculus.UCSC.mm10)
library(ChIPpeakAnno)

# Get command line arguments to customize processing for different inputs
bamfile <- snakemake@input[['bam']]
output_plot <- snakemake@output[['plot']]

# Specify BAM file tags that will be read
tags <- c("AS", "XN", "XM", "XO", "XG", "NM", "MD", "YS", "YT")

# Set up known gene database for mouse genome
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# Define sequence levels and lengths from the gene database
seqlevels_txdb <- seqlevels(txdb)
seqlengths_txdb <- seqlengths(txdb)

# Specify valid sequence levels for the mouse genome
valid_seqlevels <- paste0("chr", 1)  # Here, only "chr1" is considered to decrease computational cost/time
which <- GRanges(seqnames = valid_seqlevels, ranges = IRanges(start = 1, end = seqlengths_txdb[valid_seqlevels]))

# Read in BAM file and apply coordinate shifting
gal <- readBamFile(bamfile, tag=tags, which=which, asMates=TRUE, bigFile=TRUE)
gal1 <- shiftGAlignmentsList(gal)

# Obtain gene transcripts for chr1 from the known gene database
txs <- transcripts(txdb)
txs <- txs[seqnames(txs) %in% "chr1"]

# Define the reference genome
genome <- Mmusculus

# Split the reads into different nucleosomal patterns
objs <- splitGAlignmentsByCut(gal1, txs=txs, genome=genome)

# Define transcription start sites (TSS)
TSS <- promoters(txs, upstream=0, downstream=1)
TSS <- unique(TSS)

# Custom function to estimate library size from GAlignmentsList objects
estLibSizeFromGAlignmentsList <- function(galList) {
  sapply(galList, function(gal) length(gal))
}

# Compute library size for normalization
librarySize <- estLibSizeFromGAlignmentsList(objs[c("NucleosomeFree", "mononucleosome", "dinucleosome", "trinucleosome")])
librarySizeAdjusted <- librarySize / 2 

# Calculate the ATAC-seq signal around TSSs
NTILE <- 101
dws <- ups <- 1010
seqlev <- "chr1"
sigs <- enrichedFragments(
  gal = objs[c("NucleosomeFree", "mononucleosome", "dinucleosome", "trinucleosome")], 
  TSS = TSS,
  librarySize = librarySizeAdjusted,
  seqlev = seqlev,
  TSS.filter = 0.5,
  n.tile = NTILE,
  upstream = ups,
  downstream = dws
)

# Log-transform the signals for visualization
sigs.log2 <- lapply(sigs, function(.ele) log2(.ele+1))

# Turn off the default plotting device
graphics.off()

# Normalize signals for nucleosome-free and nucleosome-bound regions
out <- featureAlignedDistribution(
  sigs, 
  reCenterPeaks(TSS, width = ups+dws),
  zeroAt = .5, 
  n.tile = NTILE, 
  type = "l", 
  ylab = "Averaged coverage"
)

# Rescale signals to a 0-1 range
range01 <- function(x) {(x-min(x)) / (max(x)-min(x))}
out <- apply(out, 2, range01)

# Create a plot of the rescaled signals
pdf(output_plot)

matplot(
  out, type = "l", 
  xaxt = "n", 
  xlab = "Position (bp)", 
  ylab = "Fraction of signal"
)
axis(1, at = seq(0, 100, by = 10) + 1, labels = c("-1K", seq(-800, 800, by = 200), "1K"), las = 2)
abline(v = seq(0, 100, by = 10) + 1, lty = 2, col = "gray")
dev.off()

# Stop the sink function to release the connection to the log file
sink()