# Activating the sink function to redirect output and messages to the log file
log_file <- snakemake@output[['log']]
sink(log_file, append=TRUE, split=TRUE)

# Load the necessary library
library(ggplot2)

# Retrieve file paths from the Snakemake object
input_files_unfiltered <- snakemake@input[['unfiltered_peaks']]
input_files_filtered <- snakemake@input[['filtered_peaks']]

# Retrieve output file names from the Snakemake object
output_file_unfiltered <- snakemake@output[['histogram_unfiltered']]
output_file_filtered <- snakemake@output[['histogram_filtered']]
output_file_barplot <- snakemake@output[['barplot']]

# Function to count peaks
count_peaks <- function(files) {
  counts <- numeric(length(files))
  for (i in seq_along(files)) {
    peak_data <- read.table(files[i], header=FALSE)
    counts[i] <- nrow(peak_data)
  }
  return(counts)
}

peak_counts_unfiltered <- count_peaks(input_files_unfiltered)
peak_counts_filtered <- count_peaks(input_files_filtered)

# Function to generate more finely spaced breaks
generate_breaks <- function(data, factor = 2) {
  breaks <- pretty(data)
  seq(min(breaks), max(breaks), length.out = (length(breaks) - 1) * factor + 1)
}

# Histogram for unfiltered peaks
p_unfiltered <- ggplot() +
  geom_histogram(aes(x=peak_counts_unfiltered), breaks=generate_breaks(peak_counts_unfiltered)) +
  labs(x='Number of Peaks', y='Frequency', title='Unfiltered Peaks') +
  theme_minimal()

ggsave(filename=output_file_unfiltered, plot=p_unfiltered)

# Histogram for filtered peaks
p_filtered <- ggplot() +
  geom_histogram(aes(x=peak_counts_filtered), breaks=generate_breaks(peak_counts_filtered)) +
  labs(x='Number of Peaks', y='Frequency', title='Filtered Peaks') +
  theme_minimal()

ggsave(filename=output_file_filtered, plot=p_filtered)

# Stacked barplot
df <- data.frame(
  Samples = rep(basename(input_files_unfiltered), 2),
  Counts = c(peak_counts_unfiltered, peak_counts_filtered),
  Type = c(rep('Unfiltered', length(peak_counts_unfiltered)), 
           rep('Filtered', length(peak_counts_filtered)))
)

p_barplot <- ggplot(df, aes(x=factor(Samples), y=Counts, fill=Type)) +
  geom_bar(stat="identity", position="stack") +
  labs(x='Samples', y='Number of Peaks') +
  theme_minimal() +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())

ggsave(filename=output_file_barplot, plot=p_barplot)

# Stop the sink function to release the connection to the log file
sink()
