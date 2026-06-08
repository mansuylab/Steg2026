# Activating the sink function to redirect output and messages to the log file
log_file <- snakemake@log[['log']]
sink(log_file,append=TRUE, split=TRUE)

library(dplyr)

# Getting input arguments from Snakemake object
input_file <- snakemake@input[['peaks']]
output_file <- snakemake@output[['filtered_peaks']]
bed_file <- snakemake@output[['bed_file']]
threshold <- snakemake@params[['threshold']]

# Reading the peak data
data <- read.table(input_file)

# Sum of V8 (the peak signal)
sum_data <- sum(data$V8)

# Add a new column V11 for normalized signal per million
data_f <- as.data.frame(data)
data_f$V11 <- data_f$V8 / (sum_data / 1000000)

# Filter the data based on V11's value and create a new logical column V12
data_f$V12 <- data_f$V11 >= threshold
df <- data_f[which(data_f$V12 == TRUE),]

# Remove the V11 and V12 columns
df <- select(df, -V11, -V12)

# Write the filtered peaks data and a BED file
write.table(df, file = output_file, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(df[,c(1:3)], file = bed_file, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# Stop the sink function to release the connection to the log file
sink()
