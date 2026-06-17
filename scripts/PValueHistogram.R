# Function for plotting PValue histogram with ggplot2.

# Arguments:
# dataframe: Input dataframe that includes at least one column with logFC values and one with Pvalues (adjusted or not).
# PValues: Name of column in dataframe with (adjusted) PValues. Default "padj".
# output: Folder where to save the pdf file. Default "/Steg2026/".
# filename: Name to give to pdf file. Default "VolcanoGGPlot.pdf".
# Plotwidth: Width of saved pdf in cm. Default 20 cm.
# Plotheight: Height of saved pdf in cm. Default 10 cm.
# Print: Whether or not to print. Default FALSE.

PValueHistogram <- function(dataframe, PValues = "PValue", 
                            output = "/Steg2026/", filename = "PValueHist.pdf", plotwidth = 20, plotheight = 10,
                            print = FALSE){
  
  library(ggplot2)
  library(cowplot)
  
  plot <- ggplot(data = dataframe, aes(x = dataframe[,PValues])) +
    geom_histogram(color="#e9ecef", alpha = 0.9, breaks = seq(0, 1, by = 0.05)) +
    xlab("P Values") +
    ylab("") +
    labs(title = "Histgram DEA P Values")  +
    theme_cowplot()
  
  if(print == T){
    print(plot)
  }
  
  ggsave(filename, plot = plot, path=output, width = plotwidth, height = plotheight, units = "cm")
}