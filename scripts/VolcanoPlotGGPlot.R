# Function for plotting Volcano plot with ggplot2.

# Arguments:
# dataframe: Input dataframe that includes at least one column with logFC values and one with Pvalues (adjusted or not).
# logFC: Name of column in dataframe with logFC values. Default "logFC".
# PValues: Name of column in dataframe with (adjusted) PValues. Default "padj".
# logFC.thres: Threshold of logFC to draw lines and highlight genes. Default 0.5.
# PVal.thres: Threshold of logFC to draw lines and highlight genes. Default 0.05.
# xlimits: Limits for x axis in plot. If none given, they are calcuated automatically based on max(abs(logFC))
# ylimits: Limits for y axis in plot in -log10(PValues).
# color.highlight: Color to use for highlighted genes. Default "#29be4f".
# Padjust: Whether or not the PValues column contains adjusted pValues. TRUE = adjusted, FALSE = not adjusted. Default TRUE
# output: Folder where to save the pdf file. Default "/mnt/groupMansuy/leo/play/".
# filename: Name to give to pdf file. Default "VolcanoGGPlot.pdf".
# Plotwidth: Width of saved pdf in cm. Default 5 cm.
# Plotheight: Height of saved pdf in cm. Default 5 cm.
# Print: Whether or not to print. Default FALSE.

VolcanoPlotGGPlot <- function(dataframe, logFC = "logFC",PValues = "padj", 
                            logFC.thres = 0.5, PVal.thres = 0.05,
                            xlimits = NULL, ylimits = NULL,
                            color.highlight = "#29be4f", Padjust = TRUE,
                            output = "/mnt/groupMansuy/leo/play/", filename = "VolcanoGGPlot.pdf", plotwidth = 5, plotheight = 5,
                            print = FALSE){
  
  library(ggplot2)
  library(cowplot)
  
  dataframe$highlight <- with(dataframe, dataframe[,PValues] < PVal.thres & abs(dataframe[,logFC]) > logFC.thres)
  dataframe$log10P <- -log10(dataframe[,PValues])
  
  thres.fdr <- -log10(PVal.thres)
  
  if(is.null(xlimits)){
    xlimits <- max(abs(dataframe[,logFC]))+1
  }
  if(is.null(ylimits)){
    ylimits <- max(abs(dataframe$log10P))+1
  }
  
  plot <- ggplot(dataframe, aes(x = dea.logFC, y = log10P)) +
    # thin lines at x=5 and y=-0.5
    geom_vline(xintercept = -(logFC.thres), linewidth = 0.3, linetype = "dashed",color = "grey40") +
    geom_vline(xintercept = logFC.thres, linewidth = 0.3, linetype = "dashed",color = "grey40") +
    geom_hline(yintercept = thres.fdr, linewidth = 0.3, color = "grey40") +
    # filled points; highlight x>5 & y<0.5
    geom_point(aes(fill = highlight),
               shape = 21, size = 2, alpha = 1,
               color = "black", stroke = 0.2) +
    scale_fill_manual(
      values = c(`TRUE` = color.highlight, `FALSE` = "grey70"),
      labels = c(`TRUE` = "Differentially expressed", `FALSE` = "other"),
      name = NULL
    ) +
    xlim(-(xlimits),xlimits) + 
    ylim(0,ylimits) +
    labs(x = "logFC", y = ifelse(Padjust == T,"-log10(adjusted pvalue)","-log10(pvalue)")) +
    theme_cowplot(font_size = 8) +
    theme(legend.position = "top")
  
  if(print == T){
    print(plot)
  }
  
  ggsave(filename, plot = plot, path=output, width = plotwidth, height = plotheight, units = "cm")
}