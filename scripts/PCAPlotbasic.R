# Function for plotting PCA plot based on logcpm/count data with ggplot2.

# Arguments:
# input: Input matrix with count data, either raw, normalized or SVA corrected..
# colorVariable: Vector to color samples, i.e. Group. 
# colorVariableName: Name of Variable to color for. i.e. "Group". 
# labels: Labels for samples. Default NULL.
# svaadjusted: Whether or not the provided matrix contains SVA adjusted values.
# output: Folder where to save the pdf file. Default "/mnt/groupMansuy/leo/play/".
# filename: Name to give to pdf file. Default "VolcanoGGPlot.pdf".
# plotwidth: Width of saved pdf in cm. Default 15 cm.
# plotheight: Height of saved pdf in cm. Default 15 cm.
# print: Whether or not to print. Default FALSE.

PCAPlotbasic <- function(input, colorVariable, colorVariableName, 
                         labels = NULL,  svaadjusted = FALSE,
                         output = "/mnt/groupMansuy/leo/play/", filename = "PCAPlot.pdf", plotwidth = 15, plotheight = 15,
                         print = FALSE){
  
  library(ggplot2)
  library(cowplot)
  library(ggrepel)
  
  pc <- as.data.frame(prcomp(t(input))$x)
  varexplained <- paste(round(summary(prcomp(t(input)))$importance[2,]* 100, digits = 2),"%")
  
  if(nrow(pc) != length(colorVariable)){
    stop("colorVariable not same length as samples in input matrix")
  }
  if(nrow(pc) != length(labels)){
    stop("Labels not same length as samples in input matrix")
  }
  
  plot <- ggplot(data = pc, aes(x = PC1, y = PC2, col = colorVariable)) +
    geom_point(size = 4) +
    labs(colour = colorVariableName) +
    xlab(paste("PC1", varexplained[1])) +
    ylab(paste("PC2", varexplained[2])) +
    ggtitle(ifelse(svaadjusted, "Principal components for SVA corrected data", "Principal components for uncorrected data")) +
    theme_cowplot(8) 
    
  if(!(is.null(labels))){
    plot <- plot + geom_text_repel(aes(label = labels), 
                    box.padding = 0.35, point.padding = 0.5, 
                    segment.color = 'grey50', 
                    max.overlaps = Inf, 
                    size = 4.5) 
  } 
    
  if(print == T){
    print(plot)
  }
  
    
    
  ggsave(filename, plot = plot, path=output, width = plotwidth, height = plotheight, units = "cm")
}