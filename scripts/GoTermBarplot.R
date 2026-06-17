# Function for plotting Barplot with most significantly enriched pathways with ggplot2.

# Arguments:
# GO.res: Result object from gprofiler2::gost enrichment analysis.
# color.low: Color to use for lowest intersection_size value
# color.high: Color to use for highest intersection_size value.
# output: Folder where to save the pdf file. Default "/mnt/groupMansuy/leo/play/".
# filename: Name to give to pdf file. Default "GoTermBarplot.pdf".
# Plotwidth: Width of saved pdf in cm. Default 10 cm.
# Plotheight: Height of saved pdf in cm. Default 5 cm.
# Print: Whether or not to print. Default FALSE.


GoTermBarplot <- function(GO.res, color.low = "#E1E1E1" , color.high = "#09263e",
                          xlimits = NULL, pathwaysToShow = 15, 
                            output = "/mnt/groupMansuy/leo/play/", filename = "GoTermBarplot.pdf", plotwidth = 10, plotheight = 5,
                            print = FALSE){
  
  res_go <- GO.res$result
  
  if(is.null(xlimits)){
    xlimits <- max(-log10(res_go$p_value))+1
  }
  
  
  
  plot <- ggplot(res_go[1:pathwaysToShow,], aes(x = -log10(p_value), y = reorder(term_name, -p_value), fill = intersection_size) ) +
    geom_bar(stat = "identity") +
    scale_fill_gradient(low = color.low, high = color.high, name = "Pathway intersection") +
    labs(x = "-log10(adj pValue)", y = "Pathways", title = "Pathway Enrichment") +
    xlim(0,xlimits) +
    theme_cowplot(font_size = 7)
  
  ggsave(paste0(output, filename), 
         plot, width = plotwidth, height = plotheight, units = "cm")
  
  if(print == T){
    print(plot)
  }
}







