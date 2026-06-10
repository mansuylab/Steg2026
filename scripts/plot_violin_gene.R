plot_violin_gene <- function(se, assay = "logcpm",gene, name, outdir, ylimits = NULL) {
  
  library(tidyverse)
  library(SummarizedExperiment)
  library(cowplot)
  library(ggplot2)
  
  
  
  if (!file.exists(paste0(outdir))){
    dir.create(paste0(outdir))
  }
  
  # Extract logcpm matrix
  if (assay == "logcpm") {
    mat <- assays(se)$logcpm
  }
  if (assay == "corrected") {
    mat <- assays(se)$corrected
  }
  if (assay == "tpm") {
    mat <- assays(se)$tpm
  }
  
  
  # Find gene row
  if (!(gene %in% rownames(mat))) {
    stop("Gene not found in rownames(se).")
  }
  
  # Build a simple data frame
  df <- tibble(
    value  = as.numeric(mat[gene, ]),
    Group  = colData(se)$Group
  )
  
  # Plot
  plot <- ggplot(df, aes(x = Group, y = value)) +
    geom_violin(trim = FALSE) +
    geom_jitter(width = 0.12, size = 1, alpha = 0.8) +
    stat_summary(
      fun = mean,
      geom = "crossbar",
      width = 0.4,
      fatten = 0,
      color = "red"
    ) +
    labs(
      x = NULL,
      y = paste0("logCPM (", gene, ")")
    ) +
    theme_cowplot(font_size = 8)
  
  if (!is.null(ylim)) {
    plot <- plot + ylim(ylimits)
  }
  
  ggsave(name, plot = plot, device = NULL, path = outdir, width = 3, height = 3, units = "cm")
}