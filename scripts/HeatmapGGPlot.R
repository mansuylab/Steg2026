HeatmapGGPlot <- function(
    mat,
    name,
    outdir,
    width  = 4,          # cm
    height = 4,          # cm
    cluster = TRUE,      # cluster genes (rows)
    col.negative = "red3",
    col.positive = "forestgreen",
    center = 0,          # midpoint of color scale
    na.col  = "grey90",
    legend.title = NULL,
    plot.gene.names = TRUE    # NEW argument
) {
  
  library(ggplot2)
  
  if (!is.matrix(mat)) stop("mat must be a matrix with genes in rows and samples in columns.")
  if (missing(name))   stop("Please provide 'name'.")
  if (missing(outdir)) stop("Please provide 'outdir'.")
  
  # ensure row/col names
  if (is.null(rownames(mat))) {
    rownames(mat) <- paste0("Gene", seq_len(nrow(mat)))
  }
  if (is.null(colnames(mat))) {
    colnames(mat) <- paste0("Sample", seq_len(ncol(mat)))
  }
  
  # Scale per gene
  mat <- t(scale(t(mat)))
  
  # cluster rows if requested
  if (isTRUE(cluster) && nrow(mat) > 1) {
    hc <- hclust(dist(mat))
    gene_order <- rownames(mat)[hc$order]
  } else {
    gene_order <- rownames(mat)
  }
  
  # long format
  df_long <- as.data.frame(mat)
  df_long$Gene <- rownames(mat)
  
  df_long <- tidyr::pivot_longer(
    df_long,
    cols = -Gene,
    names_to = "Sample",
    values_to = "value"
  )
  
  df_long$Gene   <- factor(df_long$Gene,   levels = rev(gene_order))
  df_long$Sample <- factor(df_long$Sample, levels = colnames(mat))
  
  # base plot
  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = Sample, y = Gene, fill = value)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(
      low  = col.negative,
      mid  = "white",
      high = col.positive,
      midpoint = center,
      na.value = na.col,
      name = legend.title
    ) +
    ggplot2::theme_minimal(base_size = 8) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::labs(x = NULL, y = NULL)
  
  # toggle gene labels
  if (!plot.gene.names) {
    p <- p +
      ggplot2::theme(
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank()
      )
  }
  
  # ensure directory exists
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  outfile <- file.path(outdir, name)
  
  ggplot2::ggsave(
    filename = outfile,
    plot = p,
    width = width,
    height = height,
    units = "cm"
  )
  
  invisible(outfile)
}
