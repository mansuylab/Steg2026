rrho <- function(resultsfolder, analysis1, analysis2, outputdir = "/Steg2026/", 
                 multipletesstoption = "BH",stepsizeoption = 100, boundaryoption = 0.05,
                 label1 = NULL, label2 = NULL, pathwayanalysis = TRUE){
  
  
  # multipletesstoption "BH" or "BY"
  # methodoption "hyper" or "fisher"
  # stepsizeoption 100
  
  # Load necessary libraries
  library(RRHO2)
  library(SummarizedExperiment)
  library(ggplot2)
  library(cowplot)
  library(gprofiler2)
  
  # Create output directory
  outputdir_full <- paste0(gsub('.{1}$', '', outputdir), "_", multipletesstoption)
  
  if (!file.exists(paste0(outputdir_full))){
    dir.create(paste0(outputdir_full))
  }
  
  # Create labels
  
  label1 <- ifelse(!is.null(label1), label1, analysis1)
  label2 <- ifelse(!is.null(label2), label2, analysis2)


  # Load in data
  se1 <- readRDS(paste0(resultsfolder, analysis1, "/SE.rds"))
  se2 <- readRDS(paste0(resultsfolder, analysis2, "/SE.rds"))
  
  DEA_table1 <- rowData(se1)$dea
  DEA_table2 <- rowData(se2)$dea
  
  # Prepare data (create df with transcripts and their corresponding rrho values (-log10(P-value) * sign(effect))) and order descendingly.
  trx_names1 <- rownames(DEA_table1)
  rrho_values1 <- -log10(DEA_table1$PValue) * sign(DEA_table1$logFC)
  rrho_data1 <- data.frame(TRX = trx_names1, Value = rrho_values1)
  rrho_data1 <- rrho_data1[order(-rrho_data1$Value), ] 
  
  trx_names2 <- rownames(DEA_table2)
  rrho_values2 <- -log10(DEA_table2$PValue) * sign(DEA_table2$logFC)
  rrho_data2 <- data.frame(TRX = trx_names2, Value = rrho_values2)
  rrho_data2 <- rrho_data2[order(-rrho_data2$Value), ] 
  
  # Analysis with RRHO2: Step 1 - Identify and filter data for common trx
  common_trx <- intersect(rrho_data1$TRX, rrho_data2$TRX)
  rrho_data1 <- rrho_data1[rrho_data1$TRX %in% common_trx, ]
  rrho_data2 <- rrho_data2[rrho_data2$TRX %in% common_trx, ]

  # Step 2 - Initialize RRHO
  RRHO_obj <- RRHO2_initialize(rrho_data1, rrho_data2, labels = c(label1, label2), 
                               log10.ind = TRUE, multipleTesting = multipletesstoption, method = "hyper", boundary = boundaryoption)
  
  # Generate and save heatmap
  pdf(paste0(outputdir_full, "/Heatmap.pdf"))
  RRHO2_heatmap(RRHO_obj)
  dev.off()
  
  # Generate Venn Diagrams to visualize overlaps in regulation types
  types <- c("dd", "uu", "ud", "du")
  #dd: down regulation in list1 and down regulation in list 2
  #uu: up regulation in list1 and up regulation in list 2
  #du: down regulation in list1 and up regulation in list 2
  #ud: up regulation in list1 and down regulation in list 2
  
  for (type in types) {
    pdf(paste0(outputdir_full, "/Venndiagramm_", type, ".pdf"))
    RRHO2_vennDiagram(RRHO_obj, type = type)
    dev.off()
  }

  saveRDS(object = RRHO_obj, file = paste0(outputdir_full, "/RRHO.rds"))
  
  
  # Pathway analysis of overlapping and non-overlapping trx using GProfiler
  if(pathwayanalysis == TRUE){
    
    # Overlapping trx
    
    # Extract overlapping trx lists and adjust names
    trx_overlap <- c(RRHO_obj$genelist_dd$gene_list_overlap_dd, RRHO_obj$genelist_uu$gene_list_overlap_uu)
    trx_overlap <- gsub("\\..*", "",trx_overlap)
    
    go_overlap <- gost(query = trx_overlap, 
                       organism = "mmusculus", ordered_query = FALSE, 
                       multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                       measure_underrepresentation = FALSE, evcodes = TRUE, 
                       user_threshold = 0.05, correction_method = "fdr", 
                       domain_scope = "annotated", custom_bg = NULL, 
                       numeric_ns = "", sources = "GO:BP", as_short_link = FALSE, highlight = TRUE)
    
    gg_overlap_plot<- publish_gosttable(go_overlap,
                                        use_colors = TRUE, 
                                        show_columns = c("term_name", "p_value","term_size", "intersection_size"),
                                        filename = NULL, ggplot = TRUE) +
      ggtitle("Pathway analysis of overlapping transcripts")
    
    ggsave("GO.overlap.pdf", plot = gg_overlap_plot, device = NULL, path = outputdir_full, width = 20, height = length(go_overlap$result$term_name)/3+2,limitsize = F,bg = NULL)
    
    # Plotting barplot
    res_go_overlap <- go_overlap$result
    res_go_overlap$percentage_genes <- res_go_overlap$intersection_size / res_go_overlap$term_size
    
    barplot_p_overlap <- ggplot(res_go_overlap[1:20,], aes(x = -log10(p_value), y = reorder(term_name, -p_value)) ) +
      geom_bar(stat = "identity") +
      labs(x = "-log10(adj pValue)", y = "Pathways", title = "Pathway Enrichment") +
      theme_minimal()
    ggsave(paste0(outputdir_full, "/GO_overlap_bar_plot.svg"), barplot_p_overlap, width = 8, height = 6)
    ggsave(paste0(outputdir_full, "/GO_overlap_bar_plot.pdf"), barplot_p_overlap, width = 8, height = 6)
    
    saveRDS(object = go_overlap, file = paste0(outputdir_full, "/GO.overlap.rds"))
    # Non-overlapping trx
    trx_non_overlap <- c(RRHO_obj$genelist_du$gene_list_overlap_du, RRHO_obj$genelist_ud$gene_list_overlap_ud)
    trx_non_overlap <- gsub("\\..*", "",trx_non_overlap)
    
    go_non_overlap <- gost(query = trx_non_overlap, 
                           organism = "mmusculus", ordered_query = FALSE, 
                           multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                           measure_underrepresentation = FALSE, evcodes = TRUE, 
                           user_threshold = 0.05, correction_method = "fdr", 
                           domain_scope = "annotated", custom_bg = NULL, 
                           numeric_ns = "", sources = "GO:BP", as_short_link = FALSE, highlight = TRUE)
    
    gg_non_overlap_plot<- publish_gosttable(go_non_overlap,
                                            use_colors = TRUE, 
                                            show_columns = c("term_name", "p_value","term_size", "intersection_size"),
                                            filename = NULL, ggplot = TRUE) +
      ggtitle("Pathway analysis of non-overlapping transcripts")
    
    ggsave("GO.non.overlap.pdf", plot = gg_non_overlap_plot, device = NULL, path = outputdir_full, width = 20, height = length(go_non_overlap$result$term_name)/3+2, limitsize = F,bg = NULL)
    
    # Plotting barplot
    res_go_non_overlap <- go_non_overlap$result
    res_go_non_overlap$percentage_genes <- res_go_non_overlap$intersection_size / res_go_non_overlap$term_size
    
    barplot_p_non_overlap <- ggplot(res_go_non_overlap[1:20,], aes(x = -log10(p_value), y = reorder(term_name, -p_value)) ) +
      geom_bar(stat = "identity") +
      labs(x = "-log10(adj pValue)", y = "Pathways", title = "Pathway Enrichment") +
      theme_minimal()
    ggsave(paste0(outputdir_full, "/GO_overlap_bar_plot.svg"), barplot_p_non_overlap, width = 8, height = 6)
    ggsave(paste0(outputdir_full, "/GO_overlap_bar_plot.pdf"), barplot_p_non_overlap, width = 8, height = 6)
    
    
    saveRDS(object = go_non_overlap, file = paste0(outputdir_full, "/GO.non.overlap.rds"))
    
  }
  
}  