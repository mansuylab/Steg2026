# Function to run differential analysis on RNA-seq data, adjusted for Liver RNA-seq F1 data

dea.RNA.Liver.F1 <- function(pheno, strain, generation, sex, exclude = NULL, mergeBy = "Father",
                               VarExp = "Group", CovarBio = NULL, CovarTec = NULL, scalefor = NULL,
                               CtrGroup, ExpGroup, outputdir, suffix = NULL, salmondir){
  
  # Call relevant packages
  library(SummarizedExperiment)
  library(edgeR)
  library(SEtools)
  library(ggplot2)
  library(cowplot)
  library(pheatmap)
  library(dichromat)
  library(gprofiler2)
  library(limma)
  library(tximport)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(org.Mm.eg.db)  # Mouse gene annotation database (for Ensembl mapping)
  library(AnnotationDbi)
  library(readxl)
  library(dplyr)
  library(biomaRt)
  
  
  # Create output directory
  dirname <- paste0(paste0(strain, collapse = "+"),"-",generation,"-",paste0(sex, collapse = "+"), "-", suffix, "/")
  
  if (!file.exists(paste0(outputdir, "/", dirname))){
    dir.create(paste0(outputdir, "/", dirname))
  }
  
  # Subset se for relevant samples and extract n per group
  pheno <- pheno[!(rownames(pheno) %in% c(exclude)) & pheno$Strain %in% strain & pheno$Generation %in% generation & pheno$Sex %in% sex,]
  
  # Read in counts from salmon using tximport
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  # Extract all the transcript names (TXNAME) as keys
  k <- keys(txdb, keytype = "TXNAME")
  # Create a mapping (tx2gene) from TXNAME (transcript IDs) to GENEID (UCSC Entrez IDs)
  tx2gene <- select(txdb, keys = k, columns = "GENEID", keytype = "TXNAME")
  # Use org.Mm.eg.db to map UCSC Entrez IDs to Ensembl Gene IDs
  ensembl_map <- mapIds(org.Mm.eg.db, 
                        keys = tx2gene$GENEID, 
                        column = "ENSEMBL", 
                        keytype = "ENTREZID", 
                        multiVals = "first")  # Takes only the first value if multiple mappings exist
  # Add the Ensembl gene IDs to the tx2gene data frame
  tx2gene$ENSEMBL <- ensembl_map[tx2gene$GENEID]
  tx2gene_mapped <- tx2gene[!is.na(tx2gene$GENEID) & !is.na(tx2gene$ENSEMBL), ]
  # Convert the ENSEMBL column from a list to a character vector
  tx2gene_mapped$ENSEMBL <- sapply(tx2gene_mapped$ENSEMBL, `[`, 1)
  
  
  # Define the directories with Salmon quant.sf files
  samples <- rownames(pheno)
  sample_dirs <- paste0(samples, "_quant")
  # Manually map the sample names to the actual directory names with _SX suffix 
  
  files <- file.path(salmondir, sample_dirs, "quant.sf")
  names(files) <- samples
  
  if (!all(file.exists(files))){
    stop("Not all salmon files detected, please check salmondir argument.")
  }
  
  
  # Use tximport and extract counts
  txi <- tximport(files, type = "salmon", tx2gene = tx2gene_mapped[, c("TXNAME", "ENSEMBL")])
  counts <- txi$counts
  
  # Merge counts and pheno file
  if (length(sex) == 1) {
    
    # Identify litters
    MergeNames <- names(table(pheno[,mergeBy]))
    
    # Merge count matrix
    
    counts.merged <- matrix(data = NA, nrow = nrow(counts), ncol = length(MergeNames), dimnames = list(rownames(counts), MergeNames))

    for (i in 1:ncol(counts.merged)) {
      
      subset <- as.data.frame(counts[,pheno[,mergeBy] == colnames(counts.merged)[i]])
      
      counts.merged[,i] <- apply(subset, 1, mean)
    }
    
    ##### Prepare merged pheno file
    
    pheno.merged <- as.data.frame(matrix(data = NA, ncol = 2 , nrow = length(MergeNames), dimnames = list(MergeNames, c(mergeBy, "Group"))))
    pheno.merged[,mergeBy] <- MergeNames
    
    # Allocate groups to MergeBy
    for (i in 1:nrow(pheno.merged)) {
      pheno.merged$Group[i] <- pheno[pheno[,mergeBy] == MergeNames[i], "Group"][1]
      
    }
    
    # Take averages of technical values per MergeBy
    for (i in 1:nrow(pheno.merged)){
      
      subset <- pheno[pheno[,mergeBy] == rownames(pheno.merged)[i],colnames(pheno.merged)[-c(1:2)]]
      pheno.merged[i,-c(1:2)] <- apply(subset, 2, median)
      
    }
    
    for (i in 1:nrow(pheno.merged)) {
      
      pheno.merged$PupPerMerged[i] <- nrow(pheno[pheno[,mergeBy] == rownames(pheno.merged)[i],])
    }
    
    pheno.merged$Father <- gsub("\\-.*", "",pheno.merged[,mergeBy])
    
    
    #Identify father
    pheno.merged$Strain <- ifelse(grepl("W",pheno.merged$Father), "PWK","B6")
    
    # Add father scoresheet data
    f0.data <- read_excel("../Summary.F0.scoresheet.notes.xlsx")
    colnames(f0.data)[1] <- "Father"
    pheno.merged <- left_join(pheno.merged, f0.data, by = "Father")
    
    se <- SummarizedExperiment(assays = list(counts = as.matrix(counts.merged)),
                                      colData = DataFrame(pheno.merged))
    rowData(se) <- gconvert(query = gsub("\\..*","", rownames(se)), organism = "mmusculus", target = "ENSG", filter_na = FALSE)
    
    # Use gene names as row names
    rownames(se) <- rowData(se)$name
    
    nCtr <- table(colData(se)$Group)[CtrGroup]
    nExp <- table(colData(se)$Group)[ExpGroup]
    
    
  }else{
    
    stop("You are trying to merge animals of different Sexes.")
    
  }
  
  
  # Preparing model and variables for SVA + DEA
  if(!is.null(CovarTec)){
    
    Covar.sc <- CovarTec[CovarTec %in% scalefor]
    Covar.non.sc <- CovarTec[!CovarTec %in% scalefor]
    Covar.sc.newname <- paste0(Covar.sc, ".scaled")
    
    colData(se)[,Covar.sc.newname] <- NA
    
    for (i in 1:length(Covar.sc)) {
      
      colData(se)[,Covar.sc.newname[i]] <- as.numeric(scale(colData(se)[,Covar.sc[i]]))
      
    }
    
    if (!is.null(CovarBio)){
      basemodel <- paste0("~ ", paste0(VarExp, " + ", paste0(CovarBio, collapse = " + "), " + " ,paste0(Covar.sc.newname, collapse = " + ")))
      
      model.tec <- paste0("~ ",  paste0(CovarBio, collapse = " + "), " + " ,paste0(Covar.sc.newname, collapse = " + "))
    } else {
      basemodel <- paste0("~ ", paste0(VarExp, " + ",paste0(Covar.sc.newname, collapse = " + ")))
      
      model.tec <- paste0("~ ", paste0(Covar.sc.newname, collapse = " + "))
    }
    
    
    
  } else {
    
    
    
    if (!is.null(CovarBio)){
      basemodel <- paste0("~ ",  VarExp, " + ", paste0(CovarBio, collapse = " + "))
      model.tec <- paste0("~ ",  paste0(CovarBio, collapse = " + "))
    } else {
      basemodel <- paste0("~ ", VarExp)
      model.tec <- "~ 1"
    }
    
    
  }
  
  # Create DGE object and normalize
  geneExpr <- DGEList(assays(se)$counts)
  geneExpr <- calcNormFactors(geneExpr, method = "TMM")
  
  # Filter using filterbycounts and subset se by needed transcripts
  transcripts <- filterByExpr(counts <- geneExpr$counts, group = colData(se)$Group)
  se <- se[transcripts, ]
  geneExpr <- geneExpr[transcripts, ]
  ntranscripts <- sum(transcripts)
  
  # Normalize using calcNormFactors and add to se as "logcpm"
  assays(se)$logcpm <- log1p(cpm(geneExpr))
  
  # Apply SVA and calculate number of SVs
  form <- formula(basemodel)
  form0 <- formula(model.tec)
  
  se <- SEtools::svacor(SE = se, form = form, form0 = form0)
  nSV <- sum(grepl("SV",colnames(colData(se))))
  
  
  # Plot PCA before adjusting for SVA
  pc_bef <- as.data.frame(prcomp(t(assays(se)$counts))$x)
  varexplained <- paste(round(summary(prcomp(t(assays(se)$counts)))$importance[2,]* 100, digits = 2),"%")
  
  plotPCAbeforeSVA <- ggplot(data = pc_bef, aes(x = PC1, y = PC2, col = colData(se)$Group)) +
    geom_point(size = 4) +
    geom_text_repel(aes(label = sub("B6-F[01]-", "", rownames(pc_bef))), 
                    box.padding = 0.35, point.padding = 0.5, 
                    segment.color = 'grey50', 
                    max.overlaps = Inf, 
                    size = 4.5) + 
    labs(colour = "Group") + 
    xlab(paste("PC1", varexplained[1])) +
    ylab(paste("PC2", varexplained[2])) +
    ggtitle(paste0("Not adjusted for SVs")) +
    theme_cowplot(10) 
  
  pdf(paste0(outputdir, dirname, "PCAbeforeSVA.pdf"))
  print(plotPCAbeforeSVA)
  dev.off()
  
  # Plot PCA after adjusting for SVA
  pc_adj <- as.data.frame(prcomp(t(assays(se)$corrected))$x)
  varexplained <- paste(round(summary(prcomp(t(assays(se)$corrected)))$importance[2,]* 100, digits = 2),"%")
  
  plotPCAafterSVA <- ggplot(data = pc_adj, aes(x = PC1, y = PC2, col = colData(se)$Group)) +
    geom_point(size = 4) +
    geom_text_repel(aes(label = sub("B6-F[01]-", "", rownames(pc_adj))), 
                    box.padding = 0.35, point.padding = 0.5, 
                    segment.color = 'grey50', 
                    max.overlaps = Inf, 
                    size = 4.5) + 
    labs(colour = "Group") + 
    xlab(paste("PC1", varexplained[1])) +
    ylab(paste("PC2", varexplained[2])) +
    ggtitle(paste0("Adjusted PCs for ", nSV," SVs")) +
    theme_cowplot(10) 
  
  
  pdf(paste0(outputdir, dirname, "PCAafterSVA.pdf"))
  print(plotPCAafterSVA)
  dev.off()
  

  # DEA adjusted for SVs. Output adjusted pValues and add the results table as rowData to the SE.
  
  
  if (nSV > 0){
    formula <- formula(paste0(basemodel, " + ",paste0("SV", c(1:nSV), collapse = " + ")))
  } else {
    formula <- formula(basemodel)
  }
  
  design <- model.matrix(formula, data = colData(se))
  
  
  cat(paste0("\nThe applied formula is: ",paste0(basemodel, " + ",paste0("SV", c(1:nSV), collapse = " + "))))
  

  geneExpr <- estimateDisp(geneExpr, design)
  
  fit <- glmFit(geneExpr, design)
  lrt <- glmLRT(fit,coef=paste0(VarExp, ExpGroup))
  
  res <- lrt$table
  res$padj <- p.adjust(res$PValue, method = "BH")
  
  nSignP <- sum(res$PValue < 0.05)
  nSignFDR <- sum(res$padj < 0.05)
  
  rowData(se)$dea <- res
  
  # Print description
  cat(paste0("\nDEA was performed using ", ntranscripts, " genes and using ", nSV, " SVs . \n", 
             nSignP, " genes had a p Value < 0.05. ", nSignFDR, " genes had a FDR-adjusted p Value < 0.05\n"))
  
  # Visualization
  
  ## P value histogram
  plotPValueHist <- ggplot(data = rowData(se)$dea, aes(x = PValue)) +
    geom_histogram(color="#e9ecef", alpha = 0.9, breaks = seq(0, 1, by = 0.05)) +
    xlab("P Values") +
    ylab("") +
    labs(title = "Histgram DEA P Values") +
    theme_cowplot()
  
  pdf(paste0(outputdir, dirname, "PValueHist.pdf"))
  print(plotPValueHist)
  dev.off()
  

  # Save se with DEA 
  saveRDS(object = se, file = paste0(outputdir, dirname, "SE.rds"))
  
  # GO analysis
  
  ## Extract transcript IDs
  
  
  transcripts_up <- rownames(rowData(se)$dea[rowData(se)$dea$padj < 0.05 & rowData(se)$dea$logFC > 0.5,])
  transcripts_up_id <- gsub("\\..*", "",transcripts_up)
  
  transcripts_down <- rownames(rowData(se)$dea[rowData(se)$dea$padj < 0.05 & rowData(se)$dea$logFC < -0.5,])
  transcripts_down_id <- gsub("\\..*", "",transcripts_down)
  
  
  transcripts_merged_id <- c(transcripts_up_id, transcripts_down_id)
  
  # GO analysis
  
  ## Run Gost
  go_up <- gost(query = transcripts_up_id, 
                organism = "mmusculus", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                measure_underrepresentation = FALSE, evcodes = TRUE, 
                user_threshold = 0.05, correction_method = "g_SCS", 
                domain_scope = "annotated", custom_bg = NULL, 
                numeric_ns = "", sources = "GO:BP", as_short_link = FALSE, highlight = TRUE)
  
  go_down <- gost(query = transcripts_down_id, 
                  organism = "mmusculus", ordered_query = FALSE, 
                  multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                  measure_underrepresentation = FALSE, evcodes = TRUE, 
                  user_threshold = 0.05, correction_method = "g_SCS", 
                  domain_scope = "annotated", custom_bg = NULL, 
                  numeric_ns = "", sources = "GO:BP", as_short_link = FALSE, highlight = TRUE)
  
  go_merged <-  gost(query = transcripts_merged_id, 
                     organism = "mmusculus", ordered_query = FALSE, 
                     multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                     measure_underrepresentation = FALSE, evcodes = TRUE, 
                     user_threshold = 0.05, correction_method = "g_SCS", 
                     domain_scope = "annotated", custom_bg = NULL, 
                     numeric_ns = "", sources = "GO:BP", as_short_link = FALSE, highlight = TRUE)
  
  
  ## Modify term_name based on highlighted column
  # go_up$result$term_id[go_up$result$highlighted == TRUE] <- 
  #   paste0(go_up$result$term_id[go_up$result$highlighted == TRUE], " (Highlighted)")
  # go_down$result$term_id[go_down$result$highlighted == TRUE] <- 
  #   paste0(go_down$result$term_id[go_down$result$highlighted == TRUE], " (Highlighted)")
  
  # Add data to se
  metadata(se)$go_up <- go_up
  metadata(se)$go_down <- go_down
  metadata(se)$go_merged <- go_merged
  
  
  # Save se with all elements
  saveRDS(object = se, file = paste0(outputdir, dirname, "SE.rds"))
  
  
  
}