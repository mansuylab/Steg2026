# Function to run differential analysis on RNA-seq data, adjusted for Embryo RNA datasets

dea.RNA.Embryo <- function(pheno, stage, genome = NULL, exclude = NULL,
                    VarExp = "Group", CovarBio = NULL, CovarTec = NULL, scalefor = NULL,
                    CtrGroup, ExpGroup, outputdir, suffix = NULL, salmondir){
  
  # Call relevant packages
  library(SummarizedExperiment)
  library(edgeR)
  library(SEtools)
  library(ggplot2)
  library(DESeq2)
  library(cowplot)
  library(EnhancedVolcano)
  library(pheatmap)
  library(dichromat)
  library(gprofiler2)
  library(limma)
  library(tximport)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(org.Mm.eg.db)  # Mouse gene annotation database (for Ensembl mapping)
  library(biomaRt)
  library(AnnotationDbi)
  library(dplyr)
  
  
  # Create output directory
  if (is.null(exclude)){
    dirname <- paste0(stage,"-",suffix, "/")
  } else {
    dirname <- paste0(stage,"-",suffix, "-", paste0("excluded",paste0(sub("^[^_]+_", "", exclude), collapse = "+"),"/"))
  }
  

  if (!file.exists(paste0(outputdir, "/", dirname))){
    dir.create(paste0(outputdir, "/", dirname))
  }
  
  # Subset pheno for relevant samples and extract n per group
  pheno <- pheno[grep(stage, pheno$Cell_stage),]
  pheno <- pheno[!(rownames(pheno) %in% c(exclude)),]
  
  
  nCtr <- table(pheno$Group)[CtrGroup]
  nExp <- table(pheno$Group)[ExpGroup]
  
  for (i in 1:length(CovarBio)){
    if (is.factor(pheno[,CovarBio[i]])) {
      pheno[,CovarBio[i]] <- droplevels(pheno[,CovarBio[i]])
    }
  }
  
  # Read in counts from salmon using tximport
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  # Extract all the transcript names (TXNAME) as keys
  k <- keys(txdb, keytype = "TXNAME")
  # Create a mapping (tx2gene) from TXNAME (transcript IDs) to GENEID (UCSC Entrez IDs)
  tx2gene <- AnnotationDbi::select(txdb, keys = k, columns = "GENEID", keytype = "TXNAME")
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
  
  if (is.null(genome)){
    sample_dirs <- paste0(samples, "_quant")
  } else {
    sample_dirs <- paste0(samples, "_", genome, "_quant")
  }
  
  
  # Manually map the sample names to the actual directory names with _SX suffix 

  files <- file.path(salmondir, sample_dirs, "quant.sf")
  names(files) <- samples
  

  
  if (!all(file.exists(files))){
    stop("Not all salmon files detected, please check salmondir argument.")
  }
  
  
  # Use tximport and extract counts
  txi <- tximport(files, type = "salmon", tx2gene = tx2gene_mapped[, c("TXNAME", "ENSEMBL")])
  counts <- txi$counts
  
  storage.mode(counts) <- "integer"
  
  # Prepare SE
  
  if(!(identical(rownames(pheno), colnames(counts)) )){
    stop("Pheno rownames and Count table colnames do not agree.")
  }
       
  se <- SummarizedExperiment(assays = list(counts = as.matrix(counts)),
                             colData = DataFrame(pheno))
  rowData(se) <- gconvert(query = gsub("\\..*","", rownames(se)), organism = "mmusculus", target = "ENSG", filter_na = FALSE)
  
  
  
  # Use gene names as row names
  rownames(se) <- rowData(se)$name
  
  se <- se[!is.na(rownames(se)), ]
  # return(se)
  
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
  transcripts <- filterByExpr(counts <- geneExpr$counts, design = model.matrix(formula(basemodel), data = colData(se)))
  se <- se[transcripts, ]
  geneExpr <- geneExpr[transcripts, ]
  ntranscripts <- sum(transcripts)
  
  
  # Normalize using calcNormFactors and add to se as "logcpm"
  assays(se)$logcpm <- log1p(cpm(geneExpr))
  
   # Apply SVA and calculate number of SVs
  form <- formula(basemodel)
  form0 <- formula(model.tec)
  
  se <- SEtools::svacor(SE = se, form = form, form0 = form0, method = "svaseq")
  nSV <- sum(grepl("SV",colnames(colData(se))))
  
  # DEA adjusted for SVs. Output adjusted pValues and add the results table as rowData to the SE.

  
  if (nSV > 0){
    formula <- formula(paste0(basemodel, " + ",paste0("SV", c(1:nSV), collapse = " + ")))
  } else {
    formula <- formula(basemodel)
  }
  
  design <- model.matrix(formula, data = colData(se))
  

  cat(paste0("\nThe applied formula is: ",paste0(basemodel, " + ",paste0("SV", c(1:nSV), collapse = " + "))))
  
 
  
  # DEA
  geneExpr <- estimateDisp(geneExpr, design)
  
  fit <- glmFit(geneExpr, design)
  lrt <- glmLRT(fit,coef=paste0(VarExp, ExpGroup))
  
  res <- lrt$table
  res$padj <- p.adjust(res$PValue, method = "BH")
  
  nSignP <- sum(res$PValue < 0.05)
  nSignFDR <- sum(res$padj < 0.05)
  
  rowData(se)$dea <- res
  
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
  
  
  
  # Add data to se
  metadata(se)$go_up <- go_up
  metadata(se)$go_down <- go_down
  metadata(se)$go_merged <- go_merged
  
  write.csv(res, file = paste0(outputdir, dirname, "results.csv"))
  
  write.csv(rowData(se[rowData(se)$dea$padj < 0.05,]), file = paste0(outputdir, dirname, "FDR.sign.genes.csv"))
  
  # Save se with all elements
  saveRDS(object = se, file = paste0(outputdir, dirname, "SE.rds"))
  
  
  
}