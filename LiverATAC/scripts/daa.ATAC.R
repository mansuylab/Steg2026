daa.ATAC <- function(pheno, bamdir, peakdir,
                     strain, generation, sex, exclude = NULL, 
                     VarExp = "Group", CtrGroup = "CTR", ExpGroup = "LPD", 
                     CovarBio = NULL, CovarTec = NULL, minOverlap = 0.3, 
                     outputdir, suffix = NULL){
  
  # Call relevant packages
  library(SummarizedExperiment)
  library(edgeR)
  library(SEtools)
  library(ggplot2)
  library(cowplot)
  library(EnhancedVolcano)
  library(pheatmap)
  library(dichromat)
  library(gprofiler2)
  library(limma)
  library(DiffBind)
  library(dplyr)
  
  # Create output directory
  dirname <- paste0(paste0(strain, collapse = "+"),"-",generation,"-",paste0(sex, collapse = "+"), "-", suffix, "/")
  
  if (!file.exists(paste0(outputdir, "/", dirname))){
    dir.create(paste0(outputdir, "/", dirname))
  }
  
  # Prepare sample sheet
  
  ctr <- pheno$Sample[pheno$Strain %in% strain & pheno$Generation %in% generation & pheno$Sex %in% sex & pheno$Group == CtrGroup & !(pheno$Sample %in% exclude)]
  lpd <- pheno$Sample[pheno$Strain %in% strain & pheno$Generation %in% generation & pheno$Sex %in% sex & pheno$Group == ExpGroup & !(pheno$Sample %in% exclude)]
  
  samples <- data.frame(matrix(NA, nrow = sum(length(ctr), length(lpd)), ncol = 9, dimnames = list(c(), c("SampleID", "Tissue", "Factor", "Condition", "Treatment", "Replicate", "bamReads", "Peaks", "PeakCaller"))))
  
  samples$SampleID <- c(ctr, lpd)
  samples$Tissue <- "Liver"
  samples$Factor <- "ATAC"
  samples$Condition <- c(rep(CtrGroup, length(ctr)), rep(ExpGroup, length(lpd)))
  samples$Replicate <- c(c(1:length(ctr)),c(1:length(lpd)))
  samples$bamReads <- paste0(bamdir,samples$SampleID, ".bam")
  samples$Peaks <- paste0(peakdir,samples$SampleID, "_peaks.narrowPeak")
  samples$PeakCaller <- "narrow"
  
  nCtr <- length(ctr)
  nExp <- length(lpd)
  
  # Consensus peak via dba
  # Find consensus peaks with dba, quantify and normalize
  ATAC.dba <- dba(sampleSheet=samples, peakCaller = "narrow", minOverlap = minOverlap)
  ATAC.dba <- dba.count(ATAC.dba)
  ATAC.dba <- dba.normalize(ATAC.dba)
  
  # Construct se
  ATAC.ranges <- dba.peakset(ATAC.dba, bRetrieve=TRUE, DataType=DBA_DATA_GRANGES)
  
  colnames(samples)[1] <- "Sample"
  samples <- left_join(samples, pheno, by = "Sample")
  
  identical(samples$Sample, colnames(mcols(ATAC.ranges)))
  
  gr <- ATAC.ranges
  mcols(gr) <- NULL
  
  
  se <- SummarizedExperiment(assays = list(counts = as.matrix(mcols(ATAC.ranges))), rowRanges = gr,
                             colData = DataFrame(samples))
  
  # Normalize and add logcpm 
  countsLvl <- DGEList(assays(se)$counts)
  countsLvl <- calcNormFactors(countsLvl, method = "TMM")
  assays(se)$logcpm <- log1p(cpm(countsLvl))
  
  # Save se
  saveRDS(object = se, file = paste0(outputdir, dirname, "SE.rds"))
  
  # Preparing model and variables for SVA + DEA
  if(!is.null(CovarTec)){
    
    Covar.newname <- paste0(CovarTec, ".scaled")
    
    colData(se)[,Covar.newname] <- NA
    
    for (i in 1:length(CovarTec)) {
      
      colData(se)[,Covar.newname[i]] <- as.numeric(scale(colData(se)[,CovarTec[i]]))
      
    }
    
    if (!is.null(CovarBio)){
      basemodel <- paste0("~ ", paste0(VarExp, " + ", paste0(CovarBio, collapse = " + "), " + " ,paste0(Covar.newname, collapse = " + ")))
      
      model.tec <- paste0("~ ",  paste0(CovarBio, collapse = " + "), " + " ,paste0(Covar.newname, collapse = " + "))
    } else {
      basemodel <- paste0("~ ", paste0(VarExp, " + ",paste0(Covar.newname, collapse = " + ")))
      
      model.tec <- paste0("~ ", paste0(Covar.newname, collapse = " + "))
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
  
  # Apply SVA and calculate number of SVs
  form <- formula(basemodel)
  form0 <- formula(model.tec)
  
  se <- SEtools::svacor(SE = se, form = form, form0 = form0)
  nSV <- sum(grepl("SV",colnames(colData(se))))
  
  # DEA
  if (nSV > 0){
    formula <- formula(paste0(basemodel, " + ",paste0("SV", c(1:nSV), collapse = " + ")))
  } else {
    formula <- formula(basemodel)
  }
  
  design <- model.matrix(formula, data = colData(se))
  
  countsLvl <- DGEList(assays(se)$counts)
  countsLvl <- calcNormFactors(countsLvl, method = "TMM")
  
  cat(paste0("\nThe applied formula is: ",paste0(basemodel, " + ",paste0("SV", c(1:nSV), collapse = " + "))))
  
  countsLvl <- estimateDisp(countsLvl, design)
  
  fit <- glmFit(countsLvl, design)
  lrt <- glmLRT(fit,coef=paste0("Group", "LPD"))
  
  res <- lrt$table
  res$padj <- p.adjust(res$PValue, method = "BH")
  
  nSignP <- sum(res$PValue < 0.05 & abs(res$logFC) > 0.5)
  nSignFDR <- sum(res$padj < 0.05 & abs(res$logFC) > 0.5)
  
  rowData(se)$dea <- res
  
  # Save se
  saveRDS(object = se, file = paste0(outputdir, dirname, "SE.rds"))
  

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
  

  # Save se with all elements
  saveRDS(object = se, file = paste0(outputdir, dirname, "SE.rds"))
}
  
  
