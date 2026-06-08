PCA_exploration <- function(data,
                            pheno,
                            variables,
                            PCNumber = 10,
                            outputpath = "PCA/"){
  
  # Prepare traits table and convert all into numeric
  traits <- pheno[,variables]
  #traits[sapply(traits, is.factor)] <- lapply(traits[sapply(traits, is.factor)], as.numeric)
  
  id <- rownames(traits)
  traits <- apply(traits, 2, as.numeric)
  rownames(traits) <- id
  
  library(reshape2)
  
  # Prepare data table
  data.PCA <- prcomp(t(data))
  data.PCA.out <- data.PCA$x
  pcs <- data.PCA
  pheno_PCA <- merge(pheno, data.PCA.out , by = "row.names")
  
  # Create a results matrix
  res<-matrix(ncol=1, nrow=20)
  colnames(res)<-"R2"
  rownames(res)<-paste("PC",c(seq(1,20)), sep="")
  # Square the SD to get the variance
  eigs <- pcs$sdev^2
  for(i in 1:nrow(res)){
    res[i,1]<-round((eigs[i]/sum(eigs)),4)*100
  }
  sum(res[,1])
  
  
  varianceExp <- paste0(round(summary(data.PCA)$importance[2,]* 100, digits = 2),"%")
  xlabels <- paste0("PC", c(1:PCNumber),": ", varianceExp)
  
  # Function to create a heatmap of the chosen traits against PCS
  PCheatmap<-function(CorMat){ ## this is your correlation matrix 
    
    ## 'melt' trait cor - refomatting data for plot
    melted_cormat <- melt(CorMat) 
    
    ## make a correlation plot
    ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
      geom_tile(color = "white")+
      labs(y="PCs & Variance Explained", x = "Variable") +
      scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                           midpoint = 0, limit = c(-1,1), space = "Lab", 
                           name="Pearson\nCorrelation") +
      theme_minimal()+ 
      scale_y_discrete(labels = xlabels) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                       size = 12, hjust = 1))+
      coord_fixed()
    
  }
  
  # Look at the various potential confounders in the pheno file
  
  TraitCor <- cor(traits, pcs$x[,1:PCNumber], use = "pairwise.complete.obs")
  P<-as.data.frame(matrix(ncol=PCNumber, nrow=ncol(traits)))
  colnames(P)<-colnames(TraitCor)
  rownames(P)<-rownames(TraitCor)
  for (i in 1:PCNumber){
    for (j in 1:ncol(traits)){
      P[j,i]<-signif(cor.test(traits[,j], pcs$x[,i], use = "pairwise.complete.obs", method="pearson")$p.value,2)
      
    }
    
  }
  
  
  
  # Plot PCA maps with colors
  
  plotdata <- cbind(pheno, data.PCA.out)
  
  most.sign.variables <- c()
  for (i in 1:PCNumber){
    most.sign.variables[i] <-  rownames(P[which.min(P[,i]),])
  }
  
  library(viridis)
  
  Correlation_pvalues <- P
  Correlation_values <- TraitCor
  
  
  if (!dir.exists(outputpath)){
    dir.create(outputpath)
  }
  pdf(paste0(outputpath, "PCs_Variables_Correlations_Heatmap.pdf"))
  print(PCheatmap(TraitCor))
  screeplot(pcs, npcs = min(20, length(pcs$sdev)), type="lines")
  dev.off()
  write.csv(Correlation_values, file = paste0(outputpath, "Correlation_values.csv"))
  write.csv(Correlation_pvalues, file = paste0(outputpath, "Correlation_pvalues.csv"))
  
  
}


