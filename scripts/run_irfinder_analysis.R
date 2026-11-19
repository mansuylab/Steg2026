# Function for intron retention analysis of IRFinder-S results.
# Filters introns, aggregates IRratio to gene level, and compares intron retention between
# candidate genes and all other genes (background) with a two-sided Wilcoxon rank-sum test.

# Arguments:
# mat_IRratio: Numeric matrix of IRratio values (rows = introns, columns = samples).
#              Row names are unique intron IDs (e.g. Chr_Start_End_Name_Strand).
#              Values between 0 (fully spliced) and 1 (fully retained).
#              Subset to one genome (genome1 or genome2) before passing, if allele-specific.
# mat_IntronDepth: Numeric matrix of mean intron read depth with the same dimensions and dimnames as mat_IRratio.
# rowData: data.table with intron annotation in the same row order as mat_IRratio.
#          Must contain the columns "gene" (gene symbol) and "intron_type" (IRFinder classification, e.g. "clean").
# gene_list: Character vector of candidate gene symbols. All other genes with data are used as background.
# min_samples: Minimum number of samples with IntronDepth > 0 for an intron to be kept. Default 5.
# min_depth: Minimum IntronDepth for an IRratio value to be considered reliable.
#            IRratio values with IntronDepth < min_depth are set to NA. Default 3.

# Returns a list with:
# min_samples, min_depth: thresholds used
# n_introns_kept: number of introns kept after filtering
# n_candidate_genes, n_background_genes: number of candidate and background genes with data
# median_candidate, median_background: median of the mean IRratio per gene in each group
# p_value: Wilcoxon rank-sum test P value (candidate vs background)
# candidate_IR, background_IR: named vectors with mean IRratio per gene
# mat_gene: gene-level IRratio matrix (genes x samples)

# Example:
# result <- run_irfinder_analysis(mat_IRratio = mat_IRratio_input, mat_IntronDepth = mat_IntronDepth_input,
#                                 rowData = rowData_embryo, gene_list = gene_list,
#                                 min_samples = 3, min_depth = 5)

run_irfinder_analysis <- function(mat_IRratio,
                                  mat_IntronDepth,
                                  rowData,
                                  gene_list,
                                  min_samples = 5,
                                  min_depth = 3) {

  # Input checks
  stopifnot(identical(dim(mat_IRratio), dim(mat_IntronDepth)))
  stopifnot(identical(rownames(mat_IRratio), rownames(mat_IntronDepth)))
  stopifnot(nrow(mat_IRratio) == nrow(rowData))
  stopifnot(all(c("gene", "intron_type") %in% colnames(rowData)))
  stopifnot(is.character(gene_list))

  # Step 1: Filter introns by minimum sample coverage
  n_covered <- rowSums(mat_IntronDepth > 0)
  keep_introns <- n_covered >= min_samples

  mat_IR_filt <- mat_IRratio[keep_introns, ]
  mat_ID_filt <- mat_IntronDepth[keep_introns, ]
  rowData_filt <- rowData[keep_introns, ]

  # Step 2: Set IRratio to NA where IntronDepth < min_depth
  mat_IR_filt[mat_ID_filt < min_depth] <- NA

  # Step 3: Keep clean introns only
  clean_idx <- rowData_filt$intron_type == "clean"
  mat_IR_clean <- mat_IR_filt[clean_idx, ]
  rowData_clean <- rowData_filt[clean_idx, ]

  # Step 4: Aggregate to gene level (mean IRratio across introns per sample)
  genes_unique <- unique(rowData_clean$gene)
  mat_gene <- matrix(NA,
                     nrow = length(genes_unique),
                     ncol = ncol(mat_IR_clean),
                     dimnames = list(genes_unique, colnames(mat_IR_clean)))

  for(g in genes_unique) {
    idx <- rowData_clean$gene == g
    if(sum(idx) == 1) {
      mat_gene[g, ] <- mat_IR_clean[idx, ]
    } else {
      mat_gene[g, ] <- colMeans(mat_IR_clean[idx, , drop=FALSE], na.rm=TRUE)
    }
  }

  # Step 5: Mean IRratio per gene across samples
  mean_IRratio <- rowMeans(mat_gene, na.rm=TRUE)

  # Step 6: Split into candidate and background genes
  candidate_IR <- mean_IRratio[names(mean_IRratio) %in% gene_list]
  candidate_IR <- candidate_IR[!is.nan(candidate_IR)]

  background_IR <- mean_IRratio[!names(mean_IRratio) %in% gene_list]
  background_IR <- background_IR[!is.nan(background_IR)]

  # Step 7: Two-sided Wilcoxon rank-sum test
  wt <- wilcox.test(candidate_IR, background_IR)

  # Return results
  return(list(
    min_samples     = min_samples,
    min_depth       = min_depth,
    n_introns_kept  = sum(keep_introns),
    n_candidate_genes = length(candidate_IR),
    n_background_genes = length(background_IR),
    median_candidate = median(candidate_IR),
    median_background = median(background_IR),
    p_value         = wt$p.value,
    candidate_IR    = candidate_IR,
    background_IR   = background_IR,
    mat_gene        = mat_gene
  ))
}
