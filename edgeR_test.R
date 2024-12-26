### Test edgeR ###

setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/")
pkgs <- c("dplyr", "ggplot2", "edgeR", "tidyverse")
sapply(pkgs, require, character.only = TRUE)


sample_sizes <- c(10, 20, 40, 100)
gene_counts <- c(1000, 3000, 5000)
num_replicates <- 10

process_and_save_replicates <- function(replicate_num, n_samples, n_genes) {
  # Construct dynamic file paths for each replicate
  rna_path <- sprintf("CN-aware-DGE/simulations/data/replicates/%d_rna_join_%d_%d.csv", replicate_num, n_samples, n_genes)
  cnv_path <- sprintf("CN-aware-DGE/simulations/data/replicates/%d_cn_join_%d_%d.csv", replicate_num, n_samples, n_genes)
  metadata_path <- sprintf("CN-aware-DGE/simulations/data/replicates/%d_metadata_%d_%d.csv", replicate_num, n_samples, n_genes)
  res_naive_pydeseq_path <- sprintf("CN-aware-DGE/simulations/results/replicates_pydeseq/cn_naive/%d_res_CNnaive_%d_%d.csv", replicate_num, n_samples, n_genes)
  
  # Load data
  rna <- read.csv(rna_path) %>% remove_rownames() %>% column_to_rownames("X")
  cnv <- read.csv(cnv_path) %>% remove_rownames() %>% column_to_rownames("X")
  metadata <- read.csv(metadata_path) %>% remove_rownames() %>% column_to_rownames("X")
  res_naive_pydeseq <- read.csv(res_naive_pydeseq_path) %>% remove_rownames() %>% column_to_rownames("X")
  
  # CN-naive 
  
  cn_naive <- function(rna, metadata) {
    design <- model.matrix(~1+condition, data=metadata)
    edger.obj <- edgeR::DGEList(rna)
    edger.obj <- edgeR::calcNormFactors(edger.obj, method="TMM")
    edger.obj <- edgeR::estimateDisp(edger.obj, design)
    fit <- edgeR::glmFit(edger.obj, design)
    lrt <- edgeR::glmLRT(fit, coef=2)
    return(lrt)
  }
  
  lrt <- cn_naive(rna, metadata)
  res_naive_edge <- edgeR::topTags(lrt, n=Inf)$table
  
  # CN-aware 
  design <- model.matrix(~1+condition, data=metadata)
  edger.obj <- edgeR::DGEList(rna, group = metadata$condition)
  edger.obj <- edgeR::calcNormFactors(edger.obj)
  offset <- outer(rep(1,nrow(edger.obj)), getOffset(edger.obj)) + log(cnv)
  offset <- offset %>% filter_all(all_vars(!is.infinite(.))) %>% as.matrix()
  rna <- rna[ rownames(rna) %in% rownames(offset),]
  cnv <- cnv[ rownames(cnv) %in% rownames(rna),]
  
  cn_naive <- function(rna, metadata) {
    design <- model.matrix(~1+condition, data=metadata)
    edger.obj <- edgeR::DGEList(rna)
    edger.obj <- edgeR::calcNormFactors(edger.obj, method="TMM")
    edger.obj <- edgeR::estimateDisp(edger.obj, design)
    fit <- edgeR::glmFit(edger.obj, design)
    lrt <- edgeR::glmLRT(fit, coef=2)
    return(lrt)
  }
  
  lrt <- cn_naive(rna, metadata)
  res_naive_edge <- edgeR::topTags(lrt, n=Inf)$table
  
  fit_adj <- edgeR::glmFit(y=rna, design=design, offset=offset, dispersion = lrt[["dispersion"]])
  lrt_adj <- edgeR::glmLRT(fit_adj, coef=2)
  res_aware_edge <- edgeR::topTags(lrt_adj, n=Inf)$table
  
  
  # Align results
  
  rownames_idx <- match(rownames(res_naive_edge), rownames(res_aware_edge))
  res_aware_edge <- res_aware_edge[rownames_idx,]
  
  res_naive_edge <- res_naive_edge %>% dplyr::rename(padj = FDR)
  res_aware_edge <- res_aware_edge %>% dplyr::rename(padj = FDR)
  
  rownames_idx <- match(rownames(res_naive_pydeseq), rownames(res_naive_edge))
  res_aware_edge <- res_aware_edge[rownames_idx,] %>% na.omit()
  res_naive_edge <- res_naive_edge[rownames_idx,] %>% na.omit()
  
  # Save results
  saveRDS(res_naive_edge, file = sprintf("CN-aware-DGE/simulations/results/replicates_edgeR/cn_naive/%d_res_CNnaive_%d_%d.RDS", replicate_num, n_samples, n_genes))
  saveRDS(res_aware_edge, file = sprintf("CN-aware-DGE/simulations/results/replicates_edgeR/cn_aware/%d_res_CNaware_%d_%d.RDS", replicate_num, n_samples, n_genes))
  
  message(sprintf("Processed replicate %d for %d samples and %d genes.", replicate_num, n_samples, n_genes))
}


# Iterate over all combinations
for (repliate_num in 1:num_replicates) {
  for (n_samples in sample_sizes) {
    for (n_genes in gene_counts) {
      tryCatch({
        process_and_save_replicates(replicate_num, n_samples, n_genes)
      }, error = function(e) {
        message(sprintf("Error processing replicate %d for %d samples and %d genes: %s", replicate_num, n_samples, n_genes, e$message))
      })
    }
  }
}
message("All replicates processed.")





