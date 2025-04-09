# Run DESeq2 #

setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/")
pkgs <- c("dplyr", "ggplot2", "tidyr", "DESeq2")
sapply(pkgs, require, character.only = TRUE)


sample_sizes <- c(10, 20, 40, 100)
gene_counts <- c(2000)
num_replicates <- 1

process_and_save_replicates <- function(replicate_num, n_samples, n_genes) {
  
  # Construct dynamic file paths for each replicate
  rna_path <- sprintf("deconveilCaseStudies/simulations/data/simulations_2/replicates/%d_rna_cn_%d_%d.csv", replicate_num, n_samples, n_genes)
  metadata_path <- sprintf("deconveilCaseStudies/simulations/data/simulations_2/replicates/%d_metadata_%d_%d.csv", replicate_num, n_samples, n_genes)
  res_pydeseq_path <- sprintf("deconveilCaseStudies/simulations/results/simulation_2/replicates_pydeseq/%d_res_CNnaive_%d_%d.csv", replicate_num, n_samples, n_genes)
  
  # Load data
  rna <- read.csv(rna_path) %>% remove_rownames() %>% column_to_rownames("X")
  metadata <- read.csv(metadata_path) %>% remove_rownames() %>% column_to_rownames("X")
  res_pydeseq <- read.csv(res_pydeseq_path) %>% remove_rownames() %>% column_to_rownames("X")
  
  # Test DESeq2
  ds <- DESeq2::DESeqDataSetFromMatrix(countData = rna, colData = metadata, design = ~condition)
  dds <- DESeq2::DESeq(ds)
  res_deseq <- lfcShrink(dds, coef="condition_B_vs_A", type="apeglm")
  res_deseq_df <- as.data.frame(res_deseq@listData)
  
  # Save results
  saveRDS(res_deseq_df, file = sprintf("deconveilCaseStudies/simulations/results/simulation_2/replicates_deseq/%d_res_CNnaive_%d_%d.RDS", replicate_num, n_samples, n_genes))
  
  message(sprintf("Processed replicate %d for %d samples and %d genes.", replicate_num, n_samples, n_genes))
}


# Iterate over all combinations
for (replicate_num in 1:num_replicates) {
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


