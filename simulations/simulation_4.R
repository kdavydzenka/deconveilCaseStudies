### Robustness of DeConveil to CNV noise - real data ###

setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/")
pkgs <- c("tidyverse")
sapply(pkgs, require, character.only = TRUE)

simulate_real_data_with_cn_noise <- function(n_samples, n_genes, noise_levels = c(0.10, 0.15, 0.20, 0.25), replicate_id = 1) {
  # Load data
  rna <- read.csv("TCGA/BRCA/test/rna.csv") %>% remove_rownames() %>% column_to_rownames("X")
  cnv <- read.csv("TCGA/BRCA/test/cnv.csv") %>% remove_rownames() %>% column_to_rownames("X")
  metadata_full <- read.csv("TCGA/BRCA/test/metadata.csv", row.names = 1)
  
  # Check sample count
  total_samples <- 2 * n_samples
  if (ncol(rna) < total_samples || ncol(cnv) < total_samples) {
    stop("Not enough samples in input RNA or CNV data.")
  }
  
  # Assume first half are A, second half are B
  condition_A_ids <- colnames(rna)[1:(ncol(rna)/2)]
  condition_B_ids <- colnames(rna)[(ncol(rna)/2 + 1):ncol(rna)]
  
  # Sample n patients from available pairs
  max_pairs <- min(length(condition_A_ids), length(condition_B_ids))
  if (n_samples > max_pairs) stop("Not enough matched samples available.")
  
  idx <- sample(1:max_pairs, n_samples)
  sampled_A <- condition_A_ids[idx]
  sampled_B <- condition_B_ids[idx]
  sample_ids <- c(sampled_A, sampled_B)
  
  # Randomly sample genes
  all_genes <- rownames(rna)
  if (length(all_genes) < n_genes) stop("Not enough genes available to sample.")
  sampled_genes <- sample(all_genes, n_genes)
  
  # Subset RNA and CNV
  rna_subset <- rna[sampled_genes, sample_ids]
  cnv_subset <- cnv[sampled_genes, sample_ids]
  
  # Build metadata
  metadata_full$sampleID <- rownames(metadata_full)
  metadata <- metadata_full[sample_ids, ]
  metadata <- metadata %>% dplyr::select(condition)
  
  # Normalize CNV
  cnv_subset <- cnv_subset * 2
  cnv_subset[cnv_subset > 6] <- 6
  
  # Prefix for filenames: replicateID_nSamples_nGenes
  prefix <- paste0(replicate_id, "_", n_samples, "_", n_genes)
  
  # === Save clean version (no noise) ===
  clean_dir <- file.path("deconveilCaseStudies/simulations/data/simulations_4/", "no_noise")
  dir.create(clean_dir, recursive = TRUE, showWarnings = FALSE)
  
  write.csv(rna_subset, file = file.path(clean_dir, paste0(prefix, "_rna.csv")))
  write.csv(cnv_subset, file = file.path(clean_dir, paste0(prefix, "_cn.csv")))
  write.csv(metadata, file = file.path(clean_dir, paste0(prefix, "_metadata.csv")))
  
  # === Generate noisy versions ===
  for (noise_level in noise_levels) {
    cn_noisy <- cnv_subset
    tumor_cols <- sampled_B
    total_elements <- length(tumor_cols) * n_genes
    n_noisy <- round(noise_level * total_elements)
    
    row_idx <- sample(1:n_genes, n_noisy, replace = TRUE)
    col_idx <- sample(tumor_cols, n_noisy, replace = TRUE)
    
    # Define noise levels and their probabilities
    noise_values <- c(-2, -1, 0, 1, 2)
    noise_probs <- c(0.05, 0.3, 0.3, 0.3, 0.05)  
    
    for (i in seq_len(n_noisy)) {
      noise <- sample(noise_values, 1, prob = noise_probs)
      cn_noisy[row_idx[i], col_idx[i]] <- round(
        cn_noisy[row_idx[i], col_idx[i]] + noise
      )
    }
    
    cn_noisy[cn_noisy < 0] <- 0
    cn_noisy[cn_noisy > 6] <- 6
    
    noisy_dir <- file.path(
      "deconveilCaseStudies/simulations/data/simulations_4/",
      paste0("noise_", sprintf("%02d", noise_level * 100))
    )
    dir.create(noisy_dir, recursive = TRUE, showWarnings = FALSE)
    
    write.csv(rna_subset, file = file.path(noisy_dir, paste0(prefix, "_rna.csv")))
    write.csv(cn_noisy, file = file.path(noisy_dir, paste0(prefix, "_cn.csv")))
    write.csv(metadata, file = file.path(noisy_dir, paste0(prefix, "_metadata.csv")))
  }
}

for (size in c(10, 20, 40, 60)) {
  for (replicate in 1:5) {
    simulate_real_data_with_cn_noise(n_samples = size, n_genes = 5000, replicate_id = replicate)
  }
}


