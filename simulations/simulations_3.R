### Robustness of DeConveil to CNV noise - sinthetic data ###

setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/")
pkgs <- c("tidyverse", "compcodeR", "DESeq2")
sapply(pkgs, require, character.only = TRUE)

simulate_data_with_cn_noise <- function(n_samples, n_genes, replicate, noise_levels = c(0.05, 0.10, 0.15, 0.20)) {
  # Load input data
  rna <- read.csv("TCGA/BRCA/test/rna.csv") %>%
    remove_rownames() %>% column_to_rownames("X")
  metadata <- read.csv("TCGA/BRCA/test/metadata.csv") %>%
    remove_rownames() %>% column_to_rownames("X")
  cnv <- read.csv("TCGA/BRCA/test/cnv.csv") %>%
    remove_rownames() %>% column_to_rownames("X")
  
  # DESeq2 to extract real dispersion and mean values
  ds <- DESeqDataSetFromMatrix(countData = rna, colData = metadata, design = ~condition)
  dds <- DESeq(ds)
  dispersion <- environment(dds@dispersionFunction)$fit$model$`disps[good]`[1:n_genes]
  mean <- environment(dds@dispersionFunction)$fit$data$means[1:n_genes]
  
  # Simulate RNA data
  rna_counts_sim <- generateSyntheticData(
    dataset = paste0("sim_", replicate),
    n.vars = n_genes,
    samples.per.cond = n_samples,
    n.diffexp = 0.5 * n_genes,
    relmeans = mean,
    dispersions = dispersion,
    seqdepth = 1e7,
    repl.id = replicate,
    output.file = NULL
  )@count.matrix
  
  # Split normal vs tumor
  rna_normal <- rna_counts_sim[, 1:n_samples]
  rna_tumor <- rna_counts_sim[, (n_samples + 1):(2 * n_samples)]
  rna_join <- cbind(rna_normal, rna_tumor)
  
  # CN values
  cn_normal <- matrix(1, nrow = n_genes, ncol = n_samples)
  rownames(cn_normal) <- rownames(rna_normal)
  colnames(cn_normal) <- colnames(rna_normal)
  
  cn_tumor <- cnv[1:n_genes, 111:220] %>% .[, 1:n_samples]
  rownames(cn_tumor) <- rownames(rna_tumor)
  colnames(cn_tumor) <- colnames(rna_tumor)
  
  # Final CN
  cn_full <- cbind(cn_normal, cn_tumor)
  cn_full <- cn_full * 2
  cn_full[cn_full > 7] <- 7
  cn_full <- cn_full / 2
  
  rna_cn_clean <- ceiling(rna_join * cn_full)
  
  metadata_sim <- data.frame(
    patID = colnames(rna_cn_clean),
    condition = rep(c("A", "B"), each = n_samples)
  ) %>% remove_rownames() %>% column_to_rownames("patID")
  
  # === Save clean data ===
  clean_dir <- file.path("deconveilCaseStudies/simulations/data/simulations_3/", "no_noise")
  dir.create(clean_dir, recursive = TRUE, showWarnings = FALSE)
  
  prefix <- paste0(replicate, "_", n_samples, "_", n_genes)
  
  write.csv(rna_cn_clean, file = file.path(clean_dir, paste0(prefix, "_rna_cn.csv")))
  write.csv(cn_full, file = file.path(clean_dir, paste0(prefix, "_cn.csv")))
  write.csv(metadata_sim, file = file.path(clean_dir, paste0(prefix, "_metadata.csv")))
  
  # === Save noisy versions ===
  for (noise_level in noise_levels) {
    cn_noisy <- cn_full
    
    # Only apply noise to tumor columns
    tumor_cols <- (n_samples + 1):(2 * n_samples)
    total_elements <- length(tumor_cols) * n_genes
    n_noisy <- round(noise_level * total_elements)
    
    # Randomly sample gene/sample (row/col) indices from tumor only
    row_idx <- sample(1:n_genes, n_noisy, replace = TRUE)
    col_idx <- sample(tumor_cols, n_noisy, replace = TRUE)
    
    for (i in seq_len(n_noisy)) {
      cn_noisy[row_idx[i], col_idx[i]] <- round(cn_noisy[row_idx[i], col_idx[i]] + sample(c(-1, -0.5, 0.5, 1), 1))
    }
    
    cn_noisy[cn_noisy < 0] <- 0
    cn_noisy[cn_noisy > 7] <- 7
    
    rna_cn_noisy <- ceiling(rna_join * cn_full)
    
    noise_tag <- paste0("noise_", sprintf("%02d", noise_level * 100))
    noisy_dir <- file.path("deconveilCaseStudies/simulations/data/simulations_3/", noise_tag)
    dir.create(noisy_dir, recursive = TRUE, showWarnings = FALSE)
    
    write.csv(rna_cn_noisy, file = file.path(noisy_dir, paste0(prefix, "_rna_cn.csv")))
    write.csv(cn_noisy, file = file.path(noisy_dir, paste0(prefix, "_cn.csv")))
    write.csv(metadata_sim, file = file.path(noisy_dir, paste0(prefix, "_metadata.csv")))
  }
}


# Run all simulations
generate_and_save_simulations <- function(sample_sizes = c(10, 20, 40, 100), gene_settings = c(1000), num_replicates = 10) {
  for (n_samples in sample_sizes) {
    for (n_genes in gene_settings) {
      for (replicate in 1:num_replicates) {
        message(sprintf("Simulating: %d samples, %d genes, replicate %d", n_samples, n_genes, replicate))
        tryCatch({
          simulate_data_with_cn_noise(n_samples, n_genes, replicate)
        }, error = function(e) {
          message(sprintf("ERROR in replicate %d: %s", replicate, e$message))
        })
      }
    }
  }
}

generate_and_save_simulations()
