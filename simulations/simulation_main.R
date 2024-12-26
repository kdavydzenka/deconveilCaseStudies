### Simulating RNA-seq counts for genes with different types of dosage sensitivity ###

setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/")
pkgs <- c("tidyverse", "compcodeR", "DESeq2")
sapply(pkgs, require, character.only = TRUE)


# Function to simulate data with different sample and gene settings

simulate_data <- function(n_samples, n_genes) {
  
  # Load data
  rna <- read.csv("TCGA/brca/test/rna_test_all_genes.csv") %>% 
    remove_rownames %>% column_to_rownames(var = "X")
  metadata <- read.csv("TCGA/brca/test/metadata_all_genes.csv") %>% 
    remove_rownames %>% column_to_rownames(var = "X")
  cnv <- read.csv("TCGA/brca/test/cnv_test_all_genes.csv") %>% remove_rownames %>% column_to_rownames(var = "X") %>% 
    as.data.frame()
  
  # Calculate gene numbers for each category based on `n_genes`
  n_dosage_sensitive <- ceiling(0.1 * n_genes)
  n_dosage_insensitive <- ceiling(0.4 * n_genes)
  n_non_deg <- n_genes - n_dosage_sensitive - n_dosage_insensitive  
  
  # Create DESeq2 dataset and run DESeq
  ds <- DESeq2::DESeqDataSetFromMatrix(countData = rna, colData = metadata, design = ~condition)
  dds <- DESeq2::DESeq(ds)
  
  # Extract dispersion and mean values, subset based on gene count
  dispersion <- environment(dds@dispersionFunction)[["fit"]][["model"]][["disps[good]"]][1:n_genes]
  mean <- environment(dds@dispersionFunction)[["fit"]][["data"]][["means"]][1:n_genes]
  
  # Generate synthetic RNA count data
  rna_counts_sim <- compcodeR::generateSyntheticData(
    dataset = paste0("rna_counts_sim_", n_samples, "_", n_genes, "_brca"),
    n.vars = n_genes,
    samples.per.cond = n_samples, 
    n.diffexp = n_dosage_insensitive,
    repl.id = 1, 
    seqdepth = 1e7,
    fraction.upregulated = 0.5,
    between.group.diffdisp = FALSE,
    filter.threshold.total = 1,
    filter.threshold.mediancpm = 0,
    fraction.non.overdispersed = 0,
    relmeans = mean,
    dispersions = dispersion,
    random.outlier.high.prob = 50,
    random.outlier.low.prob = 50,
    output.file = paste0("CN-aware-DGE/simulations/results/replicates_rna_counts_sim/rna_counts_sim_", n_samples, "_", n_genes, "_brca.rds")
  )
  
  # Extract baseline normal and tumor data
  rna_counts_sim <- rna_counts_sim@count.matrix
  rna_normal_baseline <- rna_counts_sim[, 1:n_samples]
  rna_tumor_baseline <- rna_counts_sim[, (n_samples + 1):(2 * n_samples)]
  rna_baseline <- cbind(rna_normal_baseline, rna_tumor_baseline)
  colnames(cnv) <- colnames(rna_tumor_baseline)
  
  # Dosage-insensitive genes
  d_insensitive <- rna_baseline[1:n_dosage_insensitive,]
  cnv_ins <- matrix(1, nrow(d_insensitive), ncol(d_insensitive))
  colnames(cnv_ins) <- colnames(d_insensitive)
  rownames(cnv_ins) <- rownames(d_insensitive)
  cnv_ins <- as.data.frame(cnv_ins)
  
  # Dosage-sensitve genes
  d_sensitive_tum <- rna_tumor_baseline[(n_dosage_insensitive + 1):(n_dosage_insensitive + n_dosage_sensitive), ]
  d_sensitive_norm <- rna_normal_baseline[(n_dosage_insensitive + 1):(n_dosage_insensitive + n_dosage_sensitive), ]
  cn_tumor <- cnv[1:n_dosage_sensitive, ]
  cn_tumor <- cn_tumor[,111:220]
  cn_tumor <- cn_tumor[,1:n_samples]
  colnames(cn_tumor) <- colnames(d_sensitive_tum)
  rownames(cn_tumor) <- rownames(d_sensitive_tum)
  d_sensitive <- cbind(d_sensitive_norm, d_sensitive_tum)
  cn_normal <- matrix(1, nrow(d_sensitive_norm), n_samples)
  cnv_sens <- cbind(cn_normal, cn_tumor)
  colnames(cnv_sens) <- colnames(d_sensitive)
  rownames(cnv_sens) <- rownames(d_sensitive)
  #cnv_sens <- cnv_sens / 2
  d_sensitive <- ceiling(d_sensitive * cnv_sens)
  
  # Non-differentially expressed genes (non-DEGs)
  no_deg <- rna_baseline[(n_dosage_insensitive + n_dosage_sensitive + 1):n_genes, ]
  cnv_nodeg <- matrix(1, nrow(no_deg), ncol(no_deg))
  colnames(cnv_nodeg) <- colnames(no_deg)
  rownames(cnv_nodeg) <- rownames(no_deg)
  
  # Combine all gene categories
  rna_join <- rbind(d_insensitive, d_sensitive, no_deg)
  cn_join <- rbind(cnv_ins, cnv_sens, cnv_nodeg)
  rownames(cn_join) <- rownames(rna_join)
  
  # Generate metadata
  metadata <- data.frame(patID = colnames(rna_join),
                         condition = rep(c("A", "B"), each = n_samples))
  metadata <- metadata %>% remove_rownames %>% column_to_rownames(var = "patID")
  
  # Write to CSV
  write.csv(rna_join, file = paste0("CN-aware-DGE/simulations/data/replicates/rna_join_", n_samples, "_", n_genes, ".csv"))
  write.csv(cn_join, file = paste0("CN-aware-DGE/simulations/data/replicates/cn_join_", n_samples, "_", n_genes, ".csv"))
  write.csv(metadata, file = paste0("CN-aware-DGE/simulations/data/replicates/metadata_", n_samples, "_", n_genes, ".csv"))
}


generate_and_save_simulations <- function(n_samples, n_genes, num_replicates) {
  output_dir_rna <- "CN-aware-DGE/simulations/data/replicates/"
  if (!dir.exists(output_dir_rna)) {
    dir.create(output_dir_rna, recursive = TRUE)
  }
  
  for (replicate in 1:num_replicates) {
    simulate_data(n_samples, n_genes)
    
    file.rename(
      paste0(output_dir_rna, "rna_join_", n_samples, "_", n_genes, ".csv"),
      paste0(output_dir_rna, replicate, "_rna_join_", n_samples, "_", n_genes, ".csv")
    )
    file.rename(
      paste0(output_dir_rna, "cn_join_", n_samples, "_", n_genes, ".csv"),
      paste0(output_dir_rna, replicate, "_cn_join_", n_samples, "_", n_genes, ".csv")
    )
    file.rename(
      paste0(output_dir_rna, "metadata_", n_samples, "_", n_genes, ".csv"),
      paste0(output_dir_rna, replicate, "_metadata_", n_samples, "_", n_genes, ".csv")
    )
    
  }
}


sample_sizes <- c(10, 20, 40, 100)   
gene_settings <- c(1000, 3000, 5000) 
num_replicates <- 10    

for (n_samples in sample_sizes) {
  for (n_genes in gene_settings) {
    message(sprintf("Running simulation for %d samples and %d genes...", n_samples, n_genes))
    tryCatch({
      generate_and_save_simulations(n_samples = n_samples, n_genes = n_genes, num_replicates = num_replicates)
      message(sprintf("Simulation completed for %d samples and %d genes.", n_samples, n_genes))
    }, error = function(e) {
      message(sprintf("Error encountered for %d samples and %d genes: %s", n_samples, n_genes, e$message))
    })
  }
}



# Sampling of Copy Number data #
#generate_cnv_data <- function(n_cols, size, values, prob) {
  #sapply(1:n_cols, function(x) sample(x = values, size = size, replace = TRUE, prob = prob))
#}

#cnv_1 <- generate_cnv_data(n_cols = 20, size = 100, values = c(1, 2, 3), prob = c(0.80, 0.10, 0.10))
#cnv_2 <- generate_cnv_data(n_cols = 20, size = 50, values = c(1, 2, 3, 4), prob = c(0.10, 0.50, 0.20, 0.20))
#cnv_3 <- generate_cnv_data(n_cols = 20, size = 50, values = c(2, 3, 4), prob = c(0.10, 0.80, 0.10))
#cnv_4 <- generate_cnv_data(n_cols = 20, size = 150, values = c(2, 3, 4, 5), prob = c(0.05, 0.05, 0.80, 0.10))
#cnv_5 <- generate_cnv_data(n_cols = 20, size = 150, values = c(2, 3, 4, 5), prob = c(0.05, 0.05, 0.10, 0.80))
#cn_tumor <- rbind(cnv_1, cnv_2, cnv_3, cnv_4, cnv_5) %>% as.matrix()


