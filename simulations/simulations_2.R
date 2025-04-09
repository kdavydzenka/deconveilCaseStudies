### Simulating data for CN-Aware evaluation ###

setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/")
pkgs <- c("tidyverse", "compcodeR", "DESeq2")
sapply(pkgs, require, character.only = TRUE)


# Function to simulate data with different sample and gene settings

simulate_data <- function(n_samples, n_genes, replicate) {
  
  # Load data
  rna <- read.csv("TCGA/BRCA/test/rna.csv") %>% 
    remove_rownames %>% column_to_rownames(var = "X")
  metadata <- read.csv("TCGA/BRCA/test/metadata.csv") %>% 
    remove_rownames %>% column_to_rownames(var = "X")
  cnv <- read.csv("TCGA/BRCA/test/cnv.csv") %>% remove_rownames %>% column_to_rownames(var = "X") %>% 
    as.data.frame()
  
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
    n.diffexp = 0.4 * n_genes,
    repl.id = 1, 
    seqdepth = 1e7,
    fraction.upregulated = 0.5,
    between.group.diffdisp = FALSE,
    filter.threshold.total = 1,
    filter.threshold.mediancpm = 0,
    fraction.non.overdispersed = 0,
    relmeans = mean,
    dispersions = dispersion,
    random.outlier.high.prob = 0,
    random.outlier.low.prob = 0,
    output.file = paste0("deconveilCaseStudies/simulations/results/simulation_2/replicates_rna_counts_sim/", replicate, "_rna_counts_sim_", n_samples, "_", n_genes, "_brca.rds")
  )
  
  # Extract baseline normal and tumor data
  rna_counts_sim <- rna_counts_sim@count.matrix
  rna_normal <- rna_counts_sim[, 1:n_samples]
  rna_tumor <- rna_counts_sim[, (n_samples + 1):(2 * n_samples)]
  rna_join <- cbind(rna_normal, rna_tumor)
  colnames(cnv) <- colnames(rna_tumor)
  
  # CN for normal and tumor samples
  cn_normal <- matrix(1, nrow(rna_normal), ncol(rna_normal))
  colnames(cn_normal) <- colnames(rna_normal)
  rownames(cn_normal) <- rownames(rna_normal)
  
  cn_tumor <- cnv[1:n_genes,]
  cn_tumor <- cn_tumor[,111:220]
  cn_tumor <- cn_tumor[,1:n_samples]
  colnames(cn_tumor) <- colnames(rna_tumor)
  rownames(cn_tumor) <- rownames(rna_tumor)
  
  cn_tumor_1 <- cn_tumor[1:400,]
  cn_tumor_2 <- cn_tumor[401:1000,]
  cn_tumor_2 <- apply(cn_tumor_2, 2, function(x) ifelse(x > 1, 1, x)) 
  cn_tumor <- rbind(cn_tumor_1, cn_tumor_2)
  
  cn_join <- cbind(cn_normal, cn_tumor)
  cn_join <- cn_join * 2
  cn_join <- apply(cn_join, 2, function(x) ifelse(x > 7, 7, x)) 
  cn_join <- cn_join / 2
  rna_cn <- ceiling(rna_join * cn_join)
  
  # Generate metadata
  metadata <- data.frame(patID = colnames(rna_cn),
                         condition = rep(c("A", "B"), each = n_samples))
  metadata <- metadata %>% remove_rownames %>% column_to_rownames(var = "patID")
  
  # Write to CSV
  write.csv(rna_cn, file = paste0("deconveilCaseStudies/simulations/data/simulations_2/replicates/rna_cn_", n_samples, "_", n_genes, ".csv"))
  write.csv(cn_join, file = paste0("deconveilCaseStudies/simulations/data/simulations_2/replicates/cn_join_", n_samples, "_", n_genes, ".csv"))
  write.csv(metadata, file = paste0("deconveilCaseStudies/simulations/data/simulations_2/replicates/metadata_", n_samples, "_", n_genes, ".csv"))
}


generate_and_save_simulations <- function(n_samples, n_genes, num_replicates) {
  output_dir_rna <- "deconveilCaseStudies/simulations/data/simulations_2/replicates/"
  if (!dir.exists(output_dir_rna)) {
    dir.create(output_dir_rna, recursive = TRUE)
  }
  
  for (replicate in 1:num_replicates) {
    simulate_data(n_samples, n_genes, replicate)
    
    file.rename(
      paste0(output_dir_rna, "rna_cn_", n_samples, "_", n_genes, ".csv"),
      paste0(output_dir_rna, replicate, "_rna_cn_", n_samples, "_", n_genes, ".csv")
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
gene_settings <- c(1000) 
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

