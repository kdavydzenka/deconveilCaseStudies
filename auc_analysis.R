setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/")

pkgs <- c("tidyverse", "ggplot2", "edgeR", "VennDiagram", "pROC", "PRROC", "grid", "caret")
sapply(pkgs, require, character.only = TRUE)
source("CN-aware-DGE/R/utils.R")

### AUC values calculation for 10 replicates across methods and simulation configurations ###

# Load and process Results CN-aware & CN-naive methods

sample_sizes <- c(10, 20, 40, 100)
gene_counts <- c(1000, 3000, 5000)
n_replicates <- 10

load_rna_counts <- function(sample, gene) {
  file_path <- paste0("CN-aware-DGE/simulations/results/replicates_rna_counts_sim/rna_counts_sim_", sample, "_", gene, "_brca.rds")
  readRDS(file_path)
}

rna_counts_list <- list()
for (gene in gene_counts) {
  for (sample in sample_sizes) {
    rna_counts_list[[paste0("rna_", sample, "_", gene)]] <- load_rna_counts(sample, gene)
  }
}


load_results <- function(method, approach, sample, gene, replicate, extension = "csv") {
  # Set correct subdirectory paths
  method_dir <- switch(tolower(method),
                       "pydeseq" = "replicates_pydeseq",
                       "edger" = "replicates_edgeR",
                       stop("Invalid method"))
  
  approach_dir <- switch(tolower(approach),
                         "cnnaive" = "cn_naive",
                         "cnaware" = "cn_aware",
                         stop("Invalid approach"))
  
  # Construct file path and name
  file_name <- paste0(replicate, "_res_", toupper(approach), "_", sample, "_", gene, ".", extension)
  file_path <- file.path("CN-aware-DGE/simulations/results", method_dir, approach_dir, file_name)
  
  # Load file based on extension
  if (extension == "csv") {
    result <- read.csv(file_path) %>%
      remove_rownames() %>%
      column_to_rownames(var = "X")
  } else {
    result <- readRDS(file_path)
  }
  return(result)
}


results_list <- list()

# Iterate over all combinations of sample sizes, gene counts, and replicates
for (sample in sample_sizes) {
  for (gene in gene_counts) {
    for (replicate in 1:n_replicates) {
      tryCatch({
        # Load PyDESeq2 results (CN-naive and CN-aware)
        res_naive_pydeseq <- load_results("pydeseq", "cnnaive", sample, gene, replicate, "csv")
        res_aware_pydeseq <- load_results("pydeseq", "cnaware", sample, gene, replicate, "csv")
        
        # Load EdgeR results (CN-naive and CN-aware)
        res_naive_edge <- load_results("edger", "cnnaive", sample, gene, replicate, "RDS")
        res_aware_edge <- load_results("edger", "cnaware", sample, gene, replicate, "RDS")
        
        # Process results 
        common_genes <- intersect(rownames(res_aware_edge), rownames(res_aware_pydeseq))
        res_aware_pydeseq <- res_aware_pydeseq[common_genes, ]
        res_aware_edge <- res_aware_edge[common_genes, ]
        
        # Store processed results in the results list
        results_list[[paste(replicate, sample, gene, sep = "_")]] <- list(
          naive_pydeseq = res_naive_pydeseq,
          aware_pydeseq = res_aware_pydeseq,
          naive_edge = res_naive_edge,
          aware_edge = res_aware_edge
        )
        message(sprintf("Successfully loaded results for sample: %d, genes: %d, replicate: %d", sample, gene, replicate))
      }, error = function(e) {
        message(sprintf("Error loading results for sample: %d, genes: %d, replicate: %d - %s", sample, gene, replicate, e$message))
      })
    }
  }
}

process_results_pydeseq <- function(df) {
  df %>%
    dplyr::select(log2FoldChange, padj) %>%
    dplyr::rename(logFC = log2FoldChange) %>%
    na.omit()
}

process_results_edge <- function(df) {
  df %>%
    dplyr::select(logFC, padj) %>%
    na.omit()
}


processed_results <- list()

for (key in names(results_list)) {
  tryCatch({
    # Extract the current set of results
    result_set <- results_list[[key]]
    
    # Process PyDESeq2 results (CN-naive and CN-aware)
    res_naive_pydeseq <- process_results_pydeseq(result_set[["naive_pydeseq"]])
    res_aware_pydeseq <- process_results_pydeseq(result_set[["aware_pydeseq"]])
    
    # Process EdgeR results (CN-naive and CN-aware)
    res_naive_edge <- process_results_edge(result_set[["naive_edge"]])
    res_aware_edge <- process_results_edge(result_set[["aware_edge"]])
    
    # Align rownames between CN-aware PyDESeq2 and EdgeR results
    rownames_idx <- match(rownames(res_aware_pydeseq), rownames(res_aware_edge))
    res_aware_edge <- res_aware_edge[rownames_idx, ] %>% na.omit()
    res_naive_edge <- res_naive_edge[rownames_idx, ] %>% na.omit()
    
    # Save processed results in the list
    processed_results[[key]] <- list(
      naive_pydeseq = res_naive_pydeseq,
      aware_pydeseq = res_aware_pydeseq,
      naive_edge = res_naive_edge,
      aware_edge = res_aware_edge
    )
    
    message(sprintf("Successfully processed results for: %s", key))
  }, error = function(e) {
    message(sprintf("Error processing results for: %s - %s", key, e$message))
  })
}



# AUC calculation #

configurations <- c("rna_10_1000", "rna_20_1000", "rna_40_1000", "rna_100_1000", 
                    "rna_10_3000", "rna_20_3000", "rna_40_3000", "rna_100_3000",
                    "rna_10_5000", "rna_20_5000", "rna_40_5000", "rna_100_5000")

rna_counts_replicates <- list()

for(config in configurations) {
  for(i in 1:10) {
    # Generate new names for the replicates
    new_name <- paste(i, gsub("rna_", "", config), sep = "_")
    rna_counts_replicates[[new_name]] <- rna_counts_list[[config]] 
  }
}


all_auc_values <- list()

for (replicate in names(processed_results)) {
  
  # Extract the true labels for the current replicate
  true_labels <- rna_counts_replicates[[replicate]]@variable.annotations[["differential.expression"]]
  names(true_labels) <- rownames(rna_counts_replicates[[replicate]]@count.matrix)
  
  p_values <- list()
  
  for (method in c("naive_pydeseq", "aware_pydeseq", "naive_edge", "aware_edge")) {
    
    p <- as.data.frame(processed_results[[replicate]][[method]]$padj)
    rownames(p) <- rownames(processed_results[[replicate]][[method]])
    
    # Store the p-values in the p_values list
    p_values[[method]] <- p
  }
  
  common_genes <- Reduce(intersect, lapply(p_values, rownames))
  true_labels <- true_labels[common_genes]
  
  # Filter p-values for the common genes
  p_values_filtered <- lapply(p_values, function(p) p[common_genes, , drop = FALSE])
  
  # Combine p-values for all methods into a data frame
  p_values_df <- do.call(cbind, p_values_filtered)
  rownames(p_values_df) <- names(true_labels)
  colnames(p_values_df) <- c("PyDESeq2", "DeConveil", "EdgeR", "ABCD-DNA")
  
  # Calculate ROC and AUC for the current replicate
  roc_result <- auROC(
    truth = true_labels,
    p = p_values_df,
    sig = 0.05,
    x = "fpr",
    y = "tpr"
  )
  
  # Store the AUC values for the current replicate
  all_auc_values[[replicate]] <- roc_result$AUC_values
}


all_auc_values_by_config <- list()

for (config in c("10_1000", "20_1000", "40_1000", "100_1000",
                 "10_3000", "20_3000", "40_3000", "100_3000",
                 "10_5000", "20_5000", "40_5000", "100_5000")) {
  
  method_auc_values <- list()
  
  for (method in c("PyDESeq2", "DeConveil", "EdgeR", "ABCD-DNA")) {
    
    auc_values_for_method <- sapply(1:10, function(replicate_num) {
      replicate_name <- paste0(replicate_num, "_", config)  
      auc_list <- all_auc_values[[replicate_name]]  
      auc_list[[method]]  
    })
    
    method_auc_values[[method]] <- auc_values_for_method
  }

  all_auc_values_by_config[[config]] <- method_auc_values
}

# Create a summary data frame to store the mean and SD for each method and configuration
summary_auc_df <- data.frame(
  Configuration = character(),
  Method = character(),
  Mean_AUC = numeric(),
  SD_AUC = numeric(),
  stringsAsFactors = FALSE
)

# Calculate the mean and SD for each method within each configuration

for (config in names(all_auc_values_by_config)) {
  for (method in names(all_auc_values_by_config[[config]])) {
    
    auc_values <- all_auc_values_by_config[[config]][[method]]  
    
    # Check if auc_values is a list or not numeric
    if (is.list(auc_values)) {
      auc_values <- unlist(auc_values)  
    }
    auc_values <- as.numeric(auc_values)
    
    # Calculate mean and standard deviation
    if (all(!is.na(auc_values))) {  # Make sure there are no NA values before calculating
      mean_auc <- mean(auc_values)
      sd_auc <- sd(auc_values)
    } else {
      mean_auc <- NA
      sd_auc <- NA
    }
    
    # Append the results to the summary data frame
    summary_auc_df <- rbind(summary_auc_df, data.frame(
      Configuration = config,
      Method = method,
      Mean_AUC = mean_auc,
      SD_AUC = sd_auc
    ))
  }
}


summary_auc_df <- summary_auc_df %>%
  separate(Configuration, into = c("n_samples", "n_genes"), sep = "_") %>%
  dplyr::mutate(n_samples = as.numeric(n_samples), n_genes = as.numeric(n_genes))

summary_auc_df_filt <- summary_auc_df %>%
  dplyr::mutate(`Mean ± SD` = paste0(round(`Mean_AUC`, 3), " ± ", round(`SD_AUC`, 3))) %>% 
  select(n_samples, n_genes, Method, `Mean ± SD`)

reshaped_auc_df <- summary_auc_df_filt %>%
  pivot_wider(names_from = Method, values_from = `Mean ± SD`) %>%
  dplyr::rename(
    PyDESeq2 = `PyDESeq2`,
    DeConveil = `DeConveil`,
    EdgeR = `EdgeR`,
    `ABCD-DNA` = `ABCD-DNA`
  )

saveRDS(reshaped_auc_df, file = "CN-aware-DGE/simulations/auc_summary.RDS")



# Radar chart - plot AUC values across methods #

library(fmsb)
library(scales)  # For alpha transparency in colors

data_1000 <- summary_auc_df %>% 
  filter(n_genes == 1000) %>% 
  select(n_samples, Method, Mean_AUC)

data_1000 <- data_1000 %>%
  pivot_wider(names_from = Method, values_from = `Mean_AUC`) %>%
  dplyr::rename(
    PyDESeq2 = `PyDESeq2`,
    DeConveil = `DeConveil`,
    EdgeR = `EdgeR`,
    `ABCD-DNA` = `ABCD-DNA`
  )

data_3000 <- summary_auc_df %>% 
  filter(n_genes == 3000) %>% 
  select(n_samples, Method, Mean_AUC)

data_3000 <- data_3000 %>%
  pivot_wider(names_from = Method, values_from = `Mean_AUC`) %>%
  dplyr::rename(
    PyDESeq2 = `PyDESeq2`,
    DeConveil = `DeConveil`,
    EdgeR = `EdgeR`,
    `ABCD-DNA` = `ABCD-DNA`
  )

data_5000 <- summary_auc_df %>% 
  filter(n_genes == 5000) %>% 
  select(n_samples, Method, Mean_AUC)

data_5000 <- data_5000 %>%
  pivot_wider(names_from = Method, values_from = `Mean_AUC`) %>%
  dplyr::rename(
    PyDESeq2 = `PyDESeq2`,
    DeConveil = `DeConveil`,
    EdgeR = `EdgeR`,
    `ABCD-DNA` = `ABCD-DNA`
  )

prepare_radar_data <- function(data, min_val = 0.6, max_val = 1.0) {
  transposed <- as.data.frame(t(data[,-1]))  
  colnames(transposed) <- data$n_samples 
  
  # Create radar chart data with min and max limits
  radar_data <- as.data.frame(rbind(
    Max = rep(max_val, ncol(transposed)),  
    Min = rep(min_val, ncol(transposed)),
    transposed
  ))
  return(radar_data)
}

# Apply the function to all datasets
radar_data_1000 <- prepare_radar_data(data_1000)
radar_data_3000 <- prepare_radar_data(data_3000)
radar_data_5000 <- prepare_radar_data(data_5000)


method_colors <- c("#0000FF",  "#4500ACFF", "#00BA38",  "#EE3F3FFF")
method_labels <- c("ABCD-DNA", "EdgeR", "PyDESeq", "DeConveil")

create_radarchart <- function(data, color = method_colors, 
                              vlabels = colnames(data), vlcex = 1.6,
                              caxislabels = c(0.6, 0.7, 0.8, 0.9, 1.0), 
                              title = "AUC Plot") {
  radarchart(
    data, axistype = 1,
    pcol = color,                      
    pfcol = alpha(color, 0.1),         
    plwd = 4,                          
    plty = 1,                          
    cglcol = "black", cglty = 1,    
    cglwd = 1.0,                       
    axislabcol = "black",              
    vlcex = vlcex,                     
    vlabels = vlabels,                 
    caxislabels = caxislabels                           
  )
  title(main = title, font.main = 1, cex.main = 1.6)
}

# Sample size labels
vlabels <- c("n=10", "n=20", "n=40", "n=100")

# Plot for each radar dataset
create_radarchart(radar_data_1000, vlabels = vlabels, title = "AUC - 1000 genes")
create_radarchart(radar_data_3000, vlabels = vlabels, title = "AUC - 3000 genes")
create_radarchart(radar_data_5000, vlabels = vlabels, title = "AUC - 5000 genes")


legend("bottomright", legend = method_labels, col = method_colors, 
       lty = 1, lwd = 2, bty = "n", cex = 1.2, 
       title = "Methods")


#true_labels <- rna_counts_list[["rna_10_3000"]]@variable.annotations[["differential.expression"]]
#names(true_labels) <- rownames(rna_counts_list[["rna_10_3000"]]@count.matrix)

#p1 <- as.data.frame(processed_results[["2_10_3000"]][["naive_pydeseq"]]$padj)
#rownames(p1) <- rownames(processed_results[["2_10_3000"]][["naive_pydeseq"]])
#p2 <- as.data.frame(processed_results[["2_10_3000"]][["naive_edge"]]$padj)
#rownames(p2) <- rownames(processed_results[["2_10_3000"]][["naive_edge"]])
#p3 <- as.data.frame(processed_results[["2_10_3000"]][["aware_edge"]]$padj)
#rownames(p3) <- rownames(processed_results[["2_10_3000"]][["aware_edge"]])
#p4 <- as.data.frame(processed_results[["2_10_3000"]][["aware_pydeseq"]]$padj)
#rownames(p4) <- rownames(processed_results[["2_10_3000"]][["aware_pydeseq"]])

#common_genes <- Reduce(intersect, list(rownames(p1), rownames(p2), rownames(p3), rownames(p4)))
#true_labels <- true_labels[common_genes]

#p1 <- as.data.frame(p1[common_genes,]) 
#p2 <- as.data.frame(p2[common_genes,]) 
#p3 <- as.data.frame(p3[common_genes,]) 
#p4 <- as.data.frame(p4[common_genes,]) 

#p_values <- cbind(p1, p2, p3, p4)
#rownames(p_values) <- names(true_labels)
#colnames(p_values) <- c("PyDESeq2", "EdgeR", "ABCD-DNA", "DeConveil")

#roc <- auROC(
  #truth = true_labels, 
  #p = p_values, 
  #sig = 0.05, 
  #x = "fpr", 
  #y = "tpr")

