setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/")

pkgs <- c("tidyverse", "ggplot2", "edgeR", "VennDiagram", "pROC", "PRROC", "grid", "metaseqR2", "caret")
sapply(pkgs, require, character.only = TRUE)
source("deconveilCaseStudies/utils/utils.R")
source("deconveilCaseStudies/utils/utils_plot.R")

### Performance evaluation: Accuracy | Precision | Specificity ###

n_genes_list <- c(1000, 3000, 5000)
n_samples_list <- c(10, 20, 40, 100)

for (n_genes in n_genes_list) {
  for (n_samples in n_samples_list) {
    res <- evaluate_simulation_performance(n_samples = n_samples, n_genes = n_genes)
    file_path <- sprintf("deconveilCaseStudies/simulations/results/simulation_1/res_performance/res_%d_%d.RDS", n_samples, n_genes)
    saveRDS(res, file = file_path)
  }
}


## Load datasets ##
base_path <- "deconveilCaseStudies/simulations/results/simulation_1/res_performance"
sample_sizes <- c(10, 20, 40, 100)

load_results <- function(n_genes, sample_sizes, base_path) {
  files <- sprintf("%s/res_%d_%d.RDS", base_path, sample_sizes, n_genes)
  results <- lapply(files, readRDS)
  bind_rows(results)
}

data_1000 <- load_results(1000, sample_sizes, base_path)
data_3000 <- load_results(3000, sample_sizes, base_path)
data_5000 <- load_results(5000, sample_sizes, base_path)

data_1000 <- data_1000[data_1000$Method != "EdgeR-CN-aware", ]
data_3000 <- data_3000[data_3000$Method != "EdgeR-CN-aware", ]
data_5000 <- data_5000[data_5000$Method != "EdgeR-CN-aware", ]

rename_methods <- function(data) {
  data %>%
    mutate(Method = case_when(
      Method == "PyDESeq2-CN-naive" ~ "PyDESeq2",
      Method == "PyDESeq2-CN-aware" ~ "DeConveil",
      Method == "EdgeR-CN-naive" ~ "edgeR",
      TRUE ~ Method
    ))
}

data_1000 <- rename_methods(data_1000)
data_3000 <- rename_methods(data_3000)
data_5000 <- rename_methods(data_5000)

# Combine and summarize the data

summary_1000 <- summarize_performance(data_1000) %>% mutate(GeneSize = "1000 genes")
summary_3000 <- summarize_performance(data_3000) %>% mutate(GeneSize = "3000 genes")
summary_5000 <- summarize_performance(data_5000) %>% mutate(GeneSize = "5000 genes")

plot_df <- bind_rows(summary_1000, summary_3000, summary_5000)

# AUC #

auc_df <- readRDS("deconveilCaseStudies/simulations/auc_df_plot.RDS")
auc_df <- auc_df %>% 
  dplyr::rename(SampleSize = n_samples,
                            GeneSize = n_genes,
                            Method = Method,
                            Mean = Mean_AUC,
                            SD = SD_AUC) %>% 
  dplyr::mutate(Metric = "AUC") %>% 
  dplyr::select(Method, SampleSize, Metric, Mean, SD, GeneSize) %>% 
  dplyr::mutate(GeneSize = paste0(GeneSize, " genes"))

plot_df_final <- rbind(plot_df, auc_df)
plot_df_final$SampleSize <- as.factor(plot_df_final$SampleSize)

#write.xlsx(plot_df_final, file = "deconveilCaseStudies/results/simulation_performance/performance_simulation_metrics.xlsx")
saveRDS(plot_df_final, file = "deconveilCaseStudies/plots/main/Fig 2/rds/performance_simulation_metrics.rds")


# Performance metrics plot #

plot_df_final <- readRDS("deconveilCaseStudies/plots/main/Fig 2/rds/performance_simulation_metrics.rds")

method_colors <- c("DeConveil" = "#ED665D", "PyDESeq2" = "#67BF5C", "edgeR" = "#729ECE")
plot_df_final <- plot_df_final %>%
  mutate(Metric = factor(Metric, levels = c("Accuracy", "Precision", "Specificity", "AUC")))

line_plot <- performance_plot(plot_df_final, method_colors)
line_plot

ggsave("deconveilCaseStudies/plots/main/performance_plot.png", dpi = 500, width = 7.0, height = 6.0, plot = line_plot)    


## Load and process Results CN-aware & CN-naive methods ##

sample_sizes <- c(10, 20, 40, 100)
gene_counts <- c(1000, 3000, 5000)

load_rna_counts <- function(sample, gene) {
  file_path <- paste0("deconveilCaseStudies/simulations/results/rna_counts_sim/rna_counts_sim_", sample, "_", gene, "_brca.rds")
  readRDS(file_path)
}

load_results <- function(method, sample, gene) {
  file_dir <- if (grepl("edge", method, ignore.case = TRUE)) {
    "edgeR"
  } else {
    "pydeseq"
  }
  extension <- ifelse(file_dir == "edgeR", "RDS", "csv")
  clean_method <- sub("_edge", "", method)  
  file_name <- paste0("res_", clean_method, "_", sample, "_", gene, ".", extension)
  file_path <- paste0("deconveilCaseStudies/simulations/results/", file_dir, "/", file_name)
  if (file.exists(file_path)) {
    message(paste("Loading", method, "results from:", file_path))
    if (extension == "csv") {
      return(read.csv(file_path))
    } else {
      return(readRDS(file_path))
    }
  } else {
    warning(paste("Results file not found:", file_path))
    return(NULL)
  }
}

rna_counts_list <- list()
results_list <- list()

for (gene in gene_counts) {
  for (sample in sample_sizes) {
    
    # Load RNA Counts
    rna_key <- paste0("rna_", sample, "_", gene)
    rna_counts_list[[rna_key]] <- load_rna_counts(sample, gene)
    
    # Results (CN-naive & CN-aware)
    results_list[[paste0("pydeseq_naive_", sample, "_", gene)]] <- load_results("CNnaive", sample, gene)
    results_list[[paste0("deconveil_aware_", sample, "_", gene)]] <- load_results("CNaware", sample, gene)
    results_list[[paste0("edge_naive_", sample, "_", gene)]] <- load_results("CNnaive_edge", sample, gene)
  }
}


process_results_pydeseq <- function(df) {
  df %>%
    dplyr::select(X, log2FoldChange, padj) %>%
    remove_rownames() %>%
    column_to_rownames(var = "X") %>%
    dplyr::rename(logFC = log2FoldChange) %>%
    na.omit()
}

process_results_edge <- function(df) {
  df %>%
    dplyr::select(logFC, FDR) %>%
    dplyr::rename(padj = FDR) %>%
    na.omit()
}


res_pydeseq <- process_results_pydeseq(results_list[["pydeseq_naive_100_5000"]])
res_edge <- process_results_edge(results_list[["edge_naive_100_5000"]])
res_deconveil <- process_results_pydeseq(results_list[["deconveil_aware_100_5000"]])




