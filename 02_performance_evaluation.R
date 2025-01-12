setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/")

pkgs <- c("tidyverse", "ggplot2", "edgeR", "VennDiagram", "pROC", "PRROC", "grid", "metaseqR2", "caret")
sapply(pkgs, require, character.only = TRUE)
source("deconveilCaseStudies/utils.R")

### Performance evaluation: Accuracy | Precision | Specificity ###

n_genes_list <- c(1000, 3000, 5000)
n_samples_list <- c(10, 20, 40, 100)

for (n_genes in n_genes_list) {
  for (n_samples in n_samples_list) {
    res <- evaluate_simulation_performance(n_samples = n_samples, n_genes = n_genes)
    file_path <- sprintf("deconveilCaseStudies/simulations/results/res_performance/res_%d_%d.RDS", n_samples, n_genes)
    saveRDS(res, file = file_path)
  }
}

load_results <- function(n_genes, sample_sizes, base_path) {
  files <- sprintf("%s/res_%d_%d.RDS", base_path, sample_sizes, n_genes)
  results <- lapply(files, readRDS)
  bind_rows(results)
}

# Load datasets
base_path <- "deconveilCaseStudies/simulations/results/res_performance"
sample_sizes <- c(10, 20, 40, 100)

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
data_5000 <- rename_methods(data_5000)

# Combine and summarize the data

summary_1000 <- summarize_performance(data_1000) %>% mutate(GeneSize = "1000 genes")
summary_5000 <- summarize_performance(data_5000) %>% mutate(GeneSize = "5000 genes")

plot_df <- bind_rows(summary_1000, summary_5000)
plot_df$SampleSize <- as.factor(plot_df$SampleSize)


# Performance metrics plot #

method_colors <- c("DeConveil" = "#ED665D", "PyDESeq2" = "#67BF5C", "EdgeR" = "#729ECE")

performance_plot <- ggplot(plot_df, aes(x = SampleSize, y = Mean, color = Method, group = Method)) +
  geom_line(size = 1.4) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.15, size = 0.7) +
  ggh4x::facet_nested(factor(GeneSize) ~ factor(Metric)) +
  #facet_wrap(~ Metric, scales = "free_y") +
  labs(title = "", x = "sample size", y = "Performance metric") +
  scale_x_discrete(labels = c("10", "20", "40", "100")) +   
  scale_y_continuous(limits = c(0.7, 1), breaks = seq(0, 1, by = 0.1)) + 
  scale_color_manual(values = method_colors) +
  theme_bw()+
  theme(
    strip.text = element_text(size = 18, face = "plain", color = "black"),       
    axis.title.x = element_text(size = 16, color = "black"),                     
    axis.title.y = element_text(size = 16, color = "black"),                     
    axis.text.x = element_text(size = 14, color = "black"),                      
    axis.text.y = element_text(size = 14, color = "black"),                      
    legend.text = element_text(size = 13, color = "black"),                      
    legend.title = element_text(size = 15, color = "black"),
    legend.position = "right"
  )
performance_plot

ggsave("deconveilCaseStudies/plots/main/performance_plot.png", dpi = 400, width = 10.0, height = 6.0, plot = performance_plot)    


# Load and process Results CN-aware & CN-naive methods

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

#rownames_idx <- match(rownames(res_deconveil, rownames(res_edge)))
#res_edge <- res_edge[rownames_idx,] %>% na.omit()



## AUROC calculation | Supplementary  ##
true_labels <- rna_counts_list[["rna_100_5000"]]@variable.annotations[["differential.expression"]]
names(true_labels) <- rownames(rna_counts_list[["rna_100_5000"]]@count.matrix)

p1 <- as.data.frame(res_pydeseq$padj)
rownames(p1) <- rownames(res_pydeseq)
p2 <- as.data.frame(res_edge$padj)
rownames(p2) <- rownames(res_edge)
p3 <- as.data.frame(res_deconveil$padj)
rownames(p3) <- rownames(res_deconveil)

common_genes <- Reduce(intersect, list(rownames(p1), rownames(p2), rownames(p3)))
true_labels <- true_labels[common_genes]

p1 <- as.data.frame(p1[common_genes,]) 
p2 <- as.data.frame(p2[common_genes,]) 
p3 <- as.data.frame(p3[common_genes,]) 

p_values <- cbind(p1, p2, p3)
rownames(p_values) <- names(true_labels)
colnames(p_values) <- c("PyDESeq2", "edgeR", "DeConveil")


# AUROC plot

roc <- diagplotRoc(
  truth = true_labels, 
  p = p_values, 
  sig = 0.05, 
  x = "fpr", 
  y = "tpr", 
  output = "file", 
  line_colors = c("#00BA38", "#0000FF",  "#EE3F3FFF"),
  line_width = 7,
  plot_title = "Sample size: 100; Genes: 5000",
  axis_text_size = 2.0,
  legend_text_size = 1.6,
  font.main = 1,
  title_text_size = 2.4, 
  margin = c(6, 6, 6, 5),
  path = "deconveilCaseStudies/plots/supplementary/roc_100_5000.png"
)
roc


