setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/deconveilCaseStudies/")

pkgs <- c("tidyverse", "ggplot2", "ggpubr", "reshape2", "MASS", "ggnewscale", "patchwork")
sapply(pkgs, require, character.only = TRUE)
source("utils/utils.R")
source("utils/utils_plot.R")

### Performance evaluation: CNV as a counfounder in DGE analysis ###

sample_sizes <- c(10, 20, 40, 60)
gene_counts <- c(1000)
replicates <- 1:10
methods <- c("deconveil", "pydeseq")

# Helper functions

compute_metrics <- function(true_vals, est_vals, method) {
  cor_val <- cor(true_vals, est_vals, method = "pearson", use = "complete.obs")
  mse <- mean((true_vals - est_vals)^2, na.rm = TRUE)
  return(data.frame(Method = method, Pearson_Correlation = cor_val, MSE = mse))
}

# LFC comparison against ground truth - Iterate to generate results

all_results <- list()

for (rep in replicates) {
  for (n in sample_sizes) {
    for (g in gene_counts) {
      
      # Load ground truth
      gt_path <- sprintf("simulations/results/simulation_2/replicates_rna_counts_sim/%d_rna_counts_sim_%d_%d_brca.rds", rep, n, g)
      ground_truth <- readRDS(gt_path)
      
      # Extract true values
      true_logfc <- ground_truth@variable.annotations[["truelog2foldchanges"]]
      de_status <- ground_truth@variable.annotations[["differential.expression"]]
      gene_ids <- rownames(ground_truth@count.matrix)
      
      # Initialize dataframe
      df <- data.frame(
        geneID = gene_ids,
        true_log2FC = true_logfc,
        DE_status = de_status
      )
      
      # Load method results
      for (method in methods) {
        file_prefix <- switch(method,
                              "deconveil" = "replicates_deconveil",
                              "pydeseq"   = "replicates_pydeseq"
        )
        
        suffix <- if (method == "pydeseq") "CNnaive" else "CNaware"
        
        file_name <- sprintf("simulations/results/simulation_2/%s/%d_res_%s_%d_%d.csv",
                             file_prefix, rep, suffix, n, g)
        
        if (file.exists(file_name)) {
          res <- read.csv(file_name)
          df[[paste0(method, "_log2FC")]] <- res$log2FoldChange
          
          # If `padj` exists, use it; otherwise, fall back to `pvalue`
          if ("padj" %in% names(res)) {
            df[[paste0(method, "_p")]] <- res$padj
          } else {
            df[[paste0(method, "_p")]] <- NA
          }
          
        } else {
          warning(sprintf("Missing file: %s", file_name))
          df[[paste0(method, "_log2FC")]] <- NA
          df[[paste0(method, "_p")]] <- NA
        }
      }
      
      # Compute metrics
      for (method in methods) {
        logfc_col <- paste0(method, "_log2FC")
        pval_col <- paste0(method, "_p")
        
        if (!all(is.na(df[[logfc_col]]))) {
          # Pearson & MSE
          metrics <- compute_metrics(df$true_log2FC, df[[logfc_col]], method)
          
          metrics$SampleSize <- n
          metrics$GeneCount <- g
          metrics$Replicate <- rep
          
          all_results[[length(all_results) + 1]] <- metrics
        }
      }
    }
  }
}

# Combine all into a single dataframe
results_df <- dplyr::bind_rows(all_results)

summary_df <- results_df %>%
  group_by(Method, SampleSize, GeneCount) %>%
  summarise(
    Pearson_mean = mean(Pearson_Correlation, na.rm = TRUE),
    Pearson_sd   = sd(Pearson_Correlation, na.rm = TRUE),
    MSE_mean = mean(MSE, na.rm = TRUE),
    MSE_sd   = sd(MSE, na.rm = TRUE),
    .groups = "drop"
  )

#write.xlsx(summary_df, file = "results/simulation_performance/performance_metrics_cnvCounfounder.xlsx")
saveRDS(summary_df, file = "plots/main/Fig 2/rds/performance_metrics_effectSize.rds")

# Plot 

method_colors <- c("deconveil" = "#23b4dc", "pydeseq" = "#0b2483")
method_border_colors <- c(
  "DeConveil" = "#23b4dc", 
  "PyDESeq2" = "#0b2483"
)

results_df$BorderColor <- method_border_colors[as.character(results_df$Method)]

p1 <- ggplot(summary_df,
             aes(x = SampleSize, y = Pearson_mean, color = Method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = Pearson_mean - Pearson_sd,
        ymax = Pearson_mean + Pearson_sd),
    width = 1,
    linewidth = 0.6
  ) +
  labs(
    title = "Pearson corr (log2FC)",
    x = "# samples",
    y = "Pearson correlation"
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2)
  ) +
  scale_color_manual(values = method_colors,
                     labels = c("DeConveil", "PyDESeq2")) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 14),   
    legend.text = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, color = "black"),
    axis.title.x = element_text(size = 14, color = "black"),                     
    axis.title.y = element_text(size = 14, color = "black"), 
    axis.text.x = element_text(size = 12),                      
    axis.text.y = element_text(size = 12),                      
  )

p1

p2 <- ggplot(summary_df,
             aes(x = SampleSize, y = MSE_mean, color = Method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = MSE_mean - MSE_sd,
        ymax = MSE_mean + MSE_sd),
    width = 1,
    linewidth = 0.6
  ) +
  labs(
    title = "MSE (log2FC)",
    x = "# samples",
    y = "MSE"
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2)
  ) +
  scale_color_manual(values = method_colors,
                     labels = c("DeConveil", "PyDESeq2")) +
  theme_bw() +
  theme(
    legend.position = "",
    legend.title = element_text(size = 14),   
    legend.text = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, color = "black"),
    axis.title.x = element_text(size = 14, color = "black"),                     
    axis.title.y = element_text(size = 14, color = "black"), 
    axis.text.x = element_text(size = 12),                      
    axis.text.y = element_text(size = 12),                      
  )
p2

p = p2 / p1
p

ggsave("plots/main/Fig 2/png/perf_plot_lfc.png", dpi = 500, width = 4.0, height = 8.0, plot = p)   

## Performance metrics ##

all_res <-  read.csv("simulations/sim_results/data_to_plot/performance_data_stageR_sim2.csv", header = T, sep = ",")

strongCN_res <- all_res %>% 
  dplyr::filter(cn_signal == "strong")

metrics_long <- strongCN_res %>%
  dplyr::select(
    sample_size,
    replicate,
    method,
    Precision,
    Recall,
    F1,
    MCC
  ) %>%
  pivot_longer(
    cols = c(Precision, Recall, F1, MCC),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = recode(
      metric,
      F1 = "F1-score"
    )
  )

metrics_summary <- metrics_long %>%
  group_by(sample_size, method, metric) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    se   = sd(value, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(
    lower = mean - 1.96 * se,
    upper = mean + 1.96 * se
  )

metrics_summary$metric <- factor(
  metrics_summary$metric,
  levels = c("Precision", "Recall", "F1-score", "MCC")
)

saveRDS(metrics_summary, file = "plots/main/Fig 2/rds/performance_metrics_de.RDS")

metrics_long <- metrics_long %>%
  mutate(
    metric = factor(metric, levels = c("Precision", "Recall", "F1-score", "MCC"))
  )

# Plot

method_col <- c(
  "DeConveil" = "#56B4E9", 
  "PyDESeq2" = "#0b2483"
)

cn_col <- c(
  "weak"   = "#44AA99",
  "strong" = "#CC6677"
)


p_bxp <- plot_performance_metrics(
  df = metrics_long, 
  colors = method_col, 
  facet_cols = "metric")

ggsave("plots/main/Fig 2/png/perf_boxp_method.png", dpi = 400, width = 6.0, height = 4.0, plot = p_bxp)

# Pairwise comparison 
all_res_pairwise <-  read.csv("simulations/sim_results/data_to_plot/performance_data_stageR_pairwise.csv", header = T, sep = ",")

metrics_long <- all_res_pairwise %>%
  dplyr::select(
    sample_size,
    cn_signal,
    replicate,
    comparison,
    Precision,
    Recall,
    F1.score,
    MCC
  ) %>%
  pivot_longer(
    cols = c(Precision, Recall, F1.score, MCC),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = recode(
      metric,
      F1.score = "F1-score"
    )
  )

metrics_summary <- metrics_long %>%
  group_by(sample_size, cn_signal, comparison, metric) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    se   = sd(value, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(
    lower = mean - 1.96 * se,
    upper = mean + 1.96 * se
  )

metrics_long <- metrics_long %>%
  mutate(
    comparison = gsub("_vs_", " vs ", comparison),
    comparison = factor(comparison, levels = c("DCG vs DSG", "DCG vs DIG")),
    `CN signal` = factor(cn_signal, levels = c("strong", "weak")),
    metric = factor(metric, levels = c("Precision", "Recall", "F1-score", "MCC"))
  )


p_bxp <- plot_performance_metrics(
  df = metrics_long, 
  colors = cn_col, 
  facet_cols = "metric",
  facet_rows = "comparison")
p_bxp

saveRDS(metrics_long, file = "plots/main/Fig 2/rds/performance_metrics_classif.rds")
ggsave("plots/main/Fig 2/png/perf_boxp_pairwise.png", dpi = 400, width = 6.0, height = 4.0, plot = p_bxp)
