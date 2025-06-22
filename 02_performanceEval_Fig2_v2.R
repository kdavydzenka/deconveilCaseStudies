setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/deconveilCaseStudies/")

pkgs <- c("tidyverse", "ggplot2", "ggpubr", "reshape2", "MASS")
sapply(pkgs, require, character.only = TRUE)
source("utils/utils.R")

### Performance evaluation: CNV as a counfounder in DGE analysis ###

sample_sizes <- c(10, 20, 40, 100)
gene_counts <- c(1000)
replicates <- 1:10
methods <- c("deconveil", "pydeseq")

compute_metrics <- function(true_vals, est_vals, method) {
  cor_val <- cor(true_vals, est_vals, method = "pearson", use = "complete.obs")
  mse <- mean((true_vals - est_vals)^2, na.rm = TRUE)
  return(data.frame(Method = method, Pearson_Correlation = cor_val, MSE = mse))
}

compute_mcc <- function(true_labels, predicted_pvals, threshold = 0.05) {
  if (all(is.na(predicted_pvals))) return(NA)
  
  predicted_labels <- ifelse(predicted_pvals < threshold, 1, 0)
  true_labels <- as.integer(true_labels)
  
  # Filter out NAs
  valid_idx <- complete.cases(predicted_labels, true_labels)
  predicted_labels <- predicted_labels[valid_idx]
  true_labels <- true_labels[valid_idx]
  
  TP <- sum(predicted_labels == 1 & true_labels == 1)
  TN <- sum(predicted_labels == 0 & true_labels == 0)
  FP <- sum(predicted_labels == 1 & true_labels == 0)
  FN <- sum(predicted_labels == 0 & true_labels == 1)
  
  numerator <- as.numeric((TP * TN) - (FP * FN))  # avoid integer overflow
  denom_parts <- as.numeric(c(TP + FP, TP + FN, TN + FP, TN + FN))
  
  # Avoid invalid division
  if (any(denom_parts == 0)) return(NA)
  
  denominator <- sqrt(prod(denom_parts))
  if (denominator == 0) return(NA)
  
  return(numerator / denominator)
}



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
          
          # MCC
          if (!all(is.na(df[[pval_col]]))) {
            mcc_val <- compute_mcc(df$DE_status, df[[pval_col]])
          } else {
            mcc_val <- NA
          }
          
          metrics$MCC <- mcc_val
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
    Pearson_Correlation_Mean = mean(Pearson_Correlation, na.rm = TRUE),
    Pearson_Correlation_SD   = sd(Pearson_Correlation, na.rm = TRUE),
    MSE_Mean = mean(MSE, na.rm = TRUE),
    MSE_SD   = sd(MSE, na.rm = TRUE),
    MCC_Mean = mean(MCC, na.rm = TRUE),
    MCC_SD   = sd(MCC, na.rm = TRUE),
    .groups = "drop"
  )

write.xlsx(summary_df, file = "results/simulation_performance/performance_metrics_cnvCounfounder.xlsx")
saveRDS(summary_df, file = "plots/main/Fig 2/rds/performance_metrics_cnvVounfounder.rds")


library(ggnewscale)
method_colors <- c("deconveil" = "#F66D7A", "pydeseq" = "#27AD81FF")
method_border_colors <- c("DeConveil" = "#990000", "PyDESeq2" = "#007439")
results_df$BorderColor <- method_border_colors[as.character(results_df$Method)]

# Boxplot: Pearson Correlation
p1 <- ggplot(results_df, aes(x = factor(SampleSize), y = Pearson_Correlation, color = Method)) +
  geom_boxplot(fill=NA) +
  labs(
    title = "Pearson corr (log2FC)",
    x = "# samples",
    y = "Pearson Correlation (R2)"
  ) +
  scale_y_continuous(breaks = seq(0.6, 1, by = 0.05)) +
  theme_bw() +
  scale_color_manual(values = method_colors, labels = c("DeConveil", "PyDESeq2")) +
  theme(
    legend.position = "",
    plot.title = element_text(hjust = 0.5, color = "black"),
    axis.title.x = element_text(size = 14, color = "black"),                     
    axis.title.y = element_text(size = 14, color = "black"), 
    axis.text.x = element_text(size = 12),                      
    axis.text.y = element_text(size = 12),                      
  )

# Boxplot: MSE
p2 <- ggplot(results_df, aes(x = factor(SampleSize), y = MSE, color = Method)) +
  geom_boxplot(fill=NA)+
  labs(
    title = "MSE (log2FC)",
    x = "# samples",
    y = "MSE"
  ) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
  theme_bw() +
  scale_color_manual(values = method_colors, labels = c("DeConveil", "PyDESeq2")) +
  theme(
    legend.position = "",
    legend.title = element_text(size = 14),   
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 14, hjust = 0.5, color = "black"),
    axis.title.x = element_text(size = 14, color = "black"),                     
    axis.title.y = element_text(size = 14, color = "black"), 
    axis.text.x = element_text(size = 12),                      
    axis.text.y = element_text(size = 12),                      
  )


# Boxplot: MCC
p3 <- ggplot(results_df, aes(x = factor(SampleSize), y = MCC, color = Method)) +
  geom_boxplot(fill=NA) +
  labs(
    title = "MCC (adjusted p-value)",
    x = "# samples",
    y = "MCC"
  ) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
  theme_bw() +
  scale_color_manual(values = method_colors, labels = c("DeConveil", "PyDESeq2")) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),  
    plot.title = element_text(size = 14, color = "black", hjust = 0.5),
    axis.title.x = element_text(size = 14, color = "black"),                     
    axis.title.y = element_text(size = 14, color = "black"), 
    axis.text.x = element_text(size = 12),                      
    axis.text.y = element_text(size = 12)                      
  )

ggsave("plots/main/performance_plot_p2.png", dpi = 500, width = 4.0, height = 4.0, plot = p2)   
ggsave("plots/main/performance_plot_p3.png", dpi = 500, width = 4.0, height = 4.0, plot = p3)   

#ground_truth <- readRDS("simulations/results/simulation_2/replicates_rna_counts_sim/1_rna_counts_sim_10_1000_brca.rds")
#ground_truth_logfc <- ground_truth@variable.annotations[["truelog2foldchanges"]]
#de_status <- ground_truth@variable.annotations[["differential.expression"]]

#deconveil_res <- read.csv("simulations/results/simulation_2/replicates_deconveil/1_res_CNaware_20_2000.csv")
#pydeseq_res <- read.csv("simulations/results/simulation_2/replicates_pydeseq/1_res_CNnaive_20_2000.csv")

# Merge results into a single dataframe
#df <- data.frame(
  #geneID = rownames(ground_truth@count.matrix),
  #true_log2FC = ground_truth_logfc,
  #DE_status = de_status,
  #deconveil_log2FC = deconveil_res$log2FoldChange,
  #pydeseq_log2FC = pydeseq_res$log2FoldChange,
  #deconveil_p = deconveil_res$padj,
  #pydeseq_p = pydeseq_res$padj
#)

#metrics <- rbind(
  #compute_metrics("deconveil_log2FC"),
  #compute_metrics("pydeseq_log2FC"),
  #compute_metrics("edge_log2FC"),
  #compute_metrics("deseq_log2FC")
#)


# Scatter plot: Estimated vs. True Log2FC
#scatter_plot <- function(method_name, estimated_log2FC_column) {
  #ggplot(df, aes(x = true_log2FC, y = .data[[estimated_log2FC_column]])) +
    #geom_point(alpha = 0.3, color = "darkgray") +  # Light gray points for background
    #geom_density2d(aes(color = ..level..), size = 1) +  # Contour lines for density
    #scale_color_gradient(low = "blue", high = "red") +  # Color gradient
    #geom_abline(slope = 1, intercept = 0, color = "red", linetype = "solid") +  # Identity line
    #theme_minimal() +
    #labs(title = paste("Log2FC Comparison -", method_name),
         #x = "True Log2FC",
         #y = paste(method_name, "Estimated Log2FC")) +
    #annotate("text", x = -3, y = 3, 
             #label = paste0("R = ", round(cor(df$true_log2FC, df[[estimated_log2FC_column]], method = "pearson"), 2)), 
             #size = 5)
#}


#library(fields)
#p1 <- scatter_plot("DeConveil", "deconveil_log2FC")
#p2 <- scatter_plot("PyDESeq2", "pydeseq_log2FC")
#p3 <- scatter_plot("edgeR", "edge_log2FC")
#p4 <- scatter_plot("deseq2", "deseq_log2FC")








