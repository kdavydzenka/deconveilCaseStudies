setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/")
pkgs <- c("cluster", "factoextra", "heatmaply", "DESeq2", "tidyverse", "colorspace", 
          "ggpubr", "ggpointdensity", "ggeasy", "gridExtra", "metaseqR2", "ggalluvial", "ggridges", "ggforce", "ggparallel", "alluvial")
sapply(pkgs, require, character.only = TRUE)

source("deconveilCaseStudies/utils.R")

### Data preprocessing ###

# Input data
tumor_type <- c("LUAD")
tumor_type <- c("LUSC", "BRCA", "LIHC", "KIRC") # For Supplementary
base_dir <- "TCGA"

cnv_mean_all <- data.frame()

process_tumor_data <- function(tumor_type) {
  
  cnv_tumor_path <- file.path(base_dir, tumor_type, "cnv_tumor.RDS")
  rna_normal_path <- file.path(base_dir, tumor_type, "rna_normal.RDS")
  rna_tumor_path <- file.path(base_dir, tumor_type, "rna_tumor.RDS")
  
  cnv_tumor <- readRDS(cnv_tumor_path)
  rna_normal <- readRDS(rna_normal_path)
  rna_tumor <- readRDS(rna_tumor_path)
  
  # Gene filtering
  low_expression_threshold <- 10
  expression_summary <- data.frame(
    Gene = rownames(rna_normal),
    MeanExpression = rowMeans(rna_normal)
  )
  filtered_genes <- expression_summary %>%
    filter(MeanExpression > low_expression_threshold)
  
  rna_normal <- rna_normal[filtered_genes$Gene, ]
  rna_tumor <- rna_tumor[filtered_genes$Gene, ]
  rna <- cbind(rna_normal, rna_tumor)
  cnv_tumor <- cnv_tumor[filtered_genes$Gene, ]
  
  # Process CNV tumor data
  cnv_tumor <- apply(cnv_tumor, 2, function(x) ifelse(x > 15, 15, x))
  cnv_mean <- as.data.frame(cnv_tumor) %>%
    dplyr::mutate(cnv_mean = rowMeans(cnv_tumor)) %>%
    dplyr::mutate(tumor_type = tumor_type) %>%
    dplyr::select(cnv_mean, tumor_type)
  
  cnv_normal <- matrix(2, nrow = nrow(rna), ncol = ncol(rna_normal))
  rownames(cnv_normal) <- rownames(rna)
  cnv_normal <- as.data.frame(cnv_normal)
  cnv_tumor <- as.data.frame(cnv_tumor)
  
  cnv <- merge(cnv_normal, cnv_tumor, by = "row.names")
  cnv <- cnv %>% remove_rownames() %>% column_to_rownames("Row.names")
  colnames(cnv) <- colnames(rna)
  cnv <- cnv / 2
  
  colnames(rna) <- paste0("sample", 1:(ncol(rna)))
  colnames(cnv) <- colnames(rna)
  rownames_idx <- match(rownames(cnv), rownames(rna))
  rna <- rna[rownames_idx, ] %>% na.omit()
  
  metadata <- data.frame(
    patID = colnames(rna),
    condition = rep(c("A", "B"), each = ncol(rna_tumor))
  )
  metadata <- metadata %>% remove_rownames() %>% column_to_rownames(var = "patID")
  metadata$condition <- as.factor(metadata$condition)
  
  # Save processed data
  output_dir <- file.path(base_dir, tumor_type, "test")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  write.csv(cnv, file = file.path(output_dir, "cnv.csv"), row.names = TRUE)
  write.csv(rna, file = file.path(output_dir, "rna.csv"), row.names = TRUE)
  write.csv(metadata, file = file.path(output_dir, "metadata.csv"), row.names = TRUE)
  
  return(cnv_mean)
}

for (tumor_type in tumor_type) {
  cnv_mean_all <- rbind(cnv_mean_all, process_tumor_data(tumor_type))
}


hist <- ggplot(cnv_mean_all, aes(x = cnv_mean, fill = tumor_type)) +
  geom_histogram(binwidth = 0.4, color = "black", alpha = 0.7, position = "identity", fill = "#F39B7FB2") +
  labs(
    title = "",
    x = "CN state",
    y = "Frequency"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 18, color = "black", hjust = 1),
    axis.text.y = element_text(size = 18, color = "black"),
    axis.title.x = element_text(size = 18, face = "plain", color = "black"),
    axis.title.y = element_text(size = 18, face = "plain", color = "black"),
    strip.text = element_text(size = 18, face = "plain", color = "black")
  ) +
  facet_wrap(~factor(tumor_type, levels = c("LUSC", "BRCA", "LIHC", "KIRC")), nrow = 1)
hist

ggsave("deconveilCaseStudies/plots/main/hist_luad.png", dpi = 400, width = 4.0, height = 3.5, plot = hist)
ggsave("deconveilCaseStudies/plots/supplementary/hist_pancancer.png", dpi = 400, width = 10.0, height = 3.5, plot = hist)


### Downstream analysis ###

# Read data files for a tumor type
read_data <- function(tumor_type) {
  list(
    res_pydeseq = read.csv(paste0("deconveilCaseStudies/results/", tumor_type, "/res_CNnaive.csv")),
    res_deconveil = read.csv(paste0("deconveilCaseStudies/results/", tumor_type, "/res_CNaware.csv")),
    cnv_tumor = read.csv(paste0("deconveilCaseStudies/results/", tumor_type, "/cnv_tumor.csv")) %>%
      remove_rownames() %>%
      column_to_rownames(var = "X") * 2
  )
}

# Main analysis pipeline #

tumor_type <- c("LUAD") # Main
tumor_type <- c("LUSC", "BRCA", "LIHC", "KIRC") # Supplementary
lfc_cut = 1.0
pval_cut = 0.05 
loss_threshold = 0.25

results_list <- lapply(tumor_types, function(tumor_type) {
  
  data <- read_data(tumor_type)
  cnv_data <- process_cnv_data(data$cnv_tumor, loss_threshold)
  
  res_naive <- annotate_results(data$res_pydeseq, lfc_cut, pval_cut, "CN naive", tumor_type)
  res_aware <- annotate_results(data$res_deconveil, lfc_cut, pval_cut, "CN aware", tumor_type)
  cnv_mean <- data$cnv_tumor
  loss_labels <- cnv_data$loss_labels
  
  combined <- combine_results(res_naive, res_aware, cnv_mean, loss_labels)
  
  combined[["res_naive"]] <- combined[["res_naive"]] %>% dplyr::rename(logFC = log2FoldChange)
  combined[["res_aware"]] <- combined[["res_aware"]] %>% dplyr::rename(logFC = log2FoldChange)
  
  # Separate gene groups
  res_joint <- combined$res_naive %>%
    inner_join(combined$res_aware, by = "geneID", suffix = c("_naive", "_aware"))
  
  gene_groups <- list(
    d_sensitive = res_joint %>%
      filter(DEtype_naive == "Up-reg" & DEtype_aware == "n.s." | 
               DEtype_naive == "Down-reg" & DEtype_aware == "n.s."),
    d_insensitive = res_joint %>%
      filter(DEtype_naive == "Down-reg" & DEtype_aware == "Down-reg" |
               DEtype_naive == "Up-reg" & DEtype_aware == "Up-reg"),
    d_compensated = res_joint %>%
      filter(DEtype_naive == "n.s." & DEtype_aware == "Down-reg" | 
               DEtype_naive == "n.s." & DEtype_aware == "Up-reg"),
    non_deg = res_joint %>%
      filter(DEtype_naive == "n.s." & DEtype_aware == "n.s.")
  )
  
  cn_naive <- bind_rows(
    gene_groups$d_sensitive %>% mutate(gene_group = "DSGs"),
    gene_groups$d_insensitive %>% mutate(gene_group = "DIGs"),
    gene_groups$d_compensated %>% mutate(gene_group = "DCGs")
  ) %>%
    select(geneID, logFC_naive, padj_naive, isDE_naive, DEtype_naive, tumor_type_naive, method_naive, cnv_mean_naive, gene_group) %>% 
    rename_with(~ c("geneID", "log2FC", "padj", "isDE", "DEtype", "tumor_type", "method", "cnv_mean", "gene_group"), everything())
  
  cn_aware <- bind_rows(
    gene_groups$d_sensitive %>% mutate(gene_group = "DSGs"),
    gene_groups$d_insensitive %>% mutate(gene_group = "DIGs"),
    gene_groups$d_compensated %>% mutate(gene_group = "DCGs")
  ) %>%
    select(geneID, logFC_aware, padj_aware, isDE_aware, DEtype_aware, tumor_type_aware, method_aware, cnv_mean_aware, gene_group) %>% 
    rename_with(~ c("geneID", "log2FC", "padj", "isDE", "DEtype", "tumor_type", "method", "cnv_mean", "gene_group"), everything())
  
  v_plot_data <- bind_rows(cn_naive, cn_aware) %>%
    filter(abs(log2FC) < case_when(
      tumor_type == "LUAD" ~ 6.0,
      tumor_type == "LUSC" ~ 7.0,
      tumor_type == "BRCA" ~ 7.0,
      tumor_type == "LIHC" ~ 5.0,
      tumor_type == "KIRC" ~ 7.0
    ))
  
  if (tumor_type == "BRCA") {
    v_plot_data <- v_plot_data %>% filter(padj > 2.84e-169)
  } else if (tumor_type == "LIHC") {
    v_plot_data <- v_plot_data %>% filter(padj > 4.949e-74)
  }
  
  list(
    tumor_type = tumor_type,
    v_plot_data = v_plot_data,
    cnv_data = cnv_data,
    cn_naive = cn_naive,
    cn_aware = cn_aware
  )
})

# Combine all results
v_plot_data_all <- bind_rows(lapply(results_list, `[[`, "v_plot_data"))
cn_naive_all <- bind_rows(lapply(results_list, `[[`, "cn_naive"))
cn_aware_all <- bind_rows(lapply(results_list, `[[`, "cn_aware"))

v_plot_data_all$gene_group <- as.factor(v_plot_data_all$gene_group)


# Generate Vocano plot

gene_group_colors <- c("DIGs" = "#8F3931FF", "DSGs" = "#FFB977", "DCGs"="#FAE48BFF", "non-DEGs" = "#ADB6B6FF")  
cnv_colors <- c("loss" = "#0073C299", "neutral" = "#86868699", "gain" = "#cecb76", "amplification" = "#DC0000B2")

p_volcanos <- v_plot_data_all %>%
  ggplot(mapping = aes(x = log2FC, y = -log10(padj))) +
  geom_point(data = subset(v_plot_data_all, gene_group %in% c("DIGs", "non-DEGs")),
             aes(col = gene_group), size = 1.0, alpha = 0.3) +
  geom_point(data = subset(v_plot_data_all, gene_group %in% c("DSGs", "DCGs")),
             aes(col = gene_group), size = 2.0, alpha = 0.5) +
  scale_color_manual(values = gene_group_colors) +
  theme_bw() +
  scale_x_continuous(breaks = seq(floor(min(v_plot_data_all$log2FC)), 
                                  ceiling(max(v_plot_data_all$log2FC)), by = 2)) +
  #ggh4x::facet_nested(factor(method, levels = c("CN naive", "CN aware"))~factor(tumor_type, levels = c("LUSC", "BRCA", "LIHC", "KIRC")), scales ="free", independent = "y")+
  facet_wrap(~factor(method, levels = c("CN naive", "CN aware")), nrow = 1) +
  labs(x = expression(Log[2] ~ FC), y = expression(-log[10] ~ Pvalue), col = "Gene group") +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = 'dashed') +
  geom_hline(yintercept = -log10(pval_cut), linetype = "dashed") +
  ggplot2::theme(legend.position = 'left',
                 legend.text = element_text(size = 16, color = "black"),
                 legend.title = element_text(size = 18, color = "black"),  
                 strip.text = element_text(size = 20, face = "plain", color = "black"),
                 axis.text = element_text(size = 18, color = "black"),
                 axis.title = element_text(size = 20, color = "black"))+
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1)))
p_volcanos

ggsave("deconveilCaseStudies/plots/main/volcano_luad.png", dpi = 400, width = 9.0, height = 4.0, plot = p_volcanos)
ggsave("deconveilCaseStudies/plots/supplementary/volcano_pancancer.png", dpi = 400, width = 14.0, height = 6.0, plot = p_volcanos)


## CN barplot ##

combined_data <- bind_rows(cn_aware_all)

combined_data <- combined_data %>% 
  mutate(cnv_group = case_when(
    cnv_mean > 0.5 & cnv_mean <= 1.7  ~ "loss",
    cnv_mean > 1.7 & cnv_mean <= 2.5  ~ "neutral",
    cnv_mean > 2.5 & cnv_mean <= 4.0  ~ "gain",
    cnv_mean > 4.0 ~ "amplification"
  ))

loss_labels <- bind_rows(
  lapply(results_list, function(res) {
    res$cnv_data$loss_labels %>%
      rownames_to_column(var = "geneID") %>%
      mutate(tumor_type = res$tumor_type)  
  })
)

combined_data <- merge(combined_data, loss_labels, by = c("geneID", "tumor_type"), all.x = TRUE)
combined_data$cnv_group <- ifelse(combined_data$isCNloss == "loss", "loss", combined_data$cnv_group)

barplot_data <- combined_data %>%
  group_by(tumor_type, gene_group) %>%
  summarise(Count = n()) %>%
  mutate(total = sum(Count)) %>%
  mutate(percentage = (Count / total) * 100) %>%
  ungroup()

combined_data$gene_group <- factor(combined_data$gene_group, levels = c("DIGs", "DSGs", "DCGs"))
combined_data$tumor_type <- factor(combined_data$tumor_type, levels = c("LUSC", "BRCA", "LIHC", "KIRC"))


barplot_cnv <- ggplot2::ggplot(combined_data, aes(x = gene_group, fill = cnv_group)) +
  geom_bar(position = "stack", width = 0.6) + 
  #geom_text(aes(label = Count), position = position_stack(vjust = 0.5), size = 4) +  
  scale_fill_manual(values = cnv_colors) +  
  theme_classic2() +  
  facet_wrap(~tumor_type, scales = "free", nrow = 1)+
  labs(y = "gene counts", x = "", title = "", fill = "CN group") +  
  theme(
    axis.text.x = element_text(size = 18, color = "black"),  
    axis.text.y = element_text(size = 16, color = "black"),                         
    axis.title.x = element_text(size = 18, face = "plain", color = "black"),          
    axis.title.y = element_text(size = 18, face = "plain", color = "black"),
    strip.text = element_text(size = 16, face = "plain", color = "black"),
    legend.position = 'bottom',
    legend.text = element_text(size = 14, color = "black"),                          
    legend.title = element_text(size = 16, face = "plain", color = "black")           
  )
barplot_cnv


barplot_stat <- ggplot2::ggplot(barplot_data, aes(x = gene_group, y = percentage, fill = gene_group)) +
  geom_bar(stat = "identity", color = "black", width = 0.6, alpha = 0.7) +  
  geom_text(aes(label = paste0(Count, " (", round(percentage, 1), "%)")), 
            position = position_stack(vjust = 0.5), 
            size = 3, color = "black", fontface = "plain") +
  scale_fill_manual(values = gene_group_colors) +  
  theme_classic() +  
  labs(y = "Percentage of Genes", x = "", title = "", fill = "Gene Group") +
  facet_wrap(~tumor_type, scales = "free", nrow = 1) +
  theme(
    axis.text.x = element_text(size = 16, color = "black", angle = 45, hjust = 1),  
    axis.text.y = element_text(size = 16, color = "black"),                         
    axis.title.x = element_text(size = 16, face = "plain", color = "black"),          
    axis.title.y = element_text(size = 16, face = "plain", color = "black"),          
    legend.position = 'right',
    legend.text = element_text(size = 14, color = "black"),                          
    legend.title = element_text(size = 16, face = "plain", color = "black")
  )
barplot_stat

ggsave("deconveilCaseStudies/plots/main/barplot_cnv_luad.png", dpi = 400, width = 6.0, height = 6.0, plot = barplot_cnv)
ggsave("deconveilCaseStudies/plots/main/barplot_stats_brca.png", dpi = 400, width = 6.0, height = 5.0, plot = barplot_stat)
ggsave("deconveilCaseStudies/plots/supplementary/barplot_cnv_pancancer.png", dpi = 400, width = 14.0, height = 4.5, plot = barplot_cnv)


# Calculate how many DEGs affected by CN gain & amplification

deg_counts <- combined_data %>%
  filter(cnv_group %in% c("gain", "amplification")) %>%  
  count(gene_group, cnv_group, tumor_type, name = "Count") %>%  
  group_by(gene_group, cnv_group, tumor_type) %>%  
  summarise(deg_count = sum(Count), .groups = "drop")

deg_percentages <- deg_counts %>%
  filter(cnv_group %in% c("gain", "amplification")) %>%
  group_by(gene_group, tumor_type) %>%  
  summarise(deg_count_sum = sum(deg_count), .groups = "drop") %>%  # 
  group_by(tumor_type) %>%  
  mutate(total_degs_tumor = sum(deg_count_sum)) %>%  
  mutate(percentage = (deg_count_sum / total_degs_tumor) * 100) %>% 
  ungroup()  


# Scatter - comparison LFC | p-value

names(results_list) <- sapply(results_list, function(x) x[["tumor_type"]])

plot_lfc_combined_list <- list()

for (tumor_type in names(results_list)) {
  
  cn_naive <- results_list[[tumor_type]]$cn_naive
  cn_aware <- results_list[[tumor_type]]$cn_aware
  
  loss_labels <- results_list[[tumor_type]][["cnv_data"]][["loss_labels"]]
  
  lfc_naive <- cn_naive %>% dplyr::select(geneID, log2FC) %>% dplyr::rename(log2FC_naive = log2FC)
  lfc_aware <- cn_aware %>% dplyr::select(geneID, log2FC, cnv_mean) %>% dplyr::rename(log2FC_aware = log2FC)
  
  plot_lfc <- left_join(lfc_naive, lfc_aware)
  
  plot_lfc <- plot_lfc %>%
    dplyr::mutate(cnv_group = case_when(
      cnv_mean > 0.5 & cnv_mean <= 1.7  ~ "loss",
      cnv_mean > 1.7 & cnv_mean <= 2.5  ~ "neutral",
      cnv_mean > 2.5 & cnv_mean <= 4.0 ~ "gain",
      cnv_mean > 4.0 ~ "amplification"))
  
  plot_lfc <- plot_lfc %>% remove_rownames() %>% column_to_rownames(var = "geneID")
  plot_lfc <- merge(plot_lfc, loss_labels, by = "row.names")
  plot_lfc$cnv_group <- ifelse(plot_lfc$isCNloss == "loss", "loss", plot_lfc$cnv_group)
  
  # Calculate the effect size difference and add tumor type
  plot_lfc <- plot_lfc %>%
    dplyr::mutate(eff_size_diff = abs(log2FC_naive - log2FC_aware)) %>%
    dplyr::mutate(tumor_type = tumor_type)
  
  # Add the processed data for this tumor type to the combined list
  plot_lfc_combined_list[[tumor_type]] <- plot_lfc
}

plot_lfc_combined <- bind_rows(plot_lfc_combined_list)

plot_lfc_combined <- plot_lfc_combined %>% 
  dplyr::filter(abs(log2FC_aware) < 7.0)

# Scatter plot
comparison_lfc <- ggplot(plot_lfc_combined, aes(x = log2FC_aware, y = log2FC_naive, color = cnv_group)) + 
  geom_point(shape = 20, size = 3) +
  geom_abline() +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black", size = 0.5) +
  xlab("Effect size (log2) CN-aware") +
  ylab("Effect size (log2) CN-naive") +
  scale_x_continuous(breaks = seq(-10, 10, by = 4)) +
  scale_y_continuous(breaks = seq(-10, 10, by = 4)) +
  ggh4x::facet_nested(factor(tumor_type, levels = c("LUAD", "LUSC", "BRCA", "LIHC", "KIRC")) ~ 
                        factor(cnv_group, levels = c("loss", "neutral", "gain", "amplification"))) +
  scale_color_manual(name = "CN group", values = cnv_colors) +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.title.x = element_text(size = 16, color = "black"),  
    axis.title.y = element_text(size = 16, color = "black"),
    axis.text.x = element_text(size = 16, color = "black"),    
    axis.text.y = element_text(size = 16, color = "black"),
    legend.text = element_text(size = 16, color = "black"),
    strip.text = element_text(size = 16, face = "plain", color = "black"),
    legend.title = element_text(size = 18, face = "plain", color = "black")
  )

comparison_lfc

ggsave("deconveilCaseStudies/plots/main/scatter_lfc_luad.png", dpi = 400, width = 8.0, height = 3.5, plot = comparison_lfc)
ggsave("deconveilCaseStudies/plots/supplementary/scatter_lfc_pancancer.png", dpi = 400, width = 9.0, height = 9.0, plot = comparison_lfc)


# Effect size difference 

plot_lfc <- plot_lfc_combined %>% dplyr::filter(eff_size_diff < 3.0 ,)
plot_lfc$cnv_group <- factor(plot_lfc$cnv_group, levels = c("loss", "neutral", "gain", "amplification"))

violin <- ggplot(plot_lfc, aes(x = cnv_group, y = eff_size_diff, fill = cnv_group)) + 
  geom_violin(trim = FALSE, scale = "width", alpha = 0.7) + 
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.5) + 
  scale_fill_manual(values = cnv_colors) + 
  xlab("CNV group") + 
  ylab("Effect size difference (log2)") + 
  theme_classic2() +
  facet_wrap(~factor(tumor_type, levels = c("LUAD", "LUSC", "BRCA", "LIHC", "KIRC")), nrow = 4)+
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 18, color = "black"),
    axis.title.y = element_text(size = 18, color = "black"),
    axis.text.x = element_text(size = 18, angle = 30, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 18, color = "black"),
    strip.text = element_text(size = 16, face = "plain", color = "black")
  )
violin

ggsave("deconveilCaseStudies/plots/main/violin_luad.png", dpi = 400, width = 5.0, height = 5.0, plot = violin)
ggsave("deconveilCaseStudies/plots/supplementary/violin_pancancer.png", dpi = 400, width = 3.5, height = 14.0, plot = violin)

# Calculate median and SD
lfc_stats <- plot_lfc %>% 
  dplyr::group_by(tumor_type, cnv_group) %>% 
  dplyr::summarise(
    median_lfc_diff = median(eff_size_diff, na.rm = TRUE),
    IQR_lfc_diff = IQR(eff_size_diff, na.rm = TRUE)
  ) %>% 
  dplyr::mutate(summary = sprintf("%.2f ± %.2f", median_lfc_diff, IQR_lfc_diff))


# p-value

plot_pval_combined_list <- list()

for (tumor_type in names(results_list)) {
  
  cn_naive <- results_list[[tumor_type]]$cn_naive
  cn_aware <- results_list[[tumor_type]]$cn_aware
  
  loss_labels <- results_list[[tumor_type]][["cnv_data"]][["loss_labels"]]
  
  pval_naive <- cn_naive %>% 
    dplyr::select(geneID, padj) %>% 
    dplyr::rename(padj_naive = padj)
  
  pval_aware <- cn_aware %>% 
    dplyr::select(geneID, padj, cnv_mean) %>% 
    dplyr::rename(padj_aware = padj)
  
  plot_pval <- left_join(pval_naive, pval_aware, by = "geneID")
  
  plot_pval <- plot_pval %>%
    dplyr::mutate(cnv_group = case_when(
      cnv_mean > 0.0 & cnv_mean <= 1.7  ~ "loss",
      cnv_mean > 1.7 & cnv_mean <= 2.5  ~ "neutral",
      cnv_mean > 2.5 & cnv_mean <= 4.0 ~ "gain",
      cnv_mean > 4.0 ~ "amplification"
    ))
  
  plot_pval <- plot_pval %>% remove_rownames() %>% column_to_rownames(var = "geneID")
  plot_pval<- merge(plot_pval, loss_labels, by = "row.names")
  plot_pval$cnv_group <- ifelse(plot_pval$isCNloss == "loss", "loss", plot_pval$cnv_group)
  
  # Filter p-values (apply different thresholds for each tumor type)
  if (tumor_type == "LUAD") {
    plot_pval <- plot_pval %>% dplyr::filter(padj_aware > 1.844306e-97)
  } else if (tumor_type == "LUSC") {
    plot_pval <- plot_pval %>% dplyr::filter(padj_aware > 1.268579e-154)
  } else if (tumor_type == "BRCA") {
    plot_pval <- plot_pval %>% dplyr::filter(padj_aware > 4.478953e-136)
  } else if (tumor_type == "LIHC") {
    plot_pval <- plot_pval %>% dplyr::filter(padj_aware > 1.316967e-62)
  } else if (tumor_type == "KIRC") {
    plot_pval <- plot_pval %>% dplyr::filter(padj_aware > 1.316967e-62)
  }
  
  plot_pval <- plot_pval %>%
    dplyr::mutate(tumor_type = tumor_type)
  
  plot_pval_combined_list[[tumor_type]] <- plot_pval
}

joint_pvalue <- bind_rows(plot_pval_combined_list)

comparison_pval <- ggplot(joint_pvalue, aes(x=-log10(padj_naive), y=-log10(padj_aware), color = cnv_group)) + 
  geom_point(shape=20, size=3) +
  geom_abline()+
  xlab("FDR CN-aware") +
  ylab ("FDR CN-naive") +
  scale_x_continuous(breaks = seq(0, 140, by = 40))+
  scale_y_continuous(breaks = seq(0, 140, by = 40))+
  scale_color_manual(name = "CNV group", values = cnv_colors)+
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1)))+
  facet_wrap(~factor(cnv_group, levels = c("loss", "neutral", "gain", "amplification")), nrow = 1)+
  theme_bw()+
  ggh4x::facet_nested(factor(tumor_type, levels = c("LUSC", "BRCA", "LIHC", "KIRC"))~
  factor(cnv_group, levels = c("loss", "neutral", "gain", "amplification")))+
  theme(
    legend.position="",
    axis.title.x = element_text(size=16, color = "black"),  
    axis.title.y = element_text(size=16, color = "black"),
    axis.text.x = element_text(size=16, color = "black"),    
    axis.text.y = element_text(size=16, color = "black"),
    strip.text = element_text(size = 16, face = "plain", color = "black"),
    legend.text = element_text(size = 16, color = "black"),                          
    legend.title = element_text(size = 18, face = "plain", color = "black") 
  )
comparison_pval

ggsave("deconveilCaseStudies/plots/main/scatter_pvalue_luad.png", dpi = 400, width = 8.0, height = 3.0, plot = comparison_pval)
ggsave("deconveilCaseStudies/plots/supplementary/scatter_pval_pancancer.png", dpi = 400, width = 8.0, height = 8.0, plot = comparison_pval)


## Compare p-value for amplified and gain-affected genes ##
gain_amp_pval <- joint_pvalue %>% filter(cnv_group %in% c("gain", "amplification"))

mean_naive <- mean(gain_amp_pval$padj_naive, na.rm = TRUE)
mean_aware <- mean(gain_amp_pval$padj_aware, na.rm = TRUE)
mean_shift <- mean_aware - mean_naive
mean_increase <- mean_aware > mean_naive

# Count number of genes where p-value increased
num_increased <- sum(gain_amp_pval$padj_aware > gain_amp_pval$padj_naive, na.rm = TRUE)
total_genes <- nrow(gain_amp_pval)
percentage_increased <- (num_increased / total_genes) * 100

# Count number of genes losing significance
num_lost_significance <- sum(gain_amp_pval$padj_naive < 0.05 & gain_amp_pval$padj_aware > 0.05, na.rm = TRUE)

# Print results
cat("Mean p-value shift:", mean_shift, "\n")
cat("Mean p-value increased:", mean_increase, "\n")
cat("Number of genes with increased p-value:", num_increased, " (", round(percentage_increased, 2), "%)\n")
cat("Number of genes losing statistical significance:", num_lost_significance, "\n")


### Sankey dynamic gene groups transitions ###

tumor_type <- c("LUAD")
lfc_cut = 1.0
pval_cut = 0.05 

data <- read_data(tumor_type)

res_naive <- annotate_results(data$res_pydeseq, lfc_cut, pval_cut, "CN naive", tumor_type)
res_aware <- annotate_results(data$res_deconveil, lfc_cut, pval_cut, "CN aware", tumor_type)

res_naive <- res_naive %>% dplyr::select(DEtype) %>% dplyr::rename(CN_naive = DEtype)
res_aware <- res_aware %>% dplyr::select(DEtype) %>% dplyr::rename(CN_aware = DEtype)

res_join <- cbind(res_naive, res_aware)

data_flow <- res_join %>%
  group_by(CN_naive,CN_aware) %>%
  summarise(freq = n()) %>%
  ungroup()

data_ggforce <- data_flow  %>%
  gather_set_data(1:2) %>%        
  arrange(x,CN_naive,desc(CN_aware))

data_ggforce$CN_naive <- factor(data_ggforce$CN_naive)
data_ggforce$CN_aware <- factor(data_ggforce$CN_aware)

data_ggforce <- data_ggforce %>%
  group_by(CN_naive, CN_aware) %>%
  mutate(y_mid = freq / 2) %>% 
  na.omit()

g_group_colors <- c("Down-reg" = "#3C5488B2", "n.s." = "lightgray", "Up-reg" = "#CC7677")

sankey <- ggplot(data_ggforce, aes(x = x, id = id, split = y, value = freq)) +
  geom_parallel_sets(aes(fill = CN_naive), alpha = 0.9, axis.width = 0.2,
                     n = 4415, strength = 0.5, color = "black", linewidth = 0.3) +
  geom_parallel_sets_axes(axis.width = 0.25, fill = "gray93",
                          color = "gray", linewidth = 0.5) +  
  #geom_parallel_sets_labels(colour = 'gray35', size = 2.0, angle = 0, fontface = "plain") +
  scale_fill_manual(values = g_group_colors, name = "Gene group") +
  scale_color_manual(values = g_group_colors) +
  scale_x_continuous(breaks = 1:2, labels = c("CN-naive", "CN-aware"))+
  theme_classic() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 15, face = "plain"),
    legend.text = element_text(size = 13),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x  = element_blank()
  )
sankey

ggsave("deconveilCaseStudies/plots/main/sankey_luad.png", dpi = 400, width = 4.0, height = 5.0, plot = sankey)
ggsave("deconveilCaseStudies/plots/supplementary/sankey_kirc.png", dpi = 400, width = 4.0, height = 5.0, plot = sankey)

