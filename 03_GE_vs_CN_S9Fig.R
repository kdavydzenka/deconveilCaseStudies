setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/deconveilCaseStudies/")

pkgs <- c("tidyverse", "stringr", "DESeq2", "patchwork")
sapply(pkgs, require, character.only = TRUE)

source("utils/utils.R")

# Load RNA + CN data

cnv_path <- "TCGA/LUAD/test/cnv.csv"
rna_path <- "TCGA/LUAD/test/rna.csv"
metadata_path <- "TCGA/LUAD/test/metadata.csv"

cnv <- read.csv(cnv_path) %>% 
  remove_rownames() %>%
  column_to_rownames(var = "X") * 2

cnv[cnv > 8] <- 8

rna <- read.csv(rna_path) %>% 
  remove_rownames() %>%
  column_to_rownames(var = "X") 

meta <- read.csv(metadata_path)
colnames(meta) <- c ("sampleID", "condition")
meta <- meta %>%
  dplyr::mutate(condition = recode(condition,
                            "A" = "Normal",
                            "B" = "Tumor")) %>% 
  remove_rownames() %>%
  column_to_rownames(var = "sampleID") 

# Normalize counts

dds <- DESeqDataSetFromMatrix(
  countData = rna,
  colData   = meta,
  design    = ~ 1
)

vsd <- vst(dds, blind = TRUE)
expr_mat <- assay(vsd)

meta$sample_id <- rownames(meta)

tumor_s <- meta %>%
  filter(condition == "Tumor") %>%      
  pull(sample_id)

normal_s <- meta %>%
  filter(condition == "Normal") %>%     
  pull(sample_id)

expr_tumor  <- expr_mat[, tumor_s]
expr_normal <- expr_mat[, normal_s]

# Helper functions 

pick_representative_gene <- function(df, group_col, effect_col, group_name) {
  df %>%
    mutate(group_col = str_remove(.data[[group_col]], "s$")) %>%
    dplyr::filter(.data[[group_col]] == group_name) %>%
    mutate(abs_eff = abs(.data[[effect_col]])) %>%
    arrange(abs(abs_eff - median(abs_eff, na.rm = TRUE))) %>%
    dplyr::slice(1) %>%
    pull(geneID)
}

make_plot_df_for_gene <- function(gene_id, expr_mat, cn_mat, sample_info) {
  
  common_samples <- intersect(colnames(expr_mat), colnames(cn_mat))
  
  expr_vec <- expr_mat[gene_id, common_samples]
  cn_vec   <- cn_mat[gene_id, common_samples]
  
  tibble(
    sample_id = common_samples,
    geneID    = gene_id,
    expr      = as.numeric(expr_vec),
    cn        = as.numeric(cn_vec)
  ) %>%
    left_join(sample_info, by = "sample_id")  
}

read_data <- function(tumor_type) {
  list(
    res_pydeseq = read.csv(paste0("results_tcga/", tumor_type, "/new/res_CNnaive.csv")),
    res_deconveil = read.csv(paste0("results_tcga/", tumor_type, "/new/res_CNaware.csv")),
    cnv_tumor = read.csv(paste0("results_tcga/", tumor_type, "/old/cnv_tumor.csv")) %>%
      remove_rownames() %>%
      column_to_rownames(var = "X") * 2
  )
}

# Load deconveil/NB results

res_pp_nb_03 <- read.csv("results_tcga/LUAD/new/nb_regression/luad_res_pp_0.3.csv")

tumor_types <- c("LUAD") 
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
    dplyr::select(geneID, logFC_naive, padj_naive, padj_stageR_naive, DE_confirmed_naive, isDE_naive, DEtype_naive, tumor_type_naive, method_naive, cnv_mean_naive, gene_group) %>% 
    rename_with(~ c("geneID", "log2FC", "padj", "padj_stageR", "DE_conf", "isDE", "DEtype", "tumor_type", "method", "cnv_mean", "gene_group"), everything())
  
  cn_aware <- bind_rows(
    gene_groups$d_sensitive %>% mutate(gene_group = "DSGs"),
    gene_groups$d_insensitive %>% mutate(gene_group = "DIGs"),
    gene_groups$d_compensated %>% mutate(gene_group = "DCGs")
  ) %>%
    dplyr::select(geneID, logFC_aware, padj_aware, padj_stageR_aware, DE_confirmed_aware, isDE_aware, DEtype_aware, tumor_type_aware, method_aware, cnv_mean_aware, gene_group) %>% 
    rename_with(~ c("geneID", "log2FC", "padj", "padj_stageR", "DE_conf", "isDE", "DEtype", "tumor_type", "method", "cnv_mean", "gene_group"), everything())
  
  list(
    tumor_type = tumor_type,
    cnv_data = cnv_data,
    cn_naive = cn_naive,
    cn_aware = cn_aware
  )
})

res_deconveil <- results_list[[1]][["cn_aware"]]


# Pick representative genes

# DeConveil representative DSG/DCG based on medin log2FC
#rep_dsg <- pick_representative_gene(res_deconveil, "gene_group", "log2FC", "DSGs")
#rep_dcg <- pick_representative_gene(res_deconveil, "gene_group", "log2FC", "DCGs")

rep_dig <- "MUC6"
rep_dsg <- "MYCN"
rep_dcg <- "BRAF"
rep_hyper <- "GINS4" # GINS4


# For HYPER, use NB label
#rep_hyper <- res_pp_nb_03 %>%
  #dplyr::filter(label_nb == "HYPER", status == "ok") %>%
  #arrange(abs(shrunk_comp - median(shrunk_comp, na.rm = TRUE))) %>%
  #dplyr::slice(1) %>%
  #pull(gene)

df_dsg   <- make_plot_df_for_gene(rep_dsg, expr_mat, cnv, meta)
df_dcg   <- make_plot_df_for_gene(rep_dcg, expr_mat, cnv, meta)
df_dig <- make_plot_df_for_gene(rep_dig, expr_mat, cnv, meta)
df_hyper <- make_plot_df_for_gene(rep_hyper, expr_mat, cnv, meta)

plot_df <- bind_rows(df_dsg, df_dcg, df_dig, df_hyper)

y_limits <- range(plot_df$expr, na.rm = TRUE)

plot_expr_vs_cn <- function(plot_df,
                            gene_id,
                            panel_title = "",
                            class_label = NULL,  # e.g. "DSG", "DCG", "DIG", "HYPER"
                            use_normals_for_expected = FALSE,
                            diploid_window = 0.25,
                            expected_slope_per_copy = 0.6,  # linear slope on VST scale
                            y_limits = NULL,
                            y_break_step = 2) {
  
  # Packages (use namespaces if you don't want library())
  # library(dplyr); library(tidyr); library(ggplot2)
  
  df_g <- plot_df %>%
    dplyr::filter(geneID == gene_id) %>%
    dplyr::filter(!is.na(expr), !is.na(cn))
  
  if (nrow(df_g) < 8) {
    stop("Not enough points to plot for gene: ", gene_id)
  }
  
  has_type   <- "type" %in% colnames(df_g)
  has_normal <- has_type && any(df_g$type == "Normal", na.rm = TRUE)
  has_tumor  <- !has_type || any(df_g$type == "Tumor", na.rm = TRUE)
  
  # observed trend: fit in tumors if available
  df_obs <- if (has_type && any(df_g$type == "Tumor", na.rm = TRUE)) {
    df_g %>% dplyr::filter(type == "Tumor")
  } else {
    df_g
  }
  
  fit_obs <- stats::lm(expr ~ cn, data = df_obs)
  
  cn_range <- range(df_g$cn, na.rm = TRUE)
  new_cn <- data.frame(cn = seq(cn_range[1], cn_range[2], length.out = 100))
  new_cn$expr_obs <- stats::predict(fit_obs, newdata = new_cn)
  
  # expected trend (forced linear)
  if (use_normals_for_expected && has_normal) {
    df_norm <- df_g %>% dplyr::filter(type == "Normal")
    fit_exp <- stats::lm(expr ~ cn, data = df_norm)
    new_cn$expr_exp <- stats::predict(fit_exp, newdata = new_cn)
    
  } else if (has_tumor) {
    df_dip <- df_obs %>% dplyr::filter(abs(cn - 2) <= diploid_window)
    y2 <- if (nrow(df_dip) >= 5) {
      stats::median(df_dip$expr, na.rm = TRUE)
    } else {
      stats::median(df_obs$expr, na.rm = TRUE)
    }
    
    new_cn$expr_exp <- y2 + (new_cn$cn - 2) * expected_slope_per_copy
    
  } else {
    new_cn$expr_exp <- new_cn$expr_obs
  }
  
  new_cn_long <- new_cn %>%
    tidyr::pivot_longer(
      cols = c(expr_obs, expr_exp),
      names_to = "trend",
      values_to = "expr_trend"
    ) %>%
    dplyr::mutate(
      trend = dplyr::recode(trend, expr_obs = "Observed", expr_exp = "Expected")
    )
  
  # safe title fallback (avoid %||%)
  plot_title <- if (!is.null(class_label) && nzchar(class_label)) class_label else panel_title
  
  # build plot object first
  p <- ggplot2::ggplot(df_g, ggplot2::aes(x = cn, y = expr)) +
    ggplot2::geom_point(color = "grey50", alpha = 0.7, size = 1.7) +
    ggplot2::geom_line(
      data = new_cn_long,
      ggplot2::aes(y = expr_trend, color = trend, linetype = trend),
      linewidth = 1
    ) +
    ggplot2::scale_color_manual(
      values = c("Observed" = "darkblue", "Expected" = "orange"),
      name = NULL
    ) +
    ggplot2::scale_linetype_manual(
      values = c("Observed" = "solid", "Expected" = "dashed"),
      name = NULL
    ) +
    ggplot2::scale_x_continuous(
      breaks = sort(unique(df_g$cn)),
      minor_breaks = NULL
    ) +
    {
      if (!is.null(y_limits)) {
        scale_y_continuous(
          breaks = seq(
            from = floor(y_limits[1]),
            to   = ceiling(y_limits[2]),
            by   = y_break_step
          )
        )
      }
    } +
    ggplot2::labs(
      title = plot_title,
      subtitle = gene_id,
      x = "Copy number",
      y = "Normalized expression"
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "plain"),
      plot.subtitle = ggplot2::element_text(face = "plain"),
      axis.title = ggplot2::element_text(face = "plain"),
      legend.position = "none"
    )
  
  # add uniform y-limits if requested
  if (!is.null(y_limits)) {
    p <- p + ggplot2::coord_cartesian(ylim = y_limits)
  }
  
  p
}

pA <- plot_expr_vs_cn(plot_df, rep_dsg, "DSG", 
                      use_normals_for_expected = FALSE,
                      y_limits = y_limits,
                      y_break_step = 3)
pA

pB <- plot_expr_vs_cn(plot_df, rep_dcg, "DCG", 
                      use_normals_for_expected = FALSE,
                      y_limits = y_limits,
                      y_break_step = 3)
pB

pC <- plot_expr_vs_cn(plot_df, rep_dig, "DIG", 
                      use_normals_for_expected = FALSE,
                      y_limits = y_limits,
                      y_break_step = 3)
pC

pD <- plot_expr_vs_cn(plot_df, rep_hyper, "HYPER",
                      use_normals_for_expected = FALSE,
                      y_limits = y_limits,
                      y_break_step = 3)
pD

p <- pA | pB | pC | pD
p

ggsave("plots/supplementary/png/scatter_cn_ge_v2.png", dpi = 500, width = 10.0, height = 3.5, plot = p)
