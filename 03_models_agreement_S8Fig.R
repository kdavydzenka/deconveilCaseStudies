setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/deconveilCaseStudies/")

pkgs <- c("tidyverse", "forcats", "scales", "ggalluvial", "stringr", "purrr", "patchwork")
sapply(pkgs, require, character.only = TRUE)

source("utils/utils.R")
#source("deconveilCaseStudies/utils/utils_plot.R")

### Load results fit ###

tumor_types <- c("LUAD") 
lfc_cut = 1.0
pval_cut = 0.05 
loss_threshold = 0.25

# Read data files for a tumor type - DeConveil 
read_data <- function(tumor_type) {
  list(
    res_pydeseq = read.csv(paste0("results_tcga/", tumor_type, "/new/res_CNnaive.csv")),
    res_deconveil = read.csv(paste0("results_tcga/", tumor_type, "/new/res_CNaware.csv")),
    cnv_tumor = read.csv(paste0("results_tcga/", tumor_type, "/old/cnv_tumor.csv")) %>%
      remove_rownames() %>%
      column_to_rownames(var = "X") * 2
  )
}

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


run_concordance <- function(lfc_cut,
                            tumor_types = c("LUAD"),
                            pval_cut = 0.05,
                            loss_threshold = 0.25,
                            nb_dir = "results_tcga/LUAD/new/nb_regression") {
  
  # format like "0.3" or "1.0" for filenames
  lfc_tag <- format(lfc_cut, nsmall = 1, trim = TRUE)
  
  lfc_label = dplyr::case_when(
    lfc_cut == 0.3 ~ "Moderate effect (τ > 0.3)",
    lfc_cut == 1.0 ~ "Strong effect (τ > 1.0)",
    TRUE ~ paste0("|LFC| > ", lfc_tag)
  )
  
  # DeConveil processing 
  res_deconveil <- bind_rows(lapply(tumor_types, function(tumor_type) {
    
    data <- read_data(tumor_type)
    cnv_data <- process_cnv_data(data$cnv_tumor, loss_threshold)
    
    res_naive <- annotate_results(data$res_pydeseq, lfc_cut, pval_cut,
                                  "CN naive", tumor_type)
    res_aware <- annotate_results(data$res_deconveil, lfc_cut, pval_cut,
                                  "CN aware", tumor_type)
    
    combined <- combine_results(res_naive, res_aware,
                                data$cnv_tumor, cnv_data$loss_labels)
    
    combined$res_naive <- combined$res_naive %>%
      rename(logFC = log2FoldChange)
    
    combined$res_aware <- combined$res_aware %>%
      rename(logFC = log2FoldChange)
    
    res_joint <- combined$res_naive %>%
      inner_join(combined$res_aware, by = "geneID",
                 suffix = c("_naive", "_aware"))
    
    gene_groups <- list(
      d_sensitive = res_joint %>%
        filter(DEtype_naive %in% c("Up-reg", "Down-reg") &
                 DEtype_aware == "n.s."),
      d_insensitive = res_joint %>%
        filter(DEtype_naive == DEtype_aware &
                 DEtype_naive %in% c("Up-reg", "Down-reg")),
      d_compensated = res_joint %>%
        filter(DEtype_naive == "n.s." &
                 DEtype_aware %in% c("Up-reg", "Down-reg"))
    )
    
    bind_rows(
      gene_groups$d_sensitive   %>% mutate(gene_group = "DSGs"),
      gene_groups$d_insensitive %>% mutate(gene_group = "DIGs"),
      gene_groups$d_compensated %>% mutate(gene_group = "DCGs")
    ) %>%
      select(geneID, gene_group)
  }))
  
  # NB regression: read file matching lfc_cut 
  nb_file <- file.path(nb_dir, paste0("luad_res_pp_", lfc_tag, ".csv"))
  if (!file.exists(nb_file)) {
    stop("NB file not found: ", nb_file)
  }
  
  res_pp_nb <- read.csv(nb_file) %>%
    filter(!status %in% c("error", "skipped")) %>%
    transmute(
      geneID = gene,
      nb_group = label_nb
    )
  
  # Concordance 
  concord_df <- res_deconveil %>%
    mutate(decon_group = str_remove(gene_group, "s$")) %>%
    inner_join(res_pp_nb, by = "geneID")
  
  concord_summary <- concord_df %>%
    count(decon_group, nb_group) %>%
    group_by(decon_group) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup() %>%
    mutate(
      lfc_cut = lfc_cut,
      lfc_label  = lfc_label,
      decon_group = factor(decon_group, levels = c("DIG", "DSG", "DCG")),
      nb_group    = factor(nb_group,    levels = c("DSG", "DCG", "HYPER", "OTHER")),
      prop_label  = percent(prop, accuracy = 1)
    )
  
  concord_summary
}

# Run for both thresholds 
lfc_cuts <- c(0.3, 1.0)
concord_all <- bind_rows(lapply(lfc_cuts, run_concordance))

# Heatmap
hm <- ggplot(concord_all, aes(x = nb_group, y = decon_group, fill = prop)) +
  geom_tile(color = "white") +
  geom_text(aes(label = prop_label), size = 3) +
  scale_fill_gradientn(
    colours = c(
      "#ffffff",  # 0 – white
      "#deebf7",  # light blue
      "#9ecae1",  # mid blue
      "#6999CC"   # higher blue
    ),
    limits = c(0, 1),
    name   = "Proportion"
  ) +
  facet_wrap(~ lfc_label) +
  labs(
    x = "NB regression class",
    y = "DeConveil class",
    #title = "Concordance between DeConveil and NB regression"
  ) +
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "plain")
  )
hm

ggsave("plots/supplementary/png/heatmap_concordance.png", dpi = 500, width = 6.0, height = 3.0, plot = hm)

# Sankey plot

dcg_flows_prop <- concord_all %>%
  filter(decon_group == "DCG")

sankey <- ggplot(
  dcg_flows_prop,
  aes(
    axis1 = decon_group,
    axis2 = nb_group,
    y     = prop    # proportion within DCG
  )
) +
  geom_alluvium(aes(fill = nb_group), width = 0.2, alpha = 0.8) +
  geom_stratum(width = 0.2, color = "grey30") +
  geom_text(
    stat = "stratum",
    aes(label = after_stat(stratum)),
    size = 3
  ) +
  scale_x_discrete(limits = c("DeConveil", "NB regression"), expand = c(.1, .1)) +
  scale_fill_manual(
    values = c(
      DSG   = "#DF8F44",
      DCG   = "#DDCC77",
      HYPER = "#F8766D",
      OTHER = "#aaaaaa"
    ),
    name = "NB class"
  )+
  facet_wrap(~ lfc_label) +
  labs(
    x = NULL,
    y = "Proportion of DCGs",
    #title = "Reassignment of DeConveil DCGs to NB regression classes",
    #subtitle = "Faceted by effect-size regime"
  ) +
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 12, color = "black"),  
    axis.text.y = element_text(size = 12),                         
    axis.title.x = element_text(size = 12),          
    axis.title.y = element_text(size = 12),
    strip.text = element_text(face = "bold", size = 12),
    plot.title = element_text(hjust = 0.5, face = "plain"),
    plot.subtitle = element_text(hjust = 0.5)
  )
sankey

ggsave("plots/supplementary/png/sankey_dcg_nb.png", dpi = 500, width = 8.0, height = 3.5, plot = sankey)

# DCG dominance/entropy + subtype classification

shannon_entropy_norm <- function(p_vec) {
  # normalized Shannon entropy in [0,1]
  # p_vec must sum to 1
  eps <- 1e-12
  p <- pmax(p_vec, eps)
  p <- p / sum(p)
  H <- -sum(p * log(p))
  H / log(length(p))
}

# Simulate NB posterior for C_g and class probabilities
# tau = NB compensation threshold (0.3 or 1.0)
simulate_nb_probs <- function(df_nb, tau, n_draws = 2000, seed = 1) {
  set.seed(seed)
  
  df_nb %>%
    rowwise() %>%
    mutate(
      probs = list({
        # default NA values
        p_DCG   <- NA_real_
        p_DSG   <- NA_real_
        p_HYPER <- NA_real_
        
        # check we have valid posterior params
        if (!is.na(mean_b_scaling) && !is.na(sd_b_scaling) &&
            !is.na(mean_b_deviation) && !is.na(sd_b_deviation) &&
            sd_b_scaling > 0 && sd_b_deviation > 0) {
          
          b1_draw <- rnorm(n_draws, mean = mean_b_scaling,   sd = sd_b_scaling)
          b2_draw <- rnorm(n_draws, mean = mean_b_deviation, sd = sd_b_deviation)
          
          valid <- abs(b1_draw) > 1e-8
          if (any(valid)) {
            b1v <- b1_draw[valid]
            b2v <- b2_draw[valid]
            
            # compensation ratio with same shrinkage as your comp score
            C_draw <- (b2v / b1v) * (1 - p_value)
            
            p_DCG   <- mean(C_draw <= -tau)
            p_HYPER <- mean(C_draw >=  tau)
            p_DSG   <- mean(abs(C_draw) < tau)
            
            p_sum <- p_DCG + p_DSG + p_HYPER
            if (p_sum > 0) {
              p_DCG   <- p_DCG   / p_sum
              p_DSG   <- p_DSG   / p_sum
              p_HYPER <- p_HYPER / p_sum
            } else {
              p_DCG   <- NA_real_
              p_DSG   <- NA_real_
              p_HYPER <- NA_real_
            }
          }
        }
        
        tibble(p_DCG = p_DCG, p_DSG = p_DSG, p_HYPER = p_HYPER)
      })
    ) %>%
    unnest_wider(probs) %>%
    ungroup()
}

# DCG subtype classification 

classify_dcg_subtypes <- function(res_deconveil,
                                  res_pp_nb,
                                  tau,
                                  n_draws = 2000,
                                  seed = 1,
                                  directional_nb_classes = c("DSG", "HYPER")) {

  # DeConveil DCGs
  decon_dcg <- res_deconveil %>%
    mutate(decon_group = stringr::str_remove(gene_group, "s$")) %>%
    filter(decon_group == "DCG") %>%
    select(geneID, decon_group)

  # NB info (posterior means/sds etc)
  nb_clean <- res_pp_nb %>%
    transmute(
      geneID            = gene,
      status            = status,
      nb_group          = label_nb,
      mean_b_scaling    = mean_b_scaling,
      sd_b_scaling      = sd_b_scaling,
      mean_b_deviation  = mean_b_deviation,
      sd_b_deviation    = sd_b_deviation,
      p_value           = p_value
    )

  # Join all DeConveil DCGs with NB results
  dcg_join <- decon_dcg %>%
    left_join(nb_clean, by = "geneID") %>%
    mutate(
      nb_missing  = is.na(nb_group),
      nb_filtered = (!is.na(status) & status != "ok") | nb_missing
    )

  # Simulate NB posterior probs for DCGs with valid NB info
  dcg_for_sim <- dcg_join %>%
    filter(!nb_filtered)

  dcg_probs <- simulate_nb_probs(dcg_for_sim, tau = tau, n_draws = n_draws, seed = seed) %>%
    mutate(
      dominance = pmax(p_DCG, p_DSG, p_HYPER, na.rm = TRUE),
      dominant_class = c("DCG", "DSG", "HYPER")[max.col(
        cbind(p_DCG, p_DSG, p_HYPER),
        ties.method = "first"
      )],
      entropy = purrr::pmap_dbl(
        list(p_DCG, p_DSG, p_HYPER),
        ~ shannon_entropy_norm(c(..1, ..2, ..3))
      )
    )

  # Merge probs back to all DCGs and assign subtypes
  dcg_all <- dcg_join %>%
    left_join(
      dcg_probs %>%
        select(geneID, p_DCG, p_DSG, p_HYPER, dominance, dominant_class, entropy),
      by = "geneID"
    ) %>%
    mutate(
      # Filtered if NB missing/filtered OR probs could not be computed
      filtered_flag = nb_filtered | is.na(entropy) | is.na(dominance),

      dcg_subtype = case_when(
        # Filtered DCG
        filtered_flag ~ "Filtered DCG",

        # Pure DCG: NB=DCG OR entropy >= 0.7
        nb_group == "DCG" | entropy >= 0.7 ~ "Pure DCG",

        # Directional DCG: dominant class directional + dominance>=0.6
        dominance >= 0.6 & dominant_class %in% directional_nb_classes ~ "Directional DCG",

        # Ambiguous DCG: no dominant NB, entropy in [0.3, 0.7)
        dominance < 0.6 & entropy >= 0.3 & entropy < 0.7 ~ "Ambiguous DCG",

        # Fallback
        TRUE ~ "Ambiguous DCG"
      ),

      dcg_subtype = factor(
        dcg_subtype,
        levels = c("Pure DCG", "Directional DCG", "Ambiguous DCG", "Filtered DCG")
      )
    )

  dcg_all
}

# tau must match how you thresholded C in NB (0.3, 1.0, ...)

res_deconveil_03 <- results_list[[1]][["cn_aware"]]
res_pp_nb_03 <- read.csv("results_tcga/LUAD/new/nb_regression/luad_res_pp_0.3.csv")

res_deconveil_1 <- results_list[[1]][["cn_aware"]]
res_pp_nb_1 <- read.csv("results_tcga/LUAD/new/nb_regression/luad_res_pp_1.0.csv")

dcg_subtypes_03 <- classify_dcg_subtypes(
  res_deconveil = res_deconveil_03,
  res_pp_nb     = res_pp_nb_03,
  tau           = 0.3,
  n_draws       = 2000,
  seed          = 1
)

dcg_subtypes_1 <- classify_dcg_subtypes(
  res_deconveil = res_deconveil_1,
  res_pp_nb     = res_pp_nb_1,
  tau           = 1.0,
  n_draws       = 2000,
  seed          = 1
)

summarize_dcg_subtypes <- function(dcg_df) {
  dcg_df %>%
    count(dcg_subtype) %>%
    mutate(prop = n / sum(n))
}

summarize_dcg_subtypes(dcg_subtypes_03)
summarize_dcg_subtypes(dcg_subtypes_1)

plot_dominance_entropy <- function(dcg_df, title_suffix = "") {
  
  # enforce consistent order
  dcg_df <- dcg_df %>%
    mutate(
      dcg_subtype = factor(
        dcg_subtype,
        levels = c("Pure DCG", "Directional DCG", "Ambiguous DCG", "Filtered DCG")
      )
    )
  
  subtype_colors <- c(
    "Pure DCG"        = "#DDCC77",  
    "Directional DCG" = "#cc6a70b2",  
    "Ambiguous DCG"   = "#537293",
    "Filtered DCG"    = "#bdbdbd"   
  )
  
  # Dominance plot 
  p1 <- ggplot(
    dcg_df %>% filter(!is.na(dominance)),
    aes(x = dcg_subtype, y = dominance, fill = dcg_subtype)
  ) +
    geom_violin(trim = FALSE, alpha = 0.8, color = NA) +
    geom_boxplot(
      width = 0.15,
      outlier.size = 0.4,
      color = "black",
      alpha = 0.6
    ) +
    scale_fill_manual(values = subtype_colors, guide = "none") +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.2)
      #expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      x = NULL,
      y = "Dominance (max NB class probability)",
      title = paste("DCG dominance by subtype", title_suffix)
    ) +
    theme_bw() +
    geom_hline(yintercept = c(0.3, 0.7), linetype = "dashed", color = "grey50")+
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 12, angle = 30, hjust = 1),
      axis.text.y = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      legend.position = "right",
      plot.title = element_text(face = "bold", hjust = 0.5, size = 12)
    )
  
  # Entropy plot 
  p2 <- ggplot(
    dcg_df %>% filter(!is.na(entropy)),
    aes(x = dcg_subtype, y = entropy, fill = dcg_subtype)
  ) +
    geom_violin(trim = FALSE, alpha = 0.8, color = NA) +
    geom_boxplot(
      width = 0.15,
      outlier.size = 0.4,
      color = "black",
      alpha = 0.6
    ) +
    scale_fill_manual(values = subtype_colors, guide = "none") +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.2)
      #expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      x = NULL,
      y = "Normalized entropy",
      title = paste(title_suffix)
    ) +
    theme_bw() +
    geom_hline(yintercept = c(0.3, 0.7), linetype = "dashed", color = "grey50")+
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 12, angle = 30, hjust = 1),
      axis.text.y = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      legend.position = "right",
      plot.title = element_text(face = "plain", hjust = 0.5, size = 12)
    )
  
  list(
    dominance_plot = p1,
    entropy_plot   = p2
  )
}

plots_03 <- plot_dominance_entropy(dcg_subtypes_03, "Moderate effect")
plots_03
plots_1 <- plot_dominance_entropy(dcg_subtypes_1, "Strong effect")
plots_1

ggsave("plots/supplementary/png/violin_dcg_entropy_high.png", dpi = 500, width = 4.0, height = 4.0, plot = plots_1)

