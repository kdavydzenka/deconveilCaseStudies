setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/")
pkgs <- c("tidyverse", "colorspace", "ggpubr", "ggpointdensity", "ggeasy", "gridExtra", "patchwork",
          "ggalluvial", "ggridges", "ggforce", "ggparallel", "alluvial")
sapply(pkgs, require, character.only = TRUE)

source("deconveilCaseStudies/utils/utils_plot.R")

### Load data ###

# null simulation
all_p <- read.csv("DeConveil/tutorials/sim_results/sim_null/null_pval_all.csv", header = T, sep = ",")

# pipeline comparison
all_res_stager <- read.csv("DeConveil/tutorials/sim_results/data_to_plot/performance_data_stageR_all.csv", header = T, sep = ",")
all_res_standard <- read.csv("DeConveil/tutorials/sim_results/data_to_plot/performance_data_standard_all.csv", header = T, sep = ",")

all_res_stager <- all_res_stager %>% 
  mutate(pipeline = "stageR")

all_res_standard<- all_res_standard %>% 
  mutate(pipeline = "standard")

all_res <- rbind(all_res_standard, all_res_stager)


# pairwise comparison
all_res_pairwise <-  read.csv("DeConveil/tutorials/sim_results/data_to_plot/performance_data_stageR_pairwise.csv", header = T, sep = ",")

# simulation 2
all_res_sim2 <-  read.csv("DeConveil/tutorials/sim_results/data_to_plot/performance_data_stageR_sim2.csv", header = T, sep = ",")

# Reshape to long
df_long <- all_p %>%
  dplyr::select(p_naive, p_aware, sample_size, rep) %>%
  pivot_longer(
    cols = c(p_naive, p_aware),
    names_to = "model",
    values_to = "p_value"
  ) %>%
  mutate(
    model = recode(model,
                   p_naive = "CN-naive",
                   p_aware = "CN-aware")
  )

# Faceted p-value distribution (histograms)

p_hist <- ggplot(df_long, aes(x = p_value, fill = model)) +
  geom_histogram(
    bins = 20,
    aes(y = after_stat(density)),
    color = "white"
  ) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    linewidth = 0.8
  ) +
  facet_grid(
    rows = vars(model),
    cols = vars(sample_size)
  ) +
  scale_fill_manual(
    values = c(
      "CN-naive" = "#3C5488B2",   # blue
      "CN-aware" = "#6F99ADB2"    # red
    )
  ) +
  guides(fill = "none") +
  scale_x_continuous(limits = c(0, 1)) +
  labs(
    title = "Null p-value distributions",
    x = "p-value",
    y = "Density"
  ) +
  theme_bw(base_size = 12)

ggsave("deconveilCaseStudies/revision/plots/pval_hist.png", dpi = 400, width = 10.0, height = 5.0, plot = p_hist)

# QQ-plot

qq_p <- ggplot(df_long, aes(x = p_value, color = model)) +
  stat_ecdf(size = 1.1) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    color = "black"
  ) +
  facet_wrap(~ sample_size, nrow = 1) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(
    values = c(
      "CN-naive" = "#3C5488B2",   # blue
      "CN-aware" = "#6F99ADB2"    # muted red/gray
    )
  ) +
  labs(
    title = "Null p-value calibration by sample size",
    subtitle = "ECDF vs Uniform(0,1)",
    x = "p-value",
    y = "ECDF",
    color = "Model"
  ) +
  theme_bw(base_size = 12)+
  theme(
    legend.position = "bottom"
  )

qq_p

ggsave("deconveilCaseStudies/revision/plots/qq_pval.png", dpi = 400, width = 8.0, height = 3.5, plot = qq_p)

# Performance boxplot

strong_res <- all_res %>% 
  filter(cn_signal == "strong")

metrics_long <- strong_res %>%
  select(
    sample_size,
    replicate,
    #method,
    pipeline,
    cn_signal,
    comparison,
    Precision,
    Recall,
    F1.score,
    MCC
  ) %>%
  pivot_longer(
    cols = c(Precision, Recall, F1.score),
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
  group_by(sample_size, cn_signal, pipeline, comparison, metric) %>%
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
  levels = c("Precision", "Recall", "F1-score")
)

metrics_long <- metrics_long %>%
  mutate(
    comparison = gsub("_vs_rest", " vs others", comparison),
    #comparison = factor(comparison, levels = c("DCG vs DSG", "DCG vs DIG")),
    comparison = factor(comparison, levels = c("DCG vs others", "DSG vs others", "DIG vs others")),
    `CN signal` = factor(cn_signal, levels = c("strong", "weak")),
    metric = factor(metric, levels = c("Precision", "Recall", "F1-score"))
  )

cn_col <- c(
  "weak"   = "#5A959999",
  "strong" = "#FF634899"
)

pipeline_col <-  c(
  "standard"   = "#CC6677",
  "stageR" = "#FFC073")

method_col <- c(
  "DeConveil" = "#F66D7A", 
  "PyDESeq2" = "#27AD81FF"
)


p_perf <- function(
    df,
    colors,
    outpath = NULL
) {
  
  p <- ggplot(
    df,
    aes(
      x = factor(sample_size),
      y = value,
      color = pipeline
    )
  ) +
    geom_boxplot(
      fill = NA,
      linewidth = 0.3,
      outlier.shape = NA,
      position = position_dodge(width = 0.75)
    ) +
    facet_grid(
      #rows = vars(`CN signal`),
      rows = vars(comparison),
      cols = vars(metric)
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.2)
    ) +
    scale_color_manual(values = colors) +
    labs(
      title = "",
      x = "# Samples per condition",
      y = "Performance metric",
      color = "pipeline"
    ) +
    theme_bw(base_size = 10) +
    theme(
      legend.position = "bottom",
      strip.background = element_rect(fill = "grey90"),
      plot.title = element_text(hjust = 0.5)
    )
  
  if (!is.null(outpath)) {
    ggsave(outpath, p, dpi = 300, width = 16, height = 6)
  }
  
  return(p)
}


p_bxp <- p_perf(
  df = metrics_long,
  colors = pipeline_colors
)

ggsave("deconveilCaseStudies/revision/plots/perf_boxp_pipe.png", dpi = 400, width = 5.0, height = 5.0, plot = p_bxp)

# Confusion matrices

s20 <- read.csv("DeConveil/tutorials/sim_results/data_to_plot/merged_S20_CNstrong_R2_standard.csv", header = T, sep = ",")
s60 <- read.csv("DeConveil/tutorials/sim_results/data_to_plot/merged_S60_CNstrong_R2_standard.csv", header = T, sep = ",")

all_s <- rbind(s20, s60)

confusion_one_class <- function(df, target_class) {
  
  df_bin <- df %>%
    mutate(
      y_true = ifelse(truth_class == target_class, 1, 0),
      y_pred = ifelse(predicted_class == target_class, 1, 0)
    )
  
  tp <- sum(df_bin$y_true == 1 & df_bin$y_pred == 1)
  fn <- sum(df_bin$y_true == 1 & df_bin$y_pred == 0)
  fp <- sum(df_bin$y_true == 0 & df_bin$y_pred == 1)
  tn <- sum(df_bin$y_true == 0 & df_bin$y_pred == 0)
  
  matrix(
    c(tp, fn, fp, tn),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      c(paste("True", target_class), "True other"),
      c(paste("Pred", target_class), "Pred other")
    )
  )
}

plot_confusion_heatmap <- function(cm, title) {
  
  cm_df <- as.data.frame(cm) %>%
    tibble::rownames_to_column("Truth") %>%
    pivot_longer(
      -Truth,
      names_to = "Prediction",
      values_to = "Count"
    )
  
  ggplot(cm_df, aes(x = Prediction, y = Truth, fill = Count)) +
    geom_tile(color = "black", linewidth = 0.6) +
    geom_text(aes(label = Count), size = 5) +
    scale_fill_gradient(low = "white", high = "#4C72B0") +
    labs(title = title) +
    theme_minimal(base_size = 12) +
    theme(
      axis.title = element_blank(),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5)
    )
}

plot_confusion_all <- function(
    df,
    sample_size,
    cn_signal,
    replicate
) {
  
  classes <- c("DCG", "DSG", "DIG")
  
  plots <- lapply(classes, function(cls) {
    cm <- confusion_one_class(df, cls)
    
    plot_confusion_heatmap(
      cm,
      title = paste0(
        cls, " vs others\n",
        "S=", sample_size,
        ", CN=", cn_signal
        #", R=", replicate
      )
    )
  })
  
  wrap_plots(plots, nrow = 1)
}

p1 <- plot_confusion_all(
  df = s20,
  sample_size = 20,
  cn_signal = "strong",
  replicate = 2
)
p1

p2 <- plot_confusion_all(
  df = s60,
  sample_size = 60,
  cn_signal = "strong",
  replicate = 2
)
p2

p_confusion = p1 / p2
p_confusion
ggsave("deconveilCaseStudies/revision/plots/p_confusion_standard.png", 
       dpi = 400, width = 8.0, height = 5.0, plot = p_confusion)


### Sankey plot ###

gene_dosage_col <- c("DCG" = "#FAE48BFF",
            "DIG" = "#8F3931FF",
            "DSG" = "#FFB977",
            "NEUTRAL" = "#ADB6B6FF")



# Read data files for a tumor type
read_data <- function(sim) {
  list(
    res_pydeseq = read.csv(paste0("DeConveil/tutorials/sim_results/sim_res_fit/sim1/", "res_CNnaive_S60_G5000_CNstrong_R2.csv")),
    res_deconveil = read.csv(paste0("DeConveil/tutorials/sim_results/sim_res_fit/sim1/", "res_CNaware_S60_G5000_CNstrong_R2.csv")),
    truth_df = read.csv(paste0("DeConveil/tutorials/sim_results/sim_res_fit/sim1/", "truth_S60_G5000_CNstrong_R2.csv"))
  )
}

res_list <- read_data()

# Generic annotator: one method (naive OR aware)
annotate_de_results <- function(res_df,
                                method_label = c("naive","aware"),
                                lfc_cut = 0.2,
                                pval_cut = 0.05) {
  
  method_label <- match.arg(method_label)
  suf <- ifelse(method_label == "naive", "_naive", "_aware")
  
  # SAFE RENAME BLOCK 
  res_df <- res_df %>%
    dplyr::rename(
      geneID  = dplyr::any_of(c("X", "geneID")),
      log2FC  = log2FoldChange,
      padj_BH = padj
    )
  
  if (method_label == "naive") {
    DE_confirmed_vec <- if ("DE_confirmed" %in% colnames(res_df))
      res_df$DE_confirmed else rep(NA_real_, nrow(res_df))
  } else {
    DE_confirmed_vec <- if ("DE_confirmed" %in% colnames(res_df))
      res_df$DE_confirmed else rep(NA_real_, nrow(res_df))
  }
  
  # mutate DE calls 
  res_df <- res_df %>%
    mutate(
      DE_confirmed = DE_confirmed_vec,
      method      = ifelse(method_label == "naive", "CN-naive", "CN-aware"),
      
      # BH-based significance 
      isDE_BH = !is.na(padj_BH) & padj_BH < pval_cut,
      DEtype_BH = case_when(
        !isDE_BH       ~ "n.s.",
        log2FC > 0     ~ "Up-reg",
        log2FC < 0     ~ "Down-reg",
        TRUE           ~ "n.s."
      ),
      
      # stageR-based significance
      isDE_stageR = !is.na(DE_confirmed) & DE_confirmed == 1.0,
      DEtype_stageR = case_when(
        !isDE_stageR   ~ "n.s.",
        log2FC > 0     ~ "Up-reg",
        log2FC < 0     ~ "Down-reg",
        TRUE           ~ "n.s."
      )
    )
  
  # add suffixes to prepare naive + aware merge
  res_df <- res_df %>%
    dplyr::rename(
      !!paste0("log2FC",        suf) := log2FC,
      !!paste0("padj_BH",       suf) := padj_BH,
      !!paste0("DE_confirmed",  suf) := DE_confirmed,
      !!paste0("isDE_BH",       suf) := isDE_BH,
      !!paste0("isDE_stageR",   suf) := isDE_stageR,
      !!paste0("DEtype_BH",     suf) := DEtype_BH,
      !!paste0("DEtype_stageR", suf) := DEtype_stageR,
      !!paste0("method",        suf) := method
    )
  
  # select clean output
  res_df %>%
    dplyr::select(
      geneID,
      dplyr::starts_with("log2FC"),
      dplyr::starts_with("padj_"),
      dplyr::starts_with("DE_"),
      dplyr::starts_with("isDE_"),
      dplyr::starts_with("DEtype_"),
      dplyr::starts_with("method")
    )
}

#  Combine naive + aware + CN info

combine_de_and_cnv <- function(naive_annot, aware_annot) {
  res_joint <- naive_annot %>%
    inner_join(aware_annot, by = "geneID")

  res_joint
}

# Core classification rule for a *pair* of DE calls
classify_deconveil_pair <- function(de_naive, de_aware) {
  if (!de_naive && de_aware) {
    return("DCG")
  }
  if (de_naive && !de_aware) {
    return("DSG")
  }
  if (de_naive && de_aware) {
    return("DIG")
  }
  return("NEUTRAL")
}

# Add class_BH and class_stageR to the joint table

add_deconveil_classes <- function(joint_df) {
  joint_df %>%
    mutate(
      class_BH = mapply(
        classify_deconveil_pair,
        isDE_BH_naive,
        isDE_BH_aware
      ),
      class_stageR = mapply(
        classify_deconveil_pair,
        isDE_stageR_naive,
        isDE_stageR_aware
      )
    )
}

# Run full pipeline for ONE tumor type

run_pipeline <- function(pval_cut = 0.05) {
  
  res_list <- read_data()  # res_pydeseq, res_deconveil, truth_df
  
  res_naive_annot <- annotate_de_results(res_list$res_pydeseq,
                                         method_label = "naive",
                                         pval_cut = pval_cut)
  
  res_aware_annot <- annotate_de_results(res_list$res_deconveil,
                                         method_label = "aware",
                                         pval_cut = pval_cut)
  
  joint <- combine_de_and_cnv(res_naive_annot, res_aware_annot)
  joint_cls <- add_deconveil_classes(joint)
  
  # Summaries: class counts under BH vs stageR
  summary_classes <- joint_cls %>%
    summarize(
      n_BH_DCG  = sum(class_BH      == "DCG"),
      n_BH_DSG  = sum(class_BH      == "DSG"),
      n_BH_DIG  = sum(class_BH      == "DIG"),
      n_BH_NEUT = sum(class_BH      == "NEUTRAL"),
      n_stageR_DCG  = sum(class_stageR == "DCG"),
      n_stageR_DSG  = sum(class_stageR == "DSG"),
      n_stageR_DIG  = sum(class_stageR == "DIG"),
      n_stageR_NEUT = sum(class_stageR == "NEUTRAL")
    ) 
  
  list(
    joint_table     = joint_cls,
    summary_classes = summary_classes
  )
}

res <- run_pipeline(
  pval_cut = pval_cut
)

all_joint <- res$joint_table

all_class_summary <- res$summary_classes

# Plot: BH vs stageR class counts

data_flow <- all_joint %>%
  group_by(class_BH,class_stageR) %>%
  summarise(freq = n()) %>%
  ungroup()

data_ggforce <- data_flow  %>%
  gather_set_data(1:2) %>%        
  arrange(x,class_BH,desc(class_stageR))

data_ggforce$class_BH <- factor(data_ggforce$class_BH)
data_ggforce$class_stageR <- factor(data_ggforce$class_stageR)
data_ggforce$y <- factor(data_ggforce$y)

data_ggforce <- data_ggforce %>%
  group_by(class_BH, class_stageR) %>%
  mutate(y_mid = freq / 2) %>% 
  na.omit()

data_ggforce <- data_ggforce %>%
  mutate(
    x = factor(x, levels = c("class_BH", "class_stageR")),
    x_num = as.numeric(x)   # 1, 2
  )

p <- sankey_plot(
  data = data_ggforce,
  group_colors = gene_dosage_col,
  x_labels = c("class_BH", "class_stageR"),
)
p

ggsave("deconveilCaseStudies/revision/plots/sankey_sim1.png", dpi = 400, width = 4.0, height = 3.0, plot = p)
saveRDS(data_flow, file = "deconveilCaseStudies/revision/stageR/tr_table_sim1.RDS")
