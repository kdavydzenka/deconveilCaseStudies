setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/")
pkgs <- c("cluster", "factoextra", "heatmaply", "DESeq2", "tidyverse", "colorspace", 
          "ggpubr", "ggpointdensity", "ggeasy", "gridExtra", "metaseqR2", "ggalluvial", "ggridges", "ggforce", "ggparallel", "alluvial",
          "patchwork")
sapply(pkgs, require, character.only = TRUE)

source("deconveilCaseStudies/utils/utils.R")
source("deconveilCaseStudies/utils/utils_plots.R")

# Input data
tumor_types <- c("KIRC") # Main
tumor_types <- c("LUAD", "LUSC", "BRCA", "LIHC", "KIRC") # Supplementary
lfc_cut = 0.2
pval_cut = 0.05 
loss_threshold = 0.25

gene_group_colors <- c("DIGs" = "#8F3931FF", 
                       "DSGs" = "#FFB977", 
                       "DCGs"="#FAE48BFF", 
                       "non-DEGs" = "#ADB6B6FF")  

cnv_colors <- c("loss" = "#0073C299", 
                "neutral" = "#86868699", 
                "gain" = "#cecb76", 
                "amplification" = "#DC0000B2")


### Downstream analysis ###

# Read data files for a tumor type
read_tcga_data <- function(tumor_type) {
  list(
    res_pydeseq = read.csv(paste0("deconveilCaseStudies/results_tcga/", tumor_type, "/new/res_CNnaive.csv")),
    res_deconveil = read.csv(paste0("deconveilCaseStudies/results_tcga/", tumor_type, "/new/res_CNaware.csv")),
    cnv_tumor = read.csv(paste0("deconveilCaseStudies/results_tcga/", tumor_type, "/old/cnv_tumor.csv")) %>%
      remove_rownames() %>%
      column_to_rownames(var = "X") * 2
  )
}

process_cnv_data <- function(cnv_tumor, loss_threshold = 0.25) {
  # cnv_tumor: data.frame / tibble with rownames = geneID, columns = tumor samples
  cn_mat <- as.matrix(cnv_tumor)
  cn_mean <- rowMeans(cn_mat, na.rm = TRUE)
  
  # fraction of samples with CN < 2 (loss) and > 2 (gain)
  frac_loss <- rowMeans(cn_mat < 2, na.rm = TRUE)
  frac_gain <- rowMeans(cn_mat > 2, na.rm = TRUE)
  
  loss_label <- case_when(
    frac_loss > loss_threshold & frac_gain < loss_threshold ~ "LOSS",
    frac_gain > loss_threshold & frac_loss < loss_threshold ~ "GAIN",
    TRUE ~ "DIPLOID/MIXED"
  )
  
  tibble(
    geneID     = rownames(cnv_tumor),
    cnv_mean   = cn_mean,
    frac_loss  = frac_loss,
    frac_gain  = frac_gain,
    cnv_label  = loss_label
  )
}

# Generic annotator: one method (naive OR aware)
annotate_de_results <- function(res_df,
                                method_label = c("naive","aware"),
                                tumor_type,
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
  
  # EXTRACT CORRECT stageR COLUMN 
  # if missing, fill with NA_real_
  
  if (method_label == "naive") {
    DE_confirmed_vec <- if ("DE_confirmed" %in% colnames(res_df))
      res_df$DE_confirmed else rep(NA_real_, nrow(res_df))
  } else {
    DE_confirmed_vec <- if ("DE_confirmed" %in% colnames(res_df))
      res_df$DE_confirmed else rep(NA_real_, nrow(res_df))
  }
  
  # MUTATE WITH DE CALLS 
  res_df <- res_df %>%
    mutate(
      DE_confirmed = DE_confirmed_vec,
      tumor_type  = tumor_type,
      method      = ifelse(method_label == "naive", "CN-naive", "CN-aware"),
      
      # BH-based significance 
      isDE_BH = !is.na(padj_BH) & padj_BH < pval_cut & abs(log2FC) > lfc_cut,
      DEtype_BH = case_when(
        !isDE_BH       ~ "n.s.",
        log2FC > 0     ~ "Up-reg",
        log2FC < 0     ~ "Down-reg",
        TRUE           ~ "n.s."
      ),
      
      # stageR-based significance
      isDE_stageR = !is.na(DE_confirmed) & DE_confirmed == 1.0 & abs(log2FC) > lfc_cut,
      DEtype_stageR = case_when(
        !isDE_stageR   ~ "n.s.",
        log2FC > 0     ~ "Up-reg",
        log2FC < 0     ~ "Down-reg",
        TRUE           ~ "n.s."
      )
    )
  
  # ADD SUFFIXES TO PREPARE NAIVE + AWARE MERGE 
  res_df <- res_df %>%
    dplyr::rename(
      !!paste0("log2FC",        suf) := log2FC,
      !!paste0("padj_BH",       suf) := padj_BH,
      !!paste0("DE_confirmed",  suf) := DE_confirmed,
      !!paste0("isDE_BH",       suf) := isDE_BH,
      !!paste0("isDE_stageR",   suf) := isDE_stageR,
      !!paste0("DEtype_BH",     suf) := DEtype_BH,
      !!paste0("DEtype_stageR", suf) := DEtype_stageR,
      !!paste0("tumor_type",    suf) := tumor_type,
      !!paste0("method",        suf) := method
    )
  
  # SELECT CLEAN OUTPUT 
  res_df %>%
    dplyr::select(
      geneID,
      dplyr::starts_with("log2FC"),
      dplyr::starts_with("padj_"),
      dplyr::starts_with("DE_"),
      dplyr::starts_with("isDE_"),
      dplyr::starts_with("DEtype_"),
      dplyr::starts_with("tumor_type"),
      dplyr::starts_with("method")
    )
}

#  Combine naive + aware + CN info

combine_de_and_cnv <- function(naive_annot, aware_annot, cnv_summary) {
  res_joint <- naive_annot %>%
    inner_join(aware_annot, by = "geneID") %>%
    inner_join(cnv_summary, by = "geneID")
  
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

run_tcga_deconveil_for_tumor <- function(tumor_type,
                                         lfc_cut = 0.2,
                                         pval_cut = 0.05,
                                         loss_threshold = 0.25) {
  
  message(">>> Processing tumor type: ", tumor_type)
  
  dat <- read_tcga_data(tumor_type)  # res_naive, res_aware, cnv_tumor
  
  cnv_summary <- process_cnv_data(dat$cnv_tumor, loss_threshold = loss_threshold)
  
  res_naive_annot <- annotate_de_results(dat$res_pydeseq,
                                         method_label = "naive",
                                         tumor_type = tumor_type,
                                         lfc_cut = lfc_cut,
                                         pval_cut = pval_cut)
  
  res_aware_annot <- annotate_de_results(dat$res_deconveil,
                                         method_label = "aware",
                                         tumor_type = tumor_type,
                                         lfc_cut = lfc_cut,
                                         pval_cut = pval_cut)
  
  joint <- combine_de_and_cnv(res_naive_annot, res_aware_annot, cnv_summary)
  
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
    ) %>%
    mutate(tumor_type = tumor_type) %>%
    relocate(tumor_type)
  
  list(
    tumor_type      = tumor_type,
    joint_table     = joint_cls,
    summary_classes = summary_classes
  )
}

# Run across tumor types

tcga_results <- map(tumor_types, ~ run_tcga_deconveil_for_tumor(
  tumor_type      = .x,
  lfc_cut         = lfc_cut,
  pval_cut        = pval_cut,
  loss_threshold  = loss_threshold
))

# Combined per-gene table (all tumor types)
all_joint <- tcga_results %>%
  map("joint_table") %>%
  bind_rows()

# Combined class-summary table
all_class_summary <- tcga_results %>%
  map("summary_classes") %>%
  bind_rows()


# Plot: BH vs stageR class counts

# choose colors for gene groups

colors <- c("DCG" = "#FAE48BFF",
            "DIG" = "#8F3931FF",
            "DSG" = "#FFB977",
            "NEUTRAL" = "#ADB6B6FF")


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
  #x_var = "x_num",
  #fill_var = "y",
  group_colors = colors,
  x_labels = c("class_BH", "class_stageR"),
  #title = "Gene-class reassignment: LUAD"
)
p

ggsave("deconveilCaseStudies/plots/sankey_kirc.png", dpi = 400, width = 4.0, height = 3.0, plot = p)
saveRDS(data_flow, file = "deconveilCaseStudies/revision/stageR/tr_table_kirc.RDS")

add_axis_labels <- function(p, data, x_labels) {
  
  axis_df <- data %>%
    group_by(x_num, class) %>%
    summarise(freq = sum(freq), .groups = "drop")
  
  p +
    geom_text(
      data = axis_df,
      aes(
        x     = x_num - 0.12,   # shift toward left of axis
        y     = freq / 2,       # middle of block
        label = freq
      ),
      size = 3.5,
      fontface = "bold"
    )
}

p_with_counts <- add_axis_labels(p, data_ggforce, x_labels = c("BH FDR", "stageR"))
p_with_counts
