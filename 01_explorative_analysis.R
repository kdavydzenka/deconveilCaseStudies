### Explortive analysis of RNA (z-score) and  CNV  ###

setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/")
pkgs <- c("dplyr", "ggplot2", "cluster", "factoextra", "heatmaply", "DESeq2", "tidyverse", "DESeq2", "colorspace", 
          "ggpubr", "ggpointdensity", "ggeasy", "patchwork")
invisible(sapply(pkgs, require, character.only = TRUE))

source("deconveilCaseStudies/utils/utils.R")


## Input data CN | Datasets analysed LUAD | BRCA | LIHC | LUSC ##

# Define paths
cancer_types <- c("LUAD", "BRCA", "LIHC", "LUSC")
data_paths <- paste0("TCGA/", c("lung/LUAD", "BRCA", "LIHC", "LUSC"), "/cnv_tumor.RDS")
names(data_paths) <- cancer_types

cancer_specific_cols <- list(
  LUAD = c(7, 5, 13, 14, 15, 18, 22, 24, 25, 28)
)

# Cluster IDs to retain
cancer_clusters <- list(
  BRCA = c(3, 4),
  LUSC = c(1, 3),
  LIHC = c(1, 2),
  LUAD = NULL # No clustering/filtering for LUAD in this example
)

load_and_filter_cnv <- function(data_path, cols = NULL, clip_value = 15) {
  cnv_data <- readRDS(data_path) %>% as.data.frame()
  if (!is.null(cols)) cnv_data <- cnv_data[, cols]
  cnv_data <- apply(cnv_data, 2, function(x) pmin(x, clip_value))
  return(cnv_data)
}


process_cnv <- function(cancer_type) {
  message("Processing: ", cancer_type)
  path <- data_paths[[cancer_type]]
  dataset_name <- paste0(cancer_type, "_cnv")
  cols <- cancer_specific_cols[[cancer_type]]
  
  if (!is.null(cols)) {
    # Only LUAD currently has selected columns
    cnv_data <- load_and_filter_cnv(path, cols)
  } else {
    set.seed(42)
    clustering_res <- clustering_patients_cnv(dataset_name, path)
    cnv_df <- clustering_res[[1]] %>% as.data.frame()
    cluster_ids <- cancer_clusters[[cancer_type]]
    
    if (!is.null(cluster_ids)) {
      cnv_df <- cnv_df %>% filter(Cluster %in% cluster_ids)
    }
    # Transpose and clip
    cnv_data <- t(as.matrix(cnv_df))
    cnv_data <- apply(cnv_data, 2, function(x) pmin(x, 15))
  }
  
  # Compute gene-wise mean
  cnv_mean <- as.data.frame(cnv_data) %>%
    mutate(cnv_mean = rowMeans(.)) %>%
    select(cnv_mean)
  
  cnv_mean$geneID <- rownames(cnv_mean)
  
  # Return both mean and filtered CNV data
  return(list(
    cnv_mean = cnv_mean,
    cnv_filt = cnv_data
  ))
}

cnv_results <- lapply(cancer_types, process_cnv)
names(cnv_results) <- cancer_types


## Input data RNA ##

cancer_types <- c("LUAD", "BRCA", "LIHC", "LUSC")
rna_data_paths <- paste0("TCGA/", c("lung/LUAD", "BRCA", "LIHC", "LUSC"), "/rna_counts.RDS")
names(rna_data_paths) <- cancer_types

low_expression_threshold <- 20

process_rna_data <- function(cancer_type, cnv_filt = NULL, low_expression_threshold = 20) {
  message("Processing RNA for: ", cancer_type)
  
  data_path <- rna_data_paths[[cancer_type]]
  dataset_name <- paste0(cancer_type, "_rna")
  
  # Load and process RNA data 
  rna_list <- rna_processing(dataset_name, data_path, cnv_filt)
  rna_norm <- rna_list[[1]]  # Normal 
  rna_tum <- rna_list[[2]]   # Tumor 
  
  # Check for valid data
  if (is.null(rna_norm) || is.null(rna_tum)) {
    warning("Missing RNA data for ", cancer_type)
    return(NULL)
  }
  
  # Filter out low-expression genes based on normal tissue
  gene_means <- rowMeans(rna_norm)
  filtered_genes <- names(gene_means[gene_means > low_expression_threshold])
  
  if (length(filtered_genes) == 0) {
    warning("No genes passed expression threshold for ", cancer_type)
    return(NULL)
  }
  
  # Subset filtered genes
  rna_norm <- rna_norm[filtered_genes, ]
  rna_tum <- rna_tum[filtered_genes, ]
  rna_combined <- cbind(rna_norm, rna_tum)
  
  # DESeq2 normalization 
  conditions <- rep("A", ncol(rna_combined))
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = rna_combined,
    colData = DataFrame(condition = conditions),
    design = ~1
  )
  rna_vst <- DESeq2::vst(dds, blind = TRUE)
  rna_log_normalized <- assay(rna_vst)
  
  # Z-score transformation
  rna_zscore <- t(scale(t(rna_log_normalized)))
  
  # Select tumor columns 
  n_samples <- ncol(rna_zscore)
  half_index <- floor(n_samples / 2) + 1
  
  rna_zscore_tumor <- rna_zscore[, half_index:n_samples]
  
  # Compute mean z-score per gene
  zscore_summary <- as.data.frame(rna_zscore_tumor) %>%
    mutate(zscore_mean = rowMeans(.)) %>%
    select(zscore_mean)
  
  zscore_summary$geneID <- rownames(zscore_summary)
  
  # Return as list
  return(list(
    zscore_matrix = rna_zscore_tumor,
    zscore_summary = zscore_summary,
    vst_matrix = rna_log_normalized
  ))
}

rna_results <- lapply(cancer_types, function(cancer) {
  cnv_filt <- if (!is.null(cnv_results) && !is.null(cnv_results[[cancer]])) cnv_results[[cancer]][["cnv_filt"]] else NULL
  process_rna_data(cancer, cnv_filt)
})
names(rna_results) <- cancer_types

# Plot

prepare_plot_df <- function(cancer) {
  rna_summary <- rna_results[[cancer]]$zscore_summary
  cnv <- cnv_results[[cancer]]$cnv_mean
  
  # Ensure cnv is a data frame with rownames
  cnv_df <- as.data.frame(cnv)
  if (!"geneID" %in% colnames(cnv_df)) {
    cnv_df$geneID <- rownames(cnv_df)
  }
  
  # Binning CNV values
  cnv_binned <- cnv_df %>% 
    mutate(cnv = case_when(
      cnv_mean > 0.0 & cnv_mean <= 1.8 ~ "1",
      cnv_mean > 1.8 & cnv_mean <= 2.5 ~ "2",
      cnv_mean > 2.5 & cnv_mean <= 3.5 ~ "3",
      cnv_mean > 3.5 & cnv_mean <= 4.5 ~ "4",
      cnv_mean > 4.5 ~ "5"
    )) %>%
    select(geneID, cnv)
  
  # Merge with RNA z-score summary
  df <- rna_summary %>%
    rownames_to_column("gene") %>%
    merge(cnv_binned, by.x = "gene", by.y = "geneID") %>%
    na.omit() %>%
    mutate(cancer_type = cancer)
  
  return(df)
}

p_luad_t <- prepare_plot_df("LUAD")
p_brca_t <- prepare_plot_df("BRCA")
p_lusc_t <- prepare_plot_df("LUSC")
p_lihc_t <- prepare_plot_df("LIHC")

p_tumor <- bind_rows(p_brca_t, p_lusc_t, p_lihc_t)

saveRDS(p_luad_t, file = "deconveilCaseStudies/plots/main/Fig 1/rds/plot_boxplot.rds")
saveRDS(p_tumor, file = "deconveilCaseStudies/plots/supplementary/rds/plot_boxplots_tumor_types.rds")

p_tumor$cnv <- factor(p_tumor$cnv, levels = c("1", "2", "3", "4", "5"))

col <- c("1" = "dodgerblue1", "2" = "darkgray", "3" = "green4", "4" = "coral3", "5" = "hotpink3")

bxp_t <- ggplot(p_tumor, aes(x = cnv, y = zscore_mean, color = cnv)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch = F, show.legend = F)+
  #geom_smooth(method = "loess", formula = y ~ x, se=FALSE, color="darkblue", aes(group=1), linetype = 'dashed')+
  labs(x="CN group", y = "mRNA Z-score")+
  facet_wrap(~cancer_type)+
  theme(strip.text.x = element_text(size=14, color="black", face="bold.italic"))+
  ggplot2::theme(legend.position = 'none')+
  theme_bw()+
  scale_color_manual(values=col)+
  font("xy.text", size = 12, color = "black", face = "plain")+
  font("title", size = 12, color = "black")+
  font("xlab", size = 12)+
  font("ylab", size = 12)+
  ggtitle("")+
  theme(plot.title = element_text(hjust = 0.5),
        strip.text = element_text(size = 12, face = "plain", color = "black"))+
  theme(legend.position = "")
bxp_t

ggsave("deconveilCaseStudies/plots/supplementary/boxplots.png", dpi = 500, width = 7.5, height = 4.0, plot = bxp_t)

