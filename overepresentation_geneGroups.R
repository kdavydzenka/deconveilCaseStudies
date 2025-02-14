setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/")
pkgs <- c("tidyverse")
sapply(pkgs, require, character.only = TRUE)
source("deconveilCaseStudies/utils.R")

### Overrepresentation analysys of Gene categories: DSGs | DIGs | DCGs ###

# Data preprocessing

tumor_types <- c("LUAD", "LUSC", "BRCA") 

base_paths <- list(
  res_pydeseq = "deconveilCaseStudies/results/{tumor}/res_CNnaive.csv",
  res_deconveil = "deconveilCaseStudies/results/{tumor}/res_CNaware.csv"
)

# Define thresholds
lfc_cut <- 1.0
pval_cut <- 0.05

#Function to process differential expression results
process_results <- function(file_path, tumor_type, method) {
  read.csv(file_path) %>%
    mutate(
      isDE = (abs(log2FoldChange) >= lfc_cut) & (padj <= pval_cut),
      DEtype = if_else(!isDE, "n.s.", if_else(log2FoldChange > 0, "Up-reg", "Down-reg")),
      tumor_type = tumor_type,
      method = method
    ) %>%
    dplyr::select(X, log2FoldChange, padj, isDE, DEtype, tumor_type, method) 
}


# Process all tumor types
process_tumor_data <- function(tumor_types, base_paths) {
  res_pydeseq_list <- list()
  res_deconveil_list <- list()
  #cnv_list <- list()
  
  for (tumor in tumor_types) {
    # Replace placeholders in file paths
    pydeseq_path <- gsub("\\{tumor\\}", tumor, base_paths$res_pydeseq)
    deconveil_path <- gsub("\\{tumor\\}", tumor, base_paths$res_deconveil)
    
    # Process data
    res_pydeseq_list[[tumor]] <- process_results(pydeseq_path, tumor, "CN naive")
    res_deconveil_list[[tumor]] <- process_results(deconveil_path, tumor, "CN aware")
  }
  
  list(
    res_pydeseq = do.call(rbind, res_pydeseq_list),
    res_deconveil = do.call(rbind, res_deconveil_list)
  )
}

# Process data for all tumor types
all_data <- process_tumor_data(tumor_types, base_paths)
res_pydeseq_combined <- all_data$res_pydeseq
res_deconveil_combined <- all_data$res_deconveil

colnames(res_pydeseq_combined) <- c("geneID", "logFC", "padj", "isDE", "DEtype", "tumor_type", "method")
colnames(res_deconveil_combined) <- c("geneID", "logFC", "padj", "isDE", "DEtype", "tumor_type", "method")

# Combine results into a joint table
res_joint <- cbind(
  res_pydeseq_combined %>% rename_with(~ paste0(.x, "_naive")),
  res_deconveil_combined %>% rename_with(~ paste0(.x, "_aware"))
)

gene_groups <- define_gene_groups(res_joint)

prepare_cn_specific_data <- function(df, gene_group_name) {
  df %>%
    mutate(gene_group = gene_group_name) %>%
    dplyr::select(geneID = geneID_naive, log2FC = logFC_naive, padj = padj_naive, isDE = isDE_naive, DEtype = DEtype_naive, tumor_type = tumor_type_naive, method = method_naive, gene_group)
}


extract_genes <- function(gene_group_df, gene_group_name) {
  gene_group_df %>%
    mutate(gene_group = gene_group_name) %>%
    dplyr::select(geneID = geneID_aware, tumor_type = tumor_type_aware, gene_group)
}

d_sensitive <- extract_genes(gene_groups$d_sensitive, "Dosage-sensitive")
d_insensitive <- extract_genes(gene_groups$d_insensitive, "Dosage-insensitive")
d_compensated <- extract_genes(gene_groups$d_compensated, "Dosage-compensated")
non_deg <- extract_genes(gene_groups$non_deg, "non_deg")

luad_dsg <- d_sensitive[d_sensitive$tumor_type == "LUAD", ]
lusc_dsg <- d_sensitive[d_sensitive$tumor_type == "LUSC", ]
brca_dsg <- d_sensitive[d_sensitive$tumor_type == "BRCA", ]

luad_dig <- d_insensitive[d_insensitive$tumor_type == "LUAD", ]
lusc_dig <- d_insensitive[d_insensitive$tumor_type == "LUSC", ]
brca_dig <- d_insensitive[d_insensitive$tumor_type == "BRCA", ]

luad_dcg <- d_compensated[d_compensated$tumor_type == "LUAD", ]
lusc_dcg <- d_compensated[d_compensated$tumor_type == "LUSC", ]
brca_dcg <- d_compensated[d_compensated$tumor_type == "BRCA", ]


# Define background
background_data <- data.frame(
  LUAD = c(nrow(luad_dsg), nrow(luad_dig), nrow(luad_dcg)),  
  LUSC = c(nrow(lusc_dsg), nrow(lusc_dig), nrow(lusc_dcg)),
  BRCA = c(nrow(brca_dsg), nrow(brca_dig), nrow(brca_dcg)),
  row.names = c("DSGs", "DIGs", "DCGs")
)

proportions <- sweep(background_data, 2, colSums(background_data), "/")
proportions$mean <- rowMeans(proportions)
background_data$mean <- round(rowMeans(background_data))

# Shared genes
total_deg_mean <- round(mean(colSums(background_data)))

shared_genes <- data.frame(
  background_prop = proportions$mean,
  n_shared = c(33, 937, 30),
  total_geneCat = c(background_data$mean),
  row.names = c("DSGs", "DIGs", "DCGs")
)

shared_genes$shared_prop <- shared_genes$n_shared / shared_genes$total_geneCat
shared_genes$shared_prop <- shared_genes$n_shared / total_deg_mean

# Private genes
private_genes <- data.frame(
  LUAD = c(617, 1058, 766),  
  LUSC = c(694, 1895, 494),
  BRCA = c(879, 1930, 1068),
  row.names = c("DSGs", "DIGs", "DCGs")
)

proportions_private <- sweep(private_genes, 2, total_deg_mean, "/")
proportions_private$prop_mean <- rowMeans(proportions_private) 
proportions_private$background_prop <- proportions$mean

