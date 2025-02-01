setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/")
pkgs <- c("tidyverse")
sapply(pkgs, require, character.only = TRUE)
source("deconveilCaseStudies/utils.R")

### BRCA - Oncogenes & TSGs ###

# Data preprocessing

tumor_types <- c("BRCA")

base_paths <- list(
  res_pydeseq = "deconveilCaseStudies/results/{tumor}/res_CNnaive.csv",
  res_deconveil = "deconveilCaseStudies/results/{tumor}/res_CNaware.csv"
)

lfc_cut <- 1.0
pval_cut <- 0.05

process_results <- function(file_path, tumor_type, method) {
  read.csv(file_path) %>%
    mutate(
      isDE = (abs(log2FoldChange) >= lfc_cut) & (padj <= pval_cut),
      DEtype = if_else(!isDE, "n.s.", if_else(log2FoldChange > 0, "Up-reg", "Down-reg")),
      tumor_type = tumor_type,
      method = method
    ) %>%
    select(X, log2FoldChange, padj, isDE, DEtype, tumor_type, method) 
}


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


cancer_genes <- read.delim("TCGA/lung/cancerGeneList.tsv")

oncogenes <- cancer_genes %>% dplyr::filter(Is.Oncogene=="Yes") %>% 
  dplyr::select(Hugo.Symbol) %>% 
  dplyr::rename(geneID = Hugo.Symbol) %>% 
  dplyr::mutate(gene_type = "Oncogene")

tsg <- cancer_genes %>% dplyr::filter(Is.Tumor.Suppressor.Gene=="Yes") %>% 
  dplyr::select(Hugo.Symbol) %>% 
  dplyr::rename(geneID=Hugo.Symbol) %>% 
  dplyr::mutate(gene_type = "TSG")

cancer_genes_oncokb <- rbind(oncogenes, tsg)

extract_genes <- function(gene_group_df, gene_group_name) {
  gene_group_df %>%
    mutate(gene_group = gene_group_name) %>%
    select(geneID = geneID_aware, tumor_type = tumor_type_aware, gene_group)
}

cn_aware_d_sensitive <- extract_genes(gene_groups$d_sensitive, "Dosage-sensitive")
cn_aware_d_insensitive <- extract_genes(gene_groups$d_insensitive, "Dosage-insensitive")
cn_aware_d_compensated <- extract_genes(gene_groups$d_compensated, "Dosage-compensated")


# Perform intersection of DSGs, DIGs an DCGs with ONC/TSGs
dig_onc <- intersect(oncogenes$geneID, cn_aware_d_insensitive$geneID)
dsg_onc <- intersect(oncogenes$geneID, cn_aware_d_sensitive$geneID)
dcg_onc <- intersect(oncogenes$geneID, cn_aware_d_compensated$geneID)

dig_tsg <- intersect(tsg$geneID, cn_aware_d_insensitive$geneID)
dsg_tsg <- intersect(tsg$geneID, cn_aware_d_sensitive$geneID)
dcg_tsg <- intersect(tsg$geneID, cn_aware_d_compensated$geneID)

brca_onc <- data.frame(
  gene_category = c(rep("DIGs", length(dig_onc)),
                    rep("DSGs", length(dsg_onc)),
                    rep("DCGs", length(dcg_onc))),
  gene_type = "Oncogene",
  geneID = c(dig_onc, dsg_onc, dcg_onc)
)

brca_tsg <- data.frame(
  gene_category = c(rep("DIGs", length(dig_tsg)),
                    rep("DSGs", length(dsg_tsg)),
                    rep("DCGs", length(dcg_tsg))),
  gene_type = "TSG",
  geneID = c(dig_tsg, dsg_tsg, dcg_tsg)
)

onc_ts_genes <- rbind(brca_onc, brca_tsg)

distribution_summary <- onc_ts_genes %>%
  group_by(gene_category, gene_type) %>%
  summarise(gene_count = n(), .groups = "drop")

distribution_summary$gene_category <- factor(distribution_summary$gene_category, 
                                             levels = c("DSGs", "DIGs", "DCGs"))

barplot <- ggplot(distribution_summary, aes(x = gene_category, y = gene_count, fill = gene_type)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = gene_count), 
            position = position_stack(vjust = 0.5), 
            size = 4, color = "white") +  
  theme_classic() +
  labs(title = "",
       x = "Gene category",
       y = "Gene count",
       fill = "Gene type") +
  scale_fill_manual(values = c("Oncogene" = "#BB002198", "TSG" = "#00468BB2"))+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "right",
        axis.text.x = element_text(size = 16, color = "black"),  
        axis.text.y = element_text(size = 16, color = "black"),
        legend.key.size = unit(0.6, "cm"), 
        legend.text = element_text(size = 14, color = "black"),
        legend.title = element_text(size = 16, color = "black"),
        legend.spacing.y = unit(0.1, 'cm'),
        axis.text = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 14))
barplot

ggsave("deconveilCaseStudies/plots/main/barplot_onc_tsg_brca.png", dpi = 400, width = 5.5, height = 4.5, plot = barplot)

