setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/")
pkgs <- c("tidyverse", "ggvenn", "VennDiagram", "reactome.db", "fgsea", "org.Hs.eg.db", 
          "data.table", "clusterProfiler", "enrichplot", "ggpubr", "msigdbr")
sapply(pkgs, require, character.only = TRUE)

### Compare DSGs across each cancer type, identify common genes and cancer-specific genes ###

# Data preprocessing

tumor_types <- c("LUAD", "LUSC", "BRCA", "LIHC")
tumor_types <- c("BRCA")

base_paths <- list(
  res_pydeseq = "CN-aware-DGE/results/{tumor}/res_CNnaive.csv",
  res_deconveil = "CN-aware-DGE/results/{tumor}/res_CNaware.csv"
  #cnv = "TCGA/{tumor}/test/cnv_test_all_genes.csv"
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
    select(X, log2FoldChange, padj, isDE, DEtype, tumor_type, method) 
}

#read_cnv <- function(file_path) {
  #read.csv(file_path) %>%
    #remove_rownames() %>%
    #column_to_rownames(var = "X")
#}

# Process all tumor types
process_tumor_data <- function(tumor_types, base_paths) {
  res_pydeseq_list <- list()
  res_deconveil_list <- list()
  #cnv_list <- list()
  
  for (tumor in tumor_types) {
    # Replace placeholders in file paths
    pydeseq_path <- gsub("\\{tumor\\}", tumor, base_paths$res_pydeseq)
    deconveil_path <- gsub("\\{tumor\\}", tumor, base_paths$res_deconveil)
    #cnv_path <- gsub("\\{tumor\\}", tumor, base_paths$cnv)
    
    # Process data
    res_pydeseq_list[[tumor]] <- process_results(pydeseq_path, tumor, "CN naive")
    res_deconveil_list[[tumor]] <- process_results(deconveil_path, tumor, "CN aware")
    #cnv_list[[tumor]] <- read_cnv(cnv_path)
  }
  
  list(
    res_pydeseq = do.call(rbind, res_pydeseq_list),
    res_deconveil = do.call(rbind, res_deconveil_list)
    #cnv = cnv_list
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


define_gene_groups <- function(res_joint) {
  list(
    d_sensitive = res_joint %>%
      filter(
        (DEtype_naive == "Up-reg" & DEtype_aware == "n.s.") |
          (DEtype_naive == "Down-reg" & DEtype_aware == "n.s.")
      ),
    d_insensitive = res_joint %>%
      filter(
        (DEtype_naive == "Up-reg" & DEtype_aware == "Up-reg") |
          (DEtype_naive == "Down-reg" & DEtype_aware == "Down-reg")
      ),
    d_compensated = res_joint %>%
      filter(
        (DEtype_naive == "n.s." & DEtype_aware == "Up-reg") |
          (DEtype_naive == "n.s." & DEtype_aware == "Down-reg")
      ),
    non_deg = res_joint %>%
      filter(DEtype_naive == "n.s." & DEtype_aware == "n.s.")
  )
}

gene_groups <- define_gene_groups(res_joint)

prepare_cn_specific_data <- function(df, gene_group_name) {
  df %>%
    mutate(gene_group = gene_group_name) %>%
    select(geneID = geneID_naive, log2FC = logFC_naive, padj = padj_naive, isDE = isDE_naive, DEtype = DEtype_naive, tumor_type = tumor_type_naive, method = method_naive, gene_group)
}



# Extract gene categories for each tumor type

cn_aware_d_insensitive <- gene_groups$d_insensitive %>%
  mutate(gene_group = "Dosage-insensitive") %>%
  select(geneID = geneID_aware, log2FC = logFC_aware, padj = padj_aware, isDE = isDE_aware, DEtype = DEtype_aware, tumor_type = tumor_type_aware, method = method_aware, gene_group)

tumor_specific_genes <- cn_aware_d_insensitive %>%
  group_by(tumor_type) %>%
  summarise(gene_ids = list(unique(geneID)), .groups = "drop")

gene_lists <- tumor_specific_genes %>%
  with(setNames(gene_ids, tumor_type))

common_genes <- Reduce(intersect, gene_lists)

tumor_specific_only <- lapply(names(gene_lists), function(tumor) {
  setdiff(gene_lists[[tumor]], unlist(gene_lists[names(gene_lists) != tumor]))
}) %>%
  setNames(names(gene_lists))


# Plot genes overlap

ggvenn(gene_lists, 
       show_percentage = FALSE,  # Optional: Turn this off or on
       fill_color = c("#E6A9EC", "#A9CDE6", "#E6D7A9", "#A9E6B3", "#D4A9E6")) 


# Identify shared and cancer-specific oncogenes/TSGs
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

d_compensated_cancer_g <- cn_aware_d_compensated[cn_aware_d_compensated$geneID %in% cancer_genes_oncokb$geneID ,]
d_insensitive_cancer_g <- cn_aware_d_insensitive[cn_aware_d_insensitive$geneID %in% cancer_genes_oncokb$geneID ,]


cancer_gene_lists <- d_insensitive_cancer_g %>%
  group_by(tumor_type, gene_group) %>%
  summarise(geneIDs = list(geneID), .groups = "drop")

cancer_gene_lists <- cancer_gene_lists %>%
  with(setNames(geneIDs, tumor_type))

ggvenn(cancer_gene_lists, 
       show_percentage = FALSE,  #
       fill_color = c("#E6A9EC", "#A9CDE6", "#E6D7A9", "#A9E6B3")) 


tumor_specific_cancer_g <- lapply(names(cancer_gene_lists), function(tumor) {
  setdiff(cancer_gene_lists[[tumor]], unlist(cancer_gene_lists[-match(tumor, names(cancer_gene_lists))]))
})
names(tumor_specific_cancer_g) <- names(cancer_gene_lists)


common_cancer_g <- lapply(names(cancer_gene_lists), function(tumor) {
  setdiff(unlist(cancer_gene_lists), tumor_specific_cancer_g[[tumor]])
})
names(common_cancer_g) <- names(cancer_gene_lists)
common_cancer_g <- Reduce(intersect, common_cancer_g)


## Extract each gene category ##

d_sensitive_brca <- gene_groups$d_sensitive %>%
  mutate(gene_group = "Dosage-sensitive") %>%
  select(geneID = geneID_aware, log2FC = logFC_aware, padj = padj_aware, isDE = isDE_aware, DEtype = DEtype_aware, tumor_type = tumor_type_aware, method = method_aware, gene_group)

d_insensitive_brca <- gene_groups$d_insensitive %>%
  mutate(gene_group = "Dosage-insensitive") %>%
  select(geneID = geneID_aware, log2FC = logFC_aware, padj = padj_aware, isDE = isDE_aware, DEtype = DEtype_aware, tumor_type = tumor_type_aware, method = method_aware, gene_group)

d_compensated_brca <- gene_groups$d_compensated %>%
  mutate(gene_group = "Dosage-compensated") %>%
  select(geneID = geneID_aware, log2FC = logFC_aware, padj = padj_aware, isDE = isDE_aware, DEtype = DEtype_aware, tumor_type = tumor_type_aware, method = method_aware, gene_group)


# Overrepresentation GO functional enrichment analysis

d_sensitive_g <- d_sensitive_brca %>% dplyr::select(geneID, log2FC)
d_insensitive_g <- d_insensitive_brca %>% dplyr::select(geneID, log2FC)
d_compensated_g <- d_compensated_brca %>% dplyr::select(geneID, log2FC)


gene_groups <- list(
  "Dosage-sensitive" = d_sensitive_brca,
  "Dosage-compensated" = d_compensated_brca,
  "Dosage-insensitive" = d_insensitive_brca
)

perform_enrichment_GO <- function(gene_df, gene_group) {
  my_symbols <- gene_df$geneID
  gene_list <- AnnotationDbi::select(org.Hs.eg.db,
                                     keys = my_symbols,
                                     columns = c("ENTREZID", "SYMBOL"),
                                     keytype = "SYMBOL")
  gene_l <- as.vector(gene_list$ENTREZID)
  oraGO <- enrichGO(gene = gene_l, 
                    ont = "BP", 
                    OrgDb = org.Hs.eg.db, 
                    keyType = "ENTREZID", 
                    pvalueCutoff = 0.05, 
                    minGSSize = 10, 
                    maxGSSize = 350)
  
  oraGO@result %>% mutate(gene_group = gene_group)
}

perform_enrichment_H <- function(gene_df, gene_group) {
  my_symbols <- gene_df$geneID
  gene_list <- AnnotationDbi::select(org.Hs.eg.db,
                                     keys = my_symbols,
                                     columns = c("ENTREZID", "SYMBOL"),
                                     keytype = "SYMBOL")
  gene_l <- as.vector(gene_list$ENTREZID)
  m_hallmark <- msigdbr::msigdbr(species = "Homo sapiens", category = "H") 
  msig_H <- msigdbr(species = "Homo sapiens", category = "H") %>% dplyr::select(gs_name, entrez_gene)
  msig_H <- enricher(gene_l, minGSSize = 10, 
                     maxGSSize = 500,
                     pvalueCutoff = 0.05, 
                     pAdjustMethod = "BH", 
                     TERM2GENE = msig_H)
  
  msig_H@result %>% mutate(gene_group = gene_group)
}

res_ora_GO_list <- lapply(names(gene_groups), function(group_name) {
  gene_df <- gene_groups[[group_name]] %>% dplyr::select(geneID, log2FC)
  perform_enrichment_GO(gene_df, group_name)
})
names(res_ora_GO_list) <- names(gene_groups)

res_GO_d_sensitive <- res_ora_GO_list[["Dosage-sensitive"]] %>% as.data.frame()
res_GO_d_insensitive <- res_ora_GO_list[["Dosage-insensitive"]] %>% as.data.frame()
res_GO_d_compensated <- res_ora_GO_list[["Dosage-compensated"]] %>% as.data.frame()

saveRDS(res_ora_GO_list, file = "CN-aware-DGE/results/pancancer_ora_GO/lihc.RDS")


res_ora_H_list <- lapply(names(gene_groups), function(group_name) {
  gene_df <- gene_groups[[group_name]] %>% dplyr::select(geneID, log2FC)
  perform_enrichment_H(gene_df, group_name)
})
names(res_ora_H_list) <- names(gene_groups)

res_H_d_sensitive <- res_ora_H_list[["Dosage-sensitive"]] %>% as.data.frame()
res_H_d_insensitive <- res_ora_H_list[["Dosage-insensitive"]] %>% as.data.frame()
res_H_d_compensated <- res_ora_H_list[["Dosage-compensated"]] %>% as.data.frame()


# Select 5 interesting significant pathways for each gene category

GO_d_sensitive <- c("DNA replication", "positive regulation of DNA metabolic process", "protein folding", "amino acid metabolic process", 
                    "DNA recombination")

GO_d_insensitive <- c("immunoglobulin mediated immune response", "B cell mediated immunity", "mitotic nuclear division",
                      "negative regulation of cell adhesion", "positive regulation of inflammatory response")

GO_d_compensated <- c("positive regulation of T cell activation", "positive regulation of lymphocyte activation", "Ras protein signal transduction",
                      "endothelial cell migration", "T cell differentiation")

res_GO_d_sensitive <- res_GO_d_sensitive %>% filter(Description %in% GO_d_sensitive)
res_GO_d_compensated <- res_GO_d_compensated %>% filter(Description %in% GO_d_compensated)
res_GO_d_insensitive <- res_GO_d_insensitive%>% filter(Description %in% GO_d_insensitive)

res_GO_d_sensitive <- res_GO_d_sensitive %>%
  separate(GeneRatio, into = c("numerator", "denominator"), sep = "/", convert = TRUE) %>%
  mutate(GeneRatio_val = numerator / denominator)

res_GO_d_compensated <- res_GO_d_compensated %>%
  separate(GeneRatio, into = c("numerator", "denominator"), sep = "/", convert = TRUE) %>%
  mutate(GeneRatio_val = numerator / denominator)

res_GO_d_insensitive <- res_GO_d_insensitive %>%
  separate(GeneRatio, into = c("numerator", "denominator"), sep = "/", convert = TRUE) %>%
  mutate(GeneRatio_val = numerator / denominator)

res_GO_d_sensitive$log_padjust <- -log10(res_GO_d_sensitive$p.adjust)
res_GO_d_insensitive$log_padjust <- -log10(res_GO_d_insensitive$p.adjust)
res_GO_d_compensated$log_padjust <- -log10(res_GO_d_compensated$p.adjust)

p_data <- rbind(res_GO_d_sensitive, res_GO_d_compensated, res_GO_d_insensitive)

p_data$gene_group <- as.factor(p_data$gene_group)

p_data <- p_data %>%
  group_by(gene_group) %>%
  arrange(GeneRatio_val, .by_group = TRUE) %>%
  ungroup()

p_data$Description <- factor(p_data$Description, levels = unique(p_data$Description))


gene_group_colors = c("Dosage-sensitive" = "#FF7F0E", "Dosage-compensated" = "#f9a729", "Dosage-insensitive" = "#8F3931FF")

p_gse <- ggplot(p_data, aes(x = GeneRatio_val, y = Description)) +
  geom_point(aes(size = Count, color = log_padjust))+
  scale_color_gradient(low = "blue", high = "orange")+
  labs(x = "gene ratio", y = "", title = "",
       color = "padj (log10)",  
       size = "Gene count"           
  )+
  facet_wrap(~factor(gene_group, levels = c("Dosage-sensitive", "Dosage-insensitive", "Dosage-compensated")), nrow = 1)+
  theme(strip.text.x = element_text(size=10, color="black", face="bold.italic"))+
  theme_bw()+
  theme(legend.position = "right")+
  scale_x_continuous(limits = c(0, 0.07))+
  #theme(axis.text.y = element_text(color = gene_group_colors[p_data$gene_group]))+
  theme(axis.text.x = element_text(size = 16, color = "black"),  
        axis.text.y = element_text(size = 16, color = "black"),
        legend.key.size = unit(0.6, "cm"), 
        legend.text = element_text(size = 14, color = "black"),
        legend.title = element_text(size = 16, color = "black"),
        legend.spacing.y = unit(0.2, 'cm'),
        strip.text = element_text(size = 18, face = "plain", color = "black"),
        axis.text = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 16)
  )
p_gse

ggsave("CN-aware-DGE/plots/main/GO_dotplot_luad.png", dpi = 400, width = 14.0, height = 6.0, plot = p_gse)



# Find LUAD specific cancer genes in biological pathways

convert_entrez_to_symbol <- function(entrez_ids) {
  entrez_list <- unlist(strsplit(entrez_ids, split = "/"))
  gene_symbols <- mapIds(org.Hs.eg.db, entrez_list, "SYMBOL", "ENTREZID")
  return(paste(gene_symbols, collapse = ", "))
}

res_GO_d_sensitive$geneSymbol <- sapply(res_GO_d_sensitive$geneID, convert_entrez_to_symbol)

tumor_specific_luad <- c("CDKN2B", "IRS1", "CDK4", "MET", "ZNF217", "FURIN", "CD74", "FOXA1", "PAK1", "PRKDC", "ESR1", "ACTG1")
GO_pathways <- c("fibroblast proliferation", "positive regulation of DNA metabolic process", "tRNA metabolic process",
                 "protein folding", "positive regulation of DNA repair")


res_H_d_sensitive$geneSymbol <- sapply(res_H_d_sensitive$geneID, convert_entrez_to_symbol)
H_pathways_d_sensitive <- c("HALLMARK_MYC_TARGETS_V2", "HALLMARK_E2F_TARGETS", "HALLMARK_GLYCOLYSIS", "HALLMARK_G2M_CHECKPOINT")
res_H_d_sensitive <- res_H_d_sensitive %>% filter(Description %in% H_pathways_d_sensitive)

res_H_d_sensitive <- res_H_d_sensitive %>%
  separate(GeneRatio, into = c("numerator", "denominator"), sep = "/", convert = TRUE) %>%
  mutate(GeneRatio_val = numerator / denominator)

res_H_d_sensitive$log_padjust <- -log10(res_H_d_sensitive$p.adjust)
res_H_d_sensitive$Description <- gsub("HALLMARK", "H", res_H_d_sensitive$Description)
res_H_d_sensitive <- res_H_d_sensitive %>% arrange(p.adjust)
res_H_d_sensitive$Description <- factor(res_H_d_sensitive$Description, levels = unique(res_H_d_sensitive$Description))

p_hallmark <- ggplot(res_H_d_sensitive, aes(x = log_padjust, y = Description))+
  geom_bar(stat = "identity", fill = "#BF616AFF", width = 0.3) + 
  scale_x_continuous() +
  scale_y_discrete(expand = expansion(mult = c(0.2, 0.2))) +
  labs(title = "", x = "-log10(FDR)", y = "") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 16, color = "black"),  
        axis.text.y = element_text(size = 16, color = "black"),
        legend.key.size = unit(0.6, "cm"), 
        legend.text = element_text(size = 14, color = "black"),
        legend.title = element_text(size = 16, color = "black"),
        legend.spacing.y = unit(0.2, 'cm'),
        strip.text = element_text(size = 18, face = "plain", color = "black"),
        axis.text = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 16)
  )
p_hallmark
ggsave("CN-aware-DGE/plots/main/H_barchart_luad.png", dpi = 400, width = 5.0, height = 3.5, plot = p_hallmark)

  