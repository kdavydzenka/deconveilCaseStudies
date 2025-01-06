setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/")
pkgs <- c("tidyverse", "ggvenn", "VennDiagram", "reactome.db", "fgsea", "org.Hs.eg.db", 
          "data.table", "clusterProfiler", "enrichplot", "ggpubr", "msigdbr", "ggrepel")
sapply(pkgs, require, character.only = TRUE)

### Check whether DCGs/DIGs overlap significantly with the targets of known cancer specific relevant TFs ###

res_pydeseq <- read.csv("CN-aware-DGE/results/case_studies/LUSC/res_CNnaive.csv")
res_deconveil <- read.csv("CN-aware-DGE/results/case_studies/LUSC/res_CNaware.csv")
tf <- read_delim("CN-aware-DGE/hg38_disease_TFs.txt")

# Data preprocesing

tf_lung <- tf[tf$disease %in% c("LUSC", "lung|carcinoma", "LUSC-KR", "LUSC-CN", "Squamous_cell_lung_carcinoma",
                                "Non-small_cell_lung_cancer", "Lung_cancer") ,]
tf_lung <- tf_lung[!duplicated(tf_lung$symbol), ]
tf_lung <- tf_lung %>%
  filter(!grepl("^ENSG", symbol))

tf_lung_filt <- tf_lung %>% dplyr::select(symbol, gene_type) %>% 
  dplyr::rename(geneID = symbol)

tumor_types <- c("LUSC")
base_paths <- list(
  res_pydeseq = "CN-aware-DGE/results/{tumor}/res_CNnaive_all_genes.csv",
  res_deconveil = "CN-aware-DGE/results/{tumor}/res_CNaware_all_genes.csv"
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


# Extract DCGs and DIGs
cn_aware_d_insensitive <- gene_groups$d_insensitive %>%
  mutate(gene_group = "Dosage-insensitive") %>%
  select(geneID = geneID_aware, log2FC = logFC_aware, padj = padj_aware, isDE = isDE_aware, DEtype = DEtype_aware, tumor_type = tumor_type_aware, method = method_aware, gene_group)

cn_aware_d_compensated <- gene_groups$d_compensated %>%
  mutate(gene_group = "Dosage-compensated") %>%
  select(geneID = geneID_aware, log2FC = logFC_aware, padj = padj_aware, isDE = isDE_aware, DEtype = DEtype_aware, tumor_type = tumor_type_aware, method = method_aware, gene_group)

cn_aware_d_sensitive <- gene_groups$d_sensitive %>%
  mutate(gene_group = "Dosage-sensitive") %>%
  select(geneID = geneID_aware, log2FC = logFC_aware, padj = padj_aware, isDE = isDE_aware, DEtype = DEtype_aware, tumor_type = tumor_type_aware, method = method_aware, gene_group)


  
# Perform overlap with TFs lung cancer specific gene list

common_d_insensitive <- intersect(cn_aware_d_insensitive$geneID, tf_lung_filt$geneID)
common_d_compensated <- intersect(cn_aware_d_compensated$geneID, tf_lung_filt$geneID)
common_d_sensitive <- intersect(cn_aware_d_sensitive$geneID, tf_lung_filt$geneID)

cancer_gene_lists <- list(
  "Dosage-insensitive" = cn_aware_d_insensitive$geneID,
  "Dosage-compensated" = cn_aware_d_compensated$geneID,
  "TF target" = tf_lung_filt$geneID
)

gene_group_colors <- c("Dosage-insensitive" = "#8F3931FF", "Dosage-compensated"="#FAE48BFF", "TF targets" = "lightcoral")  

ggvenn(
  cancer_gene_lists,
  show_percentage = F,  
  fill_color = c("#8F3931FF", "#FAE48BFF", "lightcoral"), 
  stroke_size = 0.7,       
  set_name_size = 6       
)



# Overrepresentation GO functional enrichment analysis

d_sensitive_g <- cn_aware_d_sensitive %>% dplyr::select(geneID, log2FC)
d_insensitive_g <- cn_aware_d_insensitive %>% dplyr::select(geneID, log2FC)
d_compensated_g <- cn_aware_d_compensated %>% dplyr::select(geneID, log2FC)


gene_groups <- list(
  "Dosage-sensitive" = cn_aware_d_sensitive,
  "Dosage-compensated" = cn_aware_d_compensated,
  "Dosage-insensitive" = cn_aware_d_insensitive
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

res_ora_GO_list <- lapply(names(gene_groups), function(group_name) {
  gene_df <- gene_groups[[group_name]] %>% dplyr::select(geneID, log2FC)
  perform_enrichment_GO(gene_df, group_name)
})
names(res_ora_GO_list) <- names(gene_groups)

res_GO_d_sensitive <- res_ora_GO_list[["Dosage-sensitive"]] %>% as.data.frame()
res_GO_d_insensitive <- res_ora_GO_list[["Dosage-insensitive"]] %>% as.data.frame()
res_GO_d_compensated <- res_ora_GO_list[["Dosage-compensated"]] %>% as.data.frame()

GO_path_sensitive <- c("ribosome biogenesis", "protein-RNA complex organization", "protein processing",
                       "regulation of chromosome organization", "positive regulation of telomere maintenance")

GO_path_compensated <- c("positive regulation of lymphocyte activation", "negative regulation of type II interferon production", 
                         "negative regulation of immune effector process", "negative regulation of cytokine production involved in immune response",
                         "positive regulation of leukocyte differentiation")

GO_path_insensitive <- c("regulation of Wnt signaling pathway", "negative regulation of cell adhesion",
                         "mesenchymal cell differentiation", "platelet activation",  "T cell mediated immunity")


res_GO_sensitive <- res_GO_d_sensitive %>% filter(Description %in% GO_path_sensitive)
res_GO_compensated <- res_GO_d_compensated %>% filter(Description %in% GO_path_compensated)
res_GO_insensitive <- res_GO_d_insensitive%>% filter(Description %in% GO_path_insensitive)

res_GO_sensitive <- res_GO_sensitive %>%
  separate(GeneRatio, into = c("numerator", "denominator"), sep = "/", convert = TRUE) %>%
  mutate(GeneRatio_val = numerator / denominator)

res_GO_compensated <- res_GO_compensated %>%
  separate(GeneRatio, into = c("numerator", "denominator"), sep = "/", convert = TRUE) %>%
  mutate(GeneRatio_val = numerator / denominator)

res_GO_insensitive <- res_GO_insensitive %>%
  separate(GeneRatio, into = c("numerator", "denominator"), sep = "/", convert = TRUE) %>%
  mutate(GeneRatio_val = numerator / denominator)

res_GO_sensitive$log_padjust <- -log10(res_GO_sensitive$p.adjust)
res_GO_insensitive$log_padjust <- -log10(res_GO_insensitive$p.adjust)
res_GO_compensated$log_padjust <- -log10(res_GO_compensated$p.adjust)

convert_entrez_to_symbol <- function(entrez_ids) {
  entrez_list <- unlist(strsplit(entrez_ids, split = "/"))
  gene_symbols <- mapIds(org.Hs.eg.db, entrez_list, "SYMBOL", "ENTREZID")
  return(paste(gene_symbols, collapse = ", "))
}

res_GO_compensated$geneSymbol <- sapply(res_GO_compensated$geneID, convert_entrez_to_symbol)
common_d_compensated <- as.data.frame(common_d_compensated)


# Evidence immune-related genes withing Dosage-compensated gene category

immune_related_genes <- c("FCRL3", "IL10", "GATA3", "PCK1", "CD5", "CD59", "SMARCC1", "KLF10",
                         "BCL6", "PTPN6", "AXL", "INPP5D", "NFKBID", "NR1H4", "CD47")

gene_group_colors <- c("Dosage-insensitive" = "#8F3931FF", "Dosage-sensitive" = "#FFB977", "Dosage-compensated"="#FAE48BFF", "non-DEG" = "#ADB6B6FF")  

p_volcanos <- cn_aware_d_compensated %>%
  ggplot(mapping = aes(x = log2FC, y = -log10(padj))) +
  geom_point(aes(col = gene_group), size = 2.0, alpha = 0.7) +
  scale_color_manual(values = gene_group_colors) +
  geom_label_repel(
    data = cn_aware_d_compensated %>% filter(geneID %in% immune_related_genes),
    aes(label = geneID),
    size = 4.0,               
    fontface = "bold",        
    color = "white",          
    fill = "#1B1919B2",       
    box.padding = 0.8,        
    point.padding = 0.5,      
    max.overlaps = Inf,       
    segment.color = "black", 
    segment.size = 0.5,       
    label.padding = unit(0.15, "lines"), 
    label.r = unit(0.4, "lines"),       
    min.segment.length = 0    
  ) +
  theme_classic() +
  scale_x_continuous(breaks = seq(floor(min(cn_aware_d_compensated$log2FC)), 
                                  ceiling(max(cn_aware_d_compensated$log2FC)), by = 2)) +
  labs(x = expression(Log[2] ~ FC), y = expression(-log[10] ~ Pvalue), col = "Gene group") +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = 'dashed') +
  geom_hline(yintercept = -log10(pval_cut), linetype = "dashed") +
  ggplot2::theme(legend.position = '',
                 legend.text = element_text(size = 14, color = "black"),
                 legend.title = element_text(size = 16, color = "black"),  
                 strip.text = element_text(size = 16, face = "plain", color = "black"),
                 axis.text = element_text(size = 14, color = "black"),
                 axis.title = element_text(size = 16))+
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1)))

p_volcanos

ggsave("CN-aware-DGE/plots/main/volcano_dcg_lusc.png", dpi = 400, width = 4.0, height = 4.0, plot = p_volcanos)
