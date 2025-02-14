setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/")
pkgs <- c("tidyverse", "ggvenn", "VennDiagram", "reactome.db", "fgsea", "org.Hs.eg.db", 
          "data.table", "clusterProfiler", "enrichplot", "ggpubr", "msigdbr", "ComplexUpset", "ggVennDiagram")
sapply(pkgs, require, character.only = TRUE)
source("deconveilCaseStudies/utils.R")

### Compare DSGs across each cancer type, identify common genes and cancer-specific genes ###

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
    select(X, log2FoldChange, padj, isDE, DEtype, tumor_type, method) 
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
    select(geneID = geneID_naive, log2FC = logFC_naive, padj = padj_naive, isDE = isDE_naive, DEtype = DEtype_naive, tumor_type = tumor_type_naive, method = method_naive, gene_group)
}


# Plot genes overlap

extract_genes <- function(gene_group_df, gene_group_name) {
  gene_group_df %>%
    mutate(gene_group = gene_group_name) %>%
    select(geneID = geneID_aware, tumor_type = tumor_type_aware, gene_group)
}

cn_aware_d_sensitive <- extract_genes(gene_groups$d_sensitive, "Dosage-sensitive")
cn_aware_d_insensitive <- extract_genes(gene_groups$d_insensitive, "Dosage-insensitive")
cn_aware_d_compensated <- extract_genes(gene_groups$d_compensated, "Dosage-compensated")

gene_groups_join <- list(cn_aware_d_sensitive, cn_aware_d_insensitive, cn_aware_d_compensated)
names(gene_groups_join) <- c("d_sensitive", "d_insensitive", "d_compensated")

venn_data <- cn_aware_d_compensated %>%
  mutate(present = 1) %>% 
  pivot_wider(names_from = tumor_type, values_from = present, values_fill = 0) %>% 
  select(-gene_group) %>%
  distinct()

LUAD_genes <- venn_data %>% filter(LUAD == 1) %>% pull(geneID)
LUSC_genes <- venn_data %>% filter(LUSC == 1) %>% pull(geneID)
BRCA_genes <- venn_data %>% filter(BRCA == 1) %>% pull(geneID)

venn_list <- list(
  LUAD = LUAD_genes,
  LUSC = LUSC_genes,
  BRCA = BRCA_genes
)

venn_plot <- ggVennDiagram(venn_list) +
  scale_fill_gradient(low = "lightgray", high = "#B2474599") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    plot.background = element_rect(fill = "white"),
    legend.position = "right"
  )
venn_plot


upset_data <- cn_aware_d_compensated %>%
  mutate(present = 1) %>% 
  pivot_wider(names_from = tumor_type, values_from = present, values_fill = 0) %>% 
  select(-gene_group) %>%
  distinct()

common_genes <- upset_data %>%
  filter(LUAD == 1, LUSC == 1, BRCA == 1) %>% 
  pull(geneID)
common_genes

saveRDS(common_genes, file = "deconveilCaseStudies/results/oraGO_cancer_types/dcg_common.RDS")


## Extract each gene category ##

d_sensitive <- gene_groups$d_sensitive %>%
  mutate(gene_group = "Dosage-sensitive") %>%
  select(geneID = geneID_aware, log2FC = logFC_aware, padj = padj_aware, isDE = isDE_aware, DEtype = DEtype_aware, tumor_type = tumor_type_aware, method = method_aware, gene_group)

d_insensitive <- gene_groups$d_insensitive %>%
  mutate(gene_group = "Dosage-insensitive") %>%
  select(geneID = geneID_aware, log2FC = logFC_aware, padj = padj_aware, isDE = isDE_aware, DEtype = DEtype_aware, tumor_type = tumor_type_aware, method = method_aware, gene_group)

d_compensated <- gene_groups$d_compensated %>%
  mutate(gene_group = "Dosage-compensated") %>%
  select(geneID = geneID_aware, log2FC = logFC_aware, padj = padj_aware, isDE = isDE_aware, DEtype = DEtype_aware, tumor_type = tumor_type_aware, method = method_aware, gene_group)


# Overrepresentation GO functional enrichment analysis

d_sensitive_g <- d_sensitive %>% dplyr::select(geneID, log2FC)
d_insensitive_g <- d_insensitive %>% dplyr::select(geneID, log2FC)
d_compensated_g <- d_compensated %>% dplyr::select(geneID, log2FC)

gene_groups <- list(
  "Dosage-sensitive" = d_sensitive,
  "Dosage-compensated" = d_compensated,
  "Dosage-insensitive" = d_insensitive
)

# GO ORA
res_ora_GO_list <- lapply(names(gene_groups), function(group_name) {
  gene_df <- gene_groups[[group_name]] %>% dplyr::select(geneID, log2FC)
  perform_enrichment_GO(gene_df, group_name)
})
names(res_ora_GO_list) <- names(gene_groups)

res_GO_d_sensitive <- res_ora_GO_list[["Dosage-sensitive"]] %>% as.data.frame()
res_GO_d_insensitive <- res_ora_GO_list[["Dosage-insensitive"]] %>% as.data.frame()
res_GO_d_compensated <- res_ora_GO_list[["Dosage-compensated"]] %>% as.data.frame()

saveRDS(res_ora_GO_list, file = "deconveilCaseStudies/results/pancancer_ora_GO/lihc.RDS")


# Select  tumor related significant pathways for each gene category

common_dsg <- readRDS("deconveilCaseStudies/results/oraGO_cancer_types/dsg_common.RDS")
common_dig <- readRDS("deconveilCaseStudies/results/oraGO_cancer_types/dig_common.RDS")
common_dcg <- readRDS("deconveilCaseStudies/results/oraGO_cancer_types/dcg_common.RDS")

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


# Plot the distribution of private ONC/TSGs

private_dig_onc <- intersect(oncogenes$geneID, cn_aware_d_insensitive$geneID)
private_dsg_onc <- intersect(oncogenes$geneID, cn_aware_d_sensitive$geneID)
private_dcg_onc <- intersect(oncogenes$geneID, cn_aware_d_compensated$geneID)

private_dig_tsg <- intersect(tsg$geneID, cn_aware_d_insensitive$geneID)
private_dsg_tsg <- intersect(tsg$geneID, cn_aware_d_sensitive$geneID)
private_dcg_tsg <- intersect(tsg$geneID, cn_aware_d_compensated$geneID)


private_onc <- data.frame(
  gene_category = c(rep("DIGs", length(private_dig_onc)),
                    rep("DSGs", length(private_dsg_onc)),
                    rep("DCGs", length(private_dcg_onc))),
  gene_type = "Oncogene",
  geneID = c(private_dig_onc, private_dsg_onc, private_dcg_onc)
)

private_tsg <- data.frame(
  gene_category = c(rep("DIGs", length(private_dig_tsg)),
                    rep("DSGs", length(private_dsg_tsg)),
                    rep("DCGs", length(private_dcg_tsg))),
  gene_type = "TSG",
  geneID = c(private_dig_tsg, private_dsg_tsg, private_dcg_tsg)
)

private_genes <- rbind(private_onc, private_tsg)

distribution_summary <- private_genes %>%
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
  labs(title = "Private genes: ONC & TSGs",
       x = "Gene category",
       y = "Gene count",
       fill = "Gene type") +
  scale_fill_manual(values = c("Oncogene" = "#BB002198", "TSG" = "#00468BB2"))+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        axis.text.x = element_text(size = 16, color = "black"),  
        axis.text.y = element_text(size = 16, color = "black"),
        legend.key.size = unit(0.6, "cm"), 
        legend.text = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 14, color = "black"),
        legend.spacing.y = unit(0.2, 'cm'),
        axis.text = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 14))
barplot

ggsave("deconveilCaseStudies/plots/main/barplot_onc_tsg.png", dpi = 400, width = 4.5, height = 5.0, plot = barplot)


# Enrichment plot

ora_GO_path <- "deconveilCaseStudies/results/oraGO_cancer_types"
tumor_types <- c("luad", "lusc", "brca")

oraGO_pancancer <- list()

for (tumor in tumor_types) {
  file_path <- file.path(ora_GO_path, paste0(tumor, ".RDS"))
  oraGO_pancancer[[tumor]] <- readRDS(file_path)
}


dsg_luad <- as.data.frame(oraGO_pancancer[["luad"]][["Dosage-sensitive"]])
dsg_lusc <- as.data.frame(oraGO_pancancer[["lusc"]][["Dosage-sensitive"]])
dsg_brca <- as.data.frame(oraGO_pancancer[["brca"]][["Dosage-sensitive"]])

dig_luad <- as.data.frame(oraGO_pancancer[["luad"]][["Dosage-insensitive"]])
dig_lusc <- as.data.frame(oraGO_pancancer[["lusc"]][["Dosage-insensitive"]])
dig_brca <- as.data.frame(oraGO_pancancer[["brca"]][["Dosage-insensitive"]])

dcg_luad <- as.data.frame(oraGO_pancancer[["luad"]][["Dosage-compensated"]])
dcg_lusc <- as.data.frame(oraGO_pancancer[["lusc"]][["Dosage-compensated"]])
dcg_brca <- as.data.frame(oraGO_pancancer[["brca"]][["Dosage-compensated"]])


# Private cancer-relevant specific pathways

GO_terms <- list(
  luad = list(
    'Dosage-sensitive' = c("fibroblast proliferation",
                           "chaperone-mediated protein folding",
                           "recombinational repair"),
    
    'Dosage-insensitive' = c("immunoglobulin production",
                             "B cell mediated immunity",
                             "regulation of mitotic nuclear division"),
    
    'Dosage-compensated' = c("MHC protein complex assembly",
                             "positive regulation of T cell activation",
                             "negative regulation of kinase activity")
  ),
  lusc = list(
    'Dosage-sensitive' = c("ribosome biogenesis",
                           "epithelial tube formation",
                           "cellular response to interferon-beta"),
    
    'Dosage-insensitive' = c("mesenchymal cell differentiation",
                             "neutrophil chemotaxis",
                             "morphogenesis of a branching epithelium"),
    
    'Dosage-compensated' = c("negative regulation of cytokine production involved in immune response",
                             "positive regulation of lymphocyte differentiation")
  ),
  brca = list(
    'Dosage-sensitive' = c("regulation of calcium ion-dependent exocytosis",
                           "positive regulation of T cell mediated cytotoxicity"),
    
    'Dosage-insensitive' = c("hormone metabolic process",
                             "calcium ion transmembrane import into cytosol",
                             "striated muscle cell development"),
    
    'Dosage-compensated' = c("striated muscle tissue development",
                             "hormone metabolic process",
                             "positive regulation of insulin secretion")
  )
)


# Shared genes relevant cancer-specific pathways

GO_terms <- list(
  luad = list(
    'Dosage-sensitive' = c("positive regulation of DNA metabolic process",
                           "amino acid metabolic process",
                           "nucleotide biosynthetic process",
                           "regulation of cellular respiration",
                           "cell cycle checkpoint signaling",
                           "regulation of miRNA transcription",
                           "DNA-templated DNA replication"),
    
    'Dosage-compensated' = c("humoral immune response",
                             "cell-cell adhesion via plasma-membrane adhesion molecules",
                             "negative regulation of lipid metabolic process",
                             "leukocyte chemotaxis",
                             "canonical Wnt signaling pathway"),
    
    'Dosage-insensitive' = c("regulation of mitotic nuclear division", 
                             "cell cycle checkpoint signaling",
                             "positive regulation of kinase activity",
                             "positive regulation of ERK1 and ERK2 cascade",
                             "cell-matrix adhesion",
                             "positive regulation of interleukin-6 production",
                             "regulation of leukocyte differentiation")
  ),
  lusc = list(
    'Dosage-sensitive' = c("positive regulation of DNA metabolic process",
                           "amino acid metabolic process",
                           "nucleotide biosynthetic process",
                           "regulation of cellular respiration",
                           "cell cycle checkpoint signaling",
                           "regulation of miRNA transcription",
                           "DNA-templated DNA replication"),
    
    'Dosage-compensated' = c("humoral immune response",
                             "cell-cell adhesion via plasma-membrane adhesion molecules",
                             "negative regulation of lipid metabolic process",
                             "leukocyte chemotaxis",
                             "canonical Wnt signaling pathway"),
    
    'Dosage-insensitive' = c("regulation of mitotic nuclear division", 
                             "cell cycle checkpoint signaling",
                             "positive regulation of kinase activity",
                             "positive regulation of ERK1 and ERK2 cascade",
                             "cell-matrix adhesion",
                             "positive regulation of interleukin-6 production",
                             "regulation of leukocyte differentiation")
  ),
  brca = list(
    'Dosage-sensitive' = c("positive regulation of DNA metabolic process",
                           "amino acid metabolic process",
                           "nucleotide biosynthetic process",
                           "regulation of cellular respiration",
                           "cell cycle checkpoint signaling",
                           "regulation of miRNA transcription",
                           "DNA-templated DNA replication"),
    
    'Dosage-compensated' = c("humoral immune response",
                             "cell-cell adhesion via plasma-membrane adhesion molecules",
                             "negative regulation of lipid metabolic process",
                             "leukocyte chemotaxis",
                             "canonical Wnt signaling pathway"),
    
    'Dosage-insensitive' = c("regulation of mitotic nuclear division", 
                             "cell cycle checkpoint signaling",
                             "positive regulation of kinase activity",
                             "positive regulation of ERK1 and ERK2 cascade",
                             "cell-matrix adhesion",
                             "positive regulation of interleukin-6 production",
                             "regulation of leukocyte differentiation")
  )
)


process_pancancer_data <- function(tumor, category, oraGO_pancancer) {
  if (!is.null(oraGO_pancancer[[tumor]]) && !is.null(oraGO_pancancer[[tumor]][[category]])) {
    terms <- GO_terms[[tumor]][[category]]
    if (!is.null(terms)) {
      oraGO_pancancer[[tumor]][[category]] %>%
        filter(Description %in% terms) %>%
        separate(GeneRatio, into = c("numerator", "denominator"), sep = "/", convert = TRUE) %>%
        mutate(GeneRatio_val = numerator / denominator)
    } else {
      NULL
    }
  } else {
    NULL
  }
}

result_GO_pancancer <- lapply(names(GO_terms), function(tumor) {
  lapply(names(GO_terms[[tumor]]), function(category) {
    processed_data <- process_pancancer_data(tumor, category, oraGO_pancancer)
    if (!is.null(processed_data) && "p.adjust" %in% colnames(processed_data)) {
      processed_data$log_padjust <- -log10(processed_data$p.adjust)
    }
    processed_data
  }) %>%
    setNames(names(GO_terms[[tumor]])) 
}) %>%
  setNames(names(GO_terms)) 


dotplot_data <- bind_rows(
  lapply(names(result_GO_pancancer), function(tumor) {
    bind_rows(
      lapply(names(result_GO_pancancer[[tumor]]), function(category) {
        if (!is.null(result_GO_pancancer[[tumor]][[category]])) {
          result_GO_pancancer[[tumor]][[category]] %>%
            mutate(
              Tumor = tumor,
              Category = category
            )
        }
      }), .id = "GeneCategory"
    )
  }), .id = "TumorType"
)

dotplot_data <- dotplot_data %>%
  mutate(
    Category = recode(Category,
                      `Dosage-sensitive` = "DSGs",
                      `Dosage-insensitive` = "DIGs",
                      `Dosage-compensated` = "DCGs"),
    Tumor = recode(Tumor,
                   luad = "LUAD",
                   lusc = "LUSC",
                   brca = "BRCA")
  )

dotplot_data <- dotplot_data %>%
  group_by(Category, Description) %>%
  mutate(order = match(Tumor, c("BRCA", "LUAD", "LUSC"))) %>%
  ungroup() %>%
  arrange(Category, order) %>%
  mutate(Description = factor(Description, levels = unique(Description)))


split_in_two_lines <- function(text) {
  words <- strsplit(text, " ")[[1]]  # Split text into words
  midpoint <- ceiling(length(words) / 2)  # Find the midpoint
  paste(paste(words[1:midpoint], collapse = " "), "\n", paste(words[(midpoint + 1):length(words)], collapse = " "), sep = "")
}

dotplot_data$Description <- sapply(dotplot_data$Description, split_in_two_lines)


plot_GO <- ggplot(dotplot_data, aes(x = Tumor, y = Description)) +
  geom_point(aes(size = Count, color = log_padjust)) +
  facet_wrap(~factor(Category, levels = c("DSGs", "DIGs", "DCGs")), scales = "free", ncol = 3) +
  scale_color_gradient(low = "blue", high = "#FAA43A", name = "-log10(p.adjust)") +
  theme_bw() +
  theme(
    panel.spacing = unit(1, "lines"),
    axis.text.x = element_text(size = 16, color = "black", angle = 45, hjust = 1),  
    axis.text.y = element_text(size = 16, color = "black"),
    legend.key.size = unit(0.6, "cm"), 
    legend.text = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 14, color = "black"),
    legend.spacing.y = unit(0.2, 'cm'),
    strip.text = element_text(size = 20, face = "plain", color = "black"),
    axis.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 16)
  ) +
  labs(
    title = "",
    x = "Tumor types",
    y = "Biological Process GO term",
    size = "Gene Count"
  )
plot_GO

ggsave("deconveilCaseStudies/plots/main/GOdotplot_pancancer.png", dpi = 400, width = 18.0, height = 6.0, plot = plot_GO)
ggsave("deconveilCaseStudies/plots/supplementary/GOdotplot_pancancer_private.png", dpi = 400, width = 24.0, height = 5.0, plot = plot_GO)

