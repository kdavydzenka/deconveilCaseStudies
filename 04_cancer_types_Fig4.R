setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/")
pkgs <- c("tidyverse", "ggvenn", "VennDiagram", "reactome.db", "fgsea", "org.Hs.eg.db", 
          "data.table", "clusterProfiler", "enrichplot", "ggpubr", "msigdbr", "ggVennDiagram")
sapply(pkgs, require, character.only = TRUE)
source("deconveilCaseStudies/utils/utils.R")
source("deconveilCaseStudies/utils/utils_plot.R")

### Compare DSGs across each cancer type, identify common genes and cancer-specific genes ###

# Data preprocessing

tumor_types <- c("LUAD", "LUSC", "BRCA") 

base_paths <- list(
  res_pydeseq = "deconveilCaseStudies/results_tcga/{tumor}/res_CNnaive.csv",
  res_deconveil = "deconveilCaseStudies/results_tcga/{tumor}/res_CNaware.csv"
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


# Plot genes overlap - venn plot

extract_genes <- function(gene_group_df, gene_group_name) {
  gene_group_df %>%
    dplyr::mutate(gene_group = gene_group_name) %>%
    dplyr::select(geneID = geneID_aware, tumor_type = tumor_type_aware, gene_group)
}

venn_plots <- lapply(names(gene_groups_join), function(group_name) {
  venn_plot(gene_groups_join[[group_name]], group_name)
})

names(venn_plots) <- names(gene_groups_join)

venn_plots$d_sensitive  
venn_plots$d_insensitive  
venn_plots$d_compensated  


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


### Overrepresentation GO functional enrichment analysis ##

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

saveRDS(res_ora_GO_list, file = "deconveilCaseStudies/plots/main/Fig 4/rds/res_oraGO_list.RDS")


# Select  tumor related significant pathways for each gene category

common_dsg <- readRDS("deconveilCaseStudies/plots/main/Fig 4/rds/dsg_common.RDS")
common_dig <- readRDS("deconveilCaseStudies/plots/main/Fig 4/rds/dig_common.RDS")
common_dcg <- readRDS("deconveilCaseStudies/plots/main/Fig 4/rds/dcg_common.RDS")

cancer_genes_oncokb <- read.delim("TCGA/cancerGeneList.tsv")

cancer_genes_oncokb <- cancer_genes_oncokb %>%
  mutate(
    gene_type = case_when(
      Is.Oncogene == "Yes" & Is.Tumor.Suppressor.Gene == "Yes" ~ "Oncogene/TSG",
      Is.Oncogene == "Yes" ~ "Oncogene",
      Is.Tumor.Suppressor.Gene == "Yes" ~ "TSG",
      TRUE ~ "Other"
    )
  ) %>%
  select(geneID = Hugo.Symbol, gene_type) %>%
  distinct()

cancer_genes_oncokb <- cancer_genes_oncokb %>% 
  dplyr::filter(gene_type != "Other")
 


# Plot the distribution of private ONC/TSGs

# Subset geneIDs by exclusive gene_type
oncogene_only <- cancer_genes_oncokb %>%
  filter(gene_type == "Oncogene") %>%
  pull(geneID)

tsg_only <- cancer_genes_oncokb %>%
  filter(gene_type == "TSG") %>%
  pull(geneID)

onc_tsg <- cancer_genes_oncokb %>%
  filter(gene_type == "Oncogene/TSG") %>%
  pull(geneID)

# Get overlaps with private CN-aware gene groups
private_dig_onc <- intersect(oncogene_only, cn_aware_d_insensitive$geneID)
private_dsg_onc <- intersect(oncogene_only, cn_aware_d_sensitive$geneID)
private_dcg_onc <- intersect(oncogene_only, cn_aware_d_compensated$geneID)

private_dig_tsg <- intersect(tsg_only, cn_aware_d_insensitive$geneID)
private_dsg_tsg <- intersect(tsg_only, cn_aware_d_sensitive$geneID)
private_dcg_tsg <- intersect(tsg_only, cn_aware_d_compensated$geneID)

private_dig_onctsg <- intersect(onc_tsg, cn_aware_d_insensitive$geneID)
private_dsg_onctsg <- intersect(onc_tsg, cn_aware_d_sensitive$geneID)
private_dcg_onctsg <- intersect(onc_tsg, cn_aware_d_compensated$geneID)

# Build data frames for each type
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

private_onctsg <- data.frame(
  gene_category = c(rep("DIGs", length(private_dig_onctsg)),
                    rep("DSGs", length(private_dsg_onctsg)),
                    rep("DCGs", length(private_dcg_onctsg))),
  gene_type = "Oncogene/TSG",
  geneID = c(private_dig_onctsg, private_dsg_onctsg, private_dcg_onctsg)
)


private_genes <- bind_rows(private_onc, private_tsg, private_onctsg)

# Summarize for plotting
distribution_summary <- private_genes %>%
  count(gene_category, gene_type, name = "gene_count") %>%
  mutate(gene_category = factor(gene_category, levels = c("DSGs", "DIGs", "DCGs")))

# Plot
gene_type_colors <- c(
  "Oncogene"      = "#BB002198",
  "TSG"           = "#00468BB2",
  "Oncogene/TSG"  = "#EFC000FF"
)

barplot <- gene_barplot(distribution_summary, gene_type_colors)
barplot  

ggsave("deconveilCaseStudies/plots/main/barplot_onc_tsg.png", dpi = 500, width = 5.0, height = 4.0, plot = barplot)


## Generate Enrichment plot ##

ora_GO_path <- "deconveilCaseStudies/results_tcga/oraGO_cancer_types"
tumor_types <- c("luad", "lusc", "brca")
categories <- c("Dosage-sensitive", "Dosage-insensitive", "Dosage-compensated")
category_labels <- c("Dosage-sensitive" = "DSGs", "Dosage-insensitive" = "DIGs", "Dosage-compensated" = "DCGs")

# Load ORA results
oraGO_pancancer <- setNames(lapply(tumor_types, function(tumor) {
  readRDS(file.path(ora_GO_path, paste0(tumor, ".RDS")))
}), tumor_types)


## Shared genes relevant cancer-specific pathways ##

GO_terms_shared <- list(
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



## Private cancer-relevant specific pathways ##

GO_terms_private <- list(
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


#GO_terms_luad <- list(
  #luad = list(
    #'Dosage-sensitive' = c("DNA-templated DNA replication",
                           #"chaperone-mediated protein folding",
                           #"positive regulation of DNA metabolic process",
                           #"fibroblast proliferation",
                           #"extracellular matrix organization"),
    
    #'Dosage-insensitive' = c("B cell mediated immunity",
                             #"immunoglobulin production",
                             #"regulation of mitotic nuclear division",
                             #"cell-cell adhesion via plasma-membrane adhesion molecules",
                             #"cilium movement"),
    
    #'Dosage-compensated' = c("positive regulation of cell-cell adhesion",
                             #"positive regulation of T cell activation",
                             #"positive regulation of lymphocyte mediated immunity",
                             #"negative regulation of kinase activity",
                             #"response to hypoxia")
  #)
#)


# Function to process data

process_GO_data <- function(tumor, category, oraGO_pancancer, GO_terms) {
  terms <- GO_terms[[tumor]][[category]]
  df <- oraGO_pancancer[[tumor]][[category]]
  if (is.null(df) || is.null(terms)) return(NULL)
  
  df %>%
    filter(Description %in% terms) %>%
    separate(GeneRatio, into = c("numerator", "denominator"), sep = "/", convert = TRUE) %>%
    mutate(
      GeneRatio_val = numerator / denominator,
      log_padjust = -log10(p.adjust),
      Tumor = str_to_upper(tumor),
      Category = category_labels[category]
    )
}

# Process all combinations
dotplot_data <- map_dfr(tumor_types, function(tumor) {
  map_dfr(categories, function(category) {
    process_GO_data(tumor, category, oraGO_pancancer, GO_terms_private)
  })
})

# Split long terms into 2 lines
split_in_two_lines <- function(text) {
  words <- strsplit(text, " ")[[1]]
  midpoint <- ceiling(length(words) / 2)
  paste(paste(words[1:midpoint], collapse = " "),
        paste(words[(midpoint + 1):length(words)], collapse = " "), sep = "\n")
}

dotplot_data <- dotplot_data %>%
  mutate(
    Description = sapply(as.character(Description), split_in_two_lines),
    GeneRatio_val = round(GeneRatio_val, 2),
    Description = factor(Description, levels = unique(Description)),
    Category = factor(Category, levels = c("DSGs", "DIGs", "DCGs"))
  )


saveRDS(dotplot_data, file = "deconveilCaseStudies/plots/main/Fig 4/rds/dotplot_GO_common.RDS")
saveRDS(dotplot_data, file = "deconveilCaseStudies/plots/supplementary/rds/dotplot_GO_shared.RDS")


# Plot
plot_GO <- enrichment_dotplot(dotplot_data)
plot_GO  

ggsave("deconveilCaseStudies/plots/main/GOdotplot_pancancer_shared.png", dpi = 500, width = 13.0, height = 5.0, plot = plot_GO)
ggsave("deconveilCaseStudies/plots/supplementary/GOdotplot_pancancer_private.png", dpi = 500, width = 24.0, height = 5.0, plot = plot_GO)
ggsave("deconveilCaseStudies/plots/main/GOdotplot_luad_v2.png", dpi = 500, width = 18.0, height = 4.5, plot = plot_GO)
