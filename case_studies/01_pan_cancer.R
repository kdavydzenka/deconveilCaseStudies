setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/")
pkgs <- c("tidyverse", "ggvenn", "VennDiagram", "reactome.db", "fgsea", "org.Hs.eg.db", 
          "data.table", "clusterProfiler", "enrichplot", "ggpubr", "msigdbr", "ComplexUpset")
sapply(pkgs, require, character.only = TRUE)
source("deconveilCaseStudies/utils.R")

### Compare DSGs across each cancer type, identify common genes and cancer-specific genes ###

# Data preprocessing

tumor_types <- c("LUAD", "LUSC", "LIHC") #Pan-cancer test
#tumor_types <- c("LIHC")

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

#saveRDS(gene_groups, file = "TCGA/BRCA/case_study/gene_groups_brca.RDS")

prepare_cn_specific_data <- function(df, gene_group_name) {
  df %>%
    mutate(gene_group = gene_group_name) %>%
    select(geneID = geneID_naive, log2FC = logFC_naive, padj = padj_naive, isDE = isDE_naive, DEtype = DEtype_naive, tumor_type = tumor_type_naive, method = method_naive, gene_group)
}



# Plot genes overlap

# UpsetPlot

extract_genes <- function(gene_group_df, gene_group_name) {
  gene_group_df %>%
    mutate(gene_group = gene_group_name) %>%
    select(geneID = geneID_aware, tumor_type = tumor_type_aware, gene_group)
}

cn_aware_d_sensitive <- extract_genes(gene_groups$d_sensitive, "Dosage-sensitive")
cn_aware_d_insensitive <- extract_genes(gene_groups$d_insensitive, "Dosage-insensitive")
cn_aware_d_compensated <- extract_genes(gene_groups$d_compensated, "Dosage-compensated")

#gene_groups_join <- list(cn_aware_d_sensitive, cn_aware_d_insensitive, cn_aware_d_compensated)
#names(gene_groups_join) <- c("d_sensitive", "d_insensitive", "d_compensated")
#saveRDS(gene_groups_join, file = "TCGA/BRCA/case_study/gene_groups.RDS")

#LUAD_genes <- upset_data %>% filter(LUAD == 1) %>% pull(geneID)
#LUSC_genes <- upset_data %>% filter(LUSC == 1) %>% pull(geneID)
#LIHC_genes <- upset_data %>% filter(LIHC == 1) %>% pull(geneID)

#venn_list <- list(
  #LUAD = LUAD_genes,
  #LUSC = LUSC_genes,
  #LIHC = LIHC_genes
#)

#library(ggVennDiagram)

#venn_plot <- ggVennDiagram(venn_list) +
  #scale_fill_gradient(low = "white", high = "skyblue") +
  #theme(
    #plot.title = element_text(hjust = 0.5, size = 16),
    #plot.background = element_rect(fill = "white"),
    #legend.position = "right"
  #)
#venn_plot


upset_data <- cn_aware_d_insensitive %>%
  mutate(present = 1) %>% 
  pivot_wider(names_from = tumor_type, values_from = present, values_fill = 0) %>% 
  select(-gene_group) %>%
  distinct()

common_genes <- upset_data %>%
  filter(LUAD == 1, LUSC == 1, LIHC == 1) %>% 
  pull(geneID)
common_genes

upset_data <- upset_data %>% select(-geneID)
upset_data <- upset_data %>%
  mutate(across(all_of(tumor_types), ~ ifelse(. > 0, 1, 0)))

# UpSet plot
upset_plot <- upset(
  upset_data,
  intersect = tumor_types,
  set_sizes = upset_set_size(
    position = 'right'),
  base_annotations = list(
    'Intersection size' = intersection_size(
      counts = T,
      text_colors = c(
        on_background = 'black',
        on_bar = 'white'
      ),
      fill = "#925E9FB2"  
    ) +
      ylab('Intersection size') +
      theme(
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14)    
      )
  ),
  queries=list(
    upset_query(set='LUAD', fill='#A3BE8CFF'),
    upset_query(set='LUSC', fill='#EBCB8BFF'),
    #upset_query(set='BRCA', fill='#ff9896'),
    upset_query(set='LIHC', fill='#aec7e8'),
    upset_query(intersect=c('LUAD', 'LUSC', 'LIHC'),
      color='red', fill='red',
      only_components=c('Intersection size'))),
  width_ratio = 0.2,      
  name = 'gene sets') +
  expand_limits(y = 1000) +  
  theme_minimal()

upset_plot


# Identify shared and cancer-specific oncogenes/TSGs
#cancer_genes <- read.delim("TCGA/lung/cancerGeneList.tsv")

#oncogenes <- cancer_genes %>% dplyr::filter(Is.Oncogene=="Yes") %>% 
  #dplyr::select(Hugo.Symbol) %>% 
  #dplyr::rename(geneID = Hugo.Symbol) %>% 
  #dplyr::mutate(gene_type = "Oncogene")
  
#tsg <- cancer_genes %>% dplyr::filter(Is.Tumor.Suppressor.Gene=="Yes") %>% 
  #dplyr::select(Hugo.Symbol) %>% 
  #dplyr::rename(geneID=Hugo.Symbol) %>% 
  #dplyr::mutate(gene_type = "TSG")
  
#cancer_genes_oncokb <- rbind(oncogenes, tsg)

#d_compensated_cancer_g <- cn_aware_d_compensated[cn_aware_d_compensated$geneID %in% cancer_genes_oncokb$geneID ,]
#d_insensitive_cancer_g <- cn_aware_d_insensitive[cn_aware_d_insensitive$geneID %in% cancer_genes_oncokb$geneID ,]
#d_sensitive_cancer_g <- cn_aware_d_sensitive[cn_aware_d_sensitive$geneID %in% cancer_genes_oncokb$geneID ,]

#cancer_gene_lists <- d_sensitive_cancer_g %>%
  #group_by(tumor_type, gene_group) %>%
  #summarise(geneIDs = list(geneID), .groups = "drop")

#cancer_gene_lists <- cancer_gene_lists %>%
  #with(setNames(geneIDs, tumor_type))

#ggvenn(cancer_gene_lists, 
       #show_percentage = FALSE,  #
       #fill_color = c("#E6A9EC", "#A9CDE6", "#E6D7A9", "#A9E6B3")) 

#tumor_specific_cancer_g <- lapply(names(cancer_gene_lists), function(tumor) {
  #setdiff(cancer_gene_lists[[tumor]], unlist(cancer_gene_lists[-match(tumor, names(cancer_gene_lists))]))
#})
#names(tumor_specific_cancer_g) <- names(cancer_gene_lists)

#common_cancer_g <- lapply(names(cancer_gene_lists), function(tumor) {
  #setdiff(unlist(cancer_gene_lists), tumor_specific_cancer_g[[tumor]])
#})
#names(common_cancer_g) <- names(cancer_gene_lists)
#common_cancer_g <- Reduce(intersect, common_cancer_g)


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


# Hallmark ORA
#res_ora_H_list <- lapply(names(gene_groups), function(group_name) {
  #gene_df <- gene_groups[[group_name]] %>% dplyr::select(geneID, log2FC)
  #perform_enrichment_H(gene_df, group_name)
#})
#names(res_ora_H_list) <- names(gene_groups)

#res_H_d_sensitive <- res_ora_H_list[["Dosage-sensitive"]] %>% as.data.frame()
#res_H_d_insensitive <- res_ora_H_list[["Dosage-insensitive"]] %>% as.data.frame()
#res_H_d_compensated <- res_ora_H_list[["Dosage-compensated"]] %>% as.data.frame()

#saveRDS(res_ora_H_list, file = "deconveilCaseStudies/results/pancancer_ora_H/brca.RDS")



# Select  tumor related significant pathways for each gene category

ora_GO_path <- "deconveilCaseStudies/results/pancancer_ora_GO"
tumor_types <- c("luad", "lusc", "lihc")

oraGO_pancancer <- list()

for (tumor in tumor_types) {
  file_path <- file.path(ora_GO_path, paste0(tumor, ".RDS"))
  oraGO_pancancer[[tumor]] <- readRDS(file_path)
}

#dig_luad <- as.data.frame(oraGO_pancancer[["luad"]][["Dosage-insensitive"]])
#dig_lusc <- as.data.frame(oraGO_pancancer[["lusc"]][["Dosage-insensitive"]])
#dig_lihc <- as.data.frame(oraGO_pancancer[["lihc"]][["Dosage-insensitive"]])

#dig_luad <- dig_luad %>%
  #mutate(geneID = sapply(geneID, convert_entrez_to_symbol))

#dig_lusc <- dig_lusc %>%
  #mutate(geneID = sapply(geneID, convert_entrez_to_symbol))

#dig_lihc <- dig_lihc %>%
  #mutate(geneID = sapply(geneID, convert_entrez_to_symbol))

#common_dsg <- c("AC244033.2", "CBX8", "CCHCR1", "DCST2", "EXOSC4", "FAM189B", "FLAD1", "GBAP1", "HCN3",
                #"LAPTM4B", "MSTO1", "P4HB", "PACC1", "PYCR3", "RHPN1-AS1", "RPL22L1", "RUSC1", "SCRIB",
                #"SNRPE", "TBC1D31", "ZNF572")  

#common_dsg <- as.data.frame(common_dsg)

#common_dcg <- c("ACTB", "CFAP69", "HCK", "HMCN1", "NCALD", "RASSF5", "SLC19A2", "TECTA")


GO_terms <- list(
  luad = list(
    'Dosage-sensitive' = c("positive regulation of DNA metabolic process",
                           "rRNA metabolic process",
                           "amino acid metabolic process",
                           "nucleocytoplasmic transport",
                           "cell-cell recognition",
                           "nucleotide biosynthetic process",
                           "protein folding",
                           "RNA modification"),
    'Dosage-compensated' = c("positive regulation of lymphocyte activation",
                             "immune response-regulating cell surface receptor signaling pathway",
                             "regulation of leukocyte proliferation",
                             "cell-cell adhesion via plasma-membrane adhesion molecules",
                             "second-messenger-mediated signaling",
                             "cell-matrix adhesion"),
    'Dosage-insensitive' = c("cell cycle checkpoint signaling",
                      "mitotic nuclear division",
                      "cellular response to tumor necrosis factor",
                      "ERK1 and ERK2 cascade",
                      "negative regulation of Wnt signaling pathway",
                      "B cell mediated immunity",
                      "humoral immune response",
                      "interleukin-6 production")
  ),
  lusc = list(
    'Dosage-sensitive' = c("positive regulation of DNA metabolic process",
                           "rRNA metabolic process",
                           "amino acid metabolic process",
                           "nucleocytoplasmic transport",
                           "cell-cell recognition",
                           "nucleotide biosynthetic process",
                           "protein folding",
                           "RNA modification"),
    'Dosage-compensated' = c("positive regulation of lymphocyte activation",
                             "immune response-regulating cell surface receptor signaling pathway",
                             "regulation of leukocyte proliferation",
                             "cell-cell adhesion via plasma-membrane adhesion molecules",
                             "second-messenger-mediated signaling",
                             "cell-matrix adhesion"),
    'Dosage-insensitive' = c("cell cycle checkpoint signaling",
                             "mitotic nuclear division",
                             "cellular response to tumor necrosis factor",
                             "ERK1 and ERK2 cascade",
                             "negative regulation of Wnt signaling pathway",
                             "B cell mediated immunity",
                             "humoral immune response",
                             "interleukin-6 production")
  ),
  lihc = list(
    'Dosage-sensitive' = c("positive regulation of DNA metabolic process",
                           "rRNA metabolic process",
                           "amino acid metabolic process",
                           "nucleocytoplasmic transport",
                           "cell-cell recognition",
                           "nucleotide biosynthetic process",
                           "protein folding",
                           "RNA modification"),
    'Dosage-compensated' = c("positive regulation of lymphocyte activation",
                             "immune response-regulating cell surface receptor signaling pathway",
                             "regulation of leukocyte proliferation",
                             "cell-cell adhesion via plasma-membrane adhesion molecules",
                             "second-messenger-mediated signaling",
                             "cell-matrix adhesion"),
    'Dosage-insensitive' = c("cell cycle checkpoint signaling",
                             "mitotic nuclear division",
                             "cellular response to tumor necrosis factor",
                             "ERK1 and ERK2 cascade",
                             "negative regulation of Wnt signaling pathway",
                             "B cell mediated immunity",
                             "humoral immune response",
                             "interleukin-6 production")
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


res_GO_ds <- result_GO_pancancer[["luad"]][["Dosage-sensitive"]]
res_GO_dc <- result_GO_pancancer[["luad"]][["Dosage-compensated"]]
res_GO_dins <- result_GO_pancancer[["luad"]][["Dosage-insensitive"]]


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
                   lihc = "LIHC")
  )

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
    strip.text = element_text(size = 18, face = "plain", color = "black"),
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

ggsave("deconveilCaseStudies/plots/main/GOdotplot_pancancer_v2.png", dpi = 400, width = 22.0, height = 5.0, plot = plot_GO)

  
