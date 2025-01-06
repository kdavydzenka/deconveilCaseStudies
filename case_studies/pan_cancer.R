setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/")
pkgs <- c("tidyverse", "ggvenn", "VennDiagram", "reactome.db", "fgsea", "org.Hs.eg.db", 
          "data.table", "clusterProfiler", "enrichplot", "ggpubr", "msigdbr", "ComplexUpset")
sapply(pkgs, require, character.only = TRUE)
source("deconveilCaseStudies/utils.R")

### Compare DSGs across each cancer type, identify common genes and cancer-specific genes ###

# Data preprocessing

tumor_types <- c("LUAD", "LUSC", "BRCA", "LIHC") #Pan-cancer test
tumor_types <- c("BRCA")

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

#saveRDS(gene_groups, file = "TCGA/BRCA/case_study/gene_groups_brca.RDS")

prepare_cn_specific_data <- function(df, gene_group_name) {
  df %>%
    mutate(gene_group = gene_group_name) %>%
    select(geneID = geneID_naive, log2FC = logFC_naive, padj = padj_naive, isDE = isDE_naive, DEtype = DEtype_naive, tumor_type = tumor_type_naive, method = method_naive, gene_group)
}



# Plot genes overlap

# Venn plot
#ggvenn(gene_lists, 
       #show_percentage = FALSE,  # Optional: Turn this off or on
       #fill_color = c("#E6A9EC", "#A9CDE6", "#E6D7A9", "#A9E6B3", "#D4A9E6")) 


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

upset_data <- cn_aware_d_sensitive %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = tumor_type, values_from = present, values_fill = 0) %>%
  select(-gene_group) %>%
  distinct()

common_genes <- upset_data %>%
  filter(LUAD == 1, LUSC == 1, BRCA == 1, LIHC == 1) %>% 
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
    upset_query(set='BRCA', fill='#A3BE8CFF'),
    upset_query(set='LUAD', fill='#EBCB8BFF'),
    upset_query(set='LUSC', fill='#ff9896'),
    upset_query(set='LIHC', fill='#aec7e8'),
    upset_query(intersect=c('BRCA', 'LUAD', 'LUSC', 'LIHC'),
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
res_ora_H_list <- lapply(names(gene_groups), function(group_name) {
  gene_df <- gene_groups[[group_name]] %>% dplyr::select(geneID, log2FC)
  perform_enrichment_H(gene_df, group_name)
})
names(res_ora_H_list) <- names(gene_groups)

res_H_d_sensitive <- res_ora_H_list[["Dosage-sensitive"]] %>% as.data.frame()
res_H_d_insensitive <- res_ora_H_list[["Dosage-insensitive"]] %>% as.data.frame()
res_H_d_compensated <- res_ora_H_list[["Dosage-compensated"]] %>% as.data.frame()

saveRDS(res_ora_H_list, file = "deconveilCaseStudies/results/pancancer_ora_H/brca.RDS")



# Select  tumor related significant pathways for each gene category

ora_GO_path <- "deconveilCaseStudies/results/pancancer_ora_GO"
tumor_types <- c("luad", "lusc", "brca", "lihc")

ora_GO_pancancer <- list()

for (tumor in tumor_types) {
  file_path <- file.path(ora_GO_path, paste0(tumor, ".RDS"))
  ora_GO_pancancer[[tumor]] <- readRDS(file_path)
}

GO_terms <- list(
  luad = list(
    'Dosage-sensitive' = c("DNA replication", "DNA recombination", "positive regulation of DNA metabolic process"),
    'Dosage-compensated' = c("positive regulation of cell-cell adhesion", "positive regulation of T cell activation",
                           "positive regulation of adaptive immune response", "regulation of lymphocyte proliferation"),
    'Dosage-insensitive' = c("B cell mediated immunity", "immunoglobulin mediated immune response", "regulation of mitotic nuclear division")
  ),
  lusc = list(
    'Dosage-sensitive' = c("ribosome biogenesis", "rRNA metabolic process", "protein folding"),
    'Dosage-insensitive' = c("cell chemotaxis", "myeloid leukocyte activation", "mesenchymal cell differentiation")
  ),
  brca = list(
    'Dosage-sensitive' = c("DNA-templated DNA replication"),
    'Dosage-compensated' = c("striated muscle cell differentiation", "regulation of muscle system process"),
    'Dosage-insensitive' = c("ERK1 and ERK2 cascade", "regulation of muscle contraction", "mitotic nuclear division")
  ),
  lihc = list(
    'Dosage-compensated' = c("alpha-beta T cell activation", "leukocyte activation involved in immune response", "negative regulation of cell adhesion"),
    'Dosage-insensitive' = c("regulation of mitotic nuclear division", "nuclear chromosome segregation", "alpha-amino acid catabolic process")
  )
)


process_pancancer_data <- function(tumor, category, ora_GO_pancancer) {
  if (!is.null(ora_GO_pancancer[[tumor]]) && !is.null(ora_GO_pancancer[[tumor]][[category]])) {
    terms <- GO_terms[[tumor]][[category]]
    if (!is.null(terms)) {
      ora_GO_pancancer[[tumor]][[category]] %>%
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
    processed_data <- process_pancancer_data(tumor, category, ora_GO_pancancer)
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
                   brca = "BRCA",
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

ggsave("deconveilCaseStudies/plots/main/GO_dotplot_pancancer.png", dpi = 400, width = 22.0, height = 5.0, plot = plot_GO)

  
