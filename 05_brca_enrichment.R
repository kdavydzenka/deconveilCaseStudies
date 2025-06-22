setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/")

pkgs <- c("tidyverse", "circlize", "RColorBrewer", "org.Hs.eg.db", "ggrepel")
sapply(pkgs, require, character.only = TRUE)
source("deconveilCaseStudies/utils/utils.R")


# TCGA-BRCA #

# Visualize prognostic genes 

# Load ORA results and extract components

ora_GO <- readRDS("deconveilCaseStudies/plots/main/Fig 5/rds/brca.RDS")

go_categories <- c("Dosage-sensitive", "Dosage-insensitive", "Dosage-compensated")

convert_to_dataframe <- function(category) {
  df <- ora_GO[[category]] %>% as.data.frame()
  df$geneSymbol <- sapply(df$geneID, convert_entrez_to_symbol)
  return(df)
}

ora_GO_data <- lapply(go_categories, convert_to_dataframe)
names(ora_GO_data) <- go_categories

dcomp <- ora_GO_data[["Dosage-compensated"]]

# Define GO terms for filtering
go_terms <- list(
  "Dosage-sensitive" = c("exocytosis", "regulation of DNA replication", "production of molecular mediator of immune response",
                         "negative regulation of cell cycle process", "regulation of mitotic nuclear division"),
  "Dosage-insensitive" = c("cell-cell adhesion via plasma-membrane adhesion molecules", "humoral immune response",
                           "myeloid leukocyte activation", "immunoglobulin mediated immune response",
                           "hormone metabolic process", "positive regulation of protein kinase activity",
                           "cAMP-mediated signaling", "response to lipopolysaccharide"),
  "Dosage-compensated" = c("regulation of muscle system process", "regulation of insulin secretion", "adenylate cyclase-modulating G protein-coupled receptor signaling pathway",
                           "cellular response to salt")
)

# Filter GO results based on terms
filtered_GO_results <- lapply(names(go_terms), function(category) {
  ora_GO_data[[category]] %>% filter(Description %in% go_terms[[category]])
})
names(filtered_GO_results) <- names(go_terms)

# Prepare data for chord diagram
data_long <- filtered_GO_results[["Dosage-compensated"]] %>%
  dplyr::select(Description, geneSymbol) %>%
  tidyr::separate_rows(geneSymbol, sep = ",") %>%
  dplyr::rename(term = Description, gene = geneSymbol)
data_long$gene <- trimws(data_long$gene)

# Load prognostic genes and filter relevant data
prognostic_files <- list(
  "Dosage-sensitive" = "deconveilcaseStudies/plots/main/Fig 5/rds/lasso_dsg.RDS",
  "Dosage-insensitive" = "deconveilcaseStudies/plots/main/Fig 5/rds/lasso_dig.RDS",
  "Dosage-compensated" = "deconveilcaseStudies/plots/main/Fig 5/rds/lasso_dcg.RDS"
)

prognostic_genes <- lapply(prognostic_files, function(file) {
  readRDS(file)$Gene %>% as.vector()
})

filtered_data <- data_long %>%
  dplyr::filter(gene %in% prognostic_genes[["Dosage-compensated"]])

# Create gene-GO matrix
gene_go_matrix <- filtered_data %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = term, values_from = value, values_fill = 0) %>%
  column_to_rownames(var = "gene")
gene_go_matrix <- t(as.matrix(gene_go_matrix))

# Chord diagram
color_genes <- rep("darkgray", length(colnames(gene_go_matrix)))
#color_genes <- brewer.pal(n = 3, name = "Pastel1")  
color_go_terms <- brewer.pal(n = 7, name = "Set2")  

all_colors <- c(color_go_terms, color_genes)

chordDiagram(
  gene_go_matrix,
  grid.col = all_colors,            
  transparency = 0.5,               
  annotationTrack = c("grid"),         
  preAllocateTracks = list(track.height = 0.1)
)

go_terms <- rownames(gene_go_matrix)
gene_labels <- colnames(gene_go_matrix)

circos.track(track.index = 1, panel.fun = function(x, y) {
  # Gene labels (display labels for gene sectors only)
  if (CELL_META$sector.index %in% gene_labels) {
    gene_color <- color_genes[which(gene_labels == CELL_META$sector.index)]
    circos.text(
      CELL_META$xcenter, CELL_META$ycenter,
      labels = CELL_META$sector.index,
      facing = "bending", adj = c(0, 0.3),
      col = "black", cex = 1.0  
    )
  }
  # Skip GO term labels: do nothing for GO term sectors
  if (CELL_META$sector.index %in% go_terms) {
    return()  # No action for GO terms, so no labels will be drawn
  }
}, bg.border = NA)

# DSGs legend
legend("bottom", 
       legend = c("exocytosis", "negative regulation of cell cycle process",
                  "regulation of mitotic nuclear division", "production of molecular mediator of immune response"), 
       fill = color_go_terms, 
       cex = 0.7)

#D DIGs legend
legend("bottom", 
       legend = c("hormone metabolic process", "humoral immune response",
                  "cell-cell adhesion via plasma-membrane adhesion molecules", "positive regulation of protein kinase activity",
                  "response to lipopolysaccharide", "cAMP-mediated signaling", "myeloid leukocyte activation", "immunoglobulin mediated immune response"), 
       fill = color_go_terms, 
       cex = 0.7)

# DCGs legend
legend("bottom", 
       legend = c("regulation of muscle system process", "phagocytosis",
                  "positive regulation of epithelial cell proliferation",
                  "cellular response to salt", "mesenchymal cell differentiation",
                  "cell-matrix adhesion", "positive regulation of epithelial cell migration"), 
       fill = color_go_terms, 
       cex = 0.7)



# Volcano plot of prognostic genes

gene_groups <- readRDS("TCGA/BRCA/case_study/gene_groups_brca.RDS")

extract_genes <- function(gene_group_df, gene_group_name) {
  gene_group_df %>%
    mutate(gene_group = gene_group_name) %>%
    select(geneID = geneID_aware, tumor_type = tumor_type_aware, logFC_aware, padj_aware, gene_group)
}

d_sensitive_aware <- extract_genes(gene_groups$d_sensitive, "DSGs")
d_insensitive_aware <- extract_genes(gene_groups$d_insensitive, "DIGs")
d_compensated_aware <- extract_genes(gene_groups$d_compensated, "DCGs")

d_compensated_aware <- d_compensated_aware %>% dplyr::filter(abs(logFC_aware) < 6.0 ,)

gene_group_colors <- c("DIGs" = "#8F3931FF", "DSGs" = "#FFB977", "DCGs"="#FAE48BFF", "non-DEG" = "#ADB6B6FF")  
lfc_cut <- 1.0
pval_cut <- 0.05

DSGs <- prognostic_genes[["Dosage-sensitive"]]
DIGs <- prognostic_genes[["Dosage-insensitive"]]
DCGs <- prognostic_genes[["Dosage-compensated"]]

p_volcanos <- d_compensated_aware %>%
  ggplot(mapping = aes(x = logFC_aware, y = -log10(padj_aware))) +
  geom_point(aes(col = gene_group), size = 2.0, alpha = 0.5) +
  scale_color_manual(values = gene_group_colors) +
  geom_label_repel(
    data = d_compensated_aware %>% filter(geneID %in% DCGs),
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
  scale_x_continuous(breaks = seq(floor(min(d_compensated_aware$logFC_aware)), 
                                  ceiling(max(d_compensated_aware$logFC_aware)), by = 2)) +
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

ggsave("deconveilCaseStudies/plots/supplementary/volcano_dcg_brca.png", dpi = 400, width = 5.0, height = 5.0, plot = p_volcanos)
