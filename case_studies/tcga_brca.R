setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/")

pkgs <- c("tidyverse", "circlize", "RColorBrewer", "org.Hs.eg.db")
sapply(pkgs, require, character.only = TRUE)
source("deconveilCaseStudies/utils.R")


# TCGA-BRCA #

# Visualize prognostic genes 

# Load ORA results and extract components

ora_GO <- readRDS("deconveilCaseStudies/results/pancancer_ora_GO/brca.RDS")

go_categories <- c("Dosage-sensitive", "Dosage-insensitive", "Dosage-compensated")

convert_to_dataframe <- function(category) {
  df <- ora_GO[[category]] %>% as.data.frame()
  df$geneSymbol <- sapply(df$geneID, convert_entrez_to_symbol)
  return(df)
}

ora_GO_data <- lapply(go_categories, convert_to_dataframe)
names(ora_GO_data) <- go_categories

# Define GO terms for filtering
go_terms <- list(
  "Dosage-sensitive" = c("exocytosis", "regulation of DNA replication", "production of molecular mediator of immune response",
                         "negative regulation of cell cycle process", "regulation of mitotic nuclear division"),
  "Dosage-insensitive" = c("cell-cell adhesion via plasma-membrane adhesion molecules", "humoral immune response",
                           "myeloid leukocyte activation", "immunoglobulin mediated immune response",
                           "hormone metabolic process", "positive regulation of protein kinase activity",
                           "cAMP-mediated signaling", "response to lipopolysaccharide"),
  "Dosage-compensated" = c("regulation of muscle system process", "regulation of insulin secretion", "positive regulation of epithelial cell migration",
                           "cell-matrix adhesion", "cellular response to salt", "phagocytosis", "positive regulation of epithelial cell proliferation",
                           "mesenchymal cell differentiation")
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
  "Dosage-sensitive" = "TCGA/BRCA/case_study/prognostic_signature_ds.RDS",
  "Dosage-insensitive" = "TCGA/BRCA/case_study/prognostic_signature_dins.RDS",
  "Dosage-compensated" = "TCGA/BRCA/case_study/prognostic_signature_dcomp.RDS"
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






