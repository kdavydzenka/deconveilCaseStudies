setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/")

pkgs <- c("dplyr", "ggplot2", "cluster", "factoextra", "heatmaply", "tidyverse", "DESeq2", "colorspace", 
          "ggpubr", "ggpointdensity", "ggeasy", "ggrepel", "gridExtra", "igraph", "ggraph")
sapply(pkgs, require, character.only = TRUE)


# TCGA-BRCA #

clinical <- readRDS("TCGA/brca/clinical_filt.RDS")
rna_normal <- readRDS("TCGA/brca/rna_normal.RDS")
rna_tumor <- readRDS("TCGA/brca/rna_tumor.RDS")
cnv_tumor <- readRDS("TCGA/brca/cnv_tumor.RDS")


colnames(cnv_tumor) <- substr(colnames(cnv_tumor), 1, 12)
colnames(rna_normal) <- substr(colnames(rna_normal), 1, 12)
colnames(rna_tumor) <- substr(colnames(rna_tumor), 1, 12)

colnames(clinical) <- c("patientID", "gender", "age", "stage")

rna_tumor <- rna_tumor[,(colnames(rna_tumor) %in% clinical$patientID)]
rna_normal <- rna_normal[,(colnames(rna_normal) %in% clinical$patientID)]
cnv_tumor <- cnv_tumor[,(colnames(cnv_tumor) %in% clinical$patientID)]
clinical <- clinical[clinical$patientID %in% colnames(cnv_tumor),]

#cnv_tumor <- apply(cnv_tumor, 2, function(x) ifelse(x > 10, 10, x))
#cnv_tumor <- as.data.frame(cnv_tumor)

# Exclude genes with low expression in normal tissue #
low_expression_threshold <- 20
expression_summary <- data.frame(
  Gene = rownames(rna_normal),
  MeanExpression = rowMeans(rna_normal)
)

filtered_genes <- expression_summary %>%
  filter(MeanExpression > low_expression_threshold)

rna_normal <- rna_normal[filtered_genes$Gene, ]
rna_tumor <- rna_tumor[filtered_genes$Gene, ]
cnv_tumor <- as.data.frame(cnv_tumor)
cnv_tumor <- cnv_tumor[filtered_genes$Gene, ]

x <- colnames(rna_normal)
names(rna_normal) <- paste(x,"-11A")

x <- colnames(rna_tumor)
names(rna_tumor) <- paste(x,"-01A")

x <- colnames(cnv_tumor)
names(cnv_tumor) <- paste(x,"-01A")

rna <- cbind(rna_normal, rna_tumor)


# Gene expression variability check

gene_iqr <- apply(rna, 1, IQR)

variability_summary <- data.frame(
  Gene = rownames(rna),
  IQR = gene_iqr
)
iqr_threshold <- quantile(variability_summary$IQR, 0.25)

filtered_genes_iqr <- variability_summary %>%
  filter(IQR > iqr_threshold)

rna_filt <- rna[filtered_genes_iqr$Gene, ]
cnv_tumor <- cnv_tumor[filtered_genes_iqr$Gene, ]

cnv_tumor <- apply(cnv_tumor, 2, function(x) ifelse(x > 15, 15, x))

hist(rowMeans(cnv_tumor),
     main = "", 
     xlab = "CN state",
     ylab = "Proportion",
     col = "#E1DEFC",
     prob = TRUE,
     breaks = 15)

cnv_mean <- cnv_tumor %>% 
  as.data.frame() %>% 
  dplyr::mutate(cnv_mean = rowMeans(cnv_tumor)) %>% 
  dplyr::select(cnv_mean)

cnv_normal <- matrix(2, nrow(rna_filt), ncol(rna_normal))
rownames(cnv_normal) <- rownames(cnv_tumor)
cnv <- cbind(cnv_normal, cnv_tumor)
colnames(cnv) <- colnames(rna_filt)
cnv <- cnv/2
colnames(rna_filt) <- paste0("sample", 1:(ncol(rna_filt)))
colnames(cnv) <- colnames(rna_filt)

#Generate metadata#
metadata <- data.frame(patID = colnames(rna_filt),
                       condition = rep(c("A", "B"), each = ncol(rna_tumor)))
metadata <- metadata %>% remove_rownames %>% column_to_rownames(var = "patID") 
metadata$condition <- as.factor(metadata$condition)

write.csv(cnv, file = "TCGA/brca/case_study/cnv.csv", row.names = T)
write.csv(rna_filt, file = "TCGA/brca/case_study/rna.csv", row.names = T)
write.csv(metadata, file = "TCGA/brca/case_study/metadata.csv", row.names = T)


### Downstaream analysis ###

res_naive <- read.csv("CN-aware-DGE/results/case_studies/BRCA/res_CNnaive.csv")
res_aware <- read.csv("CN-aware-DGE/results/case_studies/BRCA/res_CNaware.csv")
cnv <- read.csv("TCGA/brca/case_study/cnv.csv") %>% remove_rownames %>% column_to_rownames(var="X")

cnv_tumor <- cnv[,110:218]
cnv_tumor <- cnv_tumor * 2

cnv_mean <- cnv_tumor %>% 
  as.data.frame() %>% 
  dplyr::mutate(cnv_mean = rowMeans(cnv_tumor)) %>% 
  dplyr::mutate(geneID = rownames(cnv_tumor)) %>% 
  dplyr::select(geneID, cnv_mean) 


# Gene groups selection: Dosage-sensitive | Dosage-insensitive | Dosage-compensated

lfc_cut <- 1.0
pval_cut <- .05

res_aware <- res_aware %>%
  dplyr::mutate(isDE = (abs(log2FoldChange) >= lfc_cut) & (padj <= pval_cut)) %>%
  dplyr::mutate(DEtype = if_else(!isDE, "n.s.", if_else(log2FoldChange > 0, "Up-reg", "Down-reg"))) %>%
  dplyr::mutate(method = "CN aware") %>%
  dplyr::select(X,log2FoldChange, padj, isDE, DEtype, method) %>% 
  remove_rownames %>% 
  column_to_rownames(var="X")

res_naive <- res_naive %>%
  dplyr::mutate(isDE = (abs(log2FoldChange) >= lfc_cut) & (padj <= pval_cut)) %>%
  dplyr::mutate(DEtype = if_else(!isDE, "n.s.", if_else(log2FoldChange > 0, "Up-reg", "Down-reg"))) %>%
  dplyr::mutate(method = "CN naive") %>% 
  dplyr::select(X,log2FoldChange, padj, isDE, DEtype, method) %>% 
  remove_rownames %>% 
  column_to_rownames(var="X")


res_aware <- cbind(res_aware, cnv_mean)
res_naive <- cbind(res_naive, cnv_mean)

res_aware <- res_aware %>% dplyr::select(log2FoldChange, padj, isDE, DEtype, method, cnv_mean) 
res_naive <- res_naive %>% dplyr::select(log2FoldChange, padj, isDE, DEtype, method, cnv_mean)

colnames(res_naive) <- c("logFC_naive", "padj_naive", "isDE_naive", "DEtype_naive", "method_naive", "cnv_mean_naive")
colnames(res_aware) <- c("logFC_aware", "padj_aware", "isDE_aware", "DEtype_aware", "method_aware", "cnv_mean_aware")

res_joint <- cbind(res_naive, res_aware)

d_sensitive <- res_joint %>%
  dplyr::filter(DEtype_naive == "Up-reg" & DEtype_aware == "n.s." | 
                  DEtype_naive == "Down-reg" & DEtype_aware == "n.s") 

d_insensitive <- res_joint %>% 
  dplyr::filter(DEtype_naive == "Down-reg" & DEtype_aware == "Down-reg" |
                  DEtype_naive == "Up-reg" & DEtype_aware == "Up-reg") 

d_compensated <- res_joint %>% 
  dplyr::filter(DEtype_naive == "n.s." & DEtype_aware == "Down-reg" | 
                  DEtype_naive == "n.s." & DEtype_aware == "Up-reg") 

non_deg <- res_joint %>% 
  dplyr::filter(DEtype_naive == "n.s." & DEtype_aware == "n.s.")
 

# CN-naive #
cn_naive_d_sensitive <- d_sensitive %>% dplyr::select(logFC_naive, padj_naive, isDE_naive, DEtype_naive, method_naive, cnv_mean_naive) %>% 
  dplyr::mutate(gene_group = "Dosage-sensitive")
colnames(cn_naive_d_sensitive) <- c("log2FC", "padj", "isDE", "DEtype", "method", "cnv_mean", "gene_group")

cn_naive_d_insensitive <- d_insensitive %>% dplyr::select(logFC_naive, padj_naive, isDE_naive, DEtype_naive, method_naive, cnv_mean_naive) %>% 
  dplyr::mutate(gene_group = "Dosage-insensitive")
colnames(cn_naive_d_insensitive) <- c("log2FC", "padj", "isDE", "DEtype", "method", "cnv_mean", "gene_group")

cn_naive_d_compensated <- d_compensated %>% dplyr::select(logFC_naive, padj_naive, isDE_naive, DEtype_naive, method_naive, cnv_mean_naive) %>% 
  dplyr::mutate(gene_group = "Dosage-compensated")
colnames(cn_naive_d_compensated) <- c("log2FC", "padj", "isDE", "DEtype", "method", "cnv_mean", "gene_group")

cn_naive_non_DE <- non_deg %>% dplyr::select(logFC_naive, padj_naive, isDE_naive, DEtype_naive, method_naive, cnv_mean_naive) %>% 
  dplyr::mutate(gene_group = "non-DEG")
colnames(cn_naive_non_DE) <- c("log2FC", "padj", "isDE", "DEtype", "method", "cnv_mean", "gene_group")

# CN-aware #
cn_aware_d_sensitive <- d_sensitive %>% dplyr::select(logFC_aware, padj_aware, isDE_aware, DEtype_aware, method_aware, cnv_mean_aware) %>% 
  dplyr::mutate(gene_group = "Dosage-sensitive")
colnames(cn_aware_d_sensitive) <- c("log2FC", "padj", "isDE", "DEtype", "method", "cnv_mean", "gene_group")

cn_aware_d_insensitive <- d_insensitive %>% dplyr::select(logFC_aware, padj_aware, isDE_aware, DEtype_aware, method_aware, cnv_mean_aware) %>% 
  dplyr::mutate(gene_group = "Dosage-insensitive")
colnames(cn_aware_d_insensitive) <- c("log2FC", "padj", "isDE", "DEtype", "method", "cnv_mean", "gene_group")

cn_aware_d_compensated <- d_compensated %>% dplyr::select(logFC_aware, padj_aware, isDE_aware, DEtype_aware, method_aware, cnv_mean_aware) %>% 
  dplyr::mutate(gene_group = "Dosage-compensated")
colnames(cn_aware_d_compensated) <- c("log2FC", "padj", "isDE", "DEtype", "method", "cnv_mean", "gene_group")

cn_aware_non_DE <- non_deg %>% dplyr::select(logFC_aware, padj_aware, isDE_aware, DEtype_aware, method_aware, cnv_mean_aware) %>% 
  dplyr::mutate(gene_group = "non-DEG")
colnames(cn_aware_non_DE) <- c("log2FC", "padj", "isDE", "DEtype", "method", "cnv_mean", "gene_group")

saveRDS(d_compensated, file = "TCGA/brca/case_study/d_compensated_genes_2.RDS")


cn_naive <- rbind(cn_naive_d_sensitive, cn_naive_d_insensitive, cn_naive_d_compensated, cn_naive_non_DE)
cn_aware <- rbind(cn_aware_d_sensitive, cn_aware_d_insensitive, cn_aware_d_compensated, cn_aware_non_DE)
#cn_naive <- cn_naive[!(row.names(cn_naive) %in% c("PYCR1")),] #LUAD

v_plot_data <- rbind(cn_naive, cn_aware)

v_plot_data <- v_plot_data %>% dplyr::filter(log2FC > -6.0 ,)

gene_group_colors <- c("Dosage-insensitive" = "#8F3931FF", "Dosage-sensitive" = "#FFB977", "Dosage-compensated"="#FAE48BFF", "non-DEG" = "#ADB6B6FF")  
cnv_colors <- c("loss" = "#0073C299", "neutral" = "#86868699", "gain" = "#cecb76", "amplification" = "#DC0000B2")

p_volcanos <- v_plot_data %>%
  ggplot(mapping = aes(x = log2FC, y = -log10(padj))) +
  geom_point(data = subset(v_plot_data, gene_group %in% c("Dosage-insensitive", "non-DEG")),
             aes(col = gene_group), size = 1.0, alpha = 0.3) +
  geom_point(data = subset(v_plot_data, gene_group %in% c("Dosage-sensitive", "Dosage-compensated")),
             aes(col = gene_group), size = 2.0, alpha = 0.5) +
  scale_color_manual(values = gene_group_colors) +
  theme_bw() +
  scale_x_continuous(breaks = seq(floor(min(v_plot_data$log2FC)), 
                                  ceiling(max(v_plot_data$log2FC)), by = 2)) +
  
  facet_wrap(~factor(method, levels = c("CN naive", "CN aware")), nrow = 1) +
  labs(x = expression(Log[2] ~ FC), y = expression(-log[10] ~ Pvalue), col = "Gene group") +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = 'dashed') +
  geom_hline(yintercept = -log10(pval_cut), linetype = "dashed") +
  ggplot2::theme(legend.position = 'left',
                 legend.text = element_text(size = 14, color = "black"),
                 legend.title = element_text(size = 16, color = "black"),  
                 strip.text = element_text(size = 16, face = "plain", color = "black"),
                 axis.text = element_text(size = 14, color = "black"),
                 axis.title = element_text(size = 16))+
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1)))
p_volcanos


# CN barplot
combined_data <- rbind(cn_aware_d_sensitive, cn_aware_d_insensitive, cn_aware_d_compensated)

classify_cn <- function(cn_value) {
  if (cn_value == 0 || cn_value == 1) {
    return("Loss")
  } else if (cn_value == 2) {
    return("Neutral")
  } else if (cn_value == 3 || cn_value == 4) {
    return("Gain")
  } else if (cn_value > 4) {
    return("Amplification")
  } else {
    return(NA)  
  }
}

cn_categories <- apply(cnv_tumor, c(1, 2), classify_cn)
loss_proportion <- apply(cn_categories, 1, function(x) mean(x == "Loss"))
loss_proportion <- as.data.frame(loss_proportion)

loss_threshold <- 0.25
loss_labels <- loss_proportion %>% 
  dplyr::mutate(isCNloss = case_when(
    loss_proportion > loss_threshold ~ "loss",
    loss_proportion < loss_threshold ~ "not loss"))

combined_data <- combined_data %>% 
  dplyr::mutate(cnv_group = case_when(
    cnv_mean > 0.0 & cnv_mean <= 1.7  ~ "loss",
    cnv_mean > 1.7 & cnv_mean <= 2.5  ~ "neutral",
    cnv_mean > 2.5 & cnv_mean <=   4.0 ~ "gain",
    cnv_mean > 4.0 ~ "amplification"))

loss_labels <- loss_labels[rownames(loss_labels) %in% rownames(combined_data), ]

rownames_idx <- match(rownames(combined_data), rownames(loss_labels))
loss_labels <- loss_labels[rownames_idx,] %>% na.omit()

combined_data <- cbind(combined_data, loss_labels)
combined_data$cnv_group <- ifelse(combined_data$isCNloss == "loss", "loss", combined_data$cnv_group)

barplot_data <- combined_data %>%
  group_by(gene_group) %>%
  summarise(Count = n()) %>%
  mutate(total = sum(Count)) %>%
  mutate(percentage = (Count / total) * 100) %>%
  ungroup()

combined_data$gene_group <- factor(combined_data$gene_group, levels = c("Dosage-insensitive", "Dosage-sensitive", "Dosage-compensated"))
barplot_data$gene_group <- factor(barplot_data$gene_group, levels = c("Dosage-insensitive", "Dosage-sensitive", "Dosage-compensated"))

barplot_stat <- ggplot2::ggplot(barplot_data, aes(x = gene_group, y = percentage, fill = gene_group)) +
  geom_bar(stat = "identity", color = "black", width = 0.6, alpha = 0.7) + 
  geom_text(aes(label = Count), position = position_stack(vjust = 0.5), size = 4, color = "black") +  
  scale_fill_manual(values = gene_group_colors) +  
  theme_classic() +  
  labs(y = "percentage of genes", x = "", title = "", fill = "Gene group") +  
  theme(
    axis.text.x = element_text(size = 16, angle = 30, hjust = 1, color = "black"),  
    axis.text.y = element_text(size = 16, color = "black"),                         
    axis.title.x = element_text(size = 16, face = "plain", color = "black"),          
    axis.title.y = element_text(size = 16, face = "plain", color = "black"),          
    legend.position = 'left',
    legend.text = element_text(size = 14, color = "black"),                          
    legend.title = element_text(size = 16, face = "plain")           
  )

barplot_stat



barplot_cnv <- ggplot2::ggplot(combined_data, aes(x = gene_group, fill = cnv_group)) +
  geom_bar(position = "stack", width = 0.6) + 
  #geom_text(aes(label = Count), position = position_stack(vjust = 0.5), size = 4) +  
  scale_fill_manual(values = cnv_colors) +  
  theme_classic() +  
  labs(y = "gene counts", x = "", title = "", fill = "CN group") +  
  theme(
    axis.text.x = element_text(size = 16, angle = 30, hjust = 1, color = "black"),  
    axis.text.y = element_text(size = 16, color = "black"),                         
    axis.title.x = element_text(size = 16, face = "plain", color = "black"),          
    axis.title.y = element_text(size = 16, face = "plain", color = "black"),          
    legend.position = 'left',
    legend.text = element_text(size = 16),                          
    legend.title = element_text(size = 18, face = "plain")           
  )
barplot_cnv


ggsave("CN-aware-DGE/case_studies/plots/brca/barplot_cnv.png", dpi = 400, width = 5.0, height = 4.0, plot = barplot_cnv)



# Enrichment analysis #

d_compensated <- d_compensated %>% dplyr::select(logFC_aware)
d_compensated$geneID <- rownames(d_compensated)

pkgs <- c("ggplot2", "dplyr","tidyr","reactome.db", "fgsea", "org.Hs.eg.db", "data.table", "clusterProfiler", "enrichplot", "ggpubr", "msigdbr")
sapply(pkgs, require, character.only = TRUE)

# Prepare data
hs <- org.Hs.eg.db
my_symbols <- d_compensated$geneID
gene_list <- AnnotationDbi::select(hs,
                                   keys = my_symbols,
                                   columns = c("ENTREZID", "SYMBOL"),
                                   keytype = "SYMBOL")
# Overrapresentation analysis #
gene_l <- as.vector(gene_list$ENTREZID)

# GO
oraGO <- enrichGO(gene = gene_l, ont = "BP", OrgDb = org.Hs.eg.db, 
                  keyType = "ENTREZID", pvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 350)

res_ora_GO_compensated <- oraGO@result %>% mutate(gene_group = "Dosage-compensated")


# KEGG
#kegg_organism = "hsa"
#oraKEGG <- enrichKEGG(gene = gene_l, organism = kegg_organism,
                      #minGSSize = 10, maxGSSize = 350, pvalueCutoff = 0.05, keyType = "ncbi-geneid")
#res_kegg_KEGG_compensated <- oraKEGG@result %>% mutate(gene_group = "Dosage-compensated")


# Visualize Gene - Biological term network

convert_entrez_to_symbol <- function(entrez_ids) {
  entrez_list <- unlist(strsplit(entrez_ids, split = "/"))
  gene_symbols <- mapIds(org.Hs.eg.db, entrez_list, "SYMBOL", "ENTREZID")
  return(paste(gene_symbols, collapse = ", "))
}

res_ora_GO_compensated$geneSymbol <- sapply(res_ora_GO_compensated$geneID, convert_entrez_to_symbol)

GO_path_compensated <- c("regulation of calcium ion transport", "regulation of Wnt signaling pathway", "regulation of small GTPase mediated signal transduction", 
                           "mesenchymal cell differentiation", "fatty acid derivative metabolic process", "phagocytosis")

GO_results <- res_ora_GO_compensated %>% filter(Description %in% GO_path_compensated)

data <- GO_results %>%
  dplyr::select(Description, geneSymbol) %>%   
  tidyr::separate_rows(geneSymbol, sep = ",") %>% 
  dplyr::rename(term = Description, gene = geneSymbol)

data$type <- ifelse(data$gene %in% data$gene, "gene", "term")

prognostic_genes <- c(" CACNB1", " DCDC2", " FGD5", " FGFR1", " OXCT1", " TUB")


g <- graph_from_data_frame(data, directed = FALSE)
V(g)$type <- ifelse(V(g)$name %in% data$gene, "gene", "term")

# Define a custom label: label only biological terms and specific genes of interest
V(g)$label <- ifelse(V(g)$type == "term" | V(g)$name %in% prognostic_genes, V(g)$name, NA)


ggraph(g, layout = "fr") +  
  geom_edge_link(aes(edge_alpha = 0.3), color = "gray", show.legend = FALSE) +  
  #geom_node_point(aes(color = type, size = ifelse(type == "term", 7, 4))) +
  geom_node_point(aes(color = type), size = 4) +  
  geom_node_text(aes(label = label, 
                     fontface = ifelse(type == "gene", "bold", "plain")),  
                 repel = TRUE, 
                 size = 5) +  
  labs(title = "", color = "Node type") +
  theme_void() +
  scale_color_manual(values = c("gene" = "#0073C299", "term" = "#E64B35B2")) +
  ggplot2::theme(legend.position = 'bottom',
                 legend.text = element_text(size = 14, color = "black"),
                 legend.title = element_text(size = 16, color = "black"))







