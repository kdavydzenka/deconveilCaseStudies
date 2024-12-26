setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/")

pkgs <- c("tidyverse", "ggplot2", "ggrepel")
sapply(pkgs, require, character.only = TRUE)

# TCGA-LUSC #

# Preprocess the data to be tested with DeConveil

rna_normal <- readRDS("TCGA/lung/LUSC/rna_normal.RDS")
rna_tumor <- readRDS("TCGA/lung/LUSC/rna_tumor.RDS")
cnv_tumor <- readRDS("TCGA/lung/LUSC/cnv_tumor.RDS")

colnames(cnv_tumor) <- substr(colnames(cnv_tumor), 1, 12)
colnames(rna_normal) <- substr(colnames(rna_normal), 1, 12)
colnames(rna_tumor) <- substr(colnames(rna_tumor), 1, 12)


# Exclude genes with low expression in normal tissue #
low_expression_threshold <- 10
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

cnv_tumor <- apply(cnv_tumor, 2, function(x) ifelse(x > 15, 15, x))

hist(rowMeans(cnv_tumor),
     main = "", 
     xlab = "CN state",
     ylab = "Frequency",
     col = "#E1DEFC",
     prob = TRUE,
     breaks = 15)

cnv_mean <- cnv_tumor %>% 
  as.data.frame() %>% 
  dplyr::mutate(cnv_mean = rowMeans(cnv_tumor)) %>% 
  dplyr::select(cnv_mean)

cnv_normal <- matrix(2, nrow(rna_normal), ncol(rna_normal))
rownames(cnv_normal) <- rownames(cnv_tumor)
cnv <- cbind(cnv_normal, cnv_tumor)
colnames(cnv) <- colnames(rna)
cnv <- cnv/2
colnames(rna) <- paste0("sample", 1:(ncol(rna)))
colnames(cnv) <- colnames(rna)

#Generate metadata#
metadata <- data.frame(patID = colnames(rna),
                       condition = rep(c("A", "B"), each = ncol(rna_tumor)))
metadata <- metadata %>% remove_rownames %>% column_to_rownames(var = "patID") 
metadata$condition <- as.factor(metadata$condition)

write.csv(cnv, file = "TCGA/lung/LUSC/cnv.csv", row.names = T)
write.csv(rna, file = "TCGA/lung/LUSC/rna.csv", row.names = T)
write.csv(metadata, file = "TCGA/lung/LUSC/metadata.csv", row.names = T)



### Dowstream analysis ###

res_naive <- read.csv("CN-aware-DGE/results/case_studies/LUSC/res_CNnaive.csv")
res_aware <- read.csv("CN-aware-DGE/results/case_studies/LUSC/res_CNaware.csv")
cnv <- read.csv("TCGA/lung/LUSC/cnv.csv")
#cancer_genes <- read.delim("TCGA/lung/cancerGeneList.tsv")
#cosmic_genes <- read.csv("TCGA/lung/Cosmic_gene_list.csv")
hk_genes <- readRDS("TCGA/lung/housekeeping_genes_lung.RDS")
hk_genes <- hk_genes %>% dplyr::select(Gene.Symbol) %>% 
  dplyr::rename(geneID = Gene.Symbol)

tf <- read_delim("CN-aware-DGE/hg38_disease_TFs.txt")

# Prepare data #

#oncogenes <- cancer_genes %>% dplyr::filter(Is.Oncogene=="Yes") %>% 
  #dplyr::select(Hugo.Symbol) %>% 
  #dplyr::rename(geneID = Hugo.Symbol) %>% 
  #dplyr::mutate(gene_type = "Oncogene")

#tsg <- cancer_genes %>% dplyr::filter(Is.Tumor.Suppressor.Gene=="Yes") %>% 
  #dplyr::select(Hugo.Symbol) %>% 
  #dplyr::rename(geneID=Hugo.Symbol) %>% 
  #dplyr::mutate(gene_type = "TSG")

#cancer_genes_oncokb <- rbind(oncogenes, tsg)

#oncogenes_cosmic <- cosmic_genes %>% dplyr::filter(Role.in.Cancer %in% c("oncogene", "oncogene, TSG", "oncogene, fusion", "oncogene, TSG, fusion"))
#tsg_cosmic <- cosmic_genes %>% dplyr::filter(Role.in.Cancer %in% c("TSG", "oncogene, TSG", "TSG, fusion", "oncogene, TSG, fusion"))

#oncogenes_cosmic <- oncogenes_cosmic %>% 
  #dplyr::select(Gene.Symbol) %>% 
  #dplyr::rename(geneID=Gene.Symbol) %>% 
  #dplyr::mutate(gene_type = "Oncogene")

#tsg_cosmic <- tsg_cosmic %>% 
  #dplyr::select(Gene.Symbol) %>% 
  #dplyr::rename(geneID=Gene.Symbol) %>% 
  #dplyr::mutate(gene_type = "TSG")

#cancer_genes_cosmic <- rbind(oncogenes_cosmic, tsg_cosmic)
#cancer_genes_joint <- rbind(cancer_genes_oncokb, cancer_genes_cosmic)

cnv_tumor <- cnv %>% 
  remove_rownames %>% 
  column_to_rownames(var="X") %>%  
  dplyr::select(52:102,)
cnv_tumor <- cnv_tumor * 2


cnv_mean <- cnv_tumor %>% 
  as.data.frame() %>% 
  dplyr::mutate(cnv_mean = rowMeans(cnv_tumor)) %>% 
  dplyr::mutate(geneID = rownames(cnv_tumor)) %>% 
  dplyr::select(geneID, cnv_mean) 

res_naive <- res_naive %>% dplyr::select(X,log2FoldChange, padj) 
res_aware <- res_aware %>% dplyr::select(X,log2FoldChange, padj) 


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


### Dosage-sensitive | Dosage-insensitive | Dosage-compensated ###

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


# Cancer genes per category
#d_sensitive_cancer_g <- d_sensitive[rownames(d_sensitive) %in% cancer_genes_joint$geneID ,]
#d_sensitive_cancer_g <- d_sensitive_cancer_g[!duplicated(rownames(d_sensitive_cancer_g)), ]

#d_compensated_cancer_g <- d_compensated[rownames(d_compensated) %in% cancer_genes_joint$geneID ,]
#d_compensated_cancer_g <- d_compensated_cancer_g[!duplicated(rownames(d_compensated_cancer_g)), ]

#d_insensitive_cancer_g <- d_insensitive[rownames(d_insensitive) %in% cancer_genes_joint$geneID ,]
#d_insensitive_cancer_g <- d_insensitive_cancer_g[!duplicated(rownames(d_insensitive_cancer_g)), ]

#d_insensitive_cosmic <- d_insensitive_cancer_g[rownames(d_insensitive_cancer_g) %in% cancer_genes_cosmic$geneID ,]

# Housekeeping genes per category
d_sensitive_hk <- d_sensitive[rownames(d_sensitive) %in% hk_genes$geneID ,]
d_sensitive_hk <- d_sensitive_hk[!duplicated(rownames(d_sensitive_hk)), ]

d_compensated_hk <- d_compensated[rownames(d_compensated) %in% hk_genes$geneID ,]
d_compensated_hk <- d_compensated_hk[!duplicated(rownames(d_compensated)), ]

d_insensitive_hk <- d_insensitive[rownames(d_insensitive) %in% hk_genes$geneID ,]
d_insensitive_hk <- d_insensitive_hk[!duplicated(rownames(d_insensitive_hk)), ]


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


# Volcano plot #
gene_group_colors <- c("Dosage-insensitive" = "#8F3931FF", "Dosage-sensitive" = "#FFB977", "Dosage-compensated"="#FAE48BFF", "non-DEG" = "#ADB6B6FF")  
cnv_colors <- c("loss" = "#0073C299", "neutral" = "#86868699", "gain" = "#cecb76", "amplification" = "#DC0000B2")

cn_naive <- rbind(cn_naive_d_sensitive, cn_naive_d_insensitive, cn_naive_d_compensated, cn_naive_non_DE)
cn_aware <- rbind(cn_aware_d_sensitive, cn_aware_d_insensitive, cn_aware_d_compensated, cn_aware_non_DE)

v_plot_data <- rbind(cn_naive, cn_aware)
v_plot_data <- v_plot_data %>% dplyr::filter(padj > 1.293780e-155 ,)

#v_plot_data <- v_plot_data %>% dplyr::filter(padj > 3.246376e-220 ,)

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


# Barplot #
#cancer_g <- c(cancer_genes$geneID)
#hk_genes <- c(hk_genes$geneID)

label_genes <- function(gene_list, cancer_genes, hk_genes) {
  ifelse(gene_list %in% cancer_genes, "Cancer genes",
         ifelse(gene_list %in% hk_genes, "Housekeeping genes", "Other genes"))
}

#cn_aware_d_sensitive$gene_subcategory <- label_genes(cn_aware_d_sensitive$geneID, cancer_genes, hk_genes)
#cn_aware_d_insensitive$gene_subcategory <- label_genes(cn_aware_d_insensitive$geneID, cancer_genes, hk_genes)
#cn_aware_d_compensated$gene_subcategory <- label_genes(cn_aware_d_compensated$geneID, cancer_genes, hk_genes)

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

loss_threshold <- 0.30
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
    axis.text.y = element_text(size = 18, color = "black"),                         
    axis.title.x = element_text(size = 16, face = "plain", color = "black"),          
    axis.title.y = element_text(size = 16, face = "plain", color = "black"),          
    legend.position = 'left',
    legend.text = element_text(size = 16),                          
    legend.title = element_text(size = 18, face = "plain")           
  )
barplot_cnv

ggsave("CN-aware-DGE/case_studies/plots/lusc/barplot_cnv.png", dpi = 400, width = 5.0, height = 4.0, plot = barplot_cnv)


#rm(barplot_data)


# Plot HK genes | Dosage-sensitive | Dosage-compensated | CN-naive | CN-aware

d_sens_naive <- d_sensitive_hk %>% dplyr::select(logFC_naive) %>% 
  dplyr::mutate(geneID = rownames(d_sensitive_hk)) %>% 
  dplyr::mutate(method = "CN naive") %>%
  dplyr::mutate(gene_group = "Dosage-sensitive") %>% 
  dplyr::rename(log2FC = logFC_naive)
rownames(d_sens_naive) <- NULL


d_insens_naive <- d_insensitive_hk %>% dplyr::select(logFC_naive) %>% 
  dplyr::mutate(geneID = rownames(d_insensitive_hk)) %>% 
  dplyr::mutate(method = "CN naive") %>%
  dplyr::mutate(gene_group = "Dosage-insensitive") %>% 
  dplyr::rename(log2FC = logFC_naive)
rownames(d_insens_naive) <- NULL

d_comp_naive <- d_compensated_hk %>% dplyr::select(logFC_naive) %>% 
  dplyr::mutate(geneID = rownames(d_compensated_hk)) %>% 
  dplyr::mutate(method = "CN naive") %>%
  dplyr::mutate(gene_group = "Dosage-compensated") %>% 
  dplyr::rename(log2FC = logFC_naive)
rownames(d_comp_naive) <- NULL

#d_sens_naive <- d_sens_naive %>% left_join(cancer_genes, by = "row.names") %>% arrange(log2FC)
#d_insens_naive <- d_insens_naive %>% left_join(cancer_genes, by = "geneID")
#d_comp_naive <- d_comp_naive %>% left_join(cancer_genes, by = "geneID") %>% arrange(log2FC)

#d_insens_naive_30 <- d_insens_naive %>%
  #dplyr::filter(log2FC < -2.13 | log2FC > 2.20) %>% 
  #arrange(log2FC)
#d_insens_naive_30 <- d_insens_naive_30[!duplicated(d_insens_naive_30$geneID), ]


d_sens_aware <- d_sensitive_hk %>% dplyr::select(logFC_aware) %>% 
  dplyr::mutate(geneID = rownames(d_sensitive_hk)) %>% 
  dplyr::mutate(method = "CN aware") %>%
  dplyr::mutate(gene_group = "Dosage-sensitive") %>% 
  dplyr::rename(log2FC = logFC_aware)
rownames(d_sens_aware) <- NULL

d_insens_aware <- d_insensitive_hk %>% dplyr::select(logFC_aware) %>% 
  dplyr::mutate(geneID = rownames(d_insensitive_hk)) %>% 
  dplyr::mutate(method = "CN aware") %>%
  dplyr::mutate(gene_group = "Dosage-insensitive") %>% 
  dplyr::rename(log2FC = logFC_aware)
rownames(d_insens_aware) <- NULL

d_comp_aware <- d_compensated_hk %>% dplyr::select(logFC_aware) %>% 
  dplyr::mutate(geneID = rownames(d_compensated_hk)) %>% 
  dplyr::mutate(method = "CN aware") %>%
  dplyr::mutate(gene_group = "Dosage-compensated") %>% 
  dplyr::rename(log2FC = logFC_aware)
rownames(d_comp_aware) <- NULL

#d_insens_aware_30 <- d_insens_aware[d_insens_aware$geneID %in% d_insens_naive_30$geneID , ]

#d_sens_aware <- d_sens_aware %>% left_join(cancer_genes, by = "geneID") %>% arrange(log2FC)
#d_insens_aware_30 <- d_insens_aware_30 %>% left_join(cancer_genes, by = "geneID") %>% arrange(log2FC)
#d_insens_aware_30 <- d_insens_aware_30[!duplicated(d_insens_aware_30$geneID), ]
#d_comp_aware <- d_comp_aware %>% left_join(cancer_genes, by = "geneID") %>% arrange(log2FC)

#d_sens_aware <- d_sens_aware[!duplicated(d_sens_aware$geneID), ]
#d_comp_aware <- d_comp_aware[!duplicated(d_comp_aware$geneID), ]
#d_sens_naive <- d_sens_naive[!duplicated(d_sens_naive$geneID), ]
#d_comp_naive <- d_comp_naive[!duplicated(d_comp_naive$geneID), ]

#d_insens_aware_30 <- d_insens_aware[d_insens_aware$geneID %in% d_insens_naive_30$geneID , ]
#d_insens_aware_30 <- d_insens_aware_30[!duplicated(d_insens_aware_30$geneID), ]

bargraph_data <- rbind(d_sens_naive, d_insens_naive, d_comp_naive, d_sens_aware, d_insens_aware, d_comp_aware)

bargraph_data <- bargraph_data %>% 
  group_by(gene_group, method) %>%      
  arrange(log2FC, .by_group = TRUE) 


bargraph_data$method <- factor(bargraph_data$method, levels = c("CN naive", "CN aware"))

# Bar graph 
bargraph <- ggplot(bargraph_data, aes(x = reorder(geneID, log2FC), y = log2FC, fill = gene_group)) +
  geom_bar(stat = "identity", color = "black", alpha = 1.0) +  
  geom_hline(yintercept = 0, linetype = "dashed") +  
  theme_classic() +  
  labs(x = "Gene symbol", y = "log2FC") +  
  scale_fill_manual(values = gene_group_colors)+
  facet_wrap(~factor(method, levels = c("CN naive", "CN aware")), nrow=2)+
  #ggh4x::facet_nested(factor(method, levels = c("CN naive", "CN aware"))~factor(gene_group, levels = c("Dosage-sensitive", "Dosage-compensated", "Dosage-insensitive")),
                      #scales = "free", independent = "y")+
  labs(y = "Effect size (log2)", x = "", title = "", fill = "Gene group") +  
  theme(
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1, color = "black"),  
    axis.text.y = element_text(size = 14, color = "black"),                         
    axis.title.x = element_text(size = 16, face = "plain", color = "black"),          
    axis.title.y = element_text(size = 16, face = "plain", color = "black"), 
    strip.text = element_text(size = 16, face = "plain", color = "black"),
    legend.position = 'left',
    legend.text = element_text(size = 14, color = "black"),                          
    legend.title = element_text(size = 16, face = "plain")           
  )
  
bargraph


### Functional Enrichment analysis ###

cn_aware_d_sensitive$geneID <- rownames(cn_aware_d_sensitive)
cn_aware_d_sensitive <- cn_aware_d_sensitive %>% dplyr::select(geneID, log2FC)

cn_aware_d_compensated$geneID <- rownames(cn_aware_d_compensated)
cn_aware_d_compensated <- cn_aware_d_compensated %>% dplyr::select(geneID, log2FC)

cn_aware_d_insensitive$geneID <- rownames(cn_aware_d_insensitive)
cn_aware_d_insensitive <- cn_aware_d_insensitive %>% dplyr::select(geneID, log2FC)

pkgs <- c("reactome.db", "fgsea", "org.Hs.eg.db", "data.table", "clusterProfiler", "enrichplot", "ggpubr", "msigdbr")
sapply(pkgs, require, character.only = TRUE)

# Prepare data
hs <- org.Hs.eg.db
my_symbols <- cn_aware_d_insensitive$geneID
gene_list <- AnnotationDbi::select(hs,
                                   keys = my_symbols,
                                   columns = c("ENTREZID", "SYMBOL"),
                                   keytype = "SYMBOL")

# Overrapresentation analysis #
gene_l <- as.vector(gene_list$ENTREZID)

# GO
oraGO <- enrichGO(gene = gene_l, ont = "BP", OrgDb = org.Hs.eg.db, 
                  keyType = "ENTREZID", pvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 350)

res_ora_GO_sensitive <- oraGO@result %>% mutate(gene_group = "Dosage-sensitive")
res_ora_GO_compensated <- oraGO@result %>% mutate(gene_group = "Dosage-compensated")
res_ora_GO_insensitive <- oraGO@result %>% mutate(gene_group = "Dosage-insensitive")


# KEGG
#kegg_organism = "hsa"
#oraKEGG <- enrichKEGG(gene = gene_l, organism = kegg_organism,
                      #minGSSize = 10, maxGSSize = 350, pvalueCutoff = 0.05, keyType = "ncbi-geneid")

#res_ora_KEGG_sensitive <- oraKEGG@result %>% mutate(gene_group = "Dosage-sensitive")
#res_ora_GO_compensated <- oraGO@result %>% mutate(gene_group = "Dosage-compensated")
#res_ora_KEGG_insensitive <- oraKEGG@result %>% mutate(gene_group = "Dosage-insensitive")

# MSigDb 
m_hallmark <- msigdbr::msigdbr(species = "Homo sapiens", category = "H") 
msig_H <- msigdbr(species = "Homo sapiens", category = "H") %>% dplyr::select(gs_name, entrez_gene)
msig_H <- enricher(gene_l, minGSSize = 10, maxGSSize = 500,
                                     pvalueCutoff = 0.05, pAdjustMethod = "BH", TERM2GENE = msig_H)

res_ora_H_sensitive <- msig_H@result %>% mutate(gene_group = "Dosage-sensitive")
res_ora_H_compensated <- msig_H@result %>% mutate(gene_group = "Dosage-compensated")
res_ora_H_insensitive <- msig_H@result %>% mutate(gene_group = "Dosage-insensitive")



# GO enrichment dotplot

GO_path_sensitive <- c("ribosome biogenesis", "protein-RNA complex organization", "protein processing",
                      "regulation of chromosome organization", "positive regulation of telomere maintenance")

GO_path_compensated <- c("positive regulation of lymphocyte activation", "negative regulation of type II interferon production", 
                         "negative regulation of immune effector process", "negative regulation of cytokine production involved in immune response",
                         "positive regulation of leukocyte differentiation")

GO_path_insensitive <- c("regulation of Wnt signaling pathway", "negative regulation of cell adhesion",
                         "mesenchymal cell differentiation", "platelet activation",  "T cell mediated immunity")


res_GO_sensitive <- res_ora_GO_sensitive %>% filter(Description %in% GO_path_sensitive)
res_GO_compensated <- res_ora_GO_compensated %>% filter(Description %in% GO_path_compensated)
res_GO_insensitive <- res_ora_GO_insensitive%>% filter(Description %in% GO_path_insensitive)

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

p_data <- rbind(res_GO_sensitive, res_GO_compensated, res_GO_insensitive)
#p_data$Description <- fct_reorder(p_data$Description, p_data$gene_group)

saveRDS(p_data, file = "CN-aware-DGE/case_studies/overrapres_gene_categories_lusc.RDS")

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

ggsave("CN-aware-DGE/case_studies/plots/lusc/GO_dotplot.png", dpi = 400, width = 17.0, height = 6.0, plot = p_gse)

# Visualize Gene - Biological term network

# Hallmark pathways
H_path_sensitive <- c("HALLMARK_MYC_TARGETS_V1", "HALLMARK_E2F_TARGETS", "HALLMARK_UNFOLDED_PROTEIN_RESPONSE")
H_path_insensitive <- c("HALLMARK_HYPOXIA", "HALLMARK_MTORC1_SIGNALING")
H_path_compensated <- c("HALLMARK_HYPOXIA")

H_path <- c("HALLMARK_MYC_TARGETS_V1", "HALLMARK_E2F_TARGETS", "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
            "HALLMARK_HYPOXIA", "HALLMARK_MTORC1_SIGNALING")

convert_entrez_to_symbol <- function(entrez_ids) {
  entrez_list <- unlist(strsplit(entrez_ids, split = "/"))
  gene_symbols <- mapIds(org.Hs.eg.db, entrez_list, "SYMBOL", "ENTREZID")
  return(paste(gene_symbols, collapse = ", "))
}

res_ora_H_sensitive$geneSymbol <- sapply(res_ora_H_sensitive$geneID, convert_entrez_to_symbol)
res_ora_H_compensated$geneSymbol <- sapply(res_ora_H_compensated$geneID, convert_entrez_to_symbol)
res_ora_H_insensitive$geneSymbol <- sapply(res_ora_H_insensitive$geneID, convert_entrez_to_symbol)
res_ora_H <- rbind(res_ora_H_sensitive, res_ora_H_compensated, res_ora_H_insensitive)

hk_genes <- c(" VAMP5", " TMBIM4", " VAMP2", " SERINC1", " WSB1", " PPP1R12C", " SNRPA", " SF3B6",
              " MLF2", " SSRP1", " CDK2AP1", " STRAP", " S100A6", " HLA-A", " PSAP", " GPI", " EEF2")



res_H_sensitive <- res_ora_H_sensitive %>% filter(Description %in% H_path_sensitive)
res_H_insensitive <- res_ora_H_insensitive %>% filter(Description %in% H_path_insensitive)
res_H_compensated <- res_ora_H_compensated %>% filter(Description %in% H_path_compensated)
res_H <- rbind(res_H_sensitive, res_H_insensitive, res_H_compensated)

data <- res_H %>%
  dplyr::select(Description, geneSymbol) %>%   
  tidyr::separate_rows(geneSymbol, sep = ",") %>% 
  dplyr::rename(term = Description, gene = geneSymbol)

data$type <- ifelse(data$gene %in% data$gene, "gene", "term")

data$term <- gsub("HALLMARK", "H", data$term)

g <- graph_from_data_frame(data, directed = FALSE)
V(g)$type <- ifelse(V(g)$name %in% data$gene, "gene", "term")

# Define a custom label: label only biological terms and specific genes of interest
V(g)$label <- ifelse(V(g)$type == "term" | V(g)$name %in% hk_genes, V(g)$name, NA)


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
  scale_color_manual(values = c("gene" = "#7AA6DC99", "term" = "#E64B35B2")) +
  ggplot2::theme(legend.position = 'bottom',
                 legend.text = element_text(size = 14),
                 legend.title = element_text(size = 16))


