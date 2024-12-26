setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/")
pkgs <- c("cluster", "factoextra", "heatmaply", "DESeq2", "tidyverse", "colorspace", 
          "ggpubr", "ggpointdensity", "ggeasy", "gridExtra", "metaseqR2", "ggalluvial", "ggridges", "ggforce", "ggparallel", "alluvial")
sapply(pkgs, require, character.only = TRUE)


### Data preprocessing ###

# Input data
#cnv_tumor <- readRDS("TCGA/lung/LUAD/data/cnv_tumor.RDS")
#rna_normal <- readRDS("TCGA/lung/LUAD/data/rna_normal.RDS")
#rna_tumor <- readRDS("TCGA/lung/LUAD/data/rna_tumor.RDS")

cnv_tumor <- readRDS("TCGA/brca/cnv_tumor.RDS")
rna_normal <- readRDS("TCGA/brca/rna_normal.RDS")
rna_tumor <- readRDS("TCGA/brca/rna_tumor.RDS")


# Gene filtering #

# Exclude genes with low expression in normal tissue 
low_expression_threshold <- 20
expression_summary <- data.frame(
  Gene = rownames(rna_normal),
  MeanExpression = rowMeans(rna_normal)
)

filtered_genes <- expression_summary %>%
  filter(MeanExpression > low_expression_threshold)

rna_normal <- rna_normal[filtered_genes$Gene, ]
rna_tumor <- rna_tumor[filtered_genes$Gene, ]
rna <- cbind(rna_normal, rna_tumor)

cnv_tumor <- cnv_tumor[filtered_genes$Gene, ]


# Gene expression variability check

#gene_iqr <- apply(rna, 1, IQR)

#variability_summary <- data.frame(
  #Gene = rownames(rna),
  #IQR = gene_iqr
#)
#iqr_threshold <- quantile(variability_summary$IQR, 0.25)

#filtered_genes_iqr <- variability_summary %>%
  #filter(IQR > iqr_threshold)

#rna_filt <- rna[filtered_genes_iqr$Gene, ]
#cnv_tumor <- cnv_tumor[filtered_genes_iqr$Gene, ]


# Plot CN data

cnv_tumor <- apply(cnv_tumor, 2, function(x) ifelse(x > 15, 15, x))

cnv_mean_coad <- cnv_tumor %>% 
  as.data.frame() %>% 
  dplyr::mutate(cnv_mean = rowMeans(cnv_tumor)) %>% 
  dplyr::mutate(tumor_type = "COAD") %>% 
  dplyr::select(cnv_mean, tumor_type) 

cnv_mean <- rbind(cnv_mean_brca, cnv_mean_lihc, cnv_mean_hnsc, cnv_mean_coad)

hist <- ggplot(cnv_mean, aes(x = cnv_mean)) +
  geom_histogram(binwidth = 0.4, fill = "#F39B7FB2", color = "black") +
  labs(
    title = "",
    x = "CN state",
    y = "Frequency"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 18, color = "black", hjust = 1),
    axis.text.y = element_text(size = 18, color = "black"),
    axis.title.x = element_text(size = 18, face = "plain", color = "black"),
    axis.title.y = element_text(size = 18, face = "plain", color = "black"),
    strip.text = element_text(size = 18, face = "plain", color = "black")
  ) 
  #facet_wrap(~factor(tumor_type, levels = c("BRCA", "LIHC", "HNSC", "COAD")), nrow = 1)
hist

ggsave("CN-aware-DGE/plots/main/hist_luad.png", dpi = 400, width = 4.0, height = 3.5, plot = hist)


cnv_normal <- matrix(2, nrow(rna), ncol(rna_normal))
rownames(cnv_normal) <- rownames(rna)

cnv_tumor <- as.data.frame(cnv_tumor)
cnv_normal <- as.data.frame(cnv_normal)

#cnv <- cbind(cnv_normal, cnv_tumor)
cnv <- merge(cnv_normal, cnv_tumor, by = "row.names")
cnv <- cnv %>% remove_rownames() %>% column_to_rownames("Row.names")
colnames(cnv) <- colnames(rna)
cnv <- cnv/2
colnames(rna) <- paste0("sample", 1:(ncol(rna)))
colnames(cnv) <- colnames(rna)

rownames_idx <- match(rownames(cnv), rownames(rna))
rna <- rna[rownames_idx,] %>% na.omit()


#Generate metadata#
metadata <- data.frame(patID = colnames(rna),
                       condition = rep(c("A", "B"), each = ncol(rna_tumor)))
metadata <- metadata %>% remove_rownames %>% column_to_rownames(var = "patID") 
metadata$condition <- as.factor(metadata$condition)


write.csv(cnv, file = "TCGA/colon/test/cnv_test_all_genes.csv", row.names = T)
write.csv(rna, file = "TCGA/colon/test/rna_test_all_genes.csv", row.names = T)
write.csv(metadata, file = "TCGA/colon/test/metadata_all_genes.csv", row.names = T)


### Downstream analysis ###

res_pydeseq <- read.csv("CN-aware-DGE/results/COAD/res_CNnaive_all_genes.csv")
res_deconveil <- read.csv("CN-aware-DGE/results/COAD/res_CNaware_all_genes.csv")
cnv <- read.csv("TCGA/colon/test/cnv_test_all_genes.csv") %>% remove_rownames %>% column_to_rownames(var="X")

cnv_tumor <- cnv[,13:24]
cnv_tumor <- cnv_tumor * 2

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

#neutral_proportion <- apply(cn_categories, 1, function(x) mean(x == "Neutral"))
#neutral_proportion <- as.data.frame(neutral_proportion)

#gain_proportion <- apply(cn_categories, 1, function(x) mean(x == "Gain"))
#gain_proportion <- as.data.frame(gain_proportion)

#amplification_proportion <- apply(cn_categories, 1, function(x) mean(x == "Amplification"))
#amplification_proportion <- as.data.frame(amplification_proportion)

#cn_proportion <- cbind(loss_proportion, neutral_proportion, gain_proportion, amplification_proportion)

# CN mean
cnv_mean <- cnv_tumor %>% 
  as.data.frame() %>% 
  dplyr::mutate(cnv_mean = rowMeans(cnv_tumor)) %>% 
  dplyr::mutate(geneID = rownames(cnv_tumor)) %>% 
  dplyr::select(geneID, cnv_mean) 

cnv_mean <- cnv_mean %>% 
  dplyr::mutate(cnv_group = case_when(
    cnv_mean > 0 & cnv_mean <= 1.7  ~ "loss",
    cnv_mean > 1.7 & cnv_mean <= 2.5  ~ "neutral",
    cnv_mean > 2.5 & cnv_mean <=   4.0 ~ "gain",
    cnv_mean > 4.0 ~ "amplification"))

cnv_mean <- cbind(cnv_mean, loss_labels)
cnv_mean$cnv_group <- ifelse(cnv_mean$isCNloss == "loss", "loss", cnv_mean$cnv_group)


lfc_cut <- 1.0
pval_cut <- .05

res_deconveil <- res_deconveil %>%
  dplyr::mutate(isDE = (abs(log2FoldChange) >= lfc_cut) & (padj <= pval_cut)) %>%
  dplyr::mutate(DEtype = if_else(!isDE, "n.s.", if_else(log2FoldChange > 0, "Up-reg", "Down-reg"))) %>%
  dplyr::mutate(tumor_type = "COAD") %>% 
  dplyr::mutate(method = "CN aware") %>%
  dplyr::select(X,log2FoldChange, padj, isDE, DEtype, tumor_type, method) %>% 
  remove_rownames %>% 
  column_to_rownames(var="X")

res_pydeseq <- res_pydeseq %>%
  dplyr::mutate(isDE = (abs(log2FoldChange) >= lfc_cut) & (padj <= pval_cut)) %>%
  dplyr::mutate(DEtype = if_else(!isDE, "n.s.", if_else(log2FoldChange > 0, "Up-reg", "Down-reg"))) %>%
  dplyr::mutate(method = "CN naive") %>% 
  dplyr::mutate(tumor_type = "COAD") %>% 
  dplyr::select(X,log2FoldChange, padj, isDE, DEtype, tumor_type, method) %>% 
  remove_rownames %>% 
  column_to_rownames(var="X")

colnames(res_pydeseq) <- c("logFC", "padj", "isDE", "DEtype", "tumor_type", "method")
colnames(res_deconveil) <- colnames(res_pydeseq) 


# Gene groups separation: Dosage-sensitive | Dosage-insensitive | Dosage-compensated

res_deconveil <- cbind(res_deconveil, cnv_mean)
res_pydeseq <- cbind(res_pydeseq, cnv_mean)

res_aware <- res_deconveil %>% dplyr::select(logFC, padj, isDE, DEtype, tumor_type, method, cnv_mean) 
res_naive <- res_pydeseq %>% dplyr::select(logFC, padj, isDE, DEtype, tumor_type, method, cnv_mean)

colnames(res_naive) <- c("logFC_naive", "padj_naive", "isDE_naive", "DEtype_naive", "tumor_type_naive", "method_naive", "cnv_mean_naive")
colnames(res_aware) <- c("logFC_aware", "padj_aware", "isDE_aware", "DEtype_aware", "tumor_type_aware", "method_aware", "cnv_mean_aware")

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
cn_naive_d_sensitive <- d_sensitive %>% dplyr::select(logFC_naive, padj_naive, isDE_naive, DEtype_naive, tumor_type_naive, method_naive, cnv_mean_naive) %>% 
  dplyr::mutate(gene_group = "Dosage-sensitive")
colnames(cn_naive_d_sensitive) <- c("log2FC", "padj", "isDE", "DEtype", "tumor_type", "method", "cnv_mean", "gene_group")

cn_naive_d_insensitive <- d_insensitive %>% dplyr::select(logFC_naive, padj_naive, isDE_naive, DEtype_naive, tumor_type_naive, method_naive, cnv_mean_naive) %>% 
  dplyr::mutate(gene_group = "Dosage-insensitive")
colnames(cn_naive_d_insensitive) <- c("log2FC", "padj", "isDE", "DEtype", "tumor_type", "method", "cnv_mean", "gene_group")

cn_naive_d_compensated <- d_compensated %>% dplyr::select(logFC_naive, padj_naive, isDE_naive, DEtype_naive, tumor_type_naive, method_naive, cnv_mean_naive) %>% 
  dplyr::mutate(gene_group = "Dosage-compensated")
colnames(cn_naive_d_compensated) <- c("log2FC", "padj", "isDE", "DEtype", "tumor_type", "method", "cnv_mean", "gene_group")

cn_naive_non_DE <- non_deg %>% dplyr::select(logFC_naive, padj_naive, isDE_naive, DEtype_naive, tumor_type_naive, method_naive, cnv_mean_naive) %>% 
  dplyr::mutate(gene_group = "non-DEG")
colnames(cn_naive_non_DE) <- c("log2FC", "padj", "isDE", "DEtype", "tumor_type", "method", "cnv_mean", "gene_group")

# CN-aware #
cn_aware_d_sensitive <- d_sensitive %>% dplyr::select(logFC_aware, padj_aware, isDE_aware, DEtype_aware, tumor_type_aware, method_aware, cnv_mean_aware) %>% 
  dplyr::mutate(gene_group = "Dosage-sensitive")
colnames(cn_aware_d_sensitive) <- c("log2FC", "padj", "isDE", "DEtype", "tumor_type", "method", "cnv_mean", "gene_group")

cn_aware_d_insensitive <- d_insensitive %>% dplyr::select(logFC_aware, padj_aware, isDE_aware, DEtype_aware, tumor_type_aware, method_aware, cnv_mean_aware) %>% 
  dplyr::mutate(gene_group = "Dosage-insensitive")
colnames(cn_aware_d_insensitive) <- c("log2FC", "padj", "isDE", "DEtype", "tumor_type", "method", "cnv_mean", "gene_group")

cn_aware_d_compensated <- d_compensated %>% dplyr::select(logFC_aware, padj_aware, isDE_aware, DEtype_aware, tumor_type_aware, method_aware, cnv_mean_aware) %>% 
  dplyr::mutate(gene_group = "Dosage-compensated")
colnames(cn_aware_d_compensated) <- c("log2FC", "padj", "isDE", "DEtype", "tumor_type", "method", "cnv_mean", "gene_group")

cn_aware_non_DE <- non_deg %>% dplyr::select(logFC_aware, padj_aware, isDE_aware, DEtype_aware, tumor_type_aware, method_aware, cnv_mean_aware) %>% 
  dplyr::mutate(gene_group = "non-DEG")
colnames(cn_aware_non_DE) <- c("log2FC", "padj", "isDE", "DEtype", "tumor_type", "method", "cnv_mean", "gene_group")

cn_naive <- rbind(cn_naive_d_sensitive, cn_naive_d_insensitive, cn_naive_d_compensated, cn_naive_non_DE)
cn_aware <- rbind(cn_aware_d_sensitive, cn_aware_d_insensitive, cn_aware_d_compensated, cn_aware_non_DE)

# Volcano plot #

#cn_naive <- cn_naive[!(row.names(cn_naive) %in% c("PYCR1")),] #LUAD

v_plot_data_luad <- rbind(cn_naive, cn_aware)
v_plot_data_luad <- v_plot_data_luad %>% dplyr::filter(log2FC > -5.0 ,)

v_plot_data_brca <- rbind(cn_naive, cn_aware)
v_plot_data_brca <- v_plot_data_brca %>% dplyr::filter(log2FC < 5.5 ,)
v_plot_data_brca <- v_plot_data_brca %>% dplyr::filter(padj > 4.929213e-191 ,)

v_plot_data_lihc <- rbind(cn_naive, cn_aware)
v_plot_data_lihc <- v_plot_data_lihc %>% dplyr::filter(log2FC < 5.0 ,)
v_plot_data_lihc <- v_plot_data_lihc %>% dplyr::filter(padj > 4.332771e-67 ,)

v_plot_data_hnsc <- rbind(cn_naive, cn_aware)
v_plot_data_hnsc <- v_plot_data_hnsc %>% dplyr::filter(log2FC > -5.0 ,)
v_plot_data_hnsc <- v_plot_data_hnsc %>% dplyr::filter(padj > 5.140434e-42 ,)

v_plot_data_coad <- rbind(cn_naive, cn_aware)
v_plot_data_coad <- v_plot_data_coad %>% dplyr::filter(log2FC < 6.0 ,)

v_plot_data <- rbind(v_plot_data_brca, v_plot_data_lihc, v_plot_data_hnsc, v_plot_data_coad)

gene_group_colors <- c("Dosage-insensitive" = "#8F3931FF", "Dosage-sensitive" = "#FFB977", "Dosage-compensated"="#FAE48BFF", "non-DEG" = "#ADB6B6FF")  
cnv_colors <- c("loss" = "#0073C299", "neutral" = "#86868699", "gain" = "#cecb76", "amplification" = "#DC0000B2")

p_volcanos <- v_plot_data_luad %>%
  ggplot(mapping = aes(x = log2FC, y = -log10(padj))) +
  geom_point(data = subset(v_plot_data_luad, gene_group %in% c("Dosage-insensitive", "non-DEG")),
             aes(col = gene_group), size = 1.0, alpha = 0.3) +
  geom_point(data = subset(v_plot_data_luad, gene_group %in% c("Dosage-sensitive", "Dosage-compensated")),
             aes(col = gene_group), size = 2.0, alpha = 0.5) +
  scale_color_manual(values = gene_group_colors) +
  theme_bw() +
  scale_x_continuous(breaks = seq(floor(min(v_plot_data_luad$log2FC)), 
                                  ceiling(max(v_plot_data_luad$log2FC)), by = 2)) +
  #ggh4x::facet_nested(factor(method, levels = c("CN naive", "CN aware"))~factor(tumor_type, levels = c("BRCA", "LIHC", "HNSC", "COAD")), scales ="free", independent = "y")+
  facet_wrap(~factor(method, levels = c("CN naive", "CN aware")), nrow = 1) +
  labs(x = expression(Log[2] ~ FC), y = expression(-log[10] ~ Pvalue), col = "Gene group") +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = 'dashed') +
  geom_hline(yintercept = -log10(pval_cut), linetype = "dashed") +
  ggplot2::theme(legend.position = 'left',
                 legend.text = element_text(size = 16, color = "black"),
                 legend.title = element_text(size = 18, color = "black"),  
                 strip.text = element_text(size = 20, face = "plain", color = "black"),
                 axis.text = element_text(size = 18, color = "black"),
                 axis.title = element_text(size = 20, color = "black"))+
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1)))
p_volcanos

ggsave("CN-aware-DGE/plots/main/volcano_luad.png", dpi = 400, width = 10.0, height = 4.0, plot = p_volcanos)


# CN barplot
combined_data <- rbind(cn_aware_d_sensitive, cn_aware_d_insensitive, cn_aware_d_compensated, cn_aware_non_DE)

combined_data <- combined_data %>% 
  dplyr::mutate(cnv_group = case_when(
    cnv_mean > 0.5 & cnv_mean <= 1.7  ~ "loss",
    cnv_mean > 1.7 & cnv_mean <= 2.5  ~ "neutral",
    cnv_mean > 2.5 & cnv_mean <=   4.0 ~ "gain",
    cnv_mean > 4.0 ~ "amplification"))

combined_data <- merge(combined_data, loss_labels, by = "row.names")
combined_data$cnv_group <- ifelse(combined_data$isCNloss == "loss", "loss", combined_data$cnv_group)

barplot_data <- combined_data %>%
  group_by(gene_group) %>%
  summarise(Count = n()) %>%
  mutate(total = sum(Count)) %>%
  mutate(percentage = (Count / total) * 100) %>%
  ungroup()

combined_data$gene_group <- factor(combined_data$gene_group, levels = c("non-DEG", "Dosage-insensitive", "Dosage-sensitive", "Dosage-compensated"))

brca_barplot <- combined_data
lihc_barplot <- combined_data
hnsc_barplot <- combined_data
coad_barplot <- combined_data
coad_barplot <- na.omit(coad_barplot)

joint_data <- rbind(brca_barplot, lihc_barplot, hnsc_barplot, coad_barplot)

joint_data$tumor_type <- factor(joint_data$tumor_type, levels = c("BRCA", "LIHC", "HNSC", "COAD"))

barplot_cnv <- ggplot2::ggplot(joint_data, aes(x = gene_group, fill = cnv_group)) +
  geom_bar(position = "stack", width = 0.6) + 
  #geom_text(aes(label = Count), position = position_stack(vjust = 0.5), size = 4) +  
  scale_fill_manual(values = cnv_colors) +  
  theme_bw() +  
  facet_wrap(~tumor_type, scales = "free", nrow = 1)+
  labs(y = "gene counts", x = "", title = "", fill = "CN group") +  
  theme(
    axis.text.x = element_text(size = 14, angle = 30, hjust = 1, color = "black"),  
    axis.text.y = element_text(size = 14, color = "black"),                         
    axis.title.x = element_text(size = 14, face = "plain", color = "black"),          
    axis.title.y = element_text(size = 14, face = "plain", color = "black"),
    strip.text = element_text(size = 16, face = "plain", color = "black"),
    legend.position = 'bottom',
    legend.text = element_text(size = 14, color = "black"),                          
    legend.title = element_text(size = 16, face = "plain", color = "black")           
  )
barplot_cnv

ggsave("CN-aware-DGE/plots/main/barplot_cnv_luad_legend.png", dpi = 400, width = 6.5, height = 6.0, plot = barplot_cnv)
ggsave("CN-aware-DGE/plots/supplementary/barplot_cnv.png", dpi = 400, width = 10.0, height = 4.5, plot = barplot_cnv)


# Scatter - comparison LFC | p-value

lfc_naive <- res_naive %>% dplyr::select(logFC_naive) %>% dplyr::rename(logFC_naive=logFC_naive)
lfc_aware <- res_aware %>% dplyr::select(logFC_aware, cnv_mean_aware) %>% dplyr::rename(logFC_aware=logFC_aware)

plot_lfc <- merge(lfc_naive, lfc_aware, by = "row.names")
plot_lfc <- plot_lfc %>% remove_rownames() %>% column_to_rownames("Row.names")

plot_lfc <- plot_lfc %>% 
  dplyr::mutate(cnv_group = case_when(
    cnv_mean_aware > 0.5 & cnv_mean_aware <= 1.7  ~ "loss",
    cnv_mean_aware > 1.7 & cnv_mean_aware <= 2.5  ~ "neutral",
    cnv_mean_aware > 2.5 & cnv_mean_aware <=   4.0 ~ "gain",
    cnv_mean_aware > 4.0 ~ "amplification"))

plot_lfc <- merge(plot_lfc, loss_labels, by = "row.names")
plot_lfc$cnv_group <- ifelse(plot_lfc$isCNloss == "loss", "loss", plot_lfc$cnv_group)

plot_lfc <- plot_lfc %>% dplyr::filter(logFC_naive > -5.0 ,)
plot_lfc <- plot_lfc %>% dplyr::filter(logFC_aware > -5.0 ,)

plot_lfc_coad <- plot_lfc %>% 
  dplyr::mutate(eff_size_diff = abs(logFC_naive - logFC_aware)) %>% 
  dplyr::mutate(tumor_type = "COAD")

plot_lfc_coad <- na.omit(plot_lfc_coad)

plot_lfc <- rbind(plot_lfc_brca, plot_lfc_lihc, plot_lfc_hnsc, plot_lfc_coad)

comparison_lfc <- ggplot(plot_lfc, aes(x=logFC_aware, y=logFC_naive, color = cnv_group)) + 
  #geom_pointdensity(shape=20) +
  #geom_point(data=subset(plot_lfc, cnv_group != "amplification"), aes(color=cnv_group), shape=20)+
  #geom_point(data=subset(plot_lfc, cnv_group == "amplification"), aes(color=cnv_group), size=2, shape=20) + 
  geom_point(shape = 20, size = 3)+
  #geom_smooth(method=lm, se=TRUE, linetype="dashed",color="darkred")+
  geom_abline()+
  geom_vline(xintercept = 0, linetype="dotted", color="black", size=0.5)+
  geom_hline(yintercept = 0, linetype="dotted", color="black", size=0.5)+
  xlab("Effect size (log2) CN-aware") +
  ylab ("Effect size (log2) CN-naive") +
  scale_x_continuous(breaks = seq(-6, 6, by = 2))+
  scale_y_continuous(breaks = seq(-6, 6, by = 2))+
  #facet_wrap(~factor(cnv_group, levels = c("loss", "neutral", "gain", "amplification")), nrow = 1)+
  ggh4x::facet_nested(factor(tumor_type, levels = c("BRCA", "LIHC", "HNSC", "COAD"))~factor(cnv_group, levels = c("loss", "neutral", "gain", "amplification")))+
  scale_color_manual(name = "CNV Group", values = cnv_colors)+
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1)))+
  theme_bw()+
  theme(
    legend.position="bottom",
    axis.title.x = element_text(size=16, color = "black"),  
    axis.title.y = element_text(size=16, color = "black"),
    axis.text.x = element_text(size=16, color = "black"),    
    axis.text.y = element_text(size=16, color = "black"),
    legend.text = element_text(size = 16, color = "black"),
    strip.text = element_text(size = 16, face = "plain", color = "black"),
    legend.title = element_text(size = 18, face = "plain", color = "black") 
  )
comparison_lfc

ggsave("CN-aware-DGE/plots/main/scatter_lfc_luad.png", dpi = 400, width = 7.0, height = 3.5, plot = comparison_lfc)
ggsave("CN-aware-DGE/plots/supplementary/scatter_lfc.png", dpi = 400, width = 9.0, height = 9.0, plot = comparison_lfc)


# Effect size difference (log2)
plot_lfc <- plot_lfc %>% dplyr::filter(eff_size_diff < 3.0 ,)

plot_lfc$cnv_group <- factor(plot_lfc$cnv_group, levels = c("loss", "neutral", "gain", "amplification"))

violin <- ggplot(plot_lfc, aes(x = cnv_group, y = eff_size_diff, fill = cnv_group)) + 
  geom_violin(trim = FALSE, scale = "width", alpha = 0.7) + 
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.5) + 
  scale_fill_manual(values = cnv_colors) + 
  xlab("CNV group") + 
  ylab("Effect size difference (log2)") + 
  theme_bw() +
  facet_wrap(~factor(tumor_type, levels = c("BRCA", "LIHC", "HNSC", "COAD")), nrow = 4)+
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 18, color = "black"),
    axis.title.y = element_text(size = 18, color = "black"),
    axis.text.x = element_text(size = 18, angle = 30, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 18, color = "black"),
    strip.text = element_text(size = 16, face = "plain", color = "black")
  )
violin

ggsave("CN-aware-DGE/plots/main/violin_luad.png", dpi = 400, width = 5.5, height = 4.5, plot = violin)
ggsave("CN-aware-DGE/plots/supplementary/violin.png", dpi = 400, width = 3.5, height = 12.0, plot = violin)


# p_value

pval_naive <- res_naive %>% dplyr::select(padj_naive) %>% dplyr::rename(padj_naive=padj_naive)
pval_aware <- res_aware %>% dplyr::select(padj_aware, cnv_mean_aware) %>% dplyr::rename(padj_aware = padj_aware)
plot_pval <- merge(pval_naive, pval_aware, by = "row.names")
plot_pval <- plot_pval %>% remove_rownames() %>% column_to_rownames("Row.names")

plot_pval <- plot_pval %>% 
  dplyr::mutate(cnv_group = case_when(
    cnv_mean_aware > 0.5 & cnv_mean_aware <= 1.7  ~ "loss",
    cnv_mean_aware > 1.7 & cnv_mean_aware <= 2.5  ~ "neutral",
    cnv_mean_aware > 2.5 & cnv_mean_aware <=   4.0 ~ "gain",
    cnv_mean_aware > 4.0 ~ "amplification"))

plot_pval<- merge(plot_pval, loss_labels, by = "row.names")
plot_pval$cnv_group <- ifelse(plot_pval$isCNloss == "loss", "loss", plot_pval$cnv_group)

brca_pval <- plot_pval %>% dplyr::mutate(tumor_type = "BRCA")
brca_pval <- brca_pval %>% dplyr::filter(padj_naive > 2.491311e-147 ,)

lihc_pval <- plot_pval %>% dplyr::mutate(tumor_type = "LIHC")
lihc_pval <- lihc_pval %>% dplyr::filter(padj_naive > 3.824830e-63 ,)

hnsc_pval <- plot_pval %>% dplyr::mutate(tumor_type = "HNSC")
hnsc_pval <- hnsc_pval %>% dplyr::filter(padj_aware >  5.140434e-42,)

coad_pval <- plot_pval %>% dplyr::mutate(tumor_type = "COAD") %>% na.omit()
coad_pval <- coad_pval %>% dplyr::filter(padj_aware > 3.173776e-71 ,)

joint_pvalue <- rbind(brca_pval, lihc_pval, hnsc_pval, coad_pval)

comparison_pval <- ggplot(hnsc_pval, aes(x=-log10(padj_naive), y=-log10(padj_aware), color = cnv_group)) + 
  geom_point(shape=20, size=3) +
  #geom_point(data=subset(plot_pval, cnv_group != "amplification"), aes(color=cnv_group), shape=20)+
  #geom_point(data=subset(plot_pval, cnv_group == "amplification"), aes(color=cnv_group), size=2, shape=20) + 
  #geom_smooth(method=lm, se=FALSE, linetype="dashed",color="darkred")+
  geom_abline()+
  #geom_vline(xintercept = 0, linetype="dotted", color="black", size=1)+
  #geom_hline(yintercept = 0, linetype="dotted", color="black", size=1)+
  #geom_abline(intercept = 0, slope = 1, linetype="dotted", color="black", size=1) +
  xlab("FDR CN-aware") +
  ylab ("FDR CN-naive") +
  scale_x_continuous(breaks = seq(0, 100, by = 20))+
  scale_y_continuous(breaks = seq(0, 100, by = 20))+
  scale_color_manual(name = "CNV group", values = cnv_colors)+
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1)))+
  #facet_wrap(~factor(cnv_group, levels = c("loss", "neutral", "gain", "amplification")), nrow = 1)+
  theme_bw()+
  ggh4x::facet_nested(factor(tumor_type, levels = c("BRCA", "LIHC", "HNSC", "COAD"))~
                        factor(cnv_group, levels = c("loss", "neutral", "gain", "amplification")))+
  theme(
    legend.position="bottom",
    axis.title.x = element_text(size=16, color = "black"),  
    axis.title.y = element_text(size=16, color = "black"),
    axis.text.x = element_text(size=16, color = "black"),    
    axis.text.y = element_text(size=16, color = "black"),
    strip.text = element_text(size = 16, face = "plain", color = "black"),
    legend.text = element_text(size = 16, color = "black"),                          
    legend.title = element_text(size = 18, face = "plain", color = "black") 
  )
comparison_pval

ggsave("CN-aware-DGE/plots/main/scatter_pvalue_luad.png", dpi = 400, width = 7.0, height = 3.5, plot = comparison_pval)
ggsave("CN-aware-DGE/plots/supplementary/pval_hnsc.png", dpi = 400, width = 8.0, height = 3.0, plot = comparison_pval)
  
#spacer <- grid::textGrob("")
#combined_plot <- gridExtra::grid.arrange(comparison_lfc, spacer, comparison_pval, nrow = 3, heights = c(1.5, 0.2, 1.5))
#ggsave("CN-aware-DGE/plots/main/scatter_lfc_pval.png", plot = grid::grobTree(combined_plot), dpi = 400, width = 6.0, height = 12, units = "in")



# Sankey dynamic gene groups transitions 

res_naive <- res_pydeseq %>% dplyr::select(DEtype) %>% dplyr::rename(CN_naive = DEtype)
res_aware <- res_deconveil %>% dplyr::select(DEtype) %>% dplyr::rename(CN_aware = DEtype)

res_join <- cbind(res_naive, res_aware)

data_flow <- res_join %>%
  group_by(CN_naive,CN_aware) %>%
  summarise(freq = n()) %>%
  ungroup()

data_ggforce <- data_flow  %>%
  gather_set_data(1:2) %>%        
  arrange(x,CN_naive,desc(CN_aware))

data_ggforce$CN_naive <- factor(data_ggforce$CN_naive)
data_ggforce$CN_aware <- factor(data_ggforce$CN_aware)

data_ggforce <- data_ggforce %>%
  group_by(CN_naive, CN_aware) %>%
  mutate(y_mid = freq / 2) %>% 
  na.omit()

g_group_colors <- c("Down-reg" = "#3C5488B2", "n.s." = "lightgray", "Up-reg" = "#CC7677")

sankey <- ggplot(data_ggforce, aes(x = x, id = id, split = y, value = freq)) +
  geom_parallel_sets(aes(fill = CN_naive), alpha = 0.9, axis.width = 0.2,
                     n = 4415, strength = 0.5, color = "black", linewidth = 0.3) +
  geom_parallel_sets_axes(axis.width = 0.25, fill = "gray93",
                          color = "gray", linewidth = 0.5) +  
  #geom_parallel_sets_labels(colour = 'gray35', size = 2.0, angle = 0, fontface = "plain") +
  scale_fill_manual(values = g_group_colors, name = "Gene group") +
  scale_color_manual(values = g_group_colors) +
  scale_x_continuous(breaks = 1:2, labels = c("CN-naive", "CN-aware"))+
  theme_classic() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 15, face = "plain"),
    legend.text = element_text(size = 13),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x  = element_blank()
  )
sankey

ggsave("CN-aware-DGE/plots/supplementary/sankey_coad.png", dpi = 400, width = 4.0, height = 5.0, plot = sankey)

