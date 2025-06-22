### Explortive analysis of RNA (z-score) and  CNV  ###

setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/")
pkgs <- c("dplyr", "ggplot2", "cluster", "factoextra", "heatmaply", "DESeq2", "tidyverse", "DESeq2", "colorspace", 
          "ggpubr", "ggpointdensity", "ggeasy", "patchwork")
invisible(sapply(pkgs, require, character.only = TRUE))

source("deconveilCaseStudies/utils/utils.R")


## Input data CN | Datasets analysed LUAD | BRCA | LIHC | LUSC ##

# LUAD
data_path <- "TCGA/lung/LUAD/cnv_tumor.RDS"
dataset_name <- "LUAD_cnv"
cnv_tumor <- readRDS(data_path)
cnv_tumor_luad <- cnv_tumor %>% as.data.frame() %>% dplyr::select(7,5,13,14,15,18,22,24,25,28) 
cnv_filt_luad <- apply(cnv_tumor_luad, 2, function(x) ifelse(x > 15, 15, x)) 

# BRCA 
data_path <- "TCGA/BRCA/cnv_tumor.RDS"
dataset_name <- "BRCA_cnv"
cnv_tumor <- readRDS(data_path)

# LIHC
data_path <- "TCGA/LIHC/cnv_tumor.RDS"
dataset_name <- "LIHC_cnv"
cnv_tumor <- readRDS(data_path)

# LUSC
data_path <- "TCGA/LUSC/cnv_tumor.RDS"
dataset_name <- "LUSC_cnv"
cnv_tumor <- readRDS(data_path)


# Clustering patients 
set.seed(42)
clustering_res <- clustering_patients_cnv(dataset_name, data_path)
clustering_res[[2]]
cnv_tumor <- clustering_res[[1]]

cancer_clusters <- list(
  BRCA = c(3,4),
  LUSC = c(1,3),
  LIHC = c(1, 2)
)

cnv_filt <- cnv_tumor %>% 
  as.data.frame() %>% 
  dplyr::filter(Cluster %in% cancer_clusters$LUSC) # LUSC | LIHC

#Select patients for boxplot
cnv_filt <- as.matrix(t(cnv_filt))
cnv_filt<- apply(cnv_filt, 2, function(x) ifelse(x > 15, 15, x)) 

hist(rowMeans(cnv_filt),
     main = "LUAD", 
     xlab = "CN state",
     ylab = "Proportion",
     col = "#E1DEFC",
     prob = TRUE,
     breaks = 15)

cnv_filt <- cnv_filt_luad
cnv_mean <- cnv_filt %>% 
  as.data.frame() %>% 
  dplyr::mutate(cnv_mean = rowMeans(cnv_filt)) %>% 
  dplyr::select(cnv_mean)

cnv_mean$geneID <- rownames(cnv_mean)


## Input data RNA ##

# LUAD
data_path <- "TCGA/lung/LUAD/rna_counts.RDS"
dataset_name <- "LUAD_rna"

# BRCA
data_path <- "TCGA/BRCA/rna_counts.RDS"
dataset_name <- "BRCA_rna"

# LUSC
data_path <- "TCGA/LUSC/rna_counts.RDS"
dataset_name <- "LUSC_rna"

# LIHC
data_path <- "TCGA/LIHC/rna_counts.RDS"
dataset_name <- "LIHC_rna"

rna <- rna_processing(dataset_name, data_path, cnv_filt)
rna_norm <- rna[[1]]
rna_tum <- rna[[2]]

# Exclude genes with low expression in normal tissue 
low_expression_threshold <- 20
expression_summary <- data.frame(
  Gene = rownames(rna_norm),
  MeanExpression = rowMeans(rna_norm)
)

filtered_genes <- expression_summary %>%
  filter(MeanExpression > low_expression_threshold)

rna_norm <- rna_norm[filtered_genes$Gene, ]
rna_tum <- rna_tum[filtered_genes$Gene, ]
rna <- cbind(rna_norm, rna_tum)

# Normalization 
dds <- DESeqDataSetFromMatrix(countData = rna, colData = DataFrame(condition = rep("A", ncol(rna))), design = ~1)
rna_vst <- vst(dds, blind = TRUE)
rna_log_normalized <- assay(rna_vst)


# Z-score trasformation
rna_zscore <- t(scale(t(rna_log_normalized)))
rna_zscore_tumor <- rna_zscore %>% as.data.frame() %>% dplyr::select(11:20) %>% as.matrix()

rna_zscore_tumor <- rna_zscore_tumor %>% 
  as.data.frame() %>% 
  dplyr::mutate(zscore_mean = rowMeans(rna_zscore_tumor)) %>% 
  dplyr::select(zscore_mean)


# CNV factorization 

cnv <- cnv_mean %>% 
  dplyr::mutate(cnv = case_when(
    cnv_mean > 0.0 & cnv_mean <= 1.8 ~ "1",
    cnv_mean > 1.8 & cnv_mean <= 2.5 ~ "2",
    cnv_mean > 2.5 & cnv_mean <= 3.5 ~ "3",
    cnv_mean > 3.5 & cnv_mean <= 4.5 ~ "4",
    cnv_mean > 4.5 ~ "5")) %>% 
  dplyr::select(cnv)

#cnv <- cnv[rownames(cnv) %in% rownames(rna_zscore_tumor),]


### Preparing data for boxplot ###

rna_zscore_tumor <- rna_zscore_tumor %>% dplyr::mutate(sample_type = "Tumor")

p_luad_t <- rna_zscore_tumor %>%
  rownames_to_column(var = "gene") %>%
  merge(cnv %>% rownames_to_column(var = "gene"), by = "gene") %>% 
  na.omit()

#p_brca_t <- p_brca_t %>% 
  #dplyr::filter(cnv=="2" & abs(zscore_mean) <= 0.5 | cnv=="1" | cnv=="3" | cnv=="4" | cnv=="5")

# Cancer types Supplementary
p_luad_t <- p_luad_t %>% dplyr::mutate(cancer_type = "LUAD")
p_brca_t <- p_brca_t %>% dplyr::mutate(cancer_type = "BRCA")
p_lusc_t <- p_lusc_t %>% dplyr::mutate(cancer_type = "LUSC")
p_lihc_t <- p_lihc_t %>% dplyr::mutate(cancer_type = "LIHC")

p_tumor <- rbind(p_brca_t, p_lusc_t, p_lihc_t)

saveRDS(p_luad_t, file = "deconveilCaseStudies/plots/main/Fig 1/rds/plot_boxplot.rds")


# Boxplot #

#col <- qualitative_hcl(5, palette = "Warm")
p_lusc_t$cnv <- factor(p_lusc_t$cnv, levels = c(1, 2, 3, 4, 5))
col <- c("1" = "dodgerblue1", "2" = "darkgray", "3" = "green4", "4" = "coral3", "5" = "hotpink3")


bxp_t <- ggplot(p_lusc_t, aes(x = cnv, y = zscore_mean, color = cnv)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch = F, show.legend = F)+
  #geom_smooth(method = "loess", formula = y ~ x, se=FALSE, color="darkblue", aes(group=1), linetype = 'dashed')+
  labs(x="CN group", y = "mRNA Z-score")+
  facet_wrap(~cancer_type)+
  theme(strip.text.x = element_text(size=14, color="black", face="bold.italic"))+
  ggplot2::theme(legend.position = 'none')+
  theme_bw()+
  scale_color_manual(values=col)+
  font("xy.text", size = 12, color = "black", face = "plain")+
  font("title", size = 12, color = "black")+
  font("xlab", size = 12)+
  font("ylab", size = 12)+
  ggtitle("")+
  theme(plot.title = element_text(hjust = 0.5),
        strip.text = element_text(size = 12, face = "plain", color = "black"))+
  theme(legend.position = "")
bxp_t

ggsave("deconveilCaseStudies/plots/supplementary/boxplots.png", dpi = 500, width = 7.5, height = 4.0, plot = bxp_t)






