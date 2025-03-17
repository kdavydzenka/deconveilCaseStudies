setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/TCGA/LUSC/")

pkgs <- c("dplyr", "tidyr", "ggplot2")
sapply(pkgs, require, character.only = TRUE)

# Load data

beta_coeff <- readRDS("beta_coeff.RDS")
gene_groups <- readRDS("gene_groups.RDS")

beta_dsg <- beta_coeff[beta_coeff$Gene %in% gene_groups[["dsg"]], ]
beta_dig <- beta_coeff[beta_coeff$Gene %in% gene_groups[["dig"]], ]
beta_dcg <- beta_coeff[beta_coeff$Gene %in% gene_groups[["dcg"]], ]

#Create binary indicators for non-zero beta coefficients
g_CN <- ifelse(beta_dcg$CN != 0 & beta_dcg$Meth == 0 & beta_dcg$pTF == 0 & beta_dcg$miRNA == 0, 1, 0)
g_Meth <- ifelse(beta_dcg$Meth != 0 & beta_dcg$CN == 0 & beta_dcg$pTF == 0 & beta_dcg$miRNA == 0, 1, 0)
g_pTF <- ifelse(beta_dcg$pTF != 0 & beta_dcg$CN == 0 & beta_dcg$Meth == 0 & beta_dcg$miRNA == 0, 1, 0)
g_miRNA <- ifelse(beta_dcg$miRNA != 0 & beta_dcg$CN == 0 & beta_dcg$Meth == 0 & beta_dcg$pTF == 0, 1, 0)

g_all_pred <- ifelse(beta_dcg$CN != 0 & beta_dcg$Meth != 0 & beta_dcg$pTF != 0 & beta_dcg$miRNA != 0, 1, 0)

g_CN_M_TF <- ifelse(beta_dcg$CN != 0 & beta_dcg$Meth != 0 & beta_dcg$pTF != 0 & beta_dcg$miRNA == 0, 1, 0)
g_CN_M_mir <- ifelse(beta_dcg$CN != 0 & beta_dcg$Meth != 0 & beta_dcg$pTF == 0 & beta_dcg$miRNA != 0, 1, 0)
g_CN_TF_mir <- ifelse(beta_dcg$CN != 0 & beta_dcg$Meth == 0 & beta_dcg$pTF != 0 & beta_dcg$miRNA != 0, 1, 0)
g_M_TF_mir <- ifelse(beta_dcg$CN == 0 & beta_dcg$Meth != 0 & beta_dcg$pTF != 0 & beta_dcg$miRNA != 0, 1, 0)

g_CN_M <- ifelse(beta_dcg$CN != 0 & beta_dcg$Meth != 0 & beta_dcg$pTF == 0 & beta_dcg$miRNA == 0, 1, 0)
g_CN_TF <- ifelse(beta_dcg$CN != 0 & beta_dcg$Meth == 0 & beta_dcg$pTF != 0 & beta_dcg$miRNA == 0, 1, 0)
g_CN_mir <- ifelse(beta_dcg$CN != 0 & beta_dcg$Meth == 0 & beta_dcg$pTF == 0 & beta_dcg$miRNA != 0, 1, 0)
g_M_TF <- ifelse(beta_dcg$CN == 0 & beta_dcg$Meth != 0 & beta_dcg$pTF != 0 & beta_dcg$miRNA == 0, 1, 0)
g_M_mir <- ifelse(beta_dcg$CN == 0 & beta_dcg$Meth != 0 & beta_dcg$pTF == 0 & beta_dcg$miRNA == 0, 1, 0)
g_TF_mir <- ifelse(beta_dcg$CN == 0 & beta_dcg$Meth == 0 & beta_dcg$pTF != 0 & beta_dcg$miRNA != 0, 1, 0)

none <- ifelse(beta_dcg$CN == 0 & beta_dcg$Meth == 0 & beta_dcg$pTF == 0 & beta_dcg$miRNA == 0, 1, 0)


# Calculate the number of genes affected by individual/combinations of predictors

g_CN <- mean(g_CN)
g_Meth <- mean(g_Meth)
g_pTF <- mean(g_pTF)
g_miRNA <- mean(g_miRNA)

g_all_pred <- mean(g_all_pred)

g_CN_M_TF <- mean(g_CN_M_TF)
g_CN_M_mir <- mean(g_CN_M_mir)
g_CN_TF_mir <- mean(g_CN_TF_mir)
g_M_TF_mir <- mean(g_M_TF_mir)

g_CN_M <- mean(g_CN_M)
g_CN_TF <- mean(g_CN_TF)
g_M_TF <- mean(g_M_TF)
g_TF_mir <- mean(g_TF_mir)
g_M_mir <- mean(g_M_mir)
g_CN_mir <- mean(g_CN_mir)
none <- mean(none)

combinations <- c("CN", "M", "TF", "miRNA", 
                  "CN+M+TF+miRNA", "CN+M+TF", "CN+M+miRNA", "CN+TF+miRNA", "M+TF+miRNA", 
                  "CN+M", "CN+TF", "M+TF", 
                  "TF+miRNA", "M+miRNA", "CN+miRNA", "None")

proportions <- c(g_CN, g_Meth, g_pTF, g_miRNA, g_all_pred,
                 g_CN_M_TF, g_CN_M_mir, g_CN_TF_mir, g_M_TF_mir,
                 g_CN_M, g_CN_TF, g_M_TF, g_TF_mir, g_M_mir, g_CN_mir, none)

# DSGs
df_res_dsg <- data.frame(
  Combination = combinations,
  Proportion = proportions
)

df_res_dsg$gene_group <- "DSGs"
sum(df_res_dsg$Proportion)

# DIGs
df_res_dig <- data.frame(
  Combination = combinations,
  Proportion = proportions
)

df_res_dig$gene_group <- "DIGs"
sum(df_res_dsg$Proportion)

# DCGs
df_res_dcg <- data.frame(
  Combination = combinations,
  Proportion = proportions
)

df_res_dcg$gene_group <- "DCGs"
sum(df_res_dcg$Proportion)

combined_data <- rbind(df_res_dsg, df_res_dig, df_res_dcg)
combined_data <- combined_data %>% filter(Proportion > 0)
combined_data$gene_group <- factor(combined_data$gene_group, levels = c("DSGs", "DIGs", "DCGs"))
combined_data <- combined_data %>% mutate(tumor_type = "LUSC")

barplot <- ggplot2::ggplot(combined_data, aes(x = gene_group, y = Proportion, fill = Combination)) +
  geom_bar(stat = "identity", position = "stack", width = 0.6) + 
  scale_fill_brewer(palette = "Set3") +  
  theme_bw() +  
  facet_wrap(~tumor_type, scales = "free", nrow = 1)+
  labs(y = "Proportion of genes", x = "", title = "", fill = "") +  
  theme(
    axis.text.x = element_text(size = 18, color = "black"),  
    axis.text.y = element_text(size = 16, color = "black"),                         
    axis.title.x = element_text(size = 18, face = "plain", color = "black"),          
    axis.title.y = element_text(size = 18, face = "plain", color = "black"),
    strip.text = element_text(size = 16, face = "plain", color = "black"),
    legend.position = 'right',
    legend.text = element_text(size = 14, color = "black"),                          
    legend.title = element_text(size = 16, face = "plain", color = "black")           
  )
barplot
