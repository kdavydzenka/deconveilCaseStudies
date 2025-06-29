### Survival prognostic model using TCGA-BRCA dataset ###

setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/")
pkgs <- c("tidyverse", "survival", "glmnet", "survminer", "survcomp", "DESeq2", "forestplot", 
          "caret", "randomForestSRC", "gridExtra", "RColorBrewer", "patchwork")
sapply(pkgs, require, character.only = TRUE)
source("deconveilCaseStudies/utils/utils_plot.R")
source("deconveilCaseStudies/utils/utils_survival.R")

# Load data #
gene_groups_CNaware <- readRDS("deconveilCaseStudies/plots/main/Fig 5/rds/gene_groups.RDS")
res_CNnaive <- read.csv("deconveilCaseStudies/results_tcga/BRCA/res_CNnaive.csv")
clinical_data <- readRDS("TCGA/BRCA/clinical_full.RDS")
rna_tumor <- readRDS("TCGA/BRCA/rna_tumor.RDS")
cnv_tumor <- readRDS("TCGA/BRCA/cnv_tumor.RDS")

# Preprocessing #
colnames(rna_tumor) <- substr(colnames(rna_tumor), 1, 12)
colnames(cnv_tumor) <- substr(colnames(cnv_tumor), 1, 12)
clinical_data <- clinical_data[clinical_data$bcr_patient_barcode %in% colnames(rna_tumor),]
rna_tumor <- rna_tumor[,(colnames(rna_tumor) %in% clinical_data$bcr_patient_barcode)]

clinical_data <- clinical_data %>%
  dplyr::select(bcr_patient_barcode, days_to_death, days_to_last_followup) %>%
  dplyr::mutate(
    event = ifelse(!is.na(days_to_death), 1, 0),
    time = ifelse(!is.na(days_to_death), days_to_death, days_to_last_followup)
  ) %>%
  dplyr::select(bcr_patient_barcode, time, event) %>%
  drop_na()

colnames(clinical_data) <- c("patientID", "time", "event")

# Gene extraction #

get_rna_matrix <- function(gene_ids, rna_matrix) {
  rna_matrix[rownames(rna_matrix) %in% gene_ids, ]
}

rna_ds     <- get_rna_matrix(gene_groups_CNaware[["d_sensitive"]]$geneID, rna_tumor)
rna_dins   <- get_rna_matrix(gene_groups_CNaware[["d_insensitive"]]$geneID, rna_tumor)
rna_dcomp  <- get_rna_matrix(gene_groups_CNaware[["d_compensated"]]$geneID, rna_tumor)

# CN-naive DE gene extraction #

res_naive <- annotate_results(res_CNnaive, lfc_cut = 1.0, pval_cut = 0.05, method = "CN naive", tumor_type = "BRCA")
res_naive_deg <- res_naive %>% filter(isDE == TRUE)
res_naive_deg$geneID <- rownames(res_naive_deg)
rna_CNnaive <- get_rna_matrix(res_naive_deg$geneID, rna_tumor)

# Normalization and data merging #

data_ds <- prepare_survival_data(rna_ds, clinical_data, cnv_tumor)
data_dins <- prepare_survival_data(rna_dins, clinical_data, cnv_tumor, apply_cnv = TRUE)
data_dcomp <- prepare_survival_data(rna_dcomp, clinical_data, cnv_tumor, apply_cnv = TRUE)  # Only this uses CN correction
data_naive <- prepare_survival_data(rna_CNnaive, clinical_data)

## Initial Gene Selection Using Cox Proportional Hazards Model ##
sig_genes_ds <- run_cox_analysis(data_ds$rna, data_ds$clinical)
sig_genes_dins <- run_cox_analysis(data_dins$rna, data_dins$clinical)
sig_genes_dcomp <- run_cox_analysis(data_dcomp$rna, data_dcomp$clinical)
sig_genes_naive <- run_cox_analysis(data_naive$rna, data_naive$clinical)

saveRDS(sig_genes_ds, file = "deconveilCaseStudies/plots/main/Fig 5/rds/cox_significant_dsg.RDS")
saveRDS(sig_genes_dins, file = "deconveilCaseStudies/plots/main/Fig 5/rds/cox_significant_dig.RDS")
saveRDS(sig_genes_dcomp, file = "deconveilCaseStudies/plots/main/Fig 5/rds/cox_significant_dcg.RDS")
saveRDS(sig_genes_naive, file = "deconveilCaseStudies/plots/main/Fig 5/rds/cox_significant_CNnaive.RDS")


# The Hazard Ratio (HR) and confidence intervals (CI) give an indication of the prognostic impact of each gene. 
# Genes with HR > 1 may be associated with poor prognosis, while HR < 1 suggests good prognosis.
# Feature Selection Using LASSO (Penalized Cox Model) #
# Build Gene Prognostic Signature Using Lasso Regression (selects the most important genes and avoids overfitting) 

## Fit LASSO and calculate prognostic score ##

sig_genes_ds <- readRDS("deconveilCaseStudies/results/brca_prognostic/cox_significant_dsg.RDS")
sig_genes_dins <- readRDS("deconveilCaseStudies/results/brca_prognostic/cox_significant_dig.RDS")
sig_genes_dcomp <- readRDS("deconveilCaseStudies/results/brca_prognostic/cox_significant_dcg.RDS")
sig_genes_naive <- readRDS("deconveilCaseStudies/results/brca_prognostic/cox_significant_CNnaive.RDS")

lasso_dsg <- fit_lasso_and_score(data_ds$rna, data_ds$clinical, sig_genes_ds)
lasso_dig <- fit_lasso_and_score(data_dins$rna, data_dins$clinical, sig_genes_dins)
lasso_dcg <- fit_lasso_and_score(data_dcomp$rna, data_dcomp$clinical, sig_genes_dcomp)
lasso_naive <- fit_lasso_and_score(data_naive$rna, data_naive$clinical, sig_genes_naive)

saveRDS(lasso_dsg, file = "deconveilCaseStudies/plots/main/Fig 5/rds/lasso_dsg.RDS")
saveRDS(lasso_dig, file = "deconveilCaseStudies/plots/main/Fig 5/rds/lasso_dig.RDS")
saveRDS(lasso_dcg, file = "deconveilCaseStudies/plots/main/Fig 5/rds/lasso_dcg.RDS")
saveRDS(lasso_naive, file = "deconveilCaseStudies/plots/main/Fig 5/rds/lasso_CNnaive.RDS")

## Forest plot ##

lasso_naive <- readRDS("deconveilCaseStudies/plots/main/Fig 5/rds/lasso_CNnaive.RDS")
lasso_dsg <- readRDS("deconveilCaseStudies/plots/main/Fig 5/rds/lasso_dsg.RDS")
lasso_dig <- readRDS("deconveilCaseStudies/plots/main/Fig 5/rds/lasso_dig.RDS")
lasso_dcg <- readRDS("deconveilCaseStudies/plots/main/Fig 5/rds/lasso_dcg.RDS")

f1 <- forest_plot(lasso_dsg, title = "DSGs")
f2 <- forest_plot(lasso_dig, title = "DIGs")
f3 <- forest_plot(lasso_dcg, title = "DCGs")
f4 <- forest_plot(lasso_naive, title = "DEGs (n = 5063)")

joint_forest <- f1 | f2 | f3 | f4
joint_forest

ggsave("deconveilCaseStudies/plots/main/Fig 5/joint_forest.png", dpi = 500, width = 14.0, height = 5.0, plot = joint_forest)



### Perform prognostic model validation using external cohort METABRIC ###

# Load data
rna_metabric <- read_delim("TCGA/brca_metabric/data_mrna_illumina_microarray.txt")
cn_metabric <- read_delim("TCGA/brca_metabric/data_cna.txt")
clinical_metabric <- read_delim("TCGA/brca_metabric/data_clinical_patient.txt")

# Preprocessing

process_expression_data <- function(df) {
  df %>%
    dplyr::filter(!duplicated(Hugo_Symbol)) %>%
    dplyr::select(-Entrez_Gene_Id) %>%
    column_to_rownames("Hugo_Symbol")
}

rna_metabric <- process_expression_data(rna_metabric)
cn_metabric  <- process_expression_data(cn_metabric)

clinical_metabric <- clinical_metabric[-c(1:4), ] %>%
  filter(Cohort == 1) %>%
  select(patientID = `#Patient Identifier`,
         time = `Overall Survival (Months)`,
         event = `Overall Survival Status`) %>%
  mutate(
    time = as.numeric(time) %>% round(),
    event = case_when(
      event == "1:DECEASED" ~ 1,
      event == "0:LIVING" ~ 0,
      TRUE ~ NA_real_
    )
  ) %>%
  drop_na()


## Match RNA/CNV/Clinical ##

common_patients <- intersect(colnames(rna_metabric), clinical_metabric$patientID)
rna_metabric <- rna_metabric[, common_patients]
cn_metabric  <- cn_metabric[, common_patients]
clinical_metabric <- clinical_metabric %>%
  filter(patientID %in% common_patients)

# Match genes across RNA and CNV
common_genes <- intersect(rownames(rna_metabric), rownames(cn_metabric))
rna_metabric <- rna_metabric[common_genes, ]
cn_metabric  <- cn_metabric[common_genes, ]

rna_metabric <- rna_metabric[, clinical_metabric$patientID]
cn_metabric  <- cn_metabric[, clinical_metabric$patientID]

cn_metabric <- cn_metabric %>%
  mutate(across(everything(), ~ recode(.,
                                       `-2` = 0, `-1` = 1, `0` = 2, `1` = 3, `2` = 4
  )))

cn_metabric[cn_metabric == 0] <- 0.001
cn_metabric <- cn_metabric / 2

saveRDS(rna_metabric, file = "deconveilCaseStudies/TCGA/brca_metabric/rna_metabric.rds")
saveRDS(cn_metabric, file = "deconveilCaseStudies/TCGA/brca_metabric/cn_metabric.rds")
saveRDS(clinical_metabric, file = "deconveilCaseStudies/TCGA/brca_metabric/clinical_metabric.rds")


# Subset METABRIC data for the same genes used in TCGA-BRCA
rna_metabric <- readRDS("deconveilCaseStudies/TCGA/brca_metabric/rna_metabric.rds")
cn_metabric <- readRDS("deconveilCaseStudies/TCGA/brca_metabric/cn_metabric.rds")
clinical_metabric <- readRDS("deconveilCaseStudies/TCGA/brca_metabric/clinical_metabric.rds")

rna_metabric_ds <- rna_metabric[(rownames(rna_metabric) %in% rownames(rna_ds)),]
cn_metabric_ds <- cn_metabric[(rownames(cn_metabric) %in% rownames(rna_ds)),]

rna_metabric_ds <- rna_metabric_ds[,(colnames(rna_metabric_ds) %in% colnames(cn_metabric_ds))]
cn_metabric_ds <- cn_metabric_ds[,(colnames(cn_metabric_ds) %in% colnames(rna_metabric_ds))]
idx <- match(colnames(rna_metabric_ds), colnames(cn_metabric_ds))
cn_metabric_ds <- cn_metabric_ds[,idx]

#rna_metabric_dins <- rna_metabric_dins * cn_metabric_dins
rna_metabric_ds <- rna_metabric_ds %>% 
  as.data.frame() %>% 
  dplyr::mutate(patientID = rownames(rna_metabric_ds))

metabric_combined <- merge(clinical_metabric, rna_metabric_ds, by = "patientID")


# CN-naive
rna_metabric_naive <- rna_metabric[(rownames(rna_metabric) %in% rownames(rna_CNnaive)),]
rna_metabric_naive <- t(rna_metabric_naive)
rna_metabric_naive <- rna_metabric_naive %>% 
  as.data.frame() %>% 
  dplyr::mutate(patientID = rownames(rna_metabric_naive))

metabric_combined <- merge(clinical_metabric, rna_metabric_naive, by = "patientID")


# Ensure Matching Gene Names

genes_in_lasso <- lasso_dsg[["selected_genes"]]
genes_in_metabric <- colnames(metabric_combined)
matching_genes <- intersect(genes_in_lasso, genes_in_metabric)
matching_genes <- intersect(colnames(rna_metabric_ds), lasso_dsg[["selected_genes"]])
metabric_filtered <- metabric_combined[, matching_genes]
lasso_coefficients <- lasso_dsg[["lasso_coefficients"]][matching_genes,]

# Calculate the risk score for METABRIC data using the LASSO model coefficients
metabric_risk_scores <- as.matrix(metabric_filtered) %*% lasso_coefficients
metabric_combined$risk_score <- as.vector(metabric_risk_scores)
metabric_combined <- metabric_combined %>% na.omit()

# Stratify METABRIC patients into high/low risk groups based on median risk score
metabric_combined$risk_group <- ifelse(metabric_combined$risk_score > median(metabric_combined$risk_score), "High", "Low")
metabric_combined$event <- as.numeric(as.character(metabric_combined$event))
metabric_survfit <- survfit(Surv(time, event) ~ risk_group, data = metabric_combined)

# Plot
pval <- surv_pvalue(metabric_survfit, data = metabric_combined)$pval
pval_label <- paste0("p = ", signif(pval, 3))

surv_plot <- ggsurvplot(
  metabric_survfit,
  data = metabric_combined,
  pval = FALSE,  
  conf.int = TRUE,
  risk.table = TRUE,
  palette = c("#F39B7FFF", "#3C5488FF"),
  title = "DSGs",
  xlab = "Time (days)",
  ylab = "Survival probability",
  font.main = c(12, "bold", "black"),
  font.x = c(12, "plain"),
  font.y = c(12, "plain"),
  font.tickslab = c(12, "plain"),
  legend = "bottom",
  legend.title = "Risk group",
  legend.labs = c("High risk", "Low risk"),
  font.legend = c(12, "plain"),
  risk.table.height = 0.25,
  risk.table.y.text = TRUE,
  risk.table.fontsize = 6,
  risk.table.title = "Number at risk",
  risk.table.col = "strata",
  ggtheme = theme_classic2() +
    theme(strip.text = element_text(size = 12, face = "plain"),
          plot.title = element_text(hjust = 0.5)),
  risk.table.title.fontface = "bold"
)

x_pos <- 0.05 * max(metabric_combined$time, na.rm = TRUE)
y_pos <- 0.2

surv_plot$plot <- surv_plot$plot +
  annotate("text", x = x_pos, y = y_pos, label = pval_label, size = 4, hjust = 0) +
  theme(legend.text = element_text(size = 12))

surv_plot$table <- surv_plot$table + 
  theme(
    axis.text.x = element_text(size = 12),
    strip.text = element_text(size = 12, face = "plain"), 
    text = element_text(size = 12)  
  )

surv_plot

#surv_plot_dcg <- surv_plot
#surv_plot_dsg <- surv_plot
#surv_plot_naive <- surv_plot
#surv_plot_dig <- surv_plot

surv_plot_dsg$plot <- surv_plot_dsg$plot + theme(legend.position = "none")
surv_plot_dig$plot <- surv_plot_dig$plot + theme(legend.position = "none")
surv_plot_naive$plot <- surv_plot_naive$plot + theme(legend.position = "none")

plot1 <- surv_plot_dsg$plot
plot2 <- surv_plot_dig$plot
plot3 <- surv_plot_dcg$plot
plot4 <- surv_plot_naive$plot

joint_survival <- plot1 | plot2 | plot3 | plot4
joint_survival

ggsave("deconveilCaseStudies/plots/main/Fig 5/joint_survival.png", dpi = 500, width = 12.0, height = 4.0, plot = joint_survival)


## Calculate Concordance Index for both datasets ##

# TCGA-BRCA #
clinical <- clinical_data
rna <- rna_log_norm  
selected_genes <- rownames(lasso_coefs_df)[lasso_coefs_df$`1` != 0]  
rna_selected <- rna[, selected_genes]
prognostic_score_tcga <- as.matrix(rna_selected) %*% lasso_coefs[selected_genes, ]
colnames(prognostic_score_tcga) <- "progn_score"
clinical <- cbind(clinical, prognostic_score_tcga)

clinical$time <- as.numeric(clinical$time)
clinical$event <- as.numeric(as.character(clinical$event))

cox_model_tcga <- coxph(Surv(time, event) ~ progn_score, data = clinical)

cindex_cox_tcga <- concordance(cox_model_tcga)
print(paste("C-index:", cindex_cox_tcga$concordance))

# METABRIC #
clinical_metabric <- clinical_metabric %>% remove_rownames() %>% column_to_rownames("patientID")
rna_metabric <- rna_metabric_dcomp

matching_genes <- intersect(selected_genes, colnames(rna_metabric))
rna_metabric_selected <- rna_metabric[,matching_genes]

rna_metabric_selected <- rna_metabric_selected[ rownames(rna_metabric_selected) %in% rownames(clinical_metabric),]
metabric_risk_scores <- as.matrix(rna_metabric_selected) %*% lasso_coefs_df[matching_genes, ]
colnames(metabric_risk_scores) <- "progn_score"
clinical_metabric <- merge(clinical_metabric, metabric_risk_scores, by = "row.names")

clinical_metabric$event <- as.numeric(as.character(clinical_metabric$event))
cox_model_metabric <- coxph(Surv(time, event) ~ progn_score, data = clinical_metabric)

cindex_cox_metabric <- concordance(cox_model_metabric)
print(paste("C-index:", cindex_cox_metabric$concordance))

