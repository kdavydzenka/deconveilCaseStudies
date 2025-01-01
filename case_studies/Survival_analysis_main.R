### Survival prognostic model using TCGA-BRCA dataset ###

setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/")
pkgs <- c("tidyverse", "survival", "glmnet", "survminer", "survcomp", "DESeq2", "forestplot", "caret", "randomForestSRC", "gridExtra", "RColorBrewer")
sapply(pkgs, require, character.only = TRUE)

# Define functions

prepare_data <- function(rna_data, clinical_data) {
  clinical_data <- clinical_data %>% 
    dplyr::select(bcr_patient_barcode, days_to_death, days_to_last_followup) %>%
    mutate(
      event = ifelse(!is.na(days_to_death), 1, 0),
      time = ifelse(!is.na(days_to_death), days_to_death, days_to_last_followup)
    ) %>%
    dplyr::select(bcr_patient_barcode, time, event) %>% 
    drop_na()
  colnames(clinical_data) <- c("patientID", "time", "event")
  
  # Ensure consistency between RNA and clinical data
  colnames(rna_data) <- substr(colnames(rna_data), 1, 12)
  clinical_data <- clinical_data[clinical_data$patientID %in% colnames(rna_data),]
  rna_data <- rna_data[, (colnames(rna_data) %in% clinical_data$patientID)]
  
  return(list(rna_data = rna_data, clinical_data = clinical_data))
}

apply_cn_correction <- function(rna, cn, apply_correction = TRUE) {
  if (apply_correction) {
    cn <- cn[(rownames(cn) %in% rownames(rna)), ]
    cn <- cn[, (colnames(cn) %in% colnames(rna))]
    cn <- apply(cn, 2, function(x) ifelse(x > 10, 10, x)) / 2
    cn[cn == 0] <- 0.001
    rna <- rna * cn
  }
  return(rna)
}

run_survival_analysis <- function(rna_data, clinical_data, output_file) {
  clinical <- clinical_data %>% remove_rownames() %>% column_to_rownames(var = "patientID")
  idx <- match(clinical$patientID, colnames(rna_data))
  rna_data <- rna_data[, idx]
  rna_data <- t(rna_data)
  data <- cbind(clinical, rna_data)
  
  # Cox proportional hazards model
  surv_object <- survival::Surv(time = clinical$time, event = clinical$event)
  cox_results <- data.frame(Gene = character(), p.value = numeric(), HR = numeric(), CI_lower = numeric(), CI_upper = numeric())
  
  for (gene in colnames(rna_data)) {
    cox_model <- coxph(surv_object ~ rna_data[, gene], data = data)
    summary_cox <- summary(cox_model)
    
    cox_results <- rbind(cox_results, data.frame(
      Gene = gene,
      p.value = summary_cox$coefficients[, "Pr(>|z|)"],
      HR = summary_cox$coefficients[, "exp(coef)"],
      CI_lower = summary_cox$conf.int[, "lower .95"],
      CI_upper = summary_cox$conf.int[, "upper .95"]
    ))
  }
  
  # Filter significant genes
  significant_genes <- cox_results %>% filter(p.value < 0.05)
  saveRDS(significant_genes, file = output_file)
  return(significant_genes)
}

build_prognostic_signature <- function(rna_data, clinical_data, significant_genes, lasso_output_file) {
  X <- as.matrix(rna_data[, significant_genes$Gene])
  y <- survival::Surv(time = clinical_data$time, event = clinical_data$event)
  
  # LASSO regularization
  cv_model <- cv.glmnet(X, y, family = "cox", alpha = 1)
  optimal_lambda <- cv_model$lambda.min
  
  lasso_model <- glmnet(X, y, family = "cox", alpha = 1, lambda = optimal_lambda)
  lasso_coefs <- coef(lasso_model, s = optimal_lambda)
  lasso_coefs_df <- as.data.frame(as.matrix(lasso_coefs))
  selected_genes <- rownames(lasso_coefs_df)[lasso_coefs_df$`1` != 0]
  
  saveRDS(selected_genes, file = lasso_output_file)
  return(selected_genes)
}

process_gene_group <- function(gene_group, rna_tumor, cn_tumor, clinical_data, apply_cn_correction_flag, output_dir) {
  rna_subset <- rna_tumor[rownames(rna_tumor) %in% gene_group$geneID, ]
  rna_corrected <- apply_cn_correction(rna_subset, cn_tumor, apply_correction = apply_cn_correction_flag)
  prepared_data <- prepare_data(rna_corrected, clinical_data)
  
  significant_genes <- run_survival_analysis(prepared_data$rna_data, prepared_data$clinical_data, 
                                             output_file = paste0(output_dir, "/significant_genes.RDS"))
  
  selected_genes <- build_prognostic_signature(prepared_data$rna_data, prepared_data$clinical_data, 
                                               significant_genes, lasso_output_file = paste0(output_dir, "/prognostic_signature.RDS"))
  
  return(selected_genes)
}

# Load data
gene_groups <- readRDS("TCGA/BRCA/case_study/gene_groups.RDS")
clinical_data <- readRDS("TCGA/brca/clinical_full.RDS")
rna_tumor <- readRDS("TCGA/BRCA/rna_tumor.RDS")
cnv_tumor <- readRDS("TCGA/BRCA/cnv_tumor.RDS")

rna_subset <- rna_tumor[rownames(rna_tumor) %in% gene_group$geneID, ]
prepared_data <- prepare_data(rna_subset, clinical_data)
significant_genes <- run_survival_analysis(prepared_data$rna_data, prepared_data$clinical_data, 
                                           output_file = paste0(output_dir, "/significant_genes.RDS"))

output_dir <- "TCGA/BRCA/case_study"

selected_genes_ds <- process_gene_group(
  gene_group = gene_groups[["d_sensitive"]],
  rna_tumor = rna_tumor,
  cn_tumor = cnv_tumor,
  clinical_data = clinical_data,
  apply_cn_correction_flag = FALSE,
  output_dir = paste0(output_dir, "/rna_ds")
)

selected_genes_dins <- process_gene_group(
  gene_group = gene_groups[["d_insensitive"]],
  rna_tumor = rna_tumor,
  cn_tumor = cnv_tumor,
  clinical_data = clinical_data,
  apply_cn_correction_flag = TRUE,
  output_dir = paste0(output_dir, "/rna_dins")
)

selected_genes_dcomp <- process_gene_group(
  gene_group = gene_groups[["d_compensated"]],
  rna_tumor = rna_tumor,
  cn_tumor = cnv_tumor,
  clinical_data = clinical_data,
  apply_cn_correction_flag = TRUE,
  output_dir = paste0(output_dir, "/rna_dcomp")
)

# Main workflow

# Load data
gene_groups <- readRDS("TCGA/BRCA/case_study/gene_groups.RDS")
clinical_data <- readRDS("TCGA/brca/clinical_full.RDS")
rna_tumor <- readRDS("TCGA/BRCA/rna_tumor.RDS")
cnv_tumor <- readRDS("TCGA/BRCA/cnv_tumor.RDS")

colnames(rna_tumor) <- substr(colnames(rna_tumor), 1, 12)
colnames(cnv_tumor) <- substr(colnames(rna_tumor), 1, 12)
clinical_data <- clinical_data[clinical_data$bcr_patient_barcode %in% colnames(rna_tumor),]
rna_tumor <- rna_tumor[,(colnames(rna_tumor) %in% clinical_data$bcr_patient_barcode)]

clinical_data <- clinical_data %>% dplyr::select(bcr_patient_barcode, days_to_death, days_to_last_followup)
clinical_data$event <- ifelse(!is.na(clinical_data$days_to_death), 1, 0)
clinical_data$time <- ifelse(!is.na(clinical_data$days_to_death), 
                            clinical_data$days_to_death, 
                            clinical_data$days_to_last_followup)

clinical_data <- clinical_data %>%
  dplyr::select(bcr_patient_barcode, time, event) %>%   
  drop_na() 
colnames(clinical_data) <- c("patientID", "time", "event")

rna_ds <- rna_tumor[(rownames(rna_tumor) %in% gene_groups[["d_sensitive"]]$geneID),]
rna_dins <- rna_tumor[(rownames(rna_tumor) %in% gene_groups[["d_insensitive"]]$geneID),]
rna_dcomp <- rna_tumor[(rownames(rna_tumor) %in% gene_groups[["d_compensated"]]$geneID),]

# GE normalization 
rna_log_norm<- rna_dcomp %>% as.matrix() %>% DESeq2::varianceStabilizingTransformation()
cnv_tumor <- cnv_tumor[(rownames(cnv_tumor) %in% rownames(rna_log_norm)) ,]
cnv_tumor <- cnv_tumor[,(colnames(cnv_tumor) %in% colnames(rna_log_norm)) ]
cnv_tumor <- apply(cnv_tumor, 2, function(x) ifelse(x > 10, 10, x))
cnv_tumor <- cnv_tumor/2
cnv_tumor[cnv_tumor == 0] <- 0.001
rna_log_norm <- rna_log_norm * cnv_tumor

# Reorder patients indexing
idx <- match(clinical_data$patientID, colnames(rna_log_norm))
rna_log_norm <- rna_log_norm[,idx]
rna_log_norm <- t(rna_log_norm)
clinical_data <- clinical_data %>% remove_rownames %>% column_to_rownames(var="patientID") 
data <- cbind(clinical_data, rna_log_norm)

rna <- data %>% select(-time, -event)
clinical <- data %>% select(time, event)


## Initial Gene Selection Using Cox Proportional Hazards Model ##

surv_object <- survival::Surv(time = clinical$time, event = clinical$event)
cox_results <- data.frame(Gene = character(), p.value = numeric(), HR = numeric(), CI_lower = numeric(), CI_upper = numeric())

for (gene in colnames(rna)) {
  cox_model <- coxph(surv_object ~ rna[, gene], data = data)
  summary_cox <- summary(cox_model)
  
  cox_results <- rbind(cox_results, data.frame(
    Gene = gene,
    p.value = summary_cox$coefficients[,"Pr(>|z|)"],
    HR = summary_cox$coefficients[,"exp(coef)"],
    CI_lower = summary_cox$conf.int[,"lower .95"],
    CI_upper = summary_cox$conf.int[,"upper .95"]
  ))
}

significant_genes <- cox_results %>% filter(p.value < 0.05)

saveRDS(significant_genes, file = "TCGA/BRCA/case_study/significant_genes_dcomp.RDS")

# The Hazard Ratio (HR) and confidence intervals (CI) give an indication of the prognostic impact of each gene. 
# Genes with HR > 1 may be associated with poor prognosis, while HR < 1 suggests good prognosis.
# Feature Selection Using LASSO (Penalized Cox Model) #
# Build Gene Prognostic Signature Using Lasso Regression (selects the most important genes and avoids overfitting) 

# LASSO for further gene selection
X <- as.matrix(rna[, significant_genes$Gene])  
y <- surv_object  

# LASSO regularization
cv_model <- cv.glmnet(X, y, family = "cox", alpha = 1)
optimal_lambda <- cv_model$lambda.min
plot(cv_model)

lasso_model <- glmnet(X, y, family = "cox", alpha = 1, lambda = optimal_lambda)

lasso_coefs <- coef(lasso_model, s = optimal_lambda)
lasso_coefs_df <- as.data.frame(as.matrix(lasso_coefs))
selected_genes <- rownames(lasso_coefs_df)[lasso_coefs_df$`1` != 0]
selected_genes


# Calculate the prognostic score for each patient
prognostic_score <- X[, selected_genes] %*% lasso_coefs[selected_genes,]
colnames(prognostic_score) <- c("progn_score")
clinical <- cbind(clinical, prognostic_score)


#Stratify Patients into risk groups
clinical$risk_group <- ifelse(clinical$progn_score > median(clinical$progn_score), 
                                    "High risk", "Low risk")
                                    
surv_fit <- survfit(Surv(time, event) ~ risk_group, data = clinical)


#ggsave("CN-aware-DGE/case_studies/plots/brca/KM_survival.png", dpi = 400, width = 5.0, height = 10.0, plot = surv_plot)

# Forest plot using ggplot

sel_genes_data <- significant_genes %>% dplyr::filter(Gene %in% c(selected_genes))

saveRDS(sel_genes_data, file = "TCGA/BRCA/case_study/prognostic_signature_dcomp.RDS")

plot_data <- sel_genes_data %>%
  mutate(
    gene = factor(Gene, levels = rev(Gene)),  # Reversing for top-to-bottom order
    is_summary = FALSE, # Summary flag for first row
    box_color = ifelse(HR > 1, "#D60C00FF", "blue3")  # Color based on HR values
  )

plot_data$is_summary[1] <- TRUE

custom_colors <- c("category1" = "#D60C00FF", "category2" = "blue3")

forest_plot <- ggplot(plot_data, aes(x = HR, y = gene)) +
  geom_pointrange(aes(xmin = CI_lower, xmax = CI_upper, color = box_color), 
                  size = 0.9, show.legend = TRUE) +
  #geom_point(data = subset(plot_data, is_summary == TRUE), aes(x = HR), 
             #shape = 23, fill = "royalblue", color = "royalblue", size = 3) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
  labs(x = "Hazard Ratio", y = NULL) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 15, color = "black"),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 15, face = "plain"),
    panel.grid.minor = element_blank(),  
    panel.grid.major.y = element_blank()) +
  scale_color_identity() +
  #scale_color_manual(values = custom_colors)+
  coord_cartesian(xlim = c(min(plot_data$CI_lower), max(plot_data$CI_upper))) +
  scale_x_continuous(trans = "identity")
forest_plot

# Create text table for gene labels, HR, and p-values (mimicking table_text)
table_text <- plot_data %>%
  transmute(
    Gene = as.character(Gene),
    `Hazard Ratio` = sprintf("%.2f", HR),
    `p-value` = sprintf("%.3f", p.value)
  )

# Create a table plot for the labels using gridExtra
table_plot <- tableGrob(table_text, rows = NULL, 
                        theme = ttheme_minimal(core = list(fg_params = list(hjust = 0, x = 0.1))))

# Combine forest plot and table plot
grid.arrange(forest_plot, table_plot, ncol = 2, widths = c(1, 2))



### Perform prognostic model validation using external cohort METABRIC ###

# Data preprocessin
rna_metabric <- read_delim("brca_metabric/data_mrna_illumina_microarray.txt")
cn_metabric <- read_delim("brca_metabric/data_cna.txt")
clinical_metabric <- read_delim("brca_metabric/data_clinical_patient.txt")

rna_metabric <- rna_metabric[!duplicated(rna_metabric$Hugo_Symbol), ]
rna_metabric <- rna_metabric %>% dplyr::select(-Entrez_Gene_Id) %>% 
  remove_rownames() %>% column_to_rownames(var = "Hugo_Symbol")

cn_metabric <- cn_metabric[!duplicated(cn_metabric$Hugo_Symbol), ]
cn_metabric <- cn_metabric %>% dplyr::select(-Entrez_Gene_Id) %>% 
  remove_rownames() %>% column_to_rownames(var = "Hugo_Symbol")

clinical_metabric <- clinical_metabric[-c(1:4),]
#clinical_metabric <- clinical_metabric[clinical_metabric$Cohort == 6,]
#clinical_metabric <- clinical_metabric[clinical_metabric$Cohort %in% c(1,4,5), ]
clinical_metabric <- clinical_metabric[clinical_metabric$Cohort %in% c(1), ] # DSGs, DIGs
clinical_metabric <- clinical_metabric %>% dplyr::select("#Patient Identifier", "Overall Survival (Months)", "Overall Survival Status")
colnames(clinical_metabric) <- c("patientID", "time", "event")
clinical_metabric$time <- as.numeric(clinical_metabric$time)
clinical_metabric$time <- round(clinical_metabric$time)

clinical_metabric <- clinical_metabric %>% na.omit() %>% 
  mutate(event = case_when(
    event == "1:DECEASED" ~ "1",
    event == "0:LIVING" ~ "0"
  ))

rna_metabric <- rna_metabric[,colnames(rna_metabric) %in% clinical_metabric$patientID]
cn_metabric <- cn_metabric[,colnames(cn_metabric) %in% clinical_metabric$patientID]
clinical_metabric <- clinical_metabric[clinical_metabric$patientID %in% colnames(rna_metabric),]

cn_metabric <- cn_metabric[rownames(cn_metabric) %in% rownames(rna_metabric),]
rna_metabric <- rna_metabric[rownames(rna_metabric) %in% rownames(cn_metabric),]

# Reorder patients indexing
idx <- match(clinical_metabric$patientID, colnames(rna_metabric))
rna_metabric <- rna_metabric[,idx]
cn_metabric <- cn_metabric[,idx]

idx <- match(rownames(rna_metabric), rownames(cn_metabric))
cn_metabric <- cn_metabric[idx,]

cn_metabric <- cn_metabric %>%
  mutate_all(~ case_when(
    . == -2 ~ 0,
    . == -1 ~ 1,
    . == 0  ~ 2,
    . == 1  ~ 3,
    . == 2  ~ 4
  ))

cn_metabric[cn_metabric == 0] <- 0.001
cn_metabric <- cn_metabric/2


# Subset METABRIC data for the same genes used in TCGA-BRCA
rna_metabric_dcomp <- rna_metabric[(rownames(rna_metabric) %in% rownames(rna_dcomp)),]
cn_metabric_dcomp <- cn_metabric[(rownames(cn_metabric) %in% rownames(rna_dcomp)),]

rna_metabric_dcomp <- rna_metabric_dcomp[,(colnames(rna_metabric_dcomp) %in% colnames(cn_metabric_dcomp))]
cn_metabric_dcomp <- cn_metabric_dcomp[,(colnames(cn_metabric_dcomp) %in% colnames(rna_metabric_dcomp))]
idx <- match(colnames(rna_metabric_dcomp), colnames(cn_metabric_dcomp))
cn_metabric_dcomp <- cn_metabric_dcomp[,idx]

rna_metabric_comp <- rna_metabric_dcomp * cn_metabric_dcomp
rna_metabric_dcomp <- t(rna_metabric_dcomp)
rna_metabric_dcomp <- rna_metabric_dcomp %>% 
  as.data.frame() %>% 
  dplyr::mutate(patientID = rownames(rna_metabric_dcomp))

metabric_combined <- merge(clinical_metabric, rna_metabric_dcomp, by = "patientID")

# Ensure Matching Gene Names
genes_in_lasso <- selected_genes
genes_in_metabric <- colnames(metabric_combined)
matching_genes <- intersect(genes_in_lasso, genes_in_metabric)
matching_genes <- intersect(colnames(rna_metabric_dcomp), rownames(coef(lasso_model)))
metabric_filtered <- metabric_combined[, matching_genes]
lasso_coefficients <- lasso_coefs_df[matching_genes,]

# Calculate the risk score for METABRIC data using the LASSO model coefficients
metabric_risk_scores <- as.matrix(metabric_filtered) %*% lasso_coefficients
metabric_combined$risk_score <- as.vector(metabric_risk_scores)
metabric_combined <- metabric_combined %>% na.omit()

# Stratify METABRIC patients into high/low risk groups based on median risk score
metabric_combined$risk_group <- ifelse(metabric_combined$risk_score > median(metabric_combined$risk_score), "High", "Low")
metabric_combined$event <- as.numeric(as.character(metabric_combined$event))
metabric_survfit <- survfit(Surv(time, event) ~ risk_group, data = metabric_combined)

ggsurvplot(metabric_survfit, data = metabric_combined, pval = TRUE, risk.table = TRUE)

surv_plot <- ggsurvplot(metabric_survfit, data = metabric_combined, pval = TRUE, 
                        conf.int = T, risk.table = TRUE, palette = c("#DD2A7B", "#515BD4"),  
                        title = "",  xlab = "Time (days)",  ylab = "Survival Probability",
                        font.main = c(16, "bold", "black"),  
                        font.x = c(16, "plain"),  
                        font.y = c(16, "plain"),  
                        font.tickslab = c(16, "plain"),  
                        legend = "bottom",  
                        legend.title = "Risk group",  
                        legend.labs = c("High risk", "Low risk"),  
                        font.legend = c(16, "plain"),
                        risk.table.height = 0.25,  
                        risk.table.y.text = TRUE,  
                        risk.table.fontsize = 6,  
                        risk.table.title = "Number at risk",  
                        risk.table.col = "strata",  
                        pval.coord = c(1000, 0.2),  pval.size = 5,  
                        ggtheme = theme_classic2()+
                          theme(
                            strip.text = element_text(size = 14, face = "plain") 
                          ),
                        risk.table.title.fontface = "bold"  
)

surv_plot$plot <- surv_plot$plot + 
  theme(
    legend.text = element_text(size = 16)) 
surv_plot$table <- surv_plot$table + 
  theme(
    text = element_text(size = 14),
    strip.text = element_text(size = 14, face = "plain"))
surv_plot

#ggsave("deconveilCaseStudies/case_studies/plots/brca/KM_ds_2.png", dpi = 400, width = 5.0, height = 10.0, plot = surv_plot)


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

