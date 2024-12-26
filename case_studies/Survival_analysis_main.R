### Survival prognostic model using TCGA-BRCA dataset ###

setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/")
pkgs <- c("tidyverse", "survival", "glmnet", "survminer", "survcomp", "DESeq2", "forestplot", "caret", "randomForestSRC", "gridExtra", "RColorBrewer")
sapply(pkgs, require, character.only = TRUE)

# Prepare data #
d_compensated <- readRDS("TCGA/brca/case_study/d_compensated_genes_2.RDS")
clinical_data <- readRDS("TCGA/brca/clinical_full.RDS")
rna <- read.csv("TCGA/brca/case_study/rna.csv")
rna_tumor <- readRDS("TCGA/brca/rna_tumor.RDS")
cnv_tumor <- readRDS("TCGA/brca/cnv_tumor.RDS")

colnames(rna_tumor) <- substr(colnames(rna_tumor), 1, 12)
colnames(cnv_tumor) <- substr(colnames(rna_tumor), 1, 12)
rna_tumor <- rna_tumor[(rownames(rna_tumor) %in% rna$X),]
cnv_tumor <- cnv_tumor[(rownames(cnv_tumor) %in% rna$X),]
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

rna_d_compensated_genes <- rna_tumor[(rownames(rna_tumor) %in% rownames(d_compensated)),]

# GE normalization 
rna_log_norm<- rna_d_compensated_genes %>% as.matrix() %>% DESeq2::varianceStabilizingTransformation()
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

saveRDS(significant_genes, file = "TCGA/brca/case_study/significant_genes_cox.RDS")

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

saveRDS(selected_genes, file = "TCGA/brca/case_study/prognostic_signature_aware.RDS")


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
clinical_metabric <- clinical_metabric[clinical_metabric$Cohort %in% c(1,4,5), ]
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
rna_metabric_d_compensated <- rna_metabric[(rownames(rna_metabric) %in% rownames(rna_d_compensated_genes)),]
cn_metabric_d_compensated <- cn_metabric[(rownames(cn_metabric) %in% rownames(rna_d_compensated_genes)),]

rna_metabric_d_compensated <- rna_metabric_d_compensated[,(colnames(rna_metabric_d_compensated) %in% colnames(cn_metabric_d_compensated))]
cn_metabric_d_compensated <- cn_metabric_d_compensated[,(colnames(cn_metabric_d_compensated) %in% colnames(rna_metabric_d_compensated))]
idx <- match(colnames(rna_metabric_d_compensated), colnames(cn_metabric_d_compensated))
cn_metabric_d_compensated <- cn_metabric_d_compensated[,idx]

rna_metabric_d_compensated <- rna_metabric_d_compensated * cn_metabric_d_compensated
rna_metabric_d_compensated <- t(rna_metabric_d_compensated)
rna_metabric_d_compensated <- rna_metabric_d_compensated %>% 
  as.data.frame() %>% 
  dplyr::mutate(patientID = rownames(rna_metabric_d_compensated))

metabric_combined <- merge(clinical_metabric, rna_metabric_d_compensated, by = "patientID")

# Ensure Matching Gene Names
genes_in_lasso <- selected_genes
genes_in_metabric <- colnames(metabric_combined)
matching_genes <- intersect(genes_in_lasso, genes_in_metabric)
matching_genes <- intersect(colnames(rna_metabric_d_compensated), rownames(coef(lasso_model)))
metabric_filtered <- metabric_combined[, matching_genes]
lasso_coefficients <- lasso_coefs_df[matching_genes,]

# Calculate the risk score for METABRIC data using the LASSO model coefficients
metabric_risk_scores <- as.matrix(metabric_filtered) %*% lasso_coefficients
metabric_combined$risk_score <- as.vector(metabric_risk_scores)

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

#ggsave("CN-aware-DGE/case_studies/plots/brca/km_TCGA.png", dpi = 400, width = 5.0, height = 10.0, plot = surv_plot)


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
print(paste("C-index:", cindex_cox$concordance))

# METABRIC #
clinical_metabric <- clinical_metabric %>% remove_rownames() %>% column_to_rownames("patientID")
rna_metabric <- rna_metabric_d_compensated 

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

