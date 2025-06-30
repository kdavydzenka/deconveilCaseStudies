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

sig_genes_ds <- readRDS("deconveilCaseStudies/plots/main/Fig 5/rds/cox_significant_dsg.RDS")
sig_genes_dins <- readRDS("deconveilCaseStudies/plots/main/Fig 5/rds/cox_significant_dig.RDS")
sig_genes_dcomp <- readRDS("deconveilCaseStudies/plots/main/Fig 5/rds/cox_significant_dcg.RDS")
sig_genes_naive <- readRDS("deconveilCaseStudies/plots/main/Fig 5/rds/cox_significant_CNnaive.RDS")

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

process_metabric <- function(rna_input, lasso_model, label = "DS") {
  message("Processing METABRIC for: ", label)
  
  # Load METABRIC RNA and CN data
  rna_metabric <- readRDS("deconveilCaseStudies/TCGA/brca_metabric/rna_metabric.rds")
  cn_metabric <- readRDS("deconveilCaseStudies/TCGA/brca_metabric/cn_metabric.rds")
  clinical_metabric <- readRDS("deconveilCaseStudies/TCGA/brca_metabric/clinical_metabric.rds")
  
  # Filter for genes in input
  gene_filter <- intersect(rownames(rna_metabric), rownames(rna_input))
  rna_filtered <- rna_metabric[gene_filter, , drop = FALSE]
  cn_filtered <- cn_metabric[gene_filter, , drop = FALSE]
  
  # Match samples
  shared_samples <- intersect(colnames(rna_filtered), colnames(cn_filtered))
  rna_filtered <- rna_filtered[, shared_samples, drop = FALSE]
  cn_filtered <- cn_filtered[, shared_samples, drop = FALSE]
  
  # Transpose and add patient ID
  rna_df <- as.data.frame(t(rna_filtered))
  rna_df$patientID <- rownames(rna_df)
  
  # Merge with clinical data
  metabric_combined <- merge(clinical_metabric, rna_df, by = "patientID")
  
  # Extract and align genes
  selected_genes <- lasso_model[["selected_genes"]]
  genes_in_data <- intersect(selected_genes, colnames(metabric_combined))
  metabric_filtered <- metabric_combined[, genes_in_data, drop = FALSE]
  
  # Get lasso coefficients
  coefs <- lasso_model[["lasso_coefficients"]]
  
  # If coefs is a data frame with 1 column, convert to named vector
  if (is.data.frame(coefs)) {
    if (ncol(coefs) == 1) {
      coefs <- coefs[[1]]
      names(coefs) <- rownames(lasso_model[["lasso_coefficients"]])
    } else {
      stop("Expected a single-column data frame for coefficients.")
    }
  }
  
  # Align coefficient vector with gene order
  coefs <- coefs[genes_in_data]
  
  # Compute risk score
  metabric_risk_scores <- as.matrix(metabric_filtered) %*% as.numeric(coefs)
  metabric_combined$risk_score <- as.vector(metabric_risk_scores)
  metabric_combined <- metabric_combined %>% na.omit()
  
  # Stratify patients
  metabric_combined$risk_group <- ifelse(
    metabric_combined$risk_score > median(metabric_combined$risk_score),
    "High", "Low"
  )
  
  metabric_combined$event <- as.numeric(as.character(metabric_combined$event))
  
  # Fit survival model
  metabric_survfit <- survfit(Surv(time, event) ~ risk_group, data = metabric_combined)
  
  return(list(
    data = metabric_combined,
    survfit = metabric_survfit
  ))
}

results_list <- list(
  ds = process_metabric(rna_ds, lasso_dsg, "DS"),
  dins = process_metabric(rna_dins, lasso_dig, "DINS"),
  dcomp = process_metabric(rna_dcomp, lasso_dcg, "DCOMP"),
  naive = process_metabric(rna_CNnaive, lasso_dsg, "CN-naive")
)

# Plot
surv_plot_dsg <- plot_survival(results_list$ds$survfit, results_list$ds$data, title = "DSGs")
surv_plot_dig <- plot_survival(results_list$dins$survfit, results_list$dins$data, title = "DIGs")
surv_plot_dcg <- plot_survival(results_list$dcomp$survfit, results_list$dcomp$data, title = "DCGs")
surv_plot_naive <- plot_survival(results_list$naive$survfit, results_list$naive$data, title = "CN-naive")

surv_plot_dsg$plot <- surv_plot_dsg$plot + theme(legend.position = "none")
surv_plot_dig$plot <- surv_plot_dig$plot + theme(legend.position = "none")
surv_plot_naive$plot <- surv_plot_naive$plot + theme(legend.position = "none")

# Combine plots using patchwork
joint_survival <- surv_plot_dsg$plot | surv_plot_dig$plot | surv_plot_dcg$plot | surv_plot_naive$plot
joint_survival

ggsave("deconveilCaseStudies/plots/main/Fig 5/joint_survival.png", dpi = 500, width = 12.0, height = 4.0, plot = joint_survival)


## Calculate Concordance Index for both datasets ##

# TCGA
#rownames(clinical_data) <- clinical_data$patientID
#tcga_results <- calculate_cindex(
  #rna = rna_ds,
  #clinical = clinical_data,
  #lasso_coefs_df = lasso_dsg[["lasso_coefficients"]],
  #lasso_coef_col = "coef",
  #dataset_name = "TCGA-BRCA",
  #transform_rna = TRUE
#)


# METABRIC

for (name in names(results_list)) {
  data <- results_list[[name]][["data"]]
  
  # Ensure the required column exists before setting rownames
  if ("patientID" %in% colnames(data)) {
    rownames(data) <- data[["patientID"]]
    
    # Remove survival and ID columns
    data <- data[, !(colnames(data) %in% c("time", "event", "patientID")) ]
    
    # Transpose the data matrix
    results_list[[name]][["data"]] <- t(data)
  } else {
    warning(paste("Skipping", name, ": 'patientID' column not found."))
  }
}

rownames(clinical_metabric) <- clinical_metabric$patientID

datasets_to_run <- c("ds", "dins", "dcomp", "naive")

metabric_results_list <- lapply(datasets_to_run, function(name) {
  if (!is.null(results_list[[name]]) && !is.null(results_list[[name]][["data"]])) {
    cat("Running calculate_cindex for:", name, "\n")
    tryCatch({
      calculate_cindex(
        rna = results_list[[name]][["data"]],
        clinical = clinical_metabric,
        lasso_coefs_df = lasso_dsg[["lasso_coefficients"]],
        lasso_coef_col = "coef",
        dataset_name = paste0("METABRIC-", name),
        transform_rna = FALSE
      )
    }, error = function(e) {
      warning(paste("Error for dataset", name, ":", e$message))
      return(NULL)
    })
  } else {
    warning(paste("Skipping", name, ": 'data' not found"))
    return(NULL)
  }
})

names(metabric_results_list) <- datasets_to_run



