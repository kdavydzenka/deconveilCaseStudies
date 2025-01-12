### Survival analysis of TCGA-BRCA dataset ###

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
  rna_data<- rna_data %>% as.matrix() %>% DESeq2::varianceStabilizingTransformation()
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
  idx <- match(rownames(clinical), colnames(rna_data))
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
  rna_data <- t(rna_data)
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
  prepared_data <- prepare_data(rna_subset, clinical_data)
  rna_corrected <- apply_cn_correction(prepared_data[["rna_data"]], cn_tumor, apply_correction = apply_cn_correction_flag)
  
  clinical_data <- prepared_data[["clinical_data"]]
  rna_data <- rna_corrected
  significant_genes <- run_survival_analysis(rna_data, clinical_data, 
                                             output_file = paste0(output_dir, "/significant_genes.RDS"))
  
  selected_genes <- build_prognostic_signature(rna_data, clinical_data, 
                                               significant_genes, lasso_output_file = paste0(output_dir, "/prognostic_signature.RDS"))
  
  return(selected_genes)
}

# Load data
gene_groups <- readRDS("TCGA/BRCA/case_study/gene_groups.RDS")
clinical_data <- readRDS("TCGA/brca/clinical_full.RDS")
rna_tumor <- readRDS("TCGA/BRCA/rna_tumor.RDS")
cnv_tumor <- readRDS("TCGA/BRCA/cnv_tumor.RDS")

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
