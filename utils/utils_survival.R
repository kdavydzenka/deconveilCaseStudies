
annotate_results <- function(results, lfc_cut, pval_cut, method, tumor_type) {
  results %>%
    mutate(
      isDE = (abs(log2FoldChange) >= lfc_cut) & (padj <= pval_cut),
      DEtype = if_else(!isDE, "n.s.", if_else(log2FoldChange > 0, "Up-reg", "Down-reg")),
      method = method,
      tumor_type = tumor_type
    ) %>%
    remove_rownames() %>%
    column_to_rownames(var = "X")
}


prepare_survival_data <- function(rna_subset, clinical_data, cnv_tumor = NULL, apply_cnv = FALSE) {
  # VST normalization
  vst_data <- DESeq2::varianceStabilizingTransformation(as.matrix(rna_subset))
  
  # Apply CNV correction if required
  if (apply_cnv && !is.null(cnv_tumor)) {
    cnv_filtered <- cnv_tumor[rownames(cnv_tumor) %in% rownames(vst_data), colnames(cnv_tumor) %in% colnames(vst_data)]
    cnv_filtered <- apply(cnv_filtered, 2, function(x) ifelse(x > 10, 10, x)) / 2
    cnv_filtered[cnv_filtered == 0] <- 0.001
    vst_data <- vst_data * cnv_filtered
  }
  
  # Reorder and merge with clinical data
  idx <- match(clinical_data$patientID, colnames(vst_data))
  vst_data <- vst_data[, idx]
  vst_data <- t(vst_data)
  
  clinical_data <- clinical_data %>%
    remove_rownames() %>%
    column_to_rownames(var = "patientID")
  
  full_data <- cbind(clinical_data, vst_data)
  list(
    rna = full_data %>% select(-time, -event),
    clinical = full_data %>% select(time, event)
  )
}


run_cox_analysis <- function(rna_data, clinical_data, pval_threshold = 0.05) {
  
  data <- cbind(clinical_data, rna_data)
  
  # Build survival object
  surv_object <- survival::Surv(time = clinical_data$time, event = clinical_data$event)
  
  # Loop over genes and compute Cox model
  cox_results <- purrr::map_dfr(colnames(rna_data), function(gene) {
    model <- coxph(surv_object ~ rna_data[, gene], data = data)
    summary_model <- summary(model)
    
    data.frame(
      Gene = gene,
      p.value = summary_model$coefficients[,"Pr(>|z|)"],
      HR = summary_model$coefficients[,"exp(coef)"],
      CI_lower = summary_model$conf.int[,"lower .95"],
      CI_upper = summary_model$conf.int[,"upper .95"]
    )
  })
  
  dplyr::filter(cox_results, p.value < pval_threshold)
}


fit_lasso_and_score <- function(rna_data, clinical_data, significant_genes, pval_threshold = 0.05, plot_cv = FALSE) {
  sig_genes <- significant_genes %>%
    dplyr::filter(p.value < pval_threshold)
  
  # Subset expression matrix
  X <- as.matrix(rna_data[, sig_genes$Gene])
  y <- survival::Surv(clinical_data$time, clinical_data$event)
  
  # Fit cross-validated LASSO model
  cv_model <- glmnet::cv.glmnet(X, y, family = "cox", alpha = 1)
  if (plot_cv) plot(cv_model)
  optimal_lambda <- cv_model$lambda.min
  
  # Final model with optimal lambda
  lasso_model <- glmnet::glmnet(X, y, family = "cox", alpha = 1, lambda = optimal_lambda)
  lasso_coefs <- coef(lasso_model, s = optimal_lambda)
  
  # Convert sparse matrix to data frame safely
  lasso_coefs_df <- as.data.frame(as.matrix(lasso_coefs))
  colnames(lasso_coefs_df) <- "coef"
  
  # Get selected genes (non-zero coefficients)
  selected_genes <- rownames(lasso_coefs_df)[lasso_coefs_df$coef != 0]
  
  if (length(selected_genes) == 0) {
    warning("No genes selected by LASSO.")
    return(NULL)
  }
  
  # Compute prognostic score
  prognostic_score <- X[, selected_genes, drop = FALSE] %*% as.matrix(lasso_coefs_df[selected_genes, , drop = FALSE])
  colnames(prognostic_score) <- "progn_score"
  
  # Merge with clinical data
  clinical_with_score <- cbind(clinical_data, prognostic_score)
  clinical_with_score$risk_group <- ifelse(
    clinical_with_score$progn_score > median(clinical_with_score$progn_score),
    "High risk", "Low risk"
  )
  
  list(
    selected_genes = selected_genes,
    lasso_coefficients = lasso_coefs_df[lasso_coefs_df$coef != 0, , drop = FALSE],
    clinical = clinical_with_score,
    sel_genes_data = dplyr::filter(significant_genes, Gene %in% selected_genes)
  )
}


calculate_cindex <- function(rna, clinical, lasso_coefs_df, 
                             lasso_coef_col = "coef", dataset_name = "Dataset", 
                             transform_rna = FALSE, 
                             genes_are_rows = TRUE,
                             clinical_id_col = "patientID",
                             time_col = "time",
                             event_col = "event") {
  
  cat("Calculating C-index for:", dataset_name, "\n")
  
  # Optional VST transform (expect rna genes x samples)
  if (transform_rna) {
    if (!requireNamespace("DESeq2", quietly = TRUE)) {
      stop("DESeq2 package is required for transform_rna = TRUE.")
    }
    rna <- DESeq2::varianceStabilizingTransformation(as.matrix(rna))
  }
  
  # If genes are rows, transpose to samples x genes (needed for matrix mult)
  if (genes_are_rows) {
    rna <- t(rna)
  }
  
  # Extract gene names from RNA and LASSO
  gene_names_rna <- colnames(rna)
  gene_names_lasso <- rownames(lasso_coefs_df)
  
  # Find matching genes and their positions (match keeps order)
  matching_idx_rna <- match(gene_names_lasso, gene_names_rna)
  matched_genes <- gene_names_lasso[!is.na(matching_idx_rna)]
  matching_idx_rna <- matching_idx_rna[!is.na(matching_idx_rna)]
  
  cat("Genes in RNA:", length(gene_names_rna), "\n")
  cat("Genes in LASSO model:", length(gene_names_lasso), "\n")
  cat("Matching genes:", length(matched_genes), "\n")
  
  if (length(matched_genes) < 2) {
    stop("Fewer than 2 matching genes found between RNA and LASSO coefficients.")
  }
  
  if (length(matched_genes) < length(gene_names_lasso)) {
    warning(sprintf("Only %d out of %d LASSO genes matched RNA data. Proceeding with matched genes.",
                    length(matched_genes), length(gene_names_lasso)))
  }
  
  # Subset RNA and coefficients by matched genes
  rna_selected <- rna[, matching_idx_rna, drop = FALSE]
  coefs_selected <- lasso_coefs_df[matched_genes, lasso_coef_col]
  
  # Calculate prognostic score (samples x 1)
  prognostic_score <- as.matrix(rna_selected) %*% as.numeric(coefs_selected)
  colnames(prognostic_score) <- "progn_score"
  
  # Prepare clinical data
  stopifnot(clinical_id_col %in% colnames(clinical),
            time_col %in% colnames(clinical),
            event_col %in% colnames(clinical))
  
  # Match samples in clinical and prognostic score
  clinical_samples <- clinical[[clinical_id_col]]
  score_samples <- rownames(prognostic_score)
  matching_samples <- intersect(clinical_samples, score_samples)
  
  if (length(matching_samples) == 0) {
    stop("No matching samples between RNA and clinical data.")
  }
  
  # Filter clinical and prognostic score by matching samples (keep order)
  clinical_idx <- match(matching_samples, clinical_samples)
  score_idx <- match(matching_samples, score_samples)
  
  clinical_sub <- clinical[clinical_idx, , drop = FALSE]
  prognostic_score_sub <- prognostic_score[score_idx, , drop = FALSE]
  
  # Add prognostic score column to clinical data
  clinical_sub$progn_score <- prognostic_score_sub[, 1]
  
  # Ensure numeric survival columns
  clinical_sub[[time_col]] <- as.numeric(clinical_sub[[time_col]])
  clinical_sub[[event_col]] <- as.numeric(as.character(clinical_sub[[event_col]]))
  
  # Fit Cox model and calculate C-index
  cox_model <- survival::coxph(survival::Surv(clinical_sub[[time_col]], clinical_sub[[event_col]]) ~ progn_score, data = clinical_sub)
  cindex <- survival::concordance(cox_model)
  
  cat(sprintf("C-index for %s: %.4f\n", dataset_name, cindex$concordance))
  
  return(list(
    model = cox_model,
    cindex = cindex$concordance,
    n_genes = length(matched_genes),
    genes = matched_genes
  ))
}




