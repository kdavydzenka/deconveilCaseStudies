
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


validate_signature_on_metabric <- function(
    rna_metabric,
    cn_metabric,
    clinical_metabric,
    gene_set,
    lasso_coefs_df,
    selected_genes,
    normalize_cnv = TRUE,
    group_name = "group"
) {
  # Subset and normalize gene expression
  common_genes <- intersect(rownames(rna_metabric), gene_set)
  if (length(common_genes) == 0) stop("No common genes between RNA and provided gene_set.")
  rna_subset <- rna_metabric[common_genes, , drop = FALSE]
  
  if (normalize_cnv && !is.null(cn_metabric)) {
    cn_subset <- cn_metabric[common_genes, colnames(rna_subset), drop = FALSE]
    rna_subset <- rna_subset * cn_subset
  }
  
  # Transpose and add patientID
  rna_df <- as.data.frame(t(rna_subset)) %>%
    dplyr::mutate(patientID = rownames(.))
  
  # Merge with clinical data
  merged_df <- merge(clinical_metabric, rna_df, by = "patientID")
  merged_df <- merged_df %>% dplyr::filter(!is.na(time), !is.na(event))
  
  if (nrow(merged_df) == 0) stop("Merged dataset has 0 rows after removing NAs in clinical data.")
  
  # Match selected genes
  matching_genes <- intersect(selected_genes, colnames(merged_df))
  if (length(matching_genes) == 0) stop("No matching genes found in merged data.")
  
  # Extract coefficients
  if (is.null(dim(lasso_coefs_df))) {
    lasso_coefs_sub <- lasso_coefs_df[matching_genes]
  } else {
    lasso_coefs_sub <- lasso_coefs_df[matching_genes, , drop = FALSE]
    if (ncol(lasso_coefs_sub) == 1) {
      lasso_coefs_sub <- as.vector(lasso_coefs_sub[, 1])
      names(lasso_coefs_sub) <- matching_genes
    } else {
      stop("lasso_coefs_df has unexpected structure with multiple columns.")
    }
  }
  
  # Calculate risk score
  risk_input <- as.matrix(merged_df[, matching_genes, drop = FALSE])
  risk_score <- as.vector(risk_input %*% lasso_coefs_sub)
  if (all(is.na(risk_score))) stop("Risk score calculation returned all NAs.")
  
  merged_df$risk_score <- risk_score
  merged_df <- merged_df[!is.na(merged_df$risk_score), ]
  
  if (nrow(merged_df) == 0) stop("All patients were removed due to NA risk scores.")
  
  # Final steps
  merged_df$risk_group <- ifelse(risk_score > median(risk_score), "High", "Low")
  merged_df$event <- as.numeric(as.character(merged_df$event))
  merged_df$time <- as.numeric(as.character(merged_df$time))
  
  surv_obj <- survival::Surv(merged_df$time, merged_df$event)
  surv_fit <- survival::survfit(surv_obj ~ risk_group, data = merged_df)
  
  return(list(
    data = merged_df,
    surv_fit = surv_fit,
    group = group_name
  ))
}

merge_cnv_rna <- function(cancer_type, cnv_df, rna_df) {
  merged <- merge(rna_df, cnv_df, by = "gene") %>%
    na.omit() %>%
    dplyr::mutate(cancer_type = cancer_type)
  return(merged)
}

# Paths to data
paths <- list(
  LUAD = list(cnv = "TCGA/lung/LUAD/cnv_tumor.RDS", rna = "TCGA/lung/LUAD/rna_counts.RDS"),
  BRCA = list(cnv = "TCGA/BRCA/cnv_tumor.RDS", rna = "TCGA/BRCA/rna_counts.RDS"),
  LUSC = list(cnv = "TCGA/LUSC/cnv_tumor.RDS", rna = "TCGA/LUSC/rna_counts.RDS"),
  LIHC = list(cnv = "TCGA/LIHC/cnv_tumor.RDS", rna = "TCGA/LIHC/rna_counts.RDS")
)

# Run analysis for selected cancer types
results <- lapply(names(paths), function(cancer) {
  message("Processing ", cancer)
  cnv <- process_cnv(cancer, paths[[cancer]]$cnv)
  rna <- process_rna(cancer, paths[[cancer]]$rna, cnv$gene)
  merge_cnv_rna(cancer, cnv, rna)
})

# Combine all results except LUAD for the final plot (LUAD often used separately)
p_all <- do.call(rbind, results[names(results) != "LUAD"])
p_luad <- results[["LUAD"]]

# Plotting
p_all$cnv <- factor(p_all$cnv, levels = c("1", "2", "3", "4", "5"))
colors <- c("1" = "dodgerblue1", "2" = "darkgray", "3" = "green4", "4" = "coral3", "5" = "hotpink3")

ggplot(p_all, aes(x = cnv, y = zscore_mean, color = cnv)) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 16, outlier.size = 2, notch = FALSE) +
  labs(x = "CNV group", y = "mRNA Z-score") +
  facet_wrap(~cancer_type) +
  scale_color_manual(values = colors) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 12, face = "plain", color = "black"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )




