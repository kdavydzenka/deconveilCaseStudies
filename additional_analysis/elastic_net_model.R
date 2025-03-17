setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/TCGA/LUSC")

pkgs <- c("dplyr", "tidyr", "caret", "glmnet", "DESeq2")
sapply(pkgs, require, character.only = TRUE)

### Build Elastic net regression model ###

# Load data
cn <- readRDS("cnv_tumor.RDS")
colnames(cn) <- gsub("-11A", "", colnames(cn))

rna <- readRDS("rna_tumor.RDS") 
colnames(rna) <- gsub("-01A", "", colnames(rna))

res_deg <- readRDS("res_CNaware_deg.RDS")
degs <- res_deg$geneID

meth <- readRDS("dna_meth_processed.RDS") %>% na.omit()
colnames(meth) <- substr(colnames(meth), 1, 12)

meth <- meth[!duplicated(meth$Gene_Name), ] %>% 
  remove_rownames %>% 
  column_to_rownames(var="Gene_Name") 

pTF <- readRDS("pTF_features.RDS")

miRNA <- readRDS("miRNA_features.RDS")

# Align gene in all data modalities

# CN
cn <- cn[rownames(cn) %in% degs ,]
#cn <- cn %>%
  #as.data.frame() %>% 
  #mutate(mean_cn = rowMeans(., na.rm = TRUE)) %>% 
  #dplyr::select(mean_cn)

# mRNA
rna <- rna[rownames(rna) %in% degs ,]
rna <- rna %>% as.matrix() %>% DESeq2::varianceStabilizingTransformation()
#rna <- rna %>%
  #as.data.frame() %>% 
  #mutate(mean_rna = rowMeans(., na.rm = TRUE)) %>% 
  #dplyr::select(mean_rna)

# DNA-meth
meth <- meth[rownames(meth) %in% degs ,]
#meth <- meth %>%
  #as.data.frame() %>% 
  #mutate(mean_meth = rowMeans(., na.rm = TRUE)) %>% 
  #dplyr::select(mean_meth)

missing_genes <- setdiff(degs, rownames(meth))
missing_data <- data.frame(matrix(0, nrow = length(missing_genes), ncol = ncol(meth)))
rownames(missing_data) <- missing_genes
colnames(missing_data) <- colnames(meth)
meth <- rbind(meth, missing_data)

# pTF
missing_genes <- setdiff(degs, rownames(pTF))
missing_data <- data.frame(matrix(0, nrow = length(missing_genes), ncol = ncol(pTF)))
rownames(missing_data) <- missing_genes
colnames(missing_data) <- colnames(pTF)
pTF <- rbind(pTF, missing_data)

# miRNA
missing_genes <- setdiff(degs, rownames(miRNA))
missing_data <- data.frame(matrix(0, nrow = length(missing_genes), ncol = ncol(miRNA)))
rownames(missing_data) <- missing_genes
colnames(missing_data) <- colnames(miRNA)
miRNA <- rbind(miRNA, missing_data)

# Align colnames
colnames(cn) <- trimws(colnames(cn))
colnames(rna) <- trimws(colnames(rna))
colnames(meth) <- trimws(colnames(meth))
colnames(pTF) <- trimws(colnames(pTF))
colnames(miRNA) <- trimws(colnames(miRNA))

common_colnames <- Reduce(intersect, list(colnames(cn), colnames(rna), colnames(meth), colnames(pTF), colnames(miRNA)))

cn <- cn[, common_colnames, drop = FALSE]
rna <- rna[, common_colnames, drop = FALSE]
meth <- meth[, common_colnames, drop = FALSE]
pTF <- pTF[, common_colnames, drop = FALSE]
miRNA <- miRNA[, common_colnames, drop = FALSE]

# Align rownames
common_genes <- Reduce(intersect, list(rownames(cn), rownames(rna), rownames(meth), rownames(pTF), rownames(miRNA)))

cn <- cn[common_genes, , drop = FALSE]
rna <- rna[common_genes, , drop = FALSE]
meth <- meth[common_genes, , drop = FALSE]
pTF <- pTF[common_genes, , drop = FALSE]
miRNA <- miRNA[common_genes, , drop = FALSE]

#cn_val <- as.vector(cn$mean_cn)
#names(cn_val) <- rownames(cn)

#meth_val <- as.vector(meth$mean_meth)
#names(meth_val) <- rownames(meth)

#pTF_val <- as.vector(pTF$mean_exp)
#names(pTF_val) <- rownames(pTF)

#miRNA_val <- as.vector(miRNA$mean_exp)
#names(miRNA_val) <- rownames(miRNA)

#rna_val <- as.vector(rna$mean_rna)
#names(rna_val) <- rownames(rna)

# Prepare the dataset
cn <- t(cn)
meth <- t(meth)
pTF <- t(pTF)
miRNA <- t(miRNA)
rna <- t(rna)

common_genes <- Reduce(intersect, list(colnames(cn), colnames(meth), colnames(pTF), colnames(miRNA), colnames(rna)))

# Create a list to store predictor matrices for each gene
predictor_list <- list()

for (gene in common_genes) {
  predictor_list[[gene]] <- cbind(
    CN = cn[, gene],
    Meth = meth[, gene],
    pTF = pTF[, gene],
    miRNA = miRNA[, gene]
  )
}

 
# Train Elastic Net Model

models <- list()

alpha <- 0.5  # Mixing parameter between Lasso and Ridge (0 = Ridge, 1 = Lasso)

# Loop over each gene in the predictor_list

for (gene_name in names(predictor_list)) {
  X <- as.matrix(predictor_list[[gene_name]]) 
  y <- rna[, gene_name]  
  
  complete_cases <- complete.cases(X, y)
  X <- X[complete_cases, ]
  y <- y[complete_cases]
  
  cv_model <- cv.glmnet(X, y, alpha = alpha, nfolds = 10)
  lambda_values <- cv_model$lambda
  min_lambda <- min(lambda_values, na.rm = TRUE)
  
  model <- glmnet(X, y, alpha = alpha, lambda = min_lambda)
  
  models[[gene_name]] <- model
  
  # Optionally, perform cross-validation to select the optimal lambda (regularization strength)
  #cv_model <- cv.glmnet(X, y, alpha = alpha)
  
  #print(paste("Gene:", gene_name, "- Best Lambda:", cv_model$lambda.min))
  
  # Optionally, make predictions (if needed)
  # y_pred <- predict(model, newx = X_test, s = cv_model$lambda.min)  # using the best lambda
}

saveRDS(models, file = "model_results.RDS")


# Extract beta coefficients

beta_coefficients_df <- data.frame()

for (gene_name in names(models)) {
  beta_coeffs <- models[[gene_name]][["beta"]]
  beta_coeffs_df <- as.data.frame(as.matrix(t(beta_coeffs)))
  beta_coeffs_df$Gene <- gene_name
  colnames(beta_coeffs_df) <- c("CN", "Meth", "pTF", "miRNA", "Gene")
  beta_coefficients_df <- rbind(beta_coefficients_df, beta_coeffs_df)
}

beta_coefficients <- beta_coefficients_df %>%
  dplyr::select(Gene, everything())

saveRDS(beta_coefficients, file = "beta_coeff.RDS")

