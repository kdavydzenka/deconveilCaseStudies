### Survival analysis of TCGA-BRCA dataset ###

setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/")
pkgs <- c("tidyverse", "survival", "glmnet", "survminer", "survcomp", "DESeq2", "forestplot", "caret", "randomForestSRC")
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
rna_log_norm <- rna_log_norm * cnv_tumor

# Reorder patients indexing
idx <- match(clinical_data$patientID, colnames(rna_log_norm))
rna_log_norm <- rna_log_norm[,idx]
rna_log_norm <- t(rna_log_norm)
clinical_data <- clinical_data %>% remove_rownames %>% column_to_rownames(var="patientID") 
data <- cbind(clinical_data, rna_log_norm)

# Split the data into training (70%) and test (30%) sets
train_index <- createDataPartition(data$event, p = 0.7, list = FALSE)
train_data <- data[train_index, ]
test_data  <- data[-train_index, ]

train_data <- data[train_index, ]
test_data  <- data[-train_index, ]

rna_train <- train_data %>% select(-time, -event)
clinical_train <- train_data %>% select(time, event)

rna_test <- test_data %>% select(-time, -event)
clinical_test <- test_data %>% select(time, event)

# Random Survival Forest
#rsf_model <- rfsrc(Surv(time, event) ~ ., data = train_data, ntree = 100)
#var.select(rsf_model)
#rsf_predictions <- predict(rsf_model, newdata = test_data)

#rsf_c_index <- rsf_model$err.rate
#print(paste("Concordance Index (RSF):", 1 - rsf_c_index))

### Cox model ###
# Initial Cox model on Training set
surv_object <- survival::Surv(time = clinical_train$time, event = clinical_train$event)
cox_results <- data.frame(Gene = character(), p.value = numeric(), HR = numeric(), CI_lower = numeric(), CI_upper = numeric())

for (gene in colnames(rna_train)) {
  cox_model <- coxph(surv_object ~ rna_train[, gene], data = train_data)
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

# LASSO for further gene selection
X <- as.matrix(rna_train[, significant_genes$Gene])  
y <- surv_object  
lasso_model <- glmnet(X, y, family = "cox", alpha = 1)
cv_model <- cv.glmnet(X, y, family = "cox", type.measure="C", alpha = 1.0, )
optimal_lambda <- cv_model$lambda.min

# Get coefficients for the best lambda (selected genes)
lasso_coefs <- coef(lasso_model, s = optimal_lambda)
lasso_coefs_df <- as.data.frame(as.matrix(lasso_coefs))
selected_genes <- rownames(lasso_coefs_df)[lasso_coefs_df$`1` != 0]

# Train final Cox model with selected genes
cox_selected <- coxph(Surv(clinical_train$time, clinical_train$event) ~ ., data = rna_train[, selected_genes])

# Validate the model on test set
cox_predictions <- predict(cox_selected, newdata = rna_test[, selected_genes], type = "lp")
c_index <- concordance(Surv(clinical_test$time, clinical_test$event) ~ cox_predictions)$concordance
cat("Concordance Index with LASSO-selected genes: ", c_index, "\n")



## Cross-validation to access the model stability ##

set.seed(123)  

k <- 10         

# Create k-fold partitions
folds <- createFolds(data$event, k = k, list = TRUE)

# Initialize a vector to store the Concordance Index for each fold
c_indices <- numeric(k)

# Cross-validation loop
for (i in 1:k) {
  # Split into training and validation folds
  train_index <- unlist(folds[-i])
  test_index <- unlist(folds[i])
  
  train_data <- data[train_index, ]
  test_data <- data[test_index, ]
  
  # Split data into RNA and clinical for training and validation
  rna_train <- train_data %>% select(-time, -event)
  clinical_train <- train_data %>% select(time, event)
  rna_test <- test_data %>% select(-time, -event)
  clinical_test <- test_data %>% select(time, event)
  
  # Initial Cox model on training set for gene selection
  surv_object <- Surv(time = clinical_train$time, event = clinical_train$event)
  cox_results <- data.frame(Gene = character(), p.value = numeric(), HR = numeric(), CI_lower = numeric(), CI_upper = numeric())
  
  for (gene in colnames(rna_train)) {
    cox_model <- coxph(surv_object ~ rna_train[, gene], data = train_data)
    summary_cox <- summary(cox_model)
    
    cox_results <- rbind(cox_results, data.frame(
      Gene = gene,
      p.value = summary_cox$coefficients[,"Pr(>|z|)"],
      HR = summary_cox$coefficients[,"exp(coef)"],
      CI_lower = summary_cox$conf.int[,"lower .95"],
      CI_upper = summary_cox$conf.int[,"upper .95"]
    ))
  }
  
  # Filter significant genes based on p-value threshold
  significant_genes <- cox_results %>% filter(p.value < 0.05)
  
  # LASSO for further gene selection on significant genes
  X <- as.matrix(rna_train[, significant_genes$Gene])
  y <- surv_object
  lasso_model <- glmnet(X, y, family = "cox", alpha = 1, maxit = 500000, standardize = TRUE, lambda.min.ratio = 0.01)
  cv_model <- cv.glmnet(X, y, family = "cox", type.measure = "C", alpha = 1.0, maxit = 500000, standardize = TRUE, lambda.min.ratio = 0.01)
  optimal_lambda <- cv_model$lambda.min
  
  # Extract selected genes based on optimal lambda
  lasso_coefs <- coef(lasso_model, s = optimal_lambda)
  lasso_coefs_df <- as.data.frame(as.matrix(lasso_coefs))
  selected_genes <- rownames(lasso_coefs_df)[lasso_coefs_df$`1` != 0]
  
  # Ensure all selected genes are in validation set, set missing to zero
  rna_train_lasso <- rna_train[, selected_genes, drop = FALSE]
  rna_test_lasso <- rna_test[, selected_genes, drop = FALSE]
  rna_test_lasso[is.na(rna_test_lasso)] <- 0
  
  # Train final Cox model with LASSO-selected genes on training set
  final_cox_model <- coxph(Surv(clinical_train$time, clinical_train$event) ~ ., data = data.frame(clinical_train, rna_train_lasso))
  
  # Validate the model on the validation set
  cox_predictions <- predict(final_cox_model, newdata = rna_test_lasso, type = "lp")
  c_index <- concordance(Surv(clinical_test$time, clinical_test$event) ~ cox_predictions)$concordance
  
  # Store the Concordance Index for the current fold
  c_indices[i] <- c_index
  cat("Fold", i, "C-index:", c_index, "\n")
}

# Calculate the average Concordance Index across all folds
mean_c_index <- mean(c_indices, na.rm = TRUE)
cat("Average Concordance Index from Cross-Validation: ", mean_c_index, "\n")

c_indices_10 <- c_indices_10 %>% as.data.frame()
saveRDS(c_indices_10, file = "CN-aware-DGE/case_studies/c_indices_10.RDS")


# Fit Random Survival Forest model (account for nonlinear relationships between features and survival)
rsf_model <- rfsrc(Surv(time, event) ~ ., data = train_data, ntree = 100)
var.select(rsf_model)

rsf_predictions <- predict(rsf_model, newdata = test_data)

rsf_c_index <- rsf_model$err.rate
print(paste("Concordance Index (RSF):", 1 - rsf_c_index))


# Box plot 
ci_cox_aware <- data.frame(C_index = c(0.70, 0.56, 0.62, 0.50, 0.57, 0.59, 0.63, 0.60, 0.52, 0.54))
median(ci_cox_aware$C_index)

ci_cox_aware <- ci_cox_aware %>% 
  dplyr::mutate(method = "CN-aware") %>% 
  dplyr::mutate(model = "Cox")

ci_cox_naive <- data.frame(C_index = c(0.38, 0.57, 0.40, 0.55, 0.41, 0.30, 0.32, 0.49, 0.43, 0.34))  
median(ci_cox_naive$C_index)
ci_cox_naive <- ci_cox_naive %>% 
  dplyr::mutate(method = "CN-naive") %>% 
  dplyr::mutate(model = "Cox")

ci_rsf_aware <- data.frame(C_index = c(0.46, 0.57, 0.56, 0.48, 0.62, 0.52, 0.45, 0.54, 0.55, 0.51))
median(ci_rsf_aware$C_index)
ci_rsf_aware <- ci_rsf_aware %>% 
  dplyr::mutate(method = "CN-aware") %>% 
  dplyr::mutate(model = "RSF")

ci_rsf_naive <- data.frame(C_index = c(0.58, 0.50, 0.52, 0.60, 0.56, 0.50, 0.50, 0.44, 0.52, 0.54))
median(ci_rsf_naive$C_index)
ci_rsf_naive <- ci_rsf_naive %>% dplyr::mutate(method = "CN-naive") %>% 
  dplyr::mutate(model = "RSF")

ci_data <- rbind(ci_cox_aware, ci_cox_naive, ci_rsf_aware, ci_rsf_naive)

method_colors <- c("#ADB17DFF", "#D49464FF")

bxp <- ggplot(ci_data, aes(x = model, y = C_index, fill = method)) + 
  geom_boxplot(position = position_dodge())+
  labs(x="model", y = "Concordance Index", title = "")+
  labs(fill = "method")+
  theme_bw()+
  ggplot2::theme(legend.position = 'bottom')+
  #facet_wrap(~cancer_type)+
  scale_fill_manual(values = method_colors) +
  font("xy.text", size = 14, color = "black", face = "plain")+
  font("title", size = 16, color = "black")+
  font("xlab", size = 14)+
  font("ylab", size = 14)+
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.text = element_text(size = 16, color = "black"),  
    legend.title = element_text(size = 18, color = "black"),
    axis.title.x = element_text(size = 18, face = "plain", color = "black"), 
    axis.title.y = element_text(size = 18, face = "plain", color = "black"),  
    axis.text.x = element_text(size = 18, color = "black", face = "plain"),  
    axis.text.y = element_text(size = 18, color = "black", face = "plain"))
bxp
