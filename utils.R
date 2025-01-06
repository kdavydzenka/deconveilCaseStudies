clustering_patients_cnv <- function(dataset_name, data_path) {
  if (dataset_name == "LUAD_cnv") {
    cnv_tumor <- readRDS(data_path)
    cnv_tumor <- as.matrix(t(cnv_tumor)) 
    cnv_data_normalized <- scale(cnv_tumor)
    pca_result <- prcomp(cnv_data_normalized, scale. = TRUE)
    pca_data <- as.data.frame(pca_result$x[, 1:10])
    #set.seed(123)  
    kmeans_result <- kmeans(pca_data, centers = 3)
    cnv_tumor <- as.data.frame(cnv_tumor)
    cnv_tumor$Cluster <- kmeans_result$cluster
    plot_clusters <- fviz_cluster(kmeans_result, data = pca_data, geom = "point", stand = FALSE) + 
      ggtitle("PCA of CNV Profiles")+
      theme_classic()
  } else if (dataset_name == "LUSC_cnv") {
    cnv_tumor <- readRDS(data_path)
    cnv_tumor <- as.matrix(t(cnv_tumor)) 
    cnv_data_normalized <- scale(cnv_tumor)
    pca_result <- prcomp(cnv_data_normalized, scale. = TRUE)
    pca_data <- as.data.frame(pca_result$x[, 1:10])
    #set.seed(123)  
    kmeans_result <- kmeans(pca_data, centers = 3)
    cnv_tumor <- as.data.frame(cnv_tumor)
    cnv_tumor$Cluster <- kmeans_result$cluster
    plot_clusters <- fviz_cluster(kmeans_result, data = pca_data, geom = "point", stand = FALSE) + 
      ggtitle("PCA of CNV Profiles")+
      theme_classic()
  } else if (dataset_name == "BRCA_cnv") {
    cnv_tumor <- readRDS(data_path)
    cnv_tumor <- as.matrix(t(cnv_tumor)) 
    cnv_data_normalized <- scale(cnv_tumor)
    pca_result <- prcomp(cnv_data_normalized, scale. = TRUE)
    pca_data <- as.data.frame(pca_result$x[, 1:10])
    #set.seed(123)  
    kmeans_result <- kmeans(pca_data, centers = 5)
    cnv_tumor <- as.data.frame(cnv_tumor)
    cnv_tumor$Cluster <- kmeans_result$cluster
    plot_clusters <- fviz_cluster(kmeans_result, data = pca_data, geom = "point", stand = FALSE) + 
      ggtitle("PCA of CNV Profiles")+
      theme_classic()
  } else if (dataset_name == "LIHC_cnv") {
    cnv_tumor <- readRDS(data_path)
    cnv_tumor <- as.matrix(t(cnv_tumor)) 
    cnv_data_normalized <- scale(cnv_tumor)
    pca_result <- prcomp(cnv_data_normalized, scale. = TRUE)
    pca_data <- as.data.frame(pca_result$x[, 1:10])
    #set.seed(123)  
    kmeans_result <- kmeans(pca_data, centers = 2)
    cnv_tumor <- as.data.frame(cnv_tumor)
    cnv_tumor$Cluster <- kmeans_result$cluster
    plot_clusters <- fviz_cluster(kmeans_result, data = pca_data, geom = "point", stand = FALSE) + 
      ggtitle("PCA of CNV Profiles")+
      theme_classic()
  } else if (dataset_name == "HNSC_cnv") {
    cnv_tumor <- readRDS(data_path)
    cnv_tumor <- as.matrix(t(cnv_tumor)) 
    cnv_data_normalized <- scale(cnv_tumor)
    non_constant_columns <- apply(cnv_data_normalized, 2, function(x) var(x) != 0)
    cnv_data_filtered <- cnv_data_normalized[, non_constant_columns]
    pca_result <- prcomp(cnv_data_filtered, scale. = TRUE)
    pca_data <- as.data.frame(pca_result$x[, 1:10])
    kmeans_result <- kmeans(pca_data, centers = 2)
    cnv_tumor <- as.data.frame(cnv_tumor)
    cnv_tumor$Cluster <- kmeans_result$cluster
    plot_clusters <- fviz_cluster(kmeans_result, data = pca_data, geom = "point", stand = FALSE) + 
      ggtitle("PCA of CNV Profiles")+
      theme_classic()
  } else if (dataset_name == "LUSC_cnv") {
    cnv_tumor <- readRDS(data_path)
    cnv_tumor <- as.matrix(t(cnv_tumor)) 
    cnv_data_normalized <- scale(cnv_tumor)
    non_constant_columns <- apply(cnv_data_normalized, 2, function(x) var(x) != 0)
    cnv_data_filtered <- cnv_data_normalized[, non_constant_columns]
    pca_result <- prcomp(cnv_data_filtered, scale. = TRUE)
    pca_data <- as.data.frame(pca_result$x[, 1:10])
    kmeans_result <- kmeans(pca_data, centers = 3)
    cnv_tumor <- as.data.frame(cnv_tumor)
    cnv_tumor$Cluster <- kmeans_result$cluster
    plot_clusters <- fviz_cluster(kmeans_result, data = pca_data, geom = "point", stand = FALSE) + 
      ggtitle("PCA of CNV Profiles")+
      theme_classic()
  }  else {
    stop("Dataset name not recognized")
  }
  return(list(cnv_tumor, plot_clusters))
}


rna_processing <- function(dataset_name, data_path, cnv_filt) {
  if (dataset_name == "LUAD_rna") {
    rna <- readRDS(data_path)
    colnames(rna) <- gsub(pattern = "\\.", replacement = "-", colnames(rna))
    rna_norm <- rna %>% dplyr::select(1:45)
    rna_tum <- rna %>% dplyr::select(46:90)
    colnames(rna_norm) <- stringr::str_sub(colnames(rna_norm),1,12)
    colnames(rna_tum) <- stringr::str_sub(colnames(rna_tum),1,12)
    cnv_filt <- as.data.frame(cnv_filt)
    rna_tum <- rna_tum[,colnames(rna_tum) %in% colnames(cnv_filt)]
    rna_norm <- rna_norm[,colnames(rna_norm) %in% colnames(cnv_filt)]
    x <- colnames(rna_norm)
    names(rna_norm) <- paste(x,"-11A")
    x <- colnames(rna_tum)
    names(rna_tum) <- paste(x,"-01A")
  } else if (dataset_name == "LUSC_rna") {
    rna <- readRDS(data_path)
    #colnames(rna) <- gsub(pattern = "\\.", replacement = "-", colnames(rna))
    rna_norm <- rna %>% select(1:51)
    rna_tum <- rna %>% select(52:102)
    colnames(rna_norm) <- stringr::str_sub(colnames(rna_norm),1,12)
    colnames(rna_tum) <- stringr::str_sub(colnames(rna_tum),1,12)
    cnv_filt <- as.data.frame(cnv_filt)
    colnames(cnv_filt) <- stringr::str_sub(colnames(cnv_filt),1,12)
    rna_tum <- rna_tum[,colnames(rna_tum) %in% colnames(cnv_filt)]
    rna_norm <- rna_norm[,colnames(rna_norm) %in% colnames(cnv_filt)]
    x <- colnames(rna_norm)
    names(rna_norm) <- paste(x,"-11A")
    x <- colnames(rna_tum)
    names(rna_tum) <- paste(x,"-01A")
  } else if (dataset_name == "BRCA_rna") {
    rna <- readRDS(data_path)
    #colnames(rna) <- gsub(pattern = "\\.", replacement = "-", colnames(rna))
    rna_norm <- rna %>% select(1:110)
    rna_tum <- rna %>% select(111:220)
    colnames(rna_norm) <- stringr::str_sub(colnames(rna_norm),1,12)
    colnames(rna_tum) <- stringr::str_sub(colnames(rna_tum),1,12)
    cnv_filt <- as.data.frame(cnv_filt)
    colnames(cnv_filt) <- stringr::str_sub(colnames(cnv_filt),1,12)
    rna_tum <- rna_tum[,colnames(rna_tum) %in% colnames(cnv_filt)]
    rna_norm <- rna_norm[,colnames(rna_norm) %in% colnames(cnv_filt)]
    x <- colnames(rna_norm)
    names(rna_norm) <- paste(x,"-11A")
    x <- colnames(rna_tum)
    names(rna_tum) <- paste(x,"-01A")
  } else if (dataset_name == "LIHC_rna") {
    rna <- readRDS(data_path)
    #colnames(rna) <- gsub(pattern = "\\.", replacement = "-", colnames(rna))
    rna_norm <- rna %>% select(1:50)
    rna_tum <- rna %>% select(51:100)
    colnames(rna_norm) <- stringr::str_sub(colnames(rna_norm),1,12)
    colnames(rna_tum) <- stringr::str_sub(colnames(rna_tum),1,12)
    cnv_filt <- as.data.frame(cnv_filt)
    colnames(cnv_filt) <- stringr::str_sub(colnames(cnv_filt),1,12)
    rna_tum <- rna_tum[,colnames(rna_tum) %in% colnames(cnv_filt)]
    rna_norm <- rna_norm[,colnames(rna_norm) %in% colnames(cnv_filt)]
    x <- colnames(rna_norm)
    names(rna_norm) <- paste(x,"-11A")
    x <- colnames(rna_tum)
    names(rna_tum) <- paste(x,"-01A")  
  } else if (dataset_name == "HNSC_rna") {
    rna <- readRDS(data_path)
    #colnames(rna) <- gsub(pattern = "\\.", replacement = "-", colnames(rna))
    rna_norm <- rna %>% select(1:21)
    rna_tum <- rna %>% select(22:42)
    colnames(rna_norm) <- stringr::str_sub(colnames(rna_norm),1,12)
    colnames(rna_tum) <- stringr::str_sub(colnames(rna_tum),1,12)
    cnv_tumor <- as.data.frame(cnv_tumor)
    colnames(cnv_tumor) <- stringr::str_sub(colnames(cnv_tumor),1,12)
    rna_tum <- rna_tum[,colnames(rna_tum) %in% colnames(cnv_tumor)]
    rna_norm <- rna_norm[,colnames(rna_norm) %in% colnames(cnv_tumor)]
    x <- colnames(rna_norm)
    names(rna_norm) <- paste(x,"-11A")
    x <- colnames(rna_tum)
    names(rna_tum) <- paste(x,"-01A")  
  } else if (dataset_name == "COLON_rna") {
    rna <- readRDS(data_path)
    #colnames(rna) <- gsub(pattern = "\\.", replacement = "-", colnames(rna))
    rna_norm <- rna %>% select(1:12)
    rna_tum <- rna %>% select(13:24)
    colnames(rna_norm) <- stringr::str_sub(colnames(rna_norm),1,12)
    colnames(rna_tum) <- stringr::str_sub(colnames(rna_tum),1,12)
    cnv_tumor <- as.data.frame(cnv_tumor)
    colnames(cnv_tumor) <- stringr::str_sub(colnames(cnv_tumor),1,12)
    rna_tum <- rna_tum[,colnames(rna_tum) %in% colnames(cnv_tumor)]
    rna_norm <- rna_norm[,colnames(rna_norm) %in% colnames(cnv_tumor)]
    x <- colnames(rna_norm)
    names(rna_norm) <- paste(x,"-11A")
    x <- colnames(rna_tum)
    names(rna_tum) <- paste(x,"-01A")  
  }  else {
    stop("Dataset name not recognized")
  }
  return(list(rna_norm, rna_tum))
}



diagplotRoc <- function(truth, p, sig = 0.05, x = "fpr", y = "tpr", output = "screen",
                        path = NULL, draw = TRUE, line_colors = NULL, line_width = 1.5, 
                        plot_title = NULL, axis_text_size = 1.2, legend_text_size = 1.0, 
                        title_text_size = 1.5, margin = c(5, 5, 5, 5), ...) {
  
  # Validate x and y arguments
  valid_metrics <- c("fpr", "fnr", "tpr", "tnr", "scrx", "sens", "spec")
  if (!(x %in% valid_metrics)) stop("Invalid x-axis metric. Choose from: ", paste(valid_metrics, collapse = ", "))
  if (!(y %in% valid_metrics)) stop("Invalid y-axis metric. Choose from: ", paste(valid_metrics, collapse = ", "))
  
  # Convert p to matrix if needed
  if (is.list(p)) {
    pmat <- do.call("cbind", p)
  } else if (is.data.frame(p)) {
    pmat <- as.matrix(p)
  } else if (is.matrix(p)) {
    pmat <- p
  }
  
  if (is.null(colnames(pmat))) colnames(pmat) <- paste("p", seq_len(ncol(pmat)), sep = "_")
  
  axName <- list(
    tpr = "True Positive Rate",
    tnr = "True Negative Rate",
    fpr = "False Positive Rate",
    fnr = "False Negative Rate",
    scrx = "Ratio of selected",
    scry = "Normalized TP/(FP+FN)",
    sens = "Sensitivity",
    spec = "1 - Specificity"
  )
  
  ROC <- vector("list", ncol(pmat))
  names(ROC) <- colnames(pmat)
  
  # Set line colors
  if (!is.null(line_colors) && length(line_colors) >= ncol(pmat)) {
    colspace <- line_colors[seq_len(ncol(pmat))]
  } else {
    colspaceUniverse <- c("red", "blue", "green", "orange", "darkgrey", "green4",
                          "black", "pink", "brown", "magenta", "yellowgreen", "pink4", "seagreen4", "darkcyan")
    colspace <- colspaceUniverse[seq_len(ncol(pmat))]
  }
  names(colspace) <- colnames(pmat)
  
  eps <- min(pmat[!is.na(pmat) & pmat > 0])
  
  AUC_values <- numeric(ncol(pmat))
  
  for (n in colnames(pmat)) {
    gg <- which(pmat[, n] <= sig)
    psample <- -log10(pmax(pmat[gg, n], eps))
    size <- seq(1, length(gg))
    cuts <- seq(-log10(sig), max(psample), length.out = length(gg))
    local.truth <- truth[gg]
    
    S <- length(size)
    TP <- FP <- FN <- TN <- FPR <- FNR <- TPR <- TNR <- SENS <- SPEC <- SCRX <- SCRY <- numeric(S)
    
    for (i in seq_len(S)) {
      TP[i] <- length(which(psample > cuts[i] & local.truth != 0))
      FP[i] <- length(which(psample > cuts[i] & local.truth == 0))
      FN[i] <- length(which(psample < cuts[i] & local.truth != 0))
      TN[i] <- length(which(psample < cuts[i] & local.truth == 0))
      SCRX[i] <- i / S
      SCRY[i] <- TP[i] / (FN[i] + FP[i])
      FPR[i] <- if (FP[i] + TN[i] == 0) 0 else FP[i] / (FP[i] + TN[i])
      FNR[i] <- FN[i] / (TP[i] + FN[i])
      TPR[i] <- TP[i] / (TP[i] + FN[i])
      TNR[i] <- if (TN[i] + FP[i] == 0) 0 else TN[i] / (TN[i] + FP[i])
      SENS[i] <- TPR[i]
      SPEC[i] <- 1 - TNR[i]
    }
    
    ROC[[n]] <- list(TP = TP, FP = FP, FN = FN, TN = TN, FPR = FPR, FNR = FNR, 
                     TPR = TPR, TNR = TNR, SCRX = SCRX, SCRY = SCRY / max(SCRY), 
                     SENS = SENS, SPEC = SPEC, AUC = NULL)
    
    auc <- 0
    for (i in 2:length(ROC[[n]][[toupper(y)]])) {
      auc <- auc + 0.5 * (ROC[[n]][[toupper(x)]][i] - ROC[[n]][[toupper(x)]][i - 1]) *
        (ROC[[n]][[toupper(y)]][i] + ROC[[n]][[toupper(y)]][i - 1])
    }
    ROC[[n]]$AUC <- abs(auc)
    AUC_values[n] <- abs(auc)
  }
  
  if (draw) {
    if (output == "file" && !is.null(path)) {
      png(filename = path, width = 800, height = 800, res = 100)
    }
    xlim <- c(0, 1)
    ylim <- c(0, 1)
    par(mar = margin, cex.axis = axis_text_size, cex.main = title_text_size, 
        cex.lab = axis_text_size, font.lab = 1, font.axis = 1, pty = "m")
    plot.new()
    plot.window(xlim, ylim)
    axis(1, at = pretty(xlim, 10))
    axis(2, at = pretty(ylim, 10))
    
    for (n in names(ROC)) {
      lines(ROC[[n]][[toupper(x)]], ROC[[n]][[toupper(y)]], col = colspace[n], lwd = line_width, ...)
    }
    
    grid()
    title(main = plot_title, xlab = axName[[x]], ylab = axName[[y]], font.main = 1)
    aucText <- vapply(ROC, function(x) round(x$AUC, digits = 3), numeric(1))
    legend("bottomright", col = colspace, lty = 1, cex = legend_text_size,
           legend = paste(names(ROC), " (AUC = ", aucText, ")", sep = ""))
    
    if (output == "file" && !is.null(path)) {
      dev.off()
    }
  }
  
  return(list(ROC = ROC, AUC_values = AUC_values, truth = truth, sigLevel = sig, xAxis = x, yAxis = y, path = path))
}


auROC <- function(truth, p, sig = 0.05, x = "fpr", y = "tpr") {
  
  # Validate x and y arguments
  valid_metrics <- c("fpr", "fnr", "tpr", "tnr", "scrx", "sens", "spec")
  if (!(x %in% valid_metrics)) stop("Invalid x-axis metric. Choose from: ", paste(valid_metrics, collapse = ", "))
  if (!(y %in% valid_metrics)) stop("Invalid y-axis metric. Choose from: ", paste(valid_metrics, collapse = ", "))
  
  # Convert p to matrix if needed
  if (is.list(p)) {
    pmat <- do.call("cbind", p)
  } else if (is.data.frame(p)) {
    pmat <- as.matrix(p)
  } else if (is.matrix(p)) {
    pmat <- p
  }
  
  if (is.null(colnames(pmat))) colnames(pmat) <- paste("p", seq_len(ncol(pmat)), sep = "_")
  
  ROC <- vector("list", ncol(pmat))
  names(ROC) <- colnames(pmat)
  
  eps <- min(pmat[!is.na(pmat) & pmat > 0])
  
  AUC_values <- numeric(ncol(pmat))
  
  for (n in colnames(pmat)) {
    gg <- which(pmat[, n] <= sig)
    psample <- -log10(pmax(pmat[gg, n], eps))
    size <- seq(1, length(gg))
    cuts <- seq(-log10(sig), max(psample), length.out = length(gg))
    local.truth <- truth[gg]
    
    S <- length(size)
    TP <- FP <- FN <- TN <- FPR <- FNR <- TPR <- TNR <- SENS <- SPEC <- SCRX <- SCRY <- numeric(S)
    
    for (i in seq_len(S)) {
      TP[i] <- length(which(psample > cuts[i] & local.truth != 0))
      FP[i] <- length(which(psample > cuts[i] & local.truth == 0))
      FN[i] <- length(which(psample < cuts[i] & local.truth != 0))
      TN[i] <- length(which(psample < cuts[i] & local.truth == 0))
      SCRX[i] <- i / S
      SCRY[i] <- TP[i] / (FN[i] + FP[i])
      FPR[i] <- if (FP[i] + TN[i] == 0) 0 else FP[i] / (FP[i] + TN[i])
      FNR[i] <- FN[i] / (TP[i] + FN[i])
      TPR[i] <- TP[i] / (TP[i] + FN[i])
      TNR[i] <- if (TN[i] + FP[i] == 0) 0 else TN[i] / (TN[i] + FP[i])
      SENS[i] <- TPR[i]
      SPEC[i] <- 1 - TNR[i]
    }
    
    ROC[[n]] <- list(TP = TP, FP = FP, FN = FN, TN = TN, FPR = FPR, FNR = FNR, 
                     TPR = TPR, TNR = TNR, SCRX = SCRX, SCRY = SCRY / max(SCRY), 
                     SENS = SENS, SPEC = SPEC, AUC = NULL)
    
    auc <- 0
    for (i in 2:length(ROC[[n]][[toupper(y)]])) {
      auc <- auc + 0.5 * (ROC[[n]][[toupper(x)]][i] - ROC[[n]][[toupper(x)]][i - 1]) *
        (ROC[[n]][[toupper(y)]][i] + ROC[[n]][[toupper(y)]][i - 1])
    }
    ROC[[n]]$AUC <- abs(auc)
    AUC_values[n] <- abs(auc)
  }
  
  return(list(ROC = ROC, AUC_values = AUC_values))
}


evaluate_simulation_performance <- function(n_samples, n_genes) {
  
  metrics_df <- data.frame()  
  
  for (replicate in 1:10) {
    # Load true labels
    true_labels_file <- paste0("deconveilCaseStudies/simulations/results/replicates_rna_counts_sim/rna_counts_sim_", n_samples, "_", n_genes, "_brca.rds")
    true_labels <- readRDS(true_labels_file)@variable.annotations$differential.expression
    
    # Load prediction results for the current replicate
    res_pydeseq <- read.csv(paste0("deconveilCaseStudies/simulations/results/replicates_pydeseq/cn_naive/", replicate, "_res_CNnaive_", n_samples, "_", n_genes, ".csv"))
    res_deconveil <- read.csv(paste0("deconveilCaseStudies/simulations/results/replicates_pydeseq/cn_aware/", replicate, "_res_CNaware_", n_samples, "_", n_genes, ".csv"))
    res_edge <- readRDS(paste0("deconveilCaseStudies/simulations/results/replicates_edgeR/cn_naive/", replicate, "_res_CNnaive_", n_samples, "_", n_genes, ".RDS"))
    
    # Process results
    res_pydeseq <- res_pydeseq %>% 
      dplyr::select(X, log2FoldChange, padj) %>% 
      remove_rownames %>% 
      column_to_rownames(var = "X") %>% 
      dplyr::rename(logFC = log2FoldChange)
    
    res_deconveil <- res_deconveil %>% 
      dplyr::select(X, log2FoldChange, padj) %>% 
      remove_rownames %>% 
      column_to_rownames(var = "X") %>% 
      dplyr::rename(logFC = log2FoldChange)
    
    res_edge <- res_edge %>% 
      dplyr::select(logFC, padj) 
    
    # Binary predictions (padj < 0.05)
    predicted_pydeseq <- ifelse(res_pydeseq$padj < 0.05, 1, 0)
    predicted_deconveil <- ifelse(res_deconveil$padj < 0.05, 1, 0)
    predicted_edge <- ifelse(res_edge$padj < 0.05, 1, 0)
    
    # Replace NA values with 0 in predictions
    predicted_edge <- ifelse(is.na(predicted_edge), 0, predicted_edge)
    
    # Function to compute performance metrics
    evaluate_performance <- function(true_labels, predicted_labels) {
      TP <- sum(true_labels == 1 & predicted_labels == 1)
      FP <- sum(true_labels == 0 & predicted_labels == 1)
      TN <- sum(true_labels == 0 & predicted_labels == 0)
      FN <- sum(true_labels == 1 & predicted_labels == 0)
      
      precision <- ifelse((TP + FP) > 0, TP / (TP + FP), NA)
      specificity <- ifelse((TN + FP) > 0, TN / (TN + FP), NA)
      accuracy <- ifelse((TP + TN + FP + FN) > 0, (TP + TN) / (TP + TN + FP + FN), 0)
      
      return(c(Precision = precision, Specificity = specificity, Accuracy = accuracy))
    }
    
    # Calculate performance metrics for each method and add to metrics_df
    metrics_df <- rbind(
      metrics_df,
      data.frame(Method = "PyDESeq2", t(evaluate_performance(true_labels, predicted_pydeseq)), SampleSize = n_samples, Replicate = replicate),
      data.frame(Method = "DeConveil", t(evaluate_performance(true_labels, predicted_deconveil)), SampleSize = n_samples, Replicate = replicate),
      data.frame(Method = "edgeR", t(evaluate_performance(true_labels, predicted_edge)), SampleSize = n_samples, Replicate = replicate)
    )
  }
  
  return(metrics_df)
}


summarize_performance <- function(data) {
  data %>%
    group_by(Method, SampleSize) %>%
    summarize(
      Precision_Mean = mean(Precision),
      Precision_SD = sd(Precision),
      Specificity_Mean = mean(Specificity),
      Specificity_SD = sd(Specificity),
      Accuracy_Mean = mean(Accuracy),
      Accuracy_SD = sd(Accuracy),
      .groups = "drop"
    ) %>%
    pivot_longer(cols = c(Precision_Mean, Precision_SD, Specificity_Mean, Specificity_SD, Accuracy_Mean, Accuracy_SD),
                 names_to = c("Metric", "Stat"),
                 names_sep = "_", values_to = "Value") %>%
    pivot_wider(names_from = Stat, values_from = Value)
}


perform_enrichment_H <- function(gene_df, gene_group) {
  my_symbols <- gene_df$geneID
  gene_list <- AnnotationDbi::select(org.Hs.eg.db,
                                     keys = my_symbols,
                                     columns = c("ENTREZID", "SYMBOL"),
                                     keytype = "SYMBOL")
  gene_l <- as.vector(gene_list$ENTREZID)
  m_hallmark <- msigdbr::msigdbr(species = "Homo sapiens", category = "H") 
  msig_H <- msigdbr(species = "Homo sapiens", category = "H") %>% dplyr::select(gs_name, entrez_gene)
  msig_H <- enricher(gene_l, minGSSize = 10, 
                     maxGSSize = 350,
                     pvalueCutoff = 0.05, 
                     pAdjustMethod = "BH", 
                     TERM2GENE = msig_H)
  
  msig_H@result %>% mutate(gene_group = gene_group)
}


perform_enrichment_GO <- function(gene_df, gene_group) {
  my_symbols <- gene_df$geneID
  gene_list <- AnnotationDbi::select(org.Hs.eg.db,
                                     keys = my_symbols,
                                     columns = c("ENTREZID", "SYMBOL"),
                                     keytype = "SYMBOL")
  gene_l <- as.vector(gene_list$ENTREZID)
  oraGO <- enrichGO(gene = gene_l, 
                    ont = "BP", 
                    OrgDb = org.Hs.eg.db, 
                    keyType = "ENTREZID", 
                    pvalueCutoff = 0.05, 
                    minGSSize = 10, 
                    maxGSSize = 350)
  
  oraGO@result %>% mutate(gene_group = gene_group)
}


convert_entrez_to_symbol <- function(entrez_ids) {
  entrez_list <- unlist(strsplit(entrez_ids, split = "/"))
  gene_symbols <- mapIds(org.Hs.eg.db, entrez_list, "SYMBOL", "ENTREZID")
  return(paste(gene_symbols, collapse = ", "))
}
