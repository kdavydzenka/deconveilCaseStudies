setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/TCGA/LUSC")

#setwd("/orfeo/LTS/CDSLab/LT_storage/kdavydzenka/lusc_test/")

pkgs <- c("dplyr", "tidyr", "data.table", "GenomicRanges", "BSgenome.Hsapiens.UCSC.hg38",
          "TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db", "miRNAtap", "TFBSTools",
          "JASPAR2022", "Biostrings", "multiMiR", "rtracklayer", "BiocParallel", "stringr",
          "DESeq2")
sapply(pkgs, require, character.only = TRUE)


### Promoter-binding TF activity ###

# TF expression levels (TFs x samples)

rna_tumor <- readRDS("rna_tumor.RDS") # extract from mRNA
colnames(rna_tumor) <- gsub("-01A", "", colnames(rna_tumor))
res_deg <- readRDS("res_CNaware_deg.RDS")
degs <- res_deg$geneID

# TF binding affinities

genome <- BSgenome.Hsapiens.UCSC.hg38
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
promoters <- promoters(genes(txdb), upstream=500, downstream=1500)
promoters_trimmed <- trim(promoters)

# Convert Entrez IDs to Gene Symbols
entrez_ids <- mcols(promoters_trimmed)$gene_id
gene_symbols <- mapIds(org.Hs.eg.db, 
                      keys=as.character(entrez_ids), 
                      column="SYMBOL", 
                      keytype="ENTREZID", 
                      multiVals="first")

geneIDs <- as.data.frame(gene_symbols)
mcols(promoters_trimmed)$gene_symbol <- geneIDs$gene_symbols

promoter_sequences <- getSeq(genome, promoters_trimmed)
names(promoter_sequences) <- mcols(promoters_trimmed)$gene_symbol
promoter_sequences <- promoter_sequences[names(promoter_sequences) %in% degs]


# Get TF Binding Motifs from JASPAR | Retrieve position weight matrices (PWMs) for TFs 
opts <- list(
  species = 9606,  # Human 
  collection = "CORE"  # Use JASPAR CORE collection
)
pwm_list <- getMatrixSet(JASPAR2022, opts)
pwm_weight_list <- lapply(pwm_list, toPWM)

tf_gene_names <- sapply(pwm_weight_list, function(pwm) pwm@name)
names(pwm_weight_list) <- tf_gene_names

# Scan promoters for TF binding sites

bp_param <- MulticoreParam(workers = 8)
scan_TFBS <- bplapply(pwm_weight_list, function(pwm) {
  scoreMatrix <- searchSeq(pwm, promoter_sequences, min.score = "80%")
}, BPPARAM = bp_param)


# Filter TF binding scores

scan_TFBS <- readRDS("scan_TFBS.RDS")

binding_scores <- lapply(scan_TFBS, function(site_set) {
  site_scores <- lapply(site_set, function(site) {
    score(site) 
  })
  return(site_scores)
})

filtered_scores <- lapply(binding_scores, function(tf_scores) {
  # For each TF, remove genes with no binding scores
  valid_scores <- tf_scores[sapply(tf_scores, function(x) length(x) > 0)]
  return(valid_scores)
})

mean_filtered_scores <- lapply(filtered_scores, function(tf_scores) {
  mean_scores <- sapply(tf_scores, function(scores) {
    if (length(scores) > 1) {
      mean(scores)  
    } else {
      scores
    }
  })
  return(mean_scores)
})


filter_binding_scores <- function(binding_scores) {
  filtered_scores <- lapply(binding_scores, function(tf_scores) {
    tf_scores <- tf_scores[tf_scores >= 5]  # Removing scores less than 4
    tf_scores <- tf_scores[tf_scores > 0]   # Removing negative scores
    return(tf_scores)
  })
  return(filtered_scores)
}

filt_scores <- filter_binding_scores(mean_filtered_scores)
filt_scores <- filt_scores[sapply(filt_scores, length) > 0]

df_long <- do.call(rbind, lapply(names(filt_scores), function(tf) {
  data.frame(TF = tf, Gene = names(filt_scores[[tf]]), Score = filt_scores[[tf]])
}))


# Results data reshaping

filt_scores <- readRDS("tf_score.RDS")

filt_scores <- filt_scores %>%
  group_by(Gene, TF) %>%
  summarise(Score = mean(Score, na.rm = TRUE), .groups = "drop") %>% 
  filter(Score >= 5.0)

filt_scores$Score <- round(filt_scores$Score, 2)

binding_matrix <- filt_scores %>%
  pivot_wider(names_from = TF, values_from = Score, values_fill = 0)

binding_matrix <- binding_matrix %>% column_to_rownames(var = "Gene")

# Function to compute p-values for each TF binding score
#filtered_sites <- lapply(pvalues, function(p) p[p < 0.01]) 

common_TFs <- colnames(binding_matrix)[colnames(binding_matrix) %in% rownames(rna_tumor)]
split_TFs <- strsplit(colnames(binding_matrix), "::")
#split_TFs <- split_TFs[-1]

# Extract TF expression data

composite_tf_expression <- sapply(split_TFs, function(tf_pair) {
  #print(paste("Processing pair:", paste(tf_pair, collapse = "::")))
  available_TFs <- tf_pair[tf_pair %in% rownames(rna_tumor)]
  print(paste("Available TFs found:", paste(available_TFs, collapse = ", ")))
  if (length(available_TFs) == 2) {
    expression_values <- rna_tumor[available_TFs, , drop = FALSE]
    mean_expression_sample_wise <- apply(expression_values, 2, mean, na.rm = TRUE)
    
    return(mean_expression_sample_wise)  # Return the mean expression for each sample
    
  } else if (length(available_TFs) == 1) {
    return(rna_tumor[available_TFs, , drop = FALSE])
    
  } else {
    return(NA)
  }
})


#binding_matrix <- binding_matrix %>% column_to_rownames(var = "Gene")
names(composite_tf_expression) <- colnames(binding_matrix)
final_exp_matrix <- do.call(rbind, composite_tf_expression)
final_exp_matrix <- na.omit(final_exp_matrix)
final_exp_matrix <- round(final_exp_matrix)

# Combine binding_matrix and final_exp_matrix to calculate sum(pTF*mu) for a given gene

binding_matrix <- binding_matrix[, rownames(final_exp_matrix), drop = FALSE]
final_exp_matrix <- final_exp_matrix %>% as.matrix() %>% DESeq2::varianceStabilizingTransformation()

#final_exp_matrix <- final_exp_matrix %>%
  #as.data.frame() %>% 
  #mutate(mean_exp = rowMeans(., na.rm = TRUE)) %>% 
  #dplyr::select(mean_exp)

common_TFs <- intersect(colnames(binding_matrix), rownames(final_exp_matrix))
binding_matrix <- binding_matrix[, common_TFs, drop = FALSE]
#exp_levels <- final_exp_matrix[common_TFs, "mean_exp", drop = FALSE]
exp_levels <- final_exp_matrix[common_TFs,]


binding_matrix <- as.matrix(binding_matrix)
exp_levels <- as.matrix(exp_levels)

pTF_features <- binding_matrix %*% exp_levels

#pTF_values <- as.vector(pTF_values)
#names(pTF_values) <- rownames(binding_matrix)

saveRDS(pTF_features, file = "pTF_features.RDS")



### miRNA-mRNA Target Mapping ###

mirna <- readRDS("mirna.RDS")
mirna <- mirna[, !grepl("^(reads_per_million|cross-mapped)", colnames(mirna))]
colnames(mirna) <- gsub(".*(TCGA-\\d{2}-\\d{4}).*", "\\1", colnames(mirna))

mirna_targets <- read.table("targetscan.txt", header = TRUE, sep = "\t")

mirna_targets <- mirna_targets %>%
  dplyr::rename(context_score = context...score) %>%
  dplyr::mutate(context_score = as.numeric(context_score)) %>% 
  dplyr::rename(geneID = Gene.Symbol) %>% 
  dplyr::select(geneID, miRNA, context_score)

mirna_targets$miRNA <- sub("-5p|-3p|-3p.1|-3p.2$", "", mirna_targets$miRNA)
mirna_targets$miRNA <- sub("miR", "mir", mirna_targets$miRNA)

mirna_targets <- mirna_targets[mirna_targets$miRNA %in% mirna$miRNA_ID, ]

# Compute binding affinity (context score * 10)
mirna_targets$context_score <- abs(mirna_targets$context_score) * 10

mirna_targets_avg <- mirna_targets %>%
  group_by(geneID, miRNA) %>%
  summarise(mean_context_score = mean(context_score))

mirna_filtered <- mirna[mirna$miRNA_ID %in% mirna_targets_avg$miRNA, ]
rownames(mirna_filtered) <- NULL
mirna_filtered <- mirna_filtered %>% as.data.frame() %>% column_to_rownames(var = "miRNA_ID")
mirna_filtered <- mirna_filtered %>% as.matrix() %>% DESeq2::varianceStabilizingTransformation()

#mirna_exp <- mirna_filtered %>%
  #as.data.frame() %>% 
  #mutate(mean_exp = rowMeans(., na.rm = TRUE)) %>% 
  #dplyr::select(mean_exp)

mirna_binding_score <- mirna_targets_avg %>%
  pivot_wider(
    names_from = miRNA,
    values_from = mean_context_score,
    values_fill = list(mean_context_score = 0)
  )

missing_genes <- setdiff(rownames(binding_matrix), mirna_binding_score$geneID)

if (length(missing_genes) > 0) {
  missing_mirna_binding_score <- data.frame(
    geneID = missing_genes,
    miRNA = rep(NA, length(missing_genes)),  
    mean_context_score = rep(0, length(missing_genes))  
  )
  
  mirna_binding_score_updated <- bind_rows(mirna_binding_score, missing_mirna_binding_score)
  
  mirna_binding_score_updated <- mirna_binding_score_updated %>%
    arrange(geneID)
  
} else {
  mirna_binding_score_updated <- mirna_binding_score
  
}

mirna_binding_score_updated <- mirna_binding_score_updated %>%
  mutate_all(~replace(., is.na(.), 0))

common_genes <- intersect(mirna_binding_score_updated$geneID, rownames(binding_matrix))
mirna_binding_score_filtered <- mirna_binding_score_updated %>%
  filter(geneID %in% common_genes)

mirna_binding <- mirna_binding_score_filtered %>% column_to_rownames(var = "geneID")

# Calculate miRNA feature for a given gene
mirna_exp <- mirna_filtered 
common_mirnas <- intersect(colnames(mirna_binding), rownames(mirna_exp))
mirna_binding <- mirna_binding[, common_mirnas, drop = FALSE]
mirna_exp<- mirna_exp[common_mirnas, , drop = FALSE]

mirna_binding <- as.matrix(mirna_binding)
mirna_exp <- as.matrix(mirna_exp)

miRNA_features <- mirna_binding %*% mirna_exp

#miRNA_values <- as.vector(miRNA_values)
#names(miRNA_values) <- rownames(mirna_binding)

saveRDS(miRNA_features, file = "miRNA_features.RDS")




# Download Enhancer-Promoter interactions from FANTOM5

#url <- "https://fantom.gsc.riken.jp/5/datafiles/latest/extra/Enhancers/human_permissive_enhancers_phase_1_and_2.bed.gz"
#download.file(url, destfile="fantom5_enhancers.bed.gz")

#enhancers <- fread("fantom5_enhancers.bed.gz", header=FALSE)
#colnames(enhancers) <- c("chr", "start", "end", "enhancer_id", "score", "strand",
#"thickStart", "thickEnd", "itemRgb", "blockCount", "expression", "activity_score")

#enhancers$strand <- ifelse(enhancers$strand %in% c(".", NA, ""), "*", enhancers$strand)
#enhancer_gr <- GRanges(seqnames=enhancers$chr,
#ranges=IRanges(start=enhancers$start, end=enhancers$end),
#strand=enhancers$strand)


