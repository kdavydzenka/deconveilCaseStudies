setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/deconveilCaseStudies/")

pkgs <- c("tidyverse", "openxlsx", "DESeq2", "forcats", "scales")
sapply(pkgs, require, character.only = TRUE)

### Load data ###

purity <- read.xlsx("TCGA/purity_TCGA.xlsx",  startRow = 3)

purity <- purity %>% 
  dplyr::filter(Cancer.type == "LUAD") %>% 
  dplyr::select(Sample.ID, CPE) %>% 
  na.omit() %>% 
  dplyr::rename(sample_id = Sample.ID,
                purity = CPE)


cnv_tumor <- readRDS("TCGA/LUAD/cnv_tumor.RDS")
rna_normal <- readRDS("TCGA/LUAD/rna_normal.RDS")
rna_tumor <- readRDS("TCGA/LUAD/rna_tumor.RDS")

# Gene filtering

low_expression_threshold <- 10

expression_summary <- data.frame(
  Gene = rownames(rna_normal),
  MeanExpression = rowMeans(rna_normal)
)
filtered_genes <- expression_summary %>%
  filter(MeanExpression > low_expression_threshold)

rna_normal <- rna_normal[filtered_genes$Gene, ]
rna_tumor <- rna_tumor[filtered_genes$Gene, ]
cnv_tumor <- cnv_tumor[filtered_genes$Gene, ]

cnv_tumor <- apply(cnv_tumor, 2, function(x) ifelse(x > 10, 10, x))
colnames(cnv_tumor) <- paste0(colnames(cnv_tumor), "-01A")

# Align genes and samples

common_rows <- intersect(rownames(rna_tumor), rownames(cnv_tumor))
counts <- rna_tumor[common_rows, , drop = FALSE]
cn <- cnv_tumor[common_rows, , drop = FALSE]

common_cols <- intersect(colnames(counts), purity$sample_id)
counts <- counts[, common_cols, drop = FALSE]
cn <- cn[, common_cols, drop = FALSE]

purity$purity <- as.numeric(as.character(purity$purity))
purity$purity <- round(purity$purity, 2)

meta <- purity %>%
  mutate(stroma = 1 - purity)


# Compute size factors

coldata <- data.frame(
  row.names = colnames(counts)
)

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = coldata,
  design    = ~ 1
)

dds <- estimateSizeFactors(dds)

sf <- sizeFactors(dds)

sf_df <- tibble(
  sample_id = names(sf),
  sf = as.numeric(sf)
)

# Keep genes with >=30% samples CN>=3 or CN<=1.5, OR mean CN away from diploid.

amp_frac <- rowMeans(cn >= 3, na.rm = TRUE)
del_frac <- rowMeans(cn <= 1.5, na.rm = TRUE)
mean_cn  <- rowMeans(cn, na.rm = TRUE)

genes_keep <- names(mean_cn)[
  (amp_frac >= 0.3) | (del_frac >= 0.3) | (mean_cn > 3.0) | (mean_cn < 1.5)
]

nonzero_n <- rowSums(counts[genes_keep, , drop = FALSE] > 0, na.rm = TRUE)
genes_keep <- genes_keep[nonzero_n >= 10]

counts <- counts[genes_keep, , drop = FALSE]
cn  <- cn[genes_keep, , drop = FALSE]

# Build long joint data (gene x sample rows)

counts_long <- as.data.frame(counts) %>%
  tibble::rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "sample_id", values_to = "expr")

cn_long <- as.data.frame(cn) %>%
  tibble::rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "sample_id", values_to = "copies")

df_long <- counts_long %>%
  inner_join(cn_long, by = c("gene", "sample_id")) %>%
  left_join(meta, by = "sample_id") %>%
  left_join(sf_df, by = "sample_id") %>%
  mutate(
    expr = as.integer(expr),
    # If you have separate tumor-cell CN (cancer_copies), use it here.
    # Otherwise, fall back to copies.
    #cancer_copies = if ("cancer_copies" %in% names(.)) cancer_copies else copies,
    eup_dev_cancer = (copies - 2) / 2,
    eup_equiv_cancer = eup_dev_cancer + 1
  )

write.csv(df_long, file = "TCGA/LUAD/test/stan_joint_long.csv", row.names = FALSE)
