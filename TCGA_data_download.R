### Download data from GDC portal ###

rm(list=ls())
setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/TCGA/")
pkgs <- c("tidyverse", "TCGAbiolinks", "SummarizedExperiment", "sesameData", "sesame")
sapply(pkgs, require, character.only = TRUE)

# Define project
project <- "TCGA-BRCA"

## Build queries ##

query_cnv <- GDCquery(
  project = project,
  data.category = "Copy Number Variation",
  data.type = "Gene Level Copy Number",
  sample.type = "Primary Tumor"
)

query_rna_tumor <- GDCquery(
  project = project,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  sample.type = "Primary Tumor",
  access = "open"
)

query_rna_normal <- GDCquery(
  project = project,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  sample.type = "Solid Tissue Normal",
  access = "open"
)


#patientID_brca <- colnames(rna_tumor)
#patientID_brca <- gsub(" -01A", "", patientID_brca)

#query_mirna <- TCGAbiolinks::GDCquery(
  #project = "TCGA-BRCA", 
  #experimental.strategy = "miRNA-Seq",
  #data.category = "Transcriptome Profiling", 
  #data.type = "miRNA Expression Quantification",
  #sample.type = "Primary Tumor",
  #barcode = patientID_brca)


#query_met <- GDCquery(
  #project= "TCGA-BRCA", 
  #data.category = "DNA Methylation", 
  #data.type = "Methylation Beta Value",
  #platform = "Illumina Human Methylation 450", 
  #sample.type = "Primary Tumor",
  #barcode = patientID_brca
#)


## Download data ##

GDCdownload(query_TCGA_cnv)
GDCdownload(query_TCGA_rna_tumor)
GDCdownload(query_TCGA_rna_normal)
#GDCdownload(query_mirna)
#GDCdownload(query_met)

#mirna <- GDCprepare(query = query_mirna)
#dna_met <- GDCprepare(query_met)
#beta <- dna_met@assays@data@listData[[1]]

#saveRDS(dna_met, file = "BRCA/dna_meth.RDS")
#saveRDS(mirna, file = "BRCA/mirna.RDS")

## Prepare data ##

prepare_expression <- function(query) {
  se <- GDCprepare(query, summarizedExperiment = TRUE)
  expr <- assay(se, "unstranded") %>% as.data.frame()
  gene_ids <- rowData(se)$gene_name
  expr <- cbind(GeneID = gene_ids, expr)
  expr <- expr[!duplicated(expr$GeneID), ] %>%
    column_to_rownames("GeneID") %>%
    na.omit()
  colnames(expr) <- substr(colnames(expr), 1, 12)
  return(expr)
}

prepare_cnv <- function(query) {
  se <- GDCprepare(query, summarizedExperiment = TRUE)
  cnv <- assay(se, "copy_number") %>% as.data.frame()
  gene_ids <- rowData(se)$gene_name
  cnv <- cbind(GeneID = gene_ids, cnv)
  cnv <- cnv[!duplicated(cnv$GeneID), ] %>%
    column_to_rownames("GeneID") %>%
    na.omit()
  colnames(cnv) <- substr(colnames(cnv), 1, 12)
  return(cnv)
}

rna_tum <- prepare_expression(query_rna_tumor)
rna_norm <- prepare_expression(query_rna_normal)
cnv_tum <- prepare_cnv(query_cnv)

## Match samples and genes ##
# Match sample columns
common_samples <- Reduce(intersect, list(colnames(rna_tum), colnames(rna_norm), colnames(cnv_tum)))
rna_tum <- rna_tum[, common_samples]
rna_norm <- rna_norm[, common_samples]
cnv_tum <- cnv_tum[, common_samples]

# Match gene rows
common_genes <- Reduce(intersect, list(rownames(rna_tum), rownames(rna_norm), rownames(cnv_tum)))
rna_tum <- rna_tum[common_genes, ]
rna_norm <- rna_norm[common_genes, ]
cnv_tum <- cnv_tum[common_genes, ]

# Rename barcodes for clarity
colnames(rna_tum) <- paste0(colnames(rna_tum), "-01A")
colnames(rna_norm) <- paste0(colnames(rna_norm), "-11A")
colnames(cnv_tum) <- paste0(colnames(cnv_tum), "-11A")


## Download Clinical data ##
colnames(rna_tum) <- substr(colnames(rna_tum), 1, 12)
patient_barcodes <- c(colnames(rna_tum))

clinical_query <- TCGAbiolinks::GDCquery(project = "TCGA-HNSC", 
                                         data.category = "Clinical",
                                         data.format = "bcr xml", 
                                         barcode = patient_barcodes)
GDCdownload(clinical_query)

clinical_data <- TCGAbiolinks::GDCprepare_clinic(clinical_query, clinical.info = "patient")

# Save data
saveRDS(rna_tum, "BRCA/rna_tumor.RDS")
saveRDS(rna_norm, "BRCA/rna_normal.RDS")
saveRDS(cnv_tum, "BRCA/cnv_tumor.RDS")
saveRDS(clinical_data, "BRCA/clinical_full.RDS")



# Formatting strings | Select samples | CN segment data
#cnv_seg$Sample <- stringr::str_sub(cnv_seg$Sample,1,12)
#colnames(cnv_tumor) <- stringr::str_sub(colnames(cnv_tumor),1,12)
#cnv_seg <- cnv_seg[cnv_seg$Sample %in% colnames(cnv_tumor),]

#cnv_seg <- cnv_seg %>%
  #dplyr::rename(
    #chr = Chromosome,
    #start = Start,
    #end = End,
    #seg_mean = Segment_Mean,
    #sample = Sample) %>%
  #dplyr::select(chr, start, end, seg_mean, sample) %>%
  #as.data.frame()

#save(lihc_cnv_seg, file = "liver_LIHC/lihc_cnv_seg.Rdata")
