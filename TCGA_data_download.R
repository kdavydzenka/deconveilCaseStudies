### Download data from GDC portal ###

#rm(list=ls())
setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/TCGA/")
pkgs <- c("tidyverse", "TCGAbiolinks", "SummarizedExperiment")
sapply(pkgs, require, character.only = TRUE)


# Get a list of projects #
gdcprojects = getGDCprojects()
getProjectSummary('TCGA-PRAD')

# build a query to retrieve data #

query_TCGA_cnv <- TCGAbiolinks::GDCquery(project = 'TCGA-HNSC',
                                         data.category = "Copy Number Variation",
                                         sample.type = c("Primary Tumor"),
                                         data.type = "Gene Level Copy Number")
                                         #data.type = "Copy Number Segment")
                                         #workflow.type = 'ASCAT3'
                                         #workflow.type = "DNAcopy")
    

query_TCGA_rna_tumor <- TCGAbiolinks::GDCquery(project = 'TCGA-HNSC',
                       data.category = 'Transcriptome Profiling',
                       experimental.strategy = 'RNA-Seq',
                       workflow.type = 'STAR - Counts',
                       data.type = "Gene Expression Quantification",
                       sample.type = "Primary Tumor",
                       access = 'open')
              

query_TCGA_rna_normal <- TCGAbiolinks::GDCquery(project = 'TCGA-HNSC',
                           data.category = 'Transcriptome Profiling',
                           experimental.strategy = 'RNA-Seq',
                           workflow.type = 'STAR - Counts',
                           data.type = "Gene Expression Quantification",
                           sample.type = "Solid Tissue Normal",
                           access = 'open')
                      
# download data
GDCdownload(query_TCGA_cnv)
GDCdownload(query_TCGA_rna_tumor)
GDCdownload(query_TCGA_rna_normal)


# Prepare data #

cnv_tumor <- TCGAbiolinks::GDCprepare(query_TCGA_cnv, summarizedExperiment = T)
cnv_tum <- assay(cnv_tumor, 'copy_number', rownames = TRUE) %>% as.data.frame()
gene_name <- as.data.frame(cnv_tumor@rowRanges@elementMetadata@listData[["gene_name"]])
colnames(gene_name)[1] <- "GeneID"
cnv_tum <- cbind(gene_name, cnv_tum)

cnv_tum <- cnv_tum[!duplicated(cnv_tum$GeneID), ] %>% 
  remove_rownames %>% 
  column_to_rownames(var="GeneID") %>%
  na.omit()

colnames(cnv_tum) <- substr(colnames(cnv_tum), 1, 12)
cnv_tum <- cnv_tum[, !duplicated(colnames(cnv_tum))]

rna_normal <- TCGAbiolinks::GDCprepare(query_TCGA_rna_normal, summarizedExperiment = TRUE)
rna_norm <- assay(rna_normal, 'unstranded', rownames = TRUE)
gene_name <- as.data.frame(rna_normal@rowRanges@elementMetadata@listData[["gene_name"]]) 
colnames(gene_name)[1] <- "GeneID"
rna_norm <- cbind(gene_name, rna_norm)

rna_norm <- rna_norm[!duplicated(rna_norm$GeneID), ] %>% 
  remove_rownames %>% 
  column_to_rownames(var="GeneID") %>%
  na.omit()


rna_tumor <- TCGAbiolinks::GDCprepare(query_TCGA_rna_tumor, summarizedExperiment = TRUE)
rna_tum <- assay(rna_tumor, 'unstranded', rownames = TRUE) %>% as.data.frame()
gene_name <- as.data.frame(rna_tumor@rowRanges@elementMetadata@listData[["gene_name"]]) 
colnames(gene_name)[1] <- "GeneID"
rna_tum <- cbind(gene_name, rna_tum)

rna_tum <- rna_tum[!duplicated(rna_tum$GeneID), ] %>% 
  remove_rownames %>% 
  column_to_rownames(var="GeneID") %>%
  na.omit()

#substring columns
colnames(rna_norm) <- substr(colnames(rna_norm), 1, 12)
colnames(rna_tum) <- substr(colnames(rna_tum), 1, 12)
colnames(cnv_tum) <- substr(colnames(cnv_tum), 1, 12)


## Clinical data ##
colnames(rna_tum) <- substr(colnames(rna_tum), 1, 12)
patient_barcodes <- c(colnames(rna_tum))

clinical_query <- TCGAbiolinks::GDCquery(project = "TCGA-LUSC", data.category = "Clinical",
                                          data.format = "bcr xml", barcode = patient_barcodes)
GDCdownload(clinical_query)

clinical_data <- TCGAbiolinks::GDCprepare_clinic(clinical_query, clinical.info = "patient")

saveRDS(clinical_data, file = "brca/clinical_full.RDS")


# Data preprocessing #
rna_tum <- rna_tum[,(colnames(rna_tum) %in% colnames(cnv_tum))] 
rna_norm <- rna_norm[,(colnames(rna_norm) %in% colnames(rna_tum))]
cnv_tum <- cnv_tum[,(colnames(cnv_tum) %in% colnames(rna_norm))]

# Reordering to match datasets #
colnames_idx <- match(colnames(rna_tum), colnames(rna_norm))
rna_norm <- rna_norm[,colnames_idx]
colnames_idx <- match(colnames(rna_tum), colnames(cnv_tum))
cnv_tum <- cnv_tum[,colnames_idx]

cnv_tum <- cnv_tum[(rownames(cnv_tum) %in% rownames(rna_tum)),] 
rna_tum <- rna_tum[(rownames(rna_tum) %in% rownames(cnv_tum)),]
rna_norm <- rna_norm[(rownames(rna_norm) %in% rownames(cnv_tum)),]

rownames_idx <- match(rownames(rna_tum), rownames(cnv_tum))
rna_tum <- rna_tum[rownames_idx,] %>% na.omit()
rna_norm <- rna_norm[rownames_idx,] %>% na.omit()
cnv_tum <- cnv_tum[rownames_idx,] %>% na.omit()


x <- colnames(rna_norm)
names(rna_norm) <- paste(x,"-11A")

x <- colnames(rna_tum)
names(rna_tum) <- paste(x,"-01A")

x <- colnames(cnv_tum)
names(cnv_tum) <- paste(x,"-11A")

saveRDS(cnv_tum, file = "HNSC/cnv_tumor.RDS")
saveRDS(rna_norm, file = "HNSC/rna_normal.RDS")
saveRDS(rna_tum, file = "HNSC/rna_tumor.RDS")



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
