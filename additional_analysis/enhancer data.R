setwd("/orfeo/LTS/CDSLab/LT_storage/kdavydzenka/TCGA/LUSC/")
pkgs <- c("tidyverse", "TCGAbiolinks", "SummarizedExperiment")
sapply(pkgs, require, character.only = TRUE)

setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/TCGA/LUSC/")

# Build a query to download BAM files

rna_tumor <- readRDS("~/Documents/PhD_AI/TCGA/LUSC/rna_tumor.RDS")
sample_ids <- colnames(rna_tumor)
sample_ids <- gsub(" -01A", "", sample_ids)
sample_ids <- sample_ids[1]

saveRDS(sample_ids, file = "~/Documents/PhD_AI/TCGA/LUSC/sample_ids.RDS")

query <- GDCquery(
  project = "TCGA-LUSC", 
  data.category = "Sequencing Reads",
  data.type = "Aligned Reads", 
  data.format = "bam",
  workflow.type = "STAR 2-Pass Transcriptome",
  sample.type = "Primary Tumor",
  barcode = sample_ids
)

GDCdownload(query, method = "client")

samples <- getResults(query, cols = "cases")
batch_size <- 2  
sample_groups <- split(samples, ceiling(seq_along(samples) / batch_size))

for (group in sample_groups) {
  query <- GDCquery(
    project = "TCGA-LUSC", 
    data.category = "Sequencing Reads",
    data.type = "Aligned Reads", 
    data.format = "bam",
    workflow.type = "STAR 2-Pass Transcriptome",
    sample.type = "Primary Tumor",
    barcode = sample_ids
  )
  GDCdownload(query)
}



