setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/deconveilCaseStudies/")

pkgs <- c("GenomicRanges", "TxDb.Hsapiens.UCSC.hg38.knownGene", "AnnotationHub", "dplyr", "GenomeInfoDb",
          "SummarizedExperiment", "org.Hs.eg.db")
sapply(pkgs, require, character.only = TRUE)

data_path <- "/Users/katsiarynadavydzenka/Documents/PhD_AI/TCGA/BRCA/"

dna_meth <- readRDS(paste0(data_path, "dna_meth.RDS"))

meth_matrix <- assay(dna_meth)
cpg_coordinates <- rowData(dna_meth)

# Define promoter regions
txdb <- makeTxDbFromUCSC(genome="hg38", tablename="refGene")
genes <- genes(txdb)
promoters <- promoters(genes, upstream=2000, downstream=2000)


# Map CpGs to genes 
cpg_coordinates <- cpg_coordinates[!is.na(cpg_coordinates$beg_A) & !is.na(cpg_coordinates$probeEnd), ]

cpg_gr <- GRanges(
  seqnames = cpg_coordinates$chrm_A,  # Chromosome names
  ranges = IRanges(
    start = cpg_coordinates$beg_A,   # Start positions
    end = cpg_coordinates$probeEnd   # End positions
  )
)

genome(cpg_gr) <- "hg38"

overlaps <- findOverlaps(cpg_gr, promoters)

promoter_cpgs <- cpg_gr[queryHits(overlaps)]
matched_genes <- genes[subjectHits(overlaps)]

filtered_cpg_coordinates <- cpg_coordinates[queryHits(overlaps), ]

promoter_cpgs <- GRanges(
  seqnames = filtered_cpg_coordinates$chrm_A,
  ranges = IRanges(
    start = filtered_cpg_coordinates$beg_A,
    end = filtered_cpg_coordinates$probeEnd
  )
)

names(promoter_cpgs) <- rownames(filtered_cpg_coordinates)

cpg_to_gene <- data.frame(
  CpG_id = names(promoter_cpgs),  
  Gene_Name = mcols(matched_genes)$gene_id  # Adjust column name if needed
)


# Combine methylation data with gene information
mapped_methylation <- merge(cpg_to_gene, meth_matrix, by.x = "CpG_id", by.y = "row.names", all.x = TRUE)

meth_mean_per_gene <- mapped_methylation %>%
  group_by(Gene_Name) %>%
  summarise(across(starts_with("TCGA"), mean, na.rm = TRUE))

entrez_ids <- meth_mean_per_gene$Gene_Name
gene_symbols <- mapIds(org.Hs.eg.db, keys = as.character(entrez_ids), column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
meth_mean_per_gene$Gene_Name<- gene_symbols

saveRDS(meth_mean_per_gene, file = "/Users/katsiarynadavydzenka/Documents/PhD_AI/TCGA/BRCA/dna_meth_processed.RDS")
