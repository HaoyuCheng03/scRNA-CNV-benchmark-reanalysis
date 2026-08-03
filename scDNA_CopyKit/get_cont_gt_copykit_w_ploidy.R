library(parallelDist)
library(copykit)
library(dplyr)
library(tibble)
library(SummarizedExperiment)
library(GenomicRanges)
library(BiocGenerics)
library(BiocParallel)
library(Seurat)

setwd("/ix1/ctseng/hac377/cnv_calling/letter")


df <- read.table("samples.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
ploidy <- c(3.76, 3.71, 3.29, 2.83, 2.60, 2.29, 2.25, 1.95, 1.94)
names(ploidy) <- df$sample

get_bulk_gene <- function(sampleID, ploidy) {
  
  load(paste0("/ix1/alee/LO_LAB/Personal/Rick/clone_review/Result/",
              sampleID, "/scDNA/copykit_raw.RData"))
  
  ## -----------------------------
  ## preprocess steps
  ## -----------------------------
  tumor <- runMetrics(tumor)
  
  # Mark euploid cells if they exist
  tumor <- findAneuploidCells(tumor, resolution = "auto")
  
  # Mark low-quality cells
  tumor <- findOutliers(tumor, resolution = 0.85)
  
  # Filter cells
  tumor <- tumor[, SummarizedExperiment::colData(tumor)$outlier == FALSE]
  tumor <- tumor[, SummarizedExperiment::colData(tumor)$is_aneuploid == TRUE]
  
  # kNN smoothing
  options(MulticoreParam = MulticoreParam(workers = 10))
  tumor <- knnSmooth(tumor)
  
  # Remove ploidy outliers, use given ploidy for each sample
  # tumor <- calcInteger(tumor, method = "scquantum", assay = "smoothed_bincounts")
  # tumor <- tumor[, colData(tumor)$ploidy_score < 0.4]
  
  tumor <- calcInteger(tumor, method = "fixed", ploidy_value = ploidy)
  
  cat("Average ploidy for sample ", sampleID, " is: ", mean(tumor@colData$ploidy), " ...\n")
  
  ## -----------------------------
  ## clustering
  ## -----------------------------
  tumor <- runUmap(tumor)
  
  # use k_range = c(2:45) for MKN45; c(40:75) for others; specific method for HGC27
  if (sampleID == "HGC27") {
    k_range <- seq(2, min(20, ncol(tumor) - 1))
    
  } else if (sampleID == "MKN45") {
    k_range = c(2:45)
    
  } else {
    k_range = c(40:75)
    
  }
  tumor <- findSuggestedK(tumor, k_range = k_range)
  tumor <- findClusters(tumor)
  
  # remove singletons outliers cluster c0
  tumor <- tumor[, colData(tumor)$subclones != "c0"]
  tumor@colData$subclones <- droplevels(tumor@colData$subclones)
  
  ## -----------------------------
  ## CNV ratio
  ## -----------------------------
  is_chrX <- seqnames(tumor) == "chrX"
  tumor <- tumor[!is_chrX, ]
  
  cnv.ratio <- SummarizedExperiment::assay(tumor, "segment_ratios")
  
  ## -----------------------------
  ## create gene metadata
  ## -----------------------------
  if (tumor@metadata$genome == "hg19") {
    genes_assembly <- hg19_genes
  } else if (tumor@metadata$genome == "hg38") {
    genes_assembly <- hg38_genes
  } else {
    stop("Unknown genome build in tumor@metadata$genome")
  }
  
  ranges <- tumor@rowRanges
  genes_features <- BiocGenerics::subset(genes_assembly)
  
  olaps <- suppressWarnings(
    GenomicRanges::findOverlaps(
      genes_features,
      ranges,
      ignore.strand = TRUE
    )
  )
  
  gene.meta <- data.frame(
    gene = as.character(data.frame(genes_features)$symbol[queryHits(olaps)]),
    chr = as.character(data.frame(genes_features)$seqnames[queryHits(olaps)]),
    start = as.character(data.frame(genes_features)$start[queryHits(olaps)]),
    end = as.character(data.frame(genes_features)$end[queryHits(olaps)]),
    bin_pos = subjectHits(olaps)
  ) %>%
    dplyr::distinct(gene, .keep_all = TRUE)
  
  gene.meta <- gene.meta[order(gene.meta$bin_pos, gene.meta$start), ]
  gene.meta <- gene.meta[!is.na(gene.meta$gene), , drop = FALSE]
  
  ## -----------------------------
  ## bulk-level continuous data
  ## -----------------------------
  cnv.ratio.bulk <- rowMeans(cnv.ratio)  # neutral ~ 1
  cnv.ratio.bulk <- matrix(cnv.ratio.bulk, ncol = 1)
  colnames(cnv.ratio.bulk) <- "bulk"
  
  bulk.gene.ratio <- cnv.ratio.bulk[gene.meta$bin_pos, , drop = FALSE]
  rownames(bulk.gene.ratio) <- gene.meta$gene
  
  out <- list(
    bulk.gene.ratio = bulk.gene.ratio,
    gene.meta = gene.meta
  )
  
  save(out,
       file = paste0("/ix1/ctseng/hac377/cnv_calling/letter/cont_gt_w_ploidy_copykit/bulk_gene/",
                     sampleID, "_bulk_gene.RData"))
  
  return(out)
}

# for test:
# out <- get_bulk_gene(sampleID = "SNU16")

for (sampleID in df$sample) {
  cat("Processing:", sampleID, "\n")
  
  tryCatch({
    out <- get_bulk_gene(sampleID = sampleID, ploidy = ploidy[sampleID])
    if (is.null(out)) next
    rm(out)
  }, error = function(e) {
    cat("Error in", sampleID, ":", e$message, "\n")
  })
  
  gc()
  
}


