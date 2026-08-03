setwd("/ix1/ctseng/hac377/cnv_calling/letter")

library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
library(BSgenome.Hsapiens.UCSC.hg38)
library(dplyr)
library(ensembldb)



chr_levels <- paste0("chr", 1:22)
bin_size <- 100000
threshold <- 0.125
edb <- EnsDb.Hsapiens.v86

make_bin_template <- function(bin_size = 100000) {
  seqlens <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)[chr_levels]
  
  bin_list <- lapply(chr_levels, function(chr) {
    chr_len <- seqlens[[chr]]
    bin_start <- seq(1, chr_len, by = bin_size)
    bin_end <- pmin(bin_start + bin_size - 1, chr_len)
    
    data.frame(
      chr = chr,
      bin_start = bin_start,
      bin_end = bin_end,
      bin_name = paste0(chr, "_", bin_start, "_", bin_end),
      stringsAsFactors = FALSE
    )
  })
  
  out <- bind_rows(bin_list)
  out$chr <- factor(out$chr, levels = chr_levels)
  out <- out %>% arrange(chr, bin_start)
  out
}

to_df <- function(x) {
  data.frame(
    gene = rownames(x),
    score = as.numeric(x[, 1]),
    stringsAsFactors = FALSE
  )
}


to_bin_df_template <- function(df_gene, gene_annot, bin_template, bin_size = 100000) {
  bin_score <- df_gene %>%
    dplyr::inner_join(gene_annot[, c("gene", "chr", "start")], by = "gene") %>%
    dplyr::mutate(
      chr = as.character(chr),
      bin_id = floor((start - 1) / bin_size),
      bin_start = bin_id * bin_size + 1
    ) %>%
    dplyr::group_by(chr, bin_start) %>%
    dplyr::summarise(score = mean(score, na.rm = TRUE), .groups = "drop")
  
  out <- bin_template %>%
    dplyr::left_join(bin_score, by = c("chr", "bin_start")) %>%
    dplyr::filter(!is.na(score))
  
  out
}


get_gene_annot_from_ensdb <- function(genes, edb = EnsDb.Hsapiens.v86) {
  gr <- genes(
    edb,
    filter = GeneNameFilter(genes),
    columns = c("gene_name")
  )

  annot <- data.frame(
    gene  = mcols(gr)$gene_name,
    chr   = paste0("chr", as.character(seqnames(gr))),
    start = start(gr),
    end   = end(gr),
    stringsAsFactors = FALSE
  )

  annot %>%
    dplyr::filter(chr %in% chr_levels) %>%
    dplyr::distinct(gene, .keep_all = TRUE) %>%
    dplyr::mutate(chr = factor(chr, levels = chr_levels)) %>%
    dplyr::arrange(chr, start)
}


df <- read.table("samples.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
ploidies <- c(3.76, 3.71, 3.29, 2.83, 2.60, 2.29, 2.25, 1.95, 1.94)
names(ploidies) <- df$sample

# create bin map
bin_template <- make_bin_template(bin_size = 100000)


##################
## CNV intensity
##################
for (sampleID in df$sample) {

  cat("Processing sample:", sampleID, "\n")
  ploidy <- ploidies[sampleID]
  
  load(paste0("/ix1/ctseng/hac377/cnv_calling/letter/cont_gt_w_ploidy_copykit/bulk_gene/", sampleID, "_bulk_gene.RData"))
  
  obj <- out$bulk.gene.ratio
  
  df_gene <- to_df(obj)
  gene_annot <- get_gene_annot_from_ensdb(df_gene$gene)
  df_gene <- df_gene %>% inner_join(gene_annot["gene"], by = "gene")
  df_gene$score <- df_gene$score * ploidy
  
  bin_df <- to_bin_df_template(df_gene, gene_annot, bin_template, bin_size = 100000)
  
  score_scaled <- (bin_df$score - ploidy) / sd(bin_df$score, na.rm = TRUE)
  score_scaled[is.na(score_scaled)] <- 0
  
  out <- data.frame(
    bulk = score_scaled,
    row.names = bin_df$bin_name,
    stringsAsFactors = FALSE
  )
  
  saveRDS(
    out,
    file = paste0("/ix1/ctseng/hac377/cnv_calling/letter/cont_gt_w_ploidy_copykit/bin100kb_scaled/", sampleID, "_bin100kb_scaled.rds")
  )
}