library(GenomeInfoDb)
library(BSgenome.Hsapiens.UCSC.hg38)
library(dplyr)
library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(AnnotationDbi)

chr_levels <- paste0("chr", 1:22)
bin_size <- 100000
# edb <- EnsDb.Hsapiens.v86

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
  if (is.matrix(x) || is.data.frame(x)) {
    data.frame(
      gene = rownames(x),
      score = as.numeric(x[, 1]),
      stringsAsFactors = FALSE
    )
  } else if (is.numeric(x) && !is.null(names(x))) {
    data.frame(
      gene = names(x),
      score = as.numeric(x),
      stringsAsFactors = FALSE
    )
  } else {
    stop("Unsupported object format in to_df(). Check class(obj) and str(obj).")
  }
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

to_bin_df_from_arm <- function(df_arm, chr_arm_df, bin_template) {
  chr_arm_df$Chrom <- paste0("chr", chr_arm_df$Chrom)
  chr_arm_df$Idf <- as.character(chr_arm_df$Idf)

  arm_score <- df_arm %>%
    dplyr::rename(Idf = arm, score = score)

  arm_annot <- chr_arm_df %>%
    dplyr::inner_join(arm_score, by = "Idf") %>%
    dplyr::select(chr = Chrom, start = Start, end = End, score)

  out <- bin_template %>%
    dplyr::mutate(chr = as.character(chr)) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      score = {
        hit <- arm_annot[
          arm_annot$chr == chr &
            bin_start >= arm_annot$start &
            bin_end <= arm_annot$end,
          "score"
        ]
        if (length(hit) == 0) NA_real_ else hit[1]
      }
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!is.na(score))

  out
}

to_df_arm <- function(x) {
  data.frame(
    arm = rownames(x),
    score = as.numeric(x[, 1]),
    stringsAsFactors = FALSE
  )
}

to_bin_df_from_seg <- function(seg_df, bin_template) {
  seg_df <- as.data.frame(seg_df)
  seg_df$CHROM <- as.character(seg_df$CHROM)
  seg_df <- seg_df %>% dplyr::filter(CHROM %in% chr_levels)

  out <- bin_template %>%
    dplyr::mutate(chr = as.character(chr)) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      score = {
        hit <- seg_df[
          seg_df$CHROM == chr &
            seg_df$seg_end >= bin_start &
            seg_df$seg_start <= bin_end,
          ,
          drop = FALSE
        ]

        if (nrow(hit) == 0) {
          NA_real_
        } else {
          overlap <- pmin(hit$seg_end, bin_end) - pmax(hit$seg_start, bin_start) + 1
          overlap[overlap < 0] <- 0
          if (all(overlap == 0)) NA_real_ else hit$score[which.max(overlap)]
        }
      }
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!is.na(score))

  out
}

setwd("/ix1/ctseng/hac377/cnv_calling/letter")

df <- read.table("samples.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
method_names <- c("inferCNV", "CaSpER", "CopyKAT", "SCEVAN", "Numbat", "CONICSmat")

# create bin map
bin_template <- make_bin_template(bin_size = 100000)


for (sampleID in df$sample) {
  cat("Processing sample:", sampleID, "\n")
  
  outdir <- paste0("/ix1/ctseng/hac377/cnv_calling/letter/cont_samples/bin100kb_scaled_ratio/", sampleID)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  for (method_name in method_names) {
    
    cat("Method:", method_name, "\n")
    
    infile <- paste0(
      "/ix1/ctseng/hac377/cnv_calling/letter/cont_samples/gene_raw_ratio/",
      sampleID,
      "/bulk_gene_",
      method_name,
      ".RData"
    )
    
    obj <- readRDS(infile)
    
    if (method_name == "CONICSmat") {
      
      df_arm <- to_df_arm(obj)
      
      chr.arm.df <- read.table(
        "/ix1/alee/LO_LAB/Personal/Rick/clone_review/Reference_file/CONICSmat_ref/chromosome_arm_positions_grch38.txt",
        header = TRUE
      )
      
      bin_df <- to_bin_df_from_arm(df_arm, chr.arm.df, bin_template)
      
    } else if (method_name %in% c("Numbat", "SCEVAN")) {
      
      bin_df <- to_bin_df_from_seg(obj, bin_template)
      
    } else {
      
      df_gene <- to_df(obj)
      gene_annot <- get_gene_annot_from_ensdb(df_gene$gene)
      df_gene <- df_gene %>% inner_join(gene_annot["gene"], by = "gene")
      
      bin_df <- to_bin_df_template(df_gene, gene_annot, bin_template, bin_size = 100000)
    }
    
    cn_logratio <- log2(bin_df$score)
    score_scaled <- cn_logratio / sd(cn_logratio, na.rm = TRUE)
    score_scaled[is.na(score_scaled)] <- 0
    
    out <- data.frame(
      bulk = score_scaled,
      row.names = bin_df$bin_name,
      stringsAsFactors = FALSE
    )
    
    saveRDS(
      out,
      file = paste0(outdir, "/", method_name, "_bin100kb_scaled.rds")
    )
  }
}


