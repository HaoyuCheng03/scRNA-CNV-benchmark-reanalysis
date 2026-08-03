######################################################################################################
## Generate discrete bin-level ground truth for partial AUC
## Output columns:
## chr | start | end | gt_state
## rownames = bin_id (chr_start_end)
##
## Two discretization methods:
## 1. given_ploidy
## 2. original_paper
######################################################################################################

library(AneuFinder)
library(GenomicRanges)
library(S4Vectors)
library(parallel)

setwd("/ix1/ctseng/hac377/cnv_calling/letter")

df <- read.table("samples.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

ploidies <- c(3.76, 3.71, 3.29, 2.83, 2.60, 2.29, 2.25, 1.95, 1.94)
names(ploidies) <- df$sample

ncores <- 8
delta <- 0.5

out_dir <- "/ix1/ctseng/hac377/cnv_calling/letter/cont_gt_raw_aneufinder/bin100kb_discrete_gt/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


########################################
## Helper function
########################################
make_gt <- function(loss_count, gain_count, keep, non_na_keep, bin_annot_keep) {
  
  loss_frac <- loss_count[keep] / non_na_keep
  gain_frac <- gain_count[keep] / non_na_keep
  
  loss_frac[is.na(loss_frac)] <- 0
  gain_frac[is.na(gain_frac)] <- 0
  
  gt_state <- rep("base", length(loss_frac))
  gt_state[loss_frac >= 0.5] <- "loss"
  gt_state[gain_frac >= 0.5] <- "gain"
  
  conflict <- loss_frac >= 0.5 & gain_frac >= 0.5
  gt_state[conflict] <- NA
  
  gt <- data.frame(
    chr = bin_annot_keep$chr,
    start = bin_annot_keep$start,
    end = bin_annot_keep$end,
    gt_state = gt_state,
    stringsAsFactors = FALSE
  )
  
  rownames(gt) <- rownames(bin_annot_keep)
  return(gt)
}


for (sampleID in df$sample) {
  
  cat("Processing sample:", sampleID, "\n")
  
  ploidy <- ploidies[sampleID]
  
  ########################################
  ## List Aneufinder model files
  ########################################
  gc_files <- list.files(
    paste0(
      "/ix1/ctseng/hac377/cnv_calling/letter/cont_gt_raw_aneufinder/cell_raw/",
      sampleID,
      "/MODELS/method-edivisive"
    ),
    full.names = TRUE
  )
  
  if (length(gc_files) == 0) {
    cat("No files found for", sampleID, "\n")
    next
  }
  
  ########################################
  ## Read one file for bin annotation
  ########################################
  e <- new.env()
  load(gc_files[1], envir = e)
  obj <- e[[ls(e)[1]]]
  bins <- obj$bins
  
  bin_annot <- data.frame(
    chr = as.character(seqnames(bins)),
    start = start(bins),
    end = end(bins),
    stringsAsFactors = FALSE
  )
  
  bin_ids <- paste0(bin_annot$chr, "_", bin_annot$start, "_", bin_annot$end)
  rownames(bin_annot) <- bin_ids
  
  rm(e, obj, bins)
  gc()
  
  ########################################
  ## Process all cells
  ########################################
  res_list <- mclapply(gc_files, function(f) {
    
    e <- new.env()
    load(f, envir = e)
    obj <- e[[ls(e)[1]]]
    bins <- obj$bins
    
    cn <- as.numeric(mcols(bins)$copy.number)
    
    if (is.null(cn)) {
      stop("No copy.number column in file: ", f)
    }
    
    non_na <- !is.na(cn)
    zero <- cn == 0
    
    ########################################
    ## Given ploidy discretization
    ########################################
    given_loss <- non_na & cn < (ploidy - delta)
    given_gain <- non_na & cn > (ploidy + delta)
    
    ########################################
    ## Original paper discretization
    ########################################
    raw_loss <- non_na & cn <= 1
    raw_gain <- non_na & cn >= 3
    
    list(
      zero = as.numeric(zero),
      non_na = as.numeric(non_na),
      
      given_loss = as.numeric(given_loss),
      given_gain = as.numeric(given_gain),
      
      raw_loss = as.numeric(raw_loss),
      raw_gain = as.numeric(raw_gain)
    )
    
  }, mc.cores = ncores)
  
  ########################################
  ## Aggregate
  ########################################
  zero_count <- Reduce(`+`, lapply(res_list, `[[`, "zero"))
  non_na_count <- Reduce(`+`, lapply(res_list, `[[`, "non_na"))
  
  given_loss_count <- Reduce(`+`, lapply(res_list, `[[`, "given_loss"))
  given_gain_count <- Reduce(`+`, lapply(res_list, `[[`, "given_gain"))
  
  raw_loss_count <- Reduce(`+`, lapply(res_list, `[[`, "raw_loss"))
  raw_gain_count <- Reduce(`+`, lapply(res_list, `[[`, "raw_gain"))
  
  ########################################
  ## Filter mostly-zero bins
  ########################################
  zero_frac <- zero_count / non_na_count
  zero_frac[is.na(zero_frac)] <- 1
  
  keep <- zero_frac < 0.85
  
  bin_annot_keep <- bin_annot[keep, , drop = FALSE]
  non_na_keep <- non_na_count[keep]

  
  ########################################
  ## Generate outputs
  ########################################
  gt_given_ploidy <- make_gt(
    loss_count = given_loss_count,
    gain_count = given_gain_count,
    keep = keep,
    non_na_keep = non_na_keep,
    bin_annot_keep = bin_annot_keep
  )
  
  gt_original_paper <- make_gt(
    loss_count = raw_loss_count,
    gain_count = raw_gain_count,
    keep = keep,
    non_na_keep = non_na_keep,
    bin_annot_keep = bin_annot_keep
  )
  
  ########################################
  ## Save RDS only
  ########################################
  saveRDS(
    gt_given_ploidy,
    file = paste0(
      out_dir,
      sampleID,
      "_bin100kb_discrete_gt_given_ploidy.rds"
    )
  )
  
  saveRDS(
    gt_original_paper,
    file = paste0(
      out_dir,
      sampleID,
      "_bin100kb_discrete_gt_original_paper.rds"
    )
  )
  
  cat("Saved:", sampleID, "\n")
  
  rm(
    res_list,
    zero_count,
    non_na_count,
    given_loss_count,
    given_gain_count,
    raw_loss_count,
    raw_gain_count,
    gt_given_ploidy,
    gt_original_paper,
    bin_annot,
    bin_annot_keep
  )
  gc()
}