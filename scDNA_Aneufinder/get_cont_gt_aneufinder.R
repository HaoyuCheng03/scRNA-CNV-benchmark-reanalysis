######################################################################################################
## given ploidy
######################################################################################################
library(AneuFinder)
library(GenomicRanges)
library(S4Vectors)
library(parallel)

setwd("/ix1/ctseng/hac377/cnv_calling/letter")

df <- read.table("samples.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
ploidies <- c(3.76, 3.71, 3.29, 2.83, 2.60, 2.29, 2.25, 1.95, 1.94)
names(ploidies) <- df$sample

ncores <- 8   # change based on your job resources

for (sampleID in df$sample) {
  
  cat("Processing sample:", sampleID, "\n")
  
  ploidy <- ploidies[sampleID]
  
  ########################################
  ## List Aneufinder model files
  ########################################
  gc_files <- list.files(
    paste0("/ix1/ctseng/hac377/cnv_calling/letter/cont_gt_raw_aneufinder/cell_raw/", 
           sampleID, "/MODELS/method-edivisive"),
    full.names = TRUE
  )
  
  if (length(gc_files) == 0) {
    cat("No files found for", sampleID, "\n")
    next
  }
  
  ########################################
  ## Read one file first to get bin annotation
  ########################################
  e <- new.env()
  load(gc_files[1], envir = e)
  obj_name <- ls(e)[1]
  obj <- e[[obj_name]]
  
  bins <- obj$bins
  
  if (is.null(mcols(bins)$copy.number)) {
    stop("Column 'copy.number' not found in: ", gc_files[1],
         "\nAvailable columns: ", paste(colnames(mcols(bins)), collapse = ", "))
  }
  
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
  ## Read all files in parallel
  ########################################
  res_list <- mclapply(gc_files, function(f) {
    e <- new.env()
    load(f, envir = e)
    obj_name <- ls(e)[1]
    obj <- e[[obj_name]]
    bins <- obj$bins
    cn <- as.numeric(mcols(bins)$copy.number)
    
    if (is.null(cn)) {
      stop("No copy.number column in file: ", f)
    }
    
    zero_vec <- (cn == 0)
    non_na_vec <- !is.na(cn)
    
    delta <- 0.5
    cn2 <- rep(2, length(cn))
    cn2[cn < (ploidy - delta)] <- 1
    cn2[cn > (ploidy + delta)] <- 3
    cn2[is.na(cn)] <- 0
    
    # cn2 <- cn
    # cn2[cn2 <= 1] <- 1
    # cn2[cn2 > 3] <- 3
    # cn2[is.na(cn2)] <- 0
    
    list(
      zero = as.numeric(zero_vec),
      non_na = as.numeric(non_na_vec),
      sum_signal = as.numeric(cn2)
    )
  }, mc.cores = ncores)
  
  ########################################
  ## Combine across cells
  ########################################
  zero_count <- Reduce(`+`, lapply(res_list, `[[`, "zero"))
  non_na_count <- Reduce(`+`, lapply(res_list, `[[`, "non_na"))
  sum_signal <- Reduce(`+`, lapply(res_list, `[[`, "sum_signal"))
  
  ########################################
  ## Filter mostly-zero bins
  ########################################
  zero_frac <- zero_count / non_na_count
  zero_frac[is.na(zero_frac)] <- 1
  keep <- zero_frac < 0.85
  
  bin_annot <- bin_annot[keep, , drop = FALSE]
  wgs_bulk_signal <- sum_signal[keep] / non_na_count[keep]
  
  ########################################
  ## Pseudobulk across all cells
  ########################################
  ground_truth <- cbind(
    bin_annot,
    raw_gc_corrected_bulk = wgs_bulk_signal
  )
  
  ########################################
  ## Standardize pseudobulk data
  ########################################
  s <- sd(wgs_bulk_signal, na.rm = TRUE)
  ground_truth$z <- (wgs_bulk_signal - 2) / s
  ground_truth$z[is.na(ground_truth$z)] <- 0
  
  saveRDS(
    ground_truth,
    file = paste0(
      "/ix1/ctseng/hac377/cnv_calling/letter/cont_gt_raw_aneufinder/bin100kb_scaled_ploidy/",
      sampleID, "_bin100kb_scaled_ploidy.rds"
    )
  )
  
  rm(res_list, zero_count, non_na_count, sum_signal, ground_truth, bin_annot, wgs_bulk_signal)
  gc()
}

######################################################################################################
## original paper
######################################################################################################
library(AneuFinder)
library(GenomicRanges)
library(S4Vectors)
library(parallel)

setwd("/ix1/ctseng/hac377/cnv_calling/letter")

df <- read.table("samples.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
ploidies <- c(3.76, 3.71, 3.29, 2.83, 2.60, 2.29, 2.25, 1.95, 1.94)
names(ploidies) <- df$sample

ncores <- 8   # change based on your job resources

for (sampleID in df$sample) {

  cat("Processing sample:", sampleID, "\n")

  ########################################
  ## List Aneufinder model files
  ########################################
  gc_files <- list.files(
    paste0("/ix1/ctseng/hac377/cnv_calling/letter/cont_gt_raw_aneufinder/cell_raw/",
           sampleID, "/MODELS/method-edivisive"),
    full.names = TRUE
  )

  if (length(gc_files) == 0) {
    cat("No files found for", sampleID, "\n")
    next
  }

  ########################################
  ## Read one file first to get bin annotation
  ########################################
  e <- new.env()
  load(gc_files[1], envir = e)
  obj_name <- ls(e)[1]
  obj <- e[[obj_name]]

  bins <- obj$bins

  if (is.null(mcols(bins)$copy.number)) {
    stop("Column 'copy.number' not found in: ", gc_files[1],
         "\nAvailable columns: ", paste(colnames(mcols(bins)), collapse = ", "))
  }

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
  ## Read all files in parallel
  ########################################
  res_list <- mclapply(gc_files, function(f) {
    e <- new.env()
    load(f, envir = e)
    obj_name <- ls(e)[1]
    obj <- e[[obj_name]]
    bins <- obj$bins
    cn <- as.numeric(mcols(bins)$copy.number)

    if (is.null(cn)) {
      stop("No copy.number column in file: ", f)
    }

    zero_vec <- (cn == 0)
    non_na_vec <- !is.na(cn)

    cn2 <- cn
    cn2[cn2 <= 1] <- 1
    cn2[cn2 > 3] <- 3
    cn2[is.na(cn2)] <- 0

    list(
      zero = as.numeric(zero_vec),
      non_na = as.numeric(non_na_vec),
      sum_signal = as.numeric(cn2)
    )
  }, mc.cores = ncores)

  ########################################
  ## Combine across cells
  ########################################
  zero_count <- Reduce(`+`, lapply(res_list, `[[`, "zero"))
  non_na_count <- Reduce(`+`, lapply(res_list, `[[`, "non_na"))
  sum_signal <- Reduce(`+`, lapply(res_list, `[[`, "sum_signal"))

  ########################################
  ## Filter mostly-zero bins
  ########################################
  zero_frac <- zero_count / non_na_count
  zero_frac[is.na(zero_frac)] <- 1
  keep <- zero_frac < 0.85

  bin_annot <- bin_annot[keep, , drop = FALSE]
  wgs_bulk_signal <- sum_signal[keep] / non_na_count[keep]

  ########################################
  ## Pseudobulk across all cells
  ########################################
  ground_truth <- cbind(
    bin_annot,
    raw_gc_corrected_bulk = wgs_bulk_signal
  )

  ########################################
  ## Standardize pseudobulk data
  ## assume neutral = 2
  ########################################
  s <- sd(wgs_bulk_signal, na.rm = TRUE)
  ground_truth$z <- (wgs_bulk_signal - 2) / s
  ground_truth$z[is.na(ground_truth$z)] <- 0

  saveRDS(
    ground_truth,
    file = paste0(
      "/ix1/ctseng/hac377/cnv_calling/letter/cont_gt_raw_aneufinder/bin100kb_scaled_raw/",
      sampleID, "_bin100kb_scaled_raw.rds"
    )
  )

  rm(res_list, zero_count, non_na_count, sum_signal, ground_truth, bin_annot, wgs_bulk_signal)
  gc()
}