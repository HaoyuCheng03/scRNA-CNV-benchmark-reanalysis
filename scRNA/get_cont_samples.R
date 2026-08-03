setwd("/ix1/ctseng/hac377/cnv_calling/letter")

library(infercnv)
library(numbat)
library(data.table)
library(CONICSmat)

df <- read.table("samples.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

for (sampleID in df$sample) {

  cat("Processing:", sampleID, "\n")

  outdir <- paste0("/ix1/ctseng/hac377/cnv_calling/letter/cont_samples/gene_raw_ratio/", sampleID, "/")
  
  if (!dir.exists(outdir)) {
    dir.create(outdir)
  }


  ####################################
  ## get consensus tumor cells
  ####################################
  # load consensus malignant cells to access tumorc cell's CNV accuracy
  consensus.tumor <- read.csv(paste0("/ix1/alee/LO_LAB/Personal/Rick/clone_review/Result/", sampleID, "/Malignant_ident/consensus_result.csv"), row.names = 1)

  consensus.tumor$res <- ifelse(consensus.tumor$consensus.score >= 2, "Tumor", "Normal")
  consensus.tumor.cell <- rownames(consensus.tumor)[consensus.tumor$res == "Tumor"]
  consensus.tumor.cell <- gsub("-1","", consensus.tumor.cell)

  
  ####################################
  ## inferCNV
  ####################################
  input.file <- paste0("/ix1/alee/LO_LAB/Personal/Rick/clone_review/Result/", sampleID, "/inferCNV_wRef/run.final.infercnv_obj")
  infercnv_obj = readRDS(input.file)
  # keep tumor cells only
  # normal ~ 1
  cell.gene.cnv.infercnv <- infercnv_obj@expr.data[,colnames(infercnv_obj@expr.data) %in% consensus.tumor.cell]
  # baseline: 1 (ratio)
  bulk.gene.cnv.infercnv <- as.matrix(rowMeans(cell.gene.cnv.infercnv))
  # bulk.gene.cnv.infercnv <- bulk.gene.cnv.infercnv + 1
  saveRDS(bulk.gene.cnv.infercnv, file = paste0(outdir, "bulk_gene_inferCNV.RData"))



  ####################################
  ## CaSpER
  ####################################
  input.file <- paste0("/ix1/alee/LO_LAB/Personal/Rick/clone_review/Result/", sampleID, "/CaSpER/rna_matrix.RData")
  # keep only consensus tumor cells
  cell.gene.cnv.CaSpER <- get(load(input.file))
  cell.gene.cnv.CaSpER <- cell.gene.cnv.CaSpER[,colnames(cell.gene.cnv.CaSpER) %in% consensus.tumor.cell]

  # remove normal cells again
  no.CNV.cell <- apply(cell.gene.cnv.CaSpER, 2, function(x) all(x == 0))
  cell.gene.cnv.CaSpER <- cell.gene.cnv.CaSpER[,!no.CNV.cell]
  # cell.gene.cnv.CaSpER <- cell.gene.cnv.CaSpER + 2 # (-1,0,1) -> (1, 2, 3)
  cell.gene.cnv.CaSpER <- (cell.gene.cnv.CaSpER + 2) / 2 # (-1,0,1) -> (0.5, 1, 1.5), then create a ratio like output
  bulk.gene.cnv.CaSpER <- rowMeans(cell.gene.cnv.CaSpER)
  saveRDS(bulk.gene.cnv.CaSpER, file = paste0(outdir, "bulk_gene_CaSpER.RData"))



  ####################################
  ## CopyKAT
  ####################################
  input.file <- paste0("/ix1/alee/LO_LAB/Personal/Rick/clone_review/Result/", sampleID, "/CopyKAT_wRef/",sampleID, "_copykat_CNA_raw_results_gene_by_cell.txt")
  copykat.raw.res <- read.table(input.file, header = T)
  cell.gene.cnv.CopyKAT <- as.matrix(copykat.raw.res[,-1:-7])
  rownames(cell.gene.cnv.CopyKAT) <- copykat.raw.res$hgnc_symbol
  # tumor identification
  tumor.prediction <- read.table(paste0("/ix1/alee/LO_LAB/Personal/Rick/clone_review/Result/", sampleID, "/CopyKAT_wRef/",sampleID, "_copykat_prediction.txt"), header = T)
  tumor.cells <- tumor.prediction[grepl("aneuploid", tumor.prediction$copykat.pred, fixed = TRUE), "cell.names"]
  is.ref <- grepl("Ref", tumor.cells)
  tumor.cells <- tumor.cells[!is.ref]
  # log-ratio like
  cell.gene.cnv.CopyKAT <- cell.gene.cnv.CopyKAT[,which(colnames(cell.gene.cnv.CopyKAT) %in% tumor.cells)]
  # create a ratio like output
  cell.gene.cnv.CopyKAT <- exp(cell.gene.cnv.CopyKAT)
  # bulks CNV intensity
  bulk.gene.cnv.CopyKAT <- as.matrix(rowMeans(cell.gene.cnv.CopyKAT))
  saveRDS(bulk.gene.cnv.CopyKAT, file = paste0(outdir, "bulk_gene_CopyKAT.RData"))



  ####################################
  ## SCEVAN
  ####################################
  input.file <- paste0("/ix1/alee/LO_LAB/Personal/Rick/clone_review/Result/", sampleID, "/SCEVAN_wRef/output/", sampleID, "_Clonal_CN.seg")
  gene.info.file <- paste0("/ix1/alee/LO_LAB/Personal/Rick/clone_review/Result/", sampleID, "/SCEVAN_wRef/output/", sampleID, "_count_mtx_annot.RData")
  bulk.CNV <- read.table(input.file)
  bulk.gene.cnv.SCEVAN <- data.frame(
    CHROM = paste0("chr", bulk.CNV$Chr),
    seg_start = bulk.CNV$Pos,
    seg_end = bulk.CNV$End,
    score = exp(bulk.CNV$segm.mean),  # segm.mean is log-ratio
    stringsAsFactors = FALSE
  )
  saveRDS(bulk.gene.cnv.SCEVAN, file = paste0(outdir, "bulk_gene_SCEVAN.RData"))



  ####################################
  ## Numbat
  ####################################
  input.file <- paste0(
    "/ix1/alee/LO_LAB/Personal/Rick/clone_review/Result/",
    sampleID,
    "/Numbat/"
  )

  nb <- Numbat$new(out_dir = input.file)

  # save seg-level output
  bulk.gene.cnv.numbat <- as.data.frame(nb$segs_consensus)[, c("CHROM", "seg_start", "seg_end", "cnv_state_post")]
  bulk.gene.cnv.numbat$CHROM <- as.character(bulk.gene.cnv.numbat$CHROM)
  bulk.gene.cnv.numbat$CHROM <- ifelse(
    grepl("^chr", bulk.gene.cnv.numbat$CHROM),
    bulk.gene.cnv.numbat$CHROM,
    paste0("chr", bulk.gene.cnv.numbat$CHROM)
  )

  # map Numbat states to simple CN values
  bulk.gene.cnv.numbat$score <- NA_real_
  bulk.gene.cnv.numbat$score[bulk.gene.cnv.numbat$cnv_state_post %in% c("bamp", "amp")] <- 3
  bulk.gene.cnv.numbat$score[bulk.gene.cnv.numbat$cnv_state_post %in% c("bdel", "del")] <- 1
  bulk.gene.cnv.numbat$score[bulk.gene.cnv.numbat$cnv_state_post %in% c("loh", "neu")] <- 2
  
  # create a ratio-like output
  bulk.gene.cnv.numbat$score <- bulk.gene.cnv.numbat$score / 2

  # keep only needed columns
  bulk.gene.cnv.numbat <- bulk.gene.cnv.numbat[, c("CHROM", "seg_start", "seg_end", "score")]

  # save seg-level using your original file name
  saveRDS(
    bulk.gene.cnv.numbat,
    file = paste0(outdir, "bulk_gene_Numbat.RData")
  )

  ####################################
  ## CONICSmat
  ####################################

  input.file <- paste0("/ix1/alee/LO_LAB/Personal/Rick/clone_review/Result/", sampleID, "/CONICSmat/final_posterior.RDS")
  tumor.prediction.file <- paste0("/ix1/alee/LO_LAB/Personal/Rick/clone_review/Result/", sampleID, "/CONICSmat/tumor_cell_pred.txt")
  CONICSmat.raw.res <- readRDS(input.file)
  tumor.prediction <- read.table(tumor.prediction.file, header = TRUE)

  is.ref <- grepl("Ref", tumor.prediction$cell.names)
  tumor.prediction$CONICSmat.pred[is.ref] <- "normal"

  normal <- which(tumor.prediction$CONICSmat.pred == "normal")
  tumor  <- which(tumor.prediction$CONICSmat.pred == "tumor")

  cell.seg.cnv.CONICSmat <- binarizeMatrix(CONICSmat.raw.res, normal, tumor, threshold = 0.8)
  cell.seg.cnv.CONICSmat[is.na(cell.seg.cnv.CONICSmat)] <- 0

  states <- sapply(colnames(cell.seg.cnv.CONICSmat), function(x) strsplit(x, "_")[[1]][1])
  states <- ifelse(states == "amp", 1, -1)

  cell.seg.cnv.CONICSmat <- t(cell.seg.cnv.CONICSmat) * states
  rownames(cell.seg.cnv.CONICSmat) <- colnames(CONICSmat.raw.res)

  cell.seg.cnv.CONICSmat <- cell.seg.cnv.CONICSmat[, tumor, drop = FALSE]
  cell.seg.cnv.CONICSmat <- cell.seg.cnv.CONICSmat + 2  # (-1,0,1) -> (1, 2, 3)
  # create ratio like
  cell.seg.cnv.CONICSmat <- cell.seg.cnv.CONICSmat / 2

  ## build full arm template and fill missing arms with 1
  chr.arm.df <- read.table(
    "/ix1/alee/LO_LAB/Personal/Rick/clone_review/Reference_file/CONICSmat_ref/chromosome_arm_positions_grch38.txt",
    header = TRUE,
    stringsAsFactors = FALSE
  )

  all.arms <- as.character(chr.arm.df$Idf)

  bulk.arm.cnv.CONICSmat <- matrix(
    1,
    nrow = length(all.arms),
    ncol = 1,
    dimnames = list(all.arms, "bulk")
  )

  arm_mean <- rowMeans(cell.seg.cnv.CONICSmat, na.rm = TRUE)
  common_arms <- intersect(names(arm_mean), all.arms)
  bulk.arm.cnv.CONICSmat[common_arms, 1] <- arm_mean[common_arms]

  ## keep your original save name
  saveRDS(
    bulk.arm.cnv.CONICSmat,
    file = paste0(outdir, "bulk_gene_CONICSmat.RData")
  )

  gc()

}

