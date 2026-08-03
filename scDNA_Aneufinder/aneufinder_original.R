library(AneuFinder)
library(BSgenome.Hsapiens.UCSC.hg38)

########################################
## Run AneuFinder on all scWGS cells
########################################
df <- read.table("samples.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

for (sampleID in df$sample) {
  
  cat("Processing:", sampleID, "\n")
  
  # sampleID <- "SNU638"
  inputfolder  <- paste0("/ix1/ctseng/hac377/cnv_calling/letter/scDNA/", sampleID, "/scDNA/cell_bam")
  outputfolder <- paste0("/ix1/ctseng/hac377/cnv_calling/letter/cont_gt_raw_aneufinder/", sampleID, "/")
  
  dir.create(outputfolder, showWarnings = FALSE, recursive = TRUE)
  
  Aneufinder(
    inputfolder = inputfolder,
    outputfolder = outputfolder,
    assembly = "hg38",
    binsizes = 1e5,
    chromosomes = c(paste0("chr", 1:22), "chrX", "chrY"),
    correction.method = "GC",
    GC.BSgenome = BSgenome.Hsapiens.UCSC.hg38,
    method = "edivisive",
    numCPU = 10,
    refine.breakpoints = FALSE
  )
  
}

