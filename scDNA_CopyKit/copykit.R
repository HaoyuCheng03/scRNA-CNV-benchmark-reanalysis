library(BiocParallel)
library(copykit)
options(MulticoreParam=MulticoreParam(workers=10))

#### preprocessing ####
setwd("/ix1/ctseng/hac377/cnv_calling/letter/")
df <- read.table("samples.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

for (sampleID in df$sample) {
  
  tumor <- runVarbin(paste0("/ix1/ctseng/hac377/cnv_calling/letter/scDNA/", sampleID, "/scDNA/cell_bam/"),
                     remove_Y = TRUE, is_paired_end = TRUE, method = "CBS") # method = "CBS" to "multipcf"
  
  save(tumor, file = paste0("/ix1/ctseng/hac377/cnv_calling/letter/", sampleID, "_copykit_raw.RData"))
}
