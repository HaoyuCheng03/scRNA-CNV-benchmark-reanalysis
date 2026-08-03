## Scripts

- **`aneufinder_original.R`** – Runs AneuFinder on single-cell WGS BAM files to generate 100-kb bin-level continuous copy number profiles.
- **`get_cont_gt_aneufinder.R`** – Generates sample-level relative continuous pseudobulk CNV profiles by aggregating cell-level AneuFinder copy number estimates using either sample-specific ploidy or the default diploid baseline.
- **`cont_to_discrete.R`** – Generates sample-level discrete CNV ground truth by classifying each cell-level AneuFinder copy-number estimate as gain, loss, or neutral using either a sample-specific ploidy threshold or the original diploid thresholds. Bin-level labels are then assigned by majority vote across cells.
