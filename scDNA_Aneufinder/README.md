## Scripts

- **`aneufinder_original.R`** – Runs AneuFinder on single-cell WGS BAM files to generate 100-kb bin-level continuous copy number profiles.
- **`get_cont_gt_aneufinder.R`** – Generates sample-level relative continuous pseudobulk CNV profiles by aggregating cell-level AneuFinder copy number estimates using either sample-specific ploidy or the default diploid baseline.
- **`cont_to_discrete.R`** – Generates sample-level discrete CNV ground truth by classifying each cell-level AneuFinder copy-number estimate as gain, loss, or neutral using either a sample-specific ploidy threshold or the original diploid thresholds. Bin-level labels are then assigned by majority vote across cells.


## Data

- **`bin100kb_cont.zip`**: Relative continuous pseudobulk CNV profiles using the default diploid baseline.
- **`bin100kb_cont_w_ploidy.zip`**: Relative continuous pseudobulk CNV profiles using sample-specific ploidy.
- **`bin100kb_discrete.zip`**: Relative discrete pseudobulk CNV ground truth (gain/loss/neutral) using the default diploid baseline.
- **`bin100kb_discrete_w_ploidy.zip`**: Relative discrete pseudobulk CNV ground truth (gain/loss/neutral) using sample-specific ploidy.
