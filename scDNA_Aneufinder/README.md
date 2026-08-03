## Scripts

- **`aneufinder_original.R`** – Runs AneuFinder on single-cell WGS BAM files to generate 100-kb bin-level continuous copy number profiles.
- **`get_cont_gt_aneufinder.R`** – Constructs sample-level relative continuous pseudobulk CNV ground truth by aggregating AneuFinder copy number estimates using either sample-specific ploidy or the default diploid threshold.
- **`semi_supervise_cont_to_discrete.R`** – Converts continuous pseudobulk CNV profiles into discrete gain/loss/neutral labels using a **semi-supervised strategy** that leverages the continuous CNV signal to preserve natural CNV segments percentage during discretization.
