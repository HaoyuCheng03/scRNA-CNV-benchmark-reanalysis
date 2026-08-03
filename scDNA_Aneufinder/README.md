## Scripts

- **`aneufinder_original.R`** – Runs AneuFinder on single-cell WGS BAM files to generate 100-kb bin-level continuous copy number profiles. :contentReference[oaicite:0]{index=0}
- **`get_cont_gt_aneufinder.R`** – Constructs sample-level continuous pseudobulk CNV ground truth by aggregating AneuFinder copy number estimates using either sample-specific ploidy or the original discretization scheme. :contentReference[oaicite:1]{index=1}
- **`semi_supervise_cont_to_discrete.R`** – Converts continuous pseudobulk CNV signals into discrete gain/loss/neutral labels using a **semi-supervised strategy**: continuous copy number estimates guide the discretization while biologically informed copy number thresholds define the final ground-truth states. :contentReference[oaicite:2]{index=2}
