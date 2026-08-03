## Scripts

- **`get_cont_samples.R`** – Extracts sample-level continuous pseudobulk CNV profiles from the outputs of six scRNA-seq CNV callers (inferCNV, CaSpER, CopyKAT, SCEVAN, Numbat, and CONICSmat).
- **`gene_to_bin100kb.R`** – Converts gene-, segment-, or chromosome arm-level CNV profiles into standardized 100-kb bin-level relative continuous CNV profiles.
- **`semi_supervise_cont_to_discrete.R`** – Converts continuous CNV profiles into discrete gain/loss/neutral states by matching the neutral-region proportion of the ground truth.

## Data

- **`bin100kb_cont/`** – 100-kb bin-level relative continuous pseudobulk CNV profiles.
- **`bin100kb_discrete/`** – 100-kb bin-level discrete CNV profiles (gain/loss/neutral).

## Note

The results of the six scRNA-seq CNV callers for the nine gastric cancer cell lines were generated using each method's default parameter settings. Therefore, the corresponding scripts and raw outputs are not included in this repository.
