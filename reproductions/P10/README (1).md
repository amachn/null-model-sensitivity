P10 - ENCODE Project Consortium (Moore et al.) 2020
Target

* Target result: Extended Data Figure 3b ~ distribution of the percent of ENCODE TF ChIP-seq peaks overlapping the cell-type-agnostic human cCRE registry
* Method: direct use of author-reported per-dataset overlap percentages (no independent peak/cCRE overlap recomputation)
* Genome build: GRCh38 / hg38

Original Source Files
Source URL: see `papers/README.md`
Code

* `human_mouse_genomes.ipynb`

Data

* `41586_2020_2493_MOESM20_ESM.xlsx` (Supplementary Table 18, sheet "Supplementary Table 18a"), downloaded directly from the publisher inside the notebook

Environment

* Operating System: Google Colab (Ubuntu-based runtime)
* Python Version: 3.x (Colab default)
* Key packages: pandas, openpyxl, matplotlib, seaborn, numpy, requests


Issues / Notes

* No independent overlap computation
   * This reproduction re-plots the authors' own precomputed per-dataset overlap percentages from Supplementary Table 18a rather than re-running peak-vs-cCRE overlap calls from raw ENCODE BED files.
   * As a result, it validates the reported summary statistic and figure shape, but does not independently test the underlying ChIP-seq/overlap methodology itself.
* Fraction vs. percent units
   * Supplementary Table 18a stores the overlap column as a fraction (0-1), not a percentage.
   * A single x100 conversion is applied and asserted (via explicit range checks) to avoid the easy mistake of double-scaling an already-percentage value.
* Reproduced vs. paper-reported median
   * Reproduced median: 90.46%
   * Paper-reported median: 90%
   * Difference: 0.46 percentage points, within the notebook's acceptance tolerance (±1 percentage point).
* Low-overlap outliers
   * 115 of 982 TF ChIP-seq datasets fall below 70% overlap; these are broken out separately in the notebook for inspection but are not treated as errors.

