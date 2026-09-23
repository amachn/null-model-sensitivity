# P03 - Chan et al. 2022

## Target
- **Target result:** Figure 5 (Supplement 4C, also listed as S18C in methods) ~ chromatin state overlap of gene annotations with genes/NADs/LADs
    - NADs are nucleolar-associated domains, and LADs are lamin-associated domains.
- **Method:** `regioneR::permTest()`
- **Permutations:** 100
- **Genome build:** hg38

## Original Source Files
**Source URL:** see [`papers/README.md`](/papers/)

### Code
- `figures/Figure_5/Figure_5.R` -> `scripts/original/`
- `figures/SuppFig_ATAC.R` -> `scripts/original/`

### Data
- `annotations/E017_25_imputed12marks_hg38lift_dense.bed` -> `data/hg38_im90_chromatin_states.bed`
- `annotations/GSE49341_LMNB1_Gro.bed` -> `data/hg38_im90_lads.bed`
- `annotations/hg38_im90_nads.bed` -> `data/`
- Gene annotations were obtained from `TxDb.Hsapiens.UCSC.hg38.knownGene` rather than a separate BED file.

## Environment
- **Operating System:** Windows 11 25H2 (n/a)
- **R Version:** 4.6.1 (4.0.3 and 3.6.2)
- **regioneR Version:** 1.44.0 (1.22.0)
- *parenthesized version indicates the version used in the original analysis

### Environment Setup
The environment is stored in a `renv.lock` file using `renv`; to recreate the environment used in this reproduction, run:

```r
install.packages("renv")
renv::restore()
```

## Results
The reproduced results generally succeeded in qualitatively reproducing the overall pattern in Figure 5-4c/S18C.
In particular:
- * list out good matches here

### Discrepancies

For `25_Quies`, repeated runs produced:
- NAD Z-score range = 157.7–187.4; mean = 175.4
- LAD Z-score range = 62.3–75.2; mean = 68.3

The paper reported approximately Z = 150 for NAD overlap and Z = 90 for LAD overlap.

Aside from `25_Quies`, most cells reproduced the published direction of enrichment/depletion. However, there were more discrepancies, including:
- * list further numerical differences in the figure

## Issues / Notes
- The paper methods refer to the figure as both Figure 5-4c and Figure S18C.
- The paper methods describe using `overlapPermTest`, while the publicly available code uses `permTest`.
  - Testing `overlapPermTest` produced drastically different Z-scores and didn't match the published heatmap.
  - The public `permTest` implementation matched the published heatmap much closer, albeit with some discrepancies.
- No random seed was used in the original analysis code, so the results do have variability.
- The figure has three columns, NADs/LADs/Genes, but the public code only tests two columns, NADs/LADs.
- The original analysis used different versions of R and Bioconductor packages, which is likely a contributing factor to the remaining quantitative differences.