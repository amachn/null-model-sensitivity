# P03 - Chan et al. 2022

## Target

- **Target result:** Figure 5 (Supplement 4C, also listed as S18C in methods) ~ chromatin state overlap of gene annotations with genes/NADs/LADs
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
Environment and input-data setup are in progress currently.

## Issues / Notes
- The paper methods refer to the figure as both Figure 5-4c and Figure S18C.