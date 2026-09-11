# P01 - Vaz et al. 2025

## Target

- **Target result:** Figure 3e ~ genomic annotation enrichment of baseline sperm sncRNA genomic locations
- **Method:** `regioneR::overlapPermTest()`
- **Genome build:** hg38

## Original Source Files
**Source URL:** see [`papers/README.md`](/papers/)

### Code
- `running-Region-miRNA.R`
- `running-Region-tRF.R`
- `running-Region-piRNA.R`

### Data
- Genomic annotation and region files were obtained from `Genomic-Annot-Enrichment.zip`

## Environment
- **Operating System:** Windows 11 25H2
- **R Version:** 4.6.1
- **regioneR Version:** 1.44.0

### Environment Setup
The environment is stored in a `renv.lock` file using `renv`; to recreate the environment used in this reproduction, run:

```r
install.packages("renv")
renv::restore()
```

## Issues / Notes
- Genome mask used in `regioneR::overlapPermTest()`
  - The original analysis script does not specify a genome mask, so the `hg38` genome in the original study may be masked or unmasked.
  - I've generated results for both a masked `hg38` genome and an unmasked one.
  - After visual inspection, it appears that the unmasked `hg38` genome looks closest to the paper.
- Extra annotation file
  - There are two annotation files for exons: `chr1-22.codingexon.bed` and `chr1-22.exonplus.bed`.
  - Only one of these appears to be used in the original figure.
  - After examining the reproduction results, it appears that the coding exons are used.