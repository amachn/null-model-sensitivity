# P01 - Vaz et al. 2025

## Target

- **Target result:** Figure 3e
- **Method:** `regioneR::overlapPermTest()`
- **Genome build:** hg38

## Original Source Files
**Source URL:** see `papers/README.md`

### Code
- `running-Region-miRNA.R`
- `running-Region-tRF.R`
- `running-Region-piRNA.R`

### Data
- Genomic annotation and region files were obtained from `Genomic-Annot-Enrichment.zip`

## Environment
- **Operating System:** Windows 11 25H2
- **R Version:** 4.5.3
- **regioneR Version:** 1.42.0

### Environment Setup
The environment is stored in a `renv.lock` file using `renv`; to recreate the environment used in this reproduction, run:

```r
install.packages("renv")
renv::restore()
```
