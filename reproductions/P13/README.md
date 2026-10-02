# P13 - Miller et al. 2022
***Reproduction Status: awaiting QC***

## Target
- **Target result:** Figure 1C ~ enrichment of R-loop peaks around predicted R-loop forming sequences (RLFS)
- **Method:** `RLSeq::analyzeRLFS` which uses `regioneR::permTest`
- **Permutations:** 100
- **Genome build:** hg38

## Original Source Files
**Source URL:** see [`papers/README.md`](/papers/)

### Code
- `figures.R` -> `scripts/original/`

### Data
- Downloaded in [script](scripts/build_figure.R) via the author's AWS bucket.

## Environment
- **Operating System:** Windows 11 25H2 (n/a)
- **R Version:** 4.6.1 (4.1.1)
- **regioneR Version:** 1.44.0 (1.25.1)
- *parenthesized version indicates the version used in the original analysis

### Environment Setup
The environment is stored in a `renv.lock` file using `renv`; to recreate the environment used in this reproduction, run:

```r
install.packages("renv")
renv::restore()
```

## Results
The reproduced results almost exactly match the original results in Figure 1C. There may be a slight quantitative
difference, but given the small number of permutations and no defined random seed, this is to be expected.

## Issues / Notes
- Package setup for `RLHub` and `RLSeq` must be performed via GitHub exactly the way described in their [README](https://github.com/Bishop-Laboratory/RLSeq).
  - This is because the packages are no longer available in Bioconductor. The `remotes` package was the easiest way
    to install them. If packages are not installed in the correct order, dependency issues will likely occur.
- When reproducing the analysis, some samples generated large numbers of `GenomicRanges` warnings during the local z-score calculation.
  - This was resolved by using a `suppressWarnings` function over the `RLSeq::analyzeRLFS` function.