P05 - Berge-Seidl et al. 2021
Reproduction Status: QC

## Target

* Target result: Table 1 / main text ~ enrichment of PD GWAS risk variants (+ LD proxies, r2 > 0.8) in open
  chromatin regions (OCRs) of neurons, across 14 brain regions from BOCA; primary reported hit is the
  superior temporal cortex (STC)
* Method: `regioneR::overlapPermTest()` (uniform `randomizeRegions` null). 
* Permutations: 1,000 (default in `run_regioner_pd.r`; original study used 10,000)
* Genome build: hg19

## Original Source Files

Source URL: https://bendlj01.u.hpc.mssm.edu/multireg/resources/boca_peaks.zip

### Code

* No original analysis code was released for this result. `run_regioner_pd.r` was written from scratch,
  based on the paper's Methods.

### Data

* BOCA ATAC-seq peak sets (28 BED files: 14 brain regions x {neuron, glia}) obtained from the paper's Data
  availability link (`boca_peaks.zip`) and stored in `data/boca_peaks/`. Only the `_neuron` files are used,
  matching the paper (glia were excluded due to cellular heterogeneity)
* PD GWAS risk variants + LD proxies (`PD_variants.bed`, 2,772 variants) obtained/derived per the Methods
  (index variants from the Nalls et al. 2019 meta-analysis, LD proxies at r2 > 0.8 via SNiPA/1000 Genomes)
* Negative-control variant sets: `IBD_variants.bed` (2,766 variants) and `PEF_variants.bed` (2,568 variants)
* Genome file: `hg19_main_chrom.sizes` (24 standard chromosomes, chr1-22/X/Y) used as the regioneR genome;
  `hg19_chrom.sizes` (93 sequences, includes alt/unplaced contigs) is also provided but not used by the script

## Environment

* Operating System: macOS Sequoia 15.7.4
* Hardware: MacBook Pro, Apple M2 Pro, 16 GB RAM
* R Version: 4.5.3

### Environment Setup

```r
install.packages("BiocManager")
BiocManager::install(c("regioneR", "GenomicRanges"))
```

Run with, e.g.:

```
Rscript run_regioner_pd.r data/boca_peaks/STC_neuron.bed data/PD_variants.bed data/hg19_main_chrom.sizes 1000
```

## Results

`run_regioner_pd.r` was run for STC_neuron vs. each of the three variant sets, 1,000 permutations, seed = 1:

| Comparison | Observed overlaps | Mean permuted | Fold enrichment | Z-score | p-value (raw) | Paper (raw, adj.) |
|---|---|---|---|---|---|---|
| PD vs STC neurons  | 55 | 28.90 | 1.90x | 3.602  | 0.003 | 0.028 / 6.94e-05 (significant) |
| IBD vs STC neurons | 29 | 28.83 | 1.01x | 0.023  | 0.489 | 1 / 1 (not significant) |
| PEF vs STC neurons | 26 | 26.73 | 0.97x | -0.107 | 0.543 | 1 / 1 (not significant) |

The qualitative pattern was reproduced: PD risk variants are enriched in STC neuron OCRs (observed overlap
~1.9x the permuted mean, p = 0.003), while neither negative control (IBD, PEF) shows enrichment — PEF's
observed overlap is in fact slightly below the permuted mean.

P-values are not directly comparable in magnitude to the paper's: this run's minimum achievable p-value is
~1/1,001 = 0.000999 (1,000 permutations; p = 0.003 corresponds to 3 of 1,000 permutations meeting or
exceeding the observed overlap), well above the paper's smallest reported adjusted p-value (6.94e-05, from
10,000 permutations). Significance at a conventional alpha (e.g. 0.05) and the direction/ranking of all three
comparisons match the paper.

## Issues / Notes

* No `-mask` applied to the genome
  * The paper's Data availability statement links a masked hg19 genome file
    (`hg19.fa.masked.gz`) alongside the BOCA peaks, implying repeat-masking may have been part of the
    original workflow, but this is not confirmed in the
    Methods text for the enrichment tests themselves. `run_regioner_pd.r` does not apply any mask;
    whether this matters has not yet been tested.
* Filename mismatch in the script's own usage example
  * The script's header comment gives the example filename `STC_neuron_clean.bed`, but the actual file in
    `boca_peaks.zip` is `STC_neuron.bed` (no `_clean` suffix, no separate cleaned version provided).
  * Likewise the usage comment writes `hg19_main.chrom.sizes` (dot before "chrom"), while the actual
    downloaded file is `hg19_main_chrom.sizes` (underscore). Cosmetic, but worth fixing in the script.
* Peak files carry extra annotation columns
  * `STC_neuron.bed` (and presumably the other region/cell-type files) are HOMER `annotatePeaks.pl`-style
    output, not a minimal 3-column BED — they include peak score, feature annotation, nearest gene, etc.
    `read_bed_gr()` only uses the first three columns (chrom/start/end), so this does not affect the test,
    but the files are not plain BEDs.
  * `STC_neuron.bed` has 76,145 lines, matching the paper's reported "No. ATAC-seq peaks" for STC in Table 1,
    confirming this is the correct, un-subsetted peak file.
* Two chrom.sizes files provided, only one used
  * `hg19_chrom.sizes` (93 sequences) is a superset of `hg19_main_chrom.sizes` (24 sequences), adding
    alt/unplaced contigs. The script uses `hg19_main_chrom.sizes`; whether the original analysis's background
    genome included alt contigs is not stated in the Methods.