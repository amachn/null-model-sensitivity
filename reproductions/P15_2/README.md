P15 - Wenz et al. 2026
Reproduction Status: awaiting QC

## Target

* Target result 1: Table S4 ~ enrichment of liver H3K27ac and H3K4me3 ChIP-seq regions (GEO GSE128072) in all
  liver ATAC-seq peaks (n = 2,518,633); published 1.86-fold (H3K27ac) and 1.83-fold (H3K4me3), permuted p < 0.001
* Target result 2: Figure S10 / Table S6 ~ enrichment of FDR5 caQTL peaks (n = 14,076) in 11 annotatr genomic
  annotation categories
* Target result 3: Figure S11 / Table S7 ~ enrichment of FDR5 caQTL peaks in 7 NIH Roadmap liver histone marks
  (H3K27ac, H3K27me3, H3K36me3, H3K4me1, H3K4me3, H3K9ac, H3K9me3)
* Method:
  * Table S4: `bedtools shuffle -chrom` (size- and chromosome-matched random regions) + `bedtools intersect`,
    as described in the Methods ("ATAC-seq/ChIP-seq enrichment")
  * Figures S10/S11: annotatr annotations with random-region nulls (bedtools shuffle and annotatr/regioneR);
    the procedure for these figures is not described in the Methods (see Issues)
* Permutations: 1,000 (same as the original study; Tables S6/S7 contain 1,000 shuffles each)
* Genome build: hg38 (GRCh38); ChIP-seq regions lifted from hg19

### Code

* The authors' GitHub repository (`bwenz91/Liver_caQTL_mapping_code`) contains only ATAC-seq processing and
  caQTL-mapping scripts, not the enrichment analyses. The following scripts were written from scratch, based on
  the paper's Methods, to reproduce them:
  * `S0_download.sh` -> `s10_s11/` — downloads Zenodo peaks/caQTLs and hg38 reference files
  * `S1_prepare_query.R` -> `s10_s11/` — maps caQTL IDs to peaks; writes caQTL-peak, all-peak and mask BEDs
  * `S2_fetch_encode_liver_marks.py`, `S2b_build_mark_variants.py` -> `s10_s11/` — retrieve the Roadmap liver
    histone marks from the ENCODE portal and build alternative mark-file inputs (Figure S11)
  * `S3_reproduce_S10_S11.R`, `S3_lib.R` -> `s10_s11/` — Figure S10/S11 reproduction under configurable null
    (`annotatr`, `bedtools`), overlap threshold and exclusion mask
  * `T0_prepare_chip.sh`, `T4_bedtools_intersect.sh` -> `s10_s11/` — Table S4 inputs and reproduction
    (three overlap definitions)
  * `null_spread_from_tables.py` -> `s10_s11/` — compares the spread of the published shuffles (Tables S6/S7)
    with independent random placement
  * `plot_reproduction_colab.py`, `plot_tableS4_colab.py` -> `s10_s11/` — figures
  * `RUN.md` -> `s10_s11/` — full command sequence

### Data

* Liver ATAC-seq peaks (`genrichAllPeaks_m10_g50_9.15.21.noBLnarrowpeak.gz`, n = 2,518,633) and FDR5 caQTLs
  (`rasqual_liver_175samples_fdr5_caqtls_4.1.24.txt.gz`, n = 14,076) were obtained from the authors' Zenodo
  record (doi:10.5281/zenodo.18329519) and stored in `data/zenodo/`
* Liver H3K27ac and H3K4me3 regions were obtained from GEO GSE128072 (`GSE128072_Pr1_*_FeatureCounts.txt.gz`;
  region coordinates converted to BED) and stored in `data/chip/`
* Roadmap liver histone-mark peaks (reference epigenome ENCSR000FYQ, 25-year-old female donor, GRCh38
  pseudoreplicated peaks) were obtained from the ENCODE portal and stored in `data/marks_F25/`; exact
  experiment and file accessions are recorded in `marks_manifest.tsv`
* ENCODE hg38 blacklist (v2), UCSC hg38 gap table, hg38 chromosome sizes and the UCSC hg19ToHg38 chain were
  obtained externally and stored in `data/ref/`
* Published values (Tables S6/S7, including all 1,000 shuffles per row) were obtained from the paper's
  supplementary material; observed and median values are stored in `s10_s11/figure_values.tsv`

## Environment

* Operating System: macOS Sequoia 15.7.4
* Hardware: MacBook Pro, Apple M2 Pro, 16 GB RAM
* Editor: VS Code
* R Version: 4.3.3 (Bioconductor 3.18; annotatr 1.28, TxDb.Hsapiens.UCSC.hg38.knownGene 3.18.0)
* Python Version: 3.11
* Key tools/packages: bedtools, CrossMap, UCSC liftOver; R: annotatr, regioneR, GenomicRanges, data.table,
  ggplot2; Python: pandas, numpy, matplotlib, openpyxl (see `install_env.R`)

### Environment Setup

The environment is a conda environment plus pinned Bioconductor 3.18 packages; to recreate it, run:

```
conda create -n s10s11 -c conda-forge -c bioconda r-base=4.3 r-recommended python=3.11 bedtools wget -y
conda activate s10s11
export PATH="$CONDA_PREFIX/bin:$PATH"
python3 -m pip install CrossMap pandas matplotlib openpyxl
Rscript s10_s11/install_env.R
```

## Results

Table S4 was reproduced. Using only the documented method, the base-pair overlap enrichment of ATAC-seq peaks
was 1.83-fold for H3K4me3 (published 1.83) and 1.75-fold for H3K27ac (published 1.86), both p < 0.001. ChIP-peak
merging, liftover tool (CrossMap vs UCSC liftOver) and number of shuffles changed the folds by < 0.005. The
published folds correspond to base-pair overlap; measured as the number of overlapping peaks, H3K4me3 reverses
to significant depletion (0.97-fold; all 1,000 shuffles exceed the observed count) and H3K27ac falls to 1.08-fold.

Figure S10 was partially reproduced. Observed counts matched Table S6 exactly for all 11 categories (with
annotatr 1.28 and a 200 bp minimum overlap). The spread of the published null matched a chromosome-matched
random shuffle (SD / binomial SD 0.84-1.00 published vs 0.86-1.00 reproduced), but its median did not: published
random-set medians were 6-60% higher for genic categories and 21% lower for intergenic than any documented null
produced. Enrichment directions matched the published values for 8-10 of 11 categories depending on the null.

Figure S11 did not reproduce. Observed and random-set counts were 76-93% lower than Table S7 under every input
tested (four ENCODE mark-file variants), so the original histone-mark input could not be recovered. Directions
matched for 4-5 of 7 marks; the published H3K9me3 enrichment (1.21-fold) was never reproduced (0.78-1.00-fold).
The published S7 null is 1.8-3.4x wider than independent random placement, unlike S6, so the two figures used
different, undocumented randomization procedures.

## Issues / Notes


* No analysis code was available
  * The linked GitHub repository contains only ATAC-seq processing and caQTL-mapping scripts.
* Figures S10/S11 method is not described in the Methods
  * The 200 bp minimum overlap was inferred from the published counts; 1,000 shuffles and the p-value definition
    were taken from Table S6.
  * The published S10 null median could not be reconstructed; annotatr/regioneR and bedtools gave identical
    results under the same null, so the difference lies in the null, not the tool.
* Figure S11 input data not specified
  * Marks are described only as Roadmap "liver samples"; ENCSR000FYQ was identified by mark coverage, but it
    mixes three donors and ~10 peak files per experiment. No tested variant reproduced the published counts.
* Inconsistent null procedures and statistics
  * Table S7 shuffles are 1.8-3.4x more variable than Table S6's, indicating a different randomization; under an
    independent-placement null, H3K27ac in Figure S11 (published p = 0.057) would be significant.
  * p-value definitions differ between Tables S6 and S7, and the Methods describe only a one-sided enrichment test.
* Unreported overlap metric for Table S4
  * Published folds match base-pair overlap; counting overlapping peaks reverses H3K4me3 to depletion, consistent
    with the paper's own Figure S7. H3K27ac remains 6% below the published fold.
* Other notes
  * GSE128072 has no peak BED files; regions were taken from its `FeatureCounts` tables.
  * The Figure S10 caption contradicts Table S6 (5'UTR not significant; intergenic enriched).
  * Bioconductor 3.18 was pinned because annotatr `randomize_regions()` is deprecated in current releases.
