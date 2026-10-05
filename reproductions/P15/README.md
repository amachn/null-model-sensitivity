# P15 – Recreation of Figure S7 (ATAC-seq) from Wenz et al. 2026, AJHG 113:260–275, doi:10.1016/j.ajhg.2026.01.003

Figure S7 compares the liver ATAC-seq peaks with size- and chromosome-matched random regions for overlap with NIH Roadmap liver histone marks (Table S3).

The authors' GitHub repo (github.com/bwenz91/Liver_caQTL_mapping_code) has **no code for S7**, only ATAC processing and caQTL mapping. The scripts here are an independent re-implementation of the paper's described method, run on the authors' data.

## Data (`data/`)
- `genrichAllPeaks_m10_g50_9.15.21.noBLnarrowpeak.gz`: authors' 2,518,633 Genrich ATAC peaks, from Zenodo 10.5281/zenodo.18329519.
- `mmc2.xlsx`: the paper's supplementary tables (Table S3 is sheet "S3 - AllLiverPeaks NIH Roadmap").
- `original_FigS7_page.png`: the published Figure S7, for comparison.
- `encode_liver_histone_used/`: 25 GRCh38 liver histone ChIP-seq peak files from the ENCODE portal (Roadmap data), listed in `files_used.tsv`.
- `hg38.autosomes.genome`: UCSC hg38 chromosome sizes, chr1–22, used by `bedtools shuffle`.

## Code (`code/`, run from `P15/`)
| script | what it does |
|---|---|
| `figS7_from_table.R` | Redraws S7 from Table S3 and re-checks the medians and enrichments from the stored shuffle counts |
| `figS7_roadmap_rerun.sh [n=1000] [par=4]` | Re-run: counts ATAC × pooled histone-mark overlaps (`bedtools intersect`), then repeats on 1,000 `bedtools shuffle -chrom` random sets |
| `figS7_plot_rerun.R` | Plots the re-run figure and writes `TableS3_rerun.csv` and `S7_rerun_vs_paper.csv` |

Requirements: R 4.x with ggplot2, readxl, tidyr and data.table, plus bedtools ≥ 2.30.

## Results (`results/`)
- `FigS7_from_TableS3.png|pdf`: figure redrawn from the authors' Table S3, an exact recreation.
- `FigS7_rerun.png|pdf`: figure re-run from the raw peaks.
- `TableS3_rerun.csv`: re-run table with observed counts, shuffled median, p-values and enrichment.
- `S7_rerun_vs_paper.csv`: paper vs re-run, side by side.
- `S7_rerun/observed.csv`, `S7_rerun/shuffled_counts.csv`: observed counts and all 1,000 shuffles.

## Reproduction status
- **Figure from Table S3: exact.**
- **Observed counts: exact for all 7 marks.**
  - The paper does not name the histone files.
  - Matching Table S3 required pooling every Roadmap liver donor's 2021 ENCODE "pseudoreplicated peaks" file per mark.
  - It also required counting all ATAC–mark overlap pairs (`bedtools intersect -a ATAC -b mark`, no `-u`).
- **Random medians: 2–4% higher than the paper**, so all re-run enrichments are slightly lower.
  - The exact shuffle settings are not given in the paper.
  - None of these variants matched: excluding gaps, hg19 sizes, no `-chrom`, adding chrX/Y/alt contigs.
- **Conclusions: 6 of 7 marks agree** (all p = 0.001).
  - H3K27ac goes from 1.010 (barely enriched) to 0.987 (barely depleted).

| Mark | Paper | Re-run |
|---|---|---|
| H3K27me3 | 1.33 | 1.28 |
| H3K4me1 | 1.21 | 1.18 |
| H3K36me3 | 1.20 | 1.16 |
| H3K9ac | 1.17 | 1.14 |
| H3K9me3 | 1.10 | 1.05 |
| H3K4me3 | 0.97 | 0.94 |
| H3K27ac | 1.01 | 0.99 |
