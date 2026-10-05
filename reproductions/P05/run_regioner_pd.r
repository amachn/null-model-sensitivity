
suppressMessages({
  library(regioneR)
  library(GenomicRanges)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript run_regioner_pd.R <ocr_peaks.bed> <variants.bed> <chrom.sizes> [ntimes]")
}

peaks_path       <- args[1]
variants_path    <- args[2]
chrom_sizes_path <- args[3]
ntimes <- ifelse(length(args) >= 4, as.integer(args[4]), 1000)

chrom_sizes <- read.table(chrom_sizes_path, col.names = c("chrom", "length"))
genome_gr <- GRanges(seqnames = chrom_sizes$chrom,
                      ranges = IRanges(start = 1, end = chrom_sizes$length))

read_bed_gr <- function(path) {
  df <- read.table(path, sep = "\t", stringsAsFactors = FALSE, quote = "")
  GRanges(seqnames = df$V1, ranges = IRanges(start = df$V2 + 1, end = df$V3))
}

peaks    <- read_bed_gr(peaks_path)
variants <- read_bed_gr(variants_path)

cat(sprintf("Loaded %d OCR peaks, %d variants\n", length(peaks), length(variants)))

set.seed(1)

pt <- overlapPermTest(A = peaks, B = variants, genome = genome_gr,
                       ntimes = ntimes, alternative = "greater",
                       non.overlapping = FALSE, verbose = FALSE)
print(pt)
cat(sprintf("\nZ-score = %.3f   p-value = %s   observed overlaps = %d   mean permuted = %.2f\n",
            pt$numOverlaps$zscore, format(pt$numOverlaps$pval, digits = 3),
            pt$numOverlaps$observed, mean(pt$numOverlaps$permuted)))

cat("\n=== COMPARE AGAINST PAPER'S REPORTED VALUES (Table 1 / Supp Table S2) ===\n")
cat("PD vs STC neurons:   GoShifter adj.p=0.028   GREGOR adj.p=6.94e-05  (SIGNIFICANT)\n")
cat("IBD vs STC neurons:  GoShifter adj.p=1        GREGOR adj.p=1         (not significant, negative control)\n")
cat("PEF vs STC neurons:  GoShifter adj.p=1        GREGOR adj.p=1         (not significant, negative control)\n")