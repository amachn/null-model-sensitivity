#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table); library(GenomicRanges); library(annotatr); library(regioneR)
})
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(here, "S3_lib.R"))

opt <- list(regions = "regions", marks = "data/marks", out = "results", n_random = 100, seed = 1,
            null = "bedtools", minoverlap = 200, blacklist = NULL, chrom_sizes = NULL,
            figvals = file.path(here, "figure_values.tsv"))
a <- commandArgs(trailingOnly = TRUE)
for (i in seq(1, length(a), by = 2)) opt[[sub("^--", "", a[i])]] <- a[i + 1]
opt$n_random <- as.integer(opt$n_random); opt$seed <- as.integer(opt$seed); opt$minoverlap <- as.integer(opt$minoverlap)
stopifnot(opt$null %in% c("bedtools", "annotatr"))
if (is.null(opt$chrom_sizes)) opt$chrom_sizes <- file.path(opt$regions, "chrom_sizes.autosomes.txt")
if (opt$null == "bedtools") {
  if (is.null(opt$blacklist)) stop("--null bedtools needs --blacklist <file> or --blacklist none")
  if (Sys.which("bedtools") == "") stop("bedtools not found on PATH")
  for (f in c(if (opt$blacklist != "none") opt$blacklist, opt$chrom_sizes)) if (!file.exists(f)) stop("missing file: ", f)
}
cat(sprintf("[info] null = %s | minoverlap = %d bp | %d random sets\n", opt$null, opt$minoverlap, opt$n_random))
dir.create(file.path(opt$out, "annotations", "categories"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(opt$out, "annotations", "marks"), recursive = TRUE, showWarnings = FALSE)
o <- function(...) file.path(opt$out, ...)
set.seed(opt$seed)

CATS <- c("hg38_genes_1to5kb", "hg38_genes_promoters", "hg38_genes_cds", "hg38_genes_5UTRs", "hg38_genes_exons",
          "hg38_genes_firstexons", "hg38_genes_introns", "hg38_genes_intronexonboundaries", "hg38_genes_3UTRs",
          "hg38_genes_intergenic", "hg38_enhancers_fantom")

bed_to_gr <- function(f) {
  d <- fread(f, header = FALSE, select = 1:3)
  gr <- GRanges(d$V1, IRanges(d$V2 + 1L, d$V3))
  gr <- keepSeqlevels(gr, unique(as.character(seqnames(gr))), pruning.mode = "coarse")
  genome(gr) <- "hg38"
  gr
}
q <- bed_to_gr(file.path(opt$regions, "caqtl_peaks.bed"))
cat(sprintf("[info] %d caQTL peaks (paper: 14,076)\n", length(q)))

cat("[info] building annotatr annotations (TxDb.Hsapiens.UCSC.hg38.knownGene",
    as.character(packageVersion("TxDb.Hsapiens.UCSC.hg38.knownGene")), ")\n")
annots <- build_annotations(genome = "hg38", annotations = CATS)
annots <- annots[seqnames(annots) %in% paste0("chr", 1:22)]

mark_files <- list.files(opt$marks, pattern = "\\.(bed|narrowPeak|broadPeak)(\\.gz)?$", full.names = TRUE)
marks <- lapply(mark_files, function(f) {
  d <- fread(f, header = FALSE, select = 1:3)
  reduce(GRanges(d$V1, IRanges(d$V2 + 1L, d$V3)))
})
names(marks) <- sub("\\.(bed|narrowPeak|broadPeak)(\\.gz)?$", "", basename(mark_files))
marks <- lapply(marks, function(m) m[seqnames(m) %in% paste0("chr", 1:22)])
cat(sprintf("[info] %d histone mark sets: %s\n", length(marks), paste(names(marks), collapse = ", ")))
if (!length(marks)) cat("[warn] no marks found in", opt$marks, "- S11 will be skipped (run S2 first)\n")

count_categories <- function(gr) {
  ar <- tryCatch(annotate_regions(gr, annots, minoverlap = opt$minoverlap, ignore.strand = TRUE, quiet = TRUE),
                 error = function(e) NULL)
  if (is.null(ar)) return(setNames(integer(length(CATS)), CATS))
  s <- summarize_annotations(ar, quiet = TRUE)
  v <- setNames(s$n, s$annot.type)
  out <- setNames(integer(length(CATS)), CATS); out[names(v)] <- v; out
}
count_marks <- function(gr) {
  gr <- unique(gr)
  vapply(marks, function(m) sum(overlapsAny(gr, m, minoverlap = opt$minoverlap, ignore.strand = TRUE)), 0L)
}
obs_cat <- count_categories(q)
obs_mark <- if (length(marks)) count_marks(q) else integer(0)

qbed <- file.path(opt$regions, "caqtl_peaks.bed")
shuffle_once <- function(i) {
  if (opt$null == "annotatr")
    return(suppressWarnings(randomize_regions(q, allow.overlaps = TRUE, per.chromosome = TRUE, quiet = TRUE)))
  tmp <- tempfile(fileext = ".bed")
  excl <- if (opt$blacklist == "none") character(0) else c("-excl", opt$blacklist)
  rc <- system2("bedtools", c("shuffle", "-i", qbed, "-g", opt$chrom_sizes, "-chrom", excl, "-seed", i), stdout = tmp)
  if (rc != 0) stop("bedtools shuffle failed")
  r <- bed_to_gr(tmp); unlink(tmp); r
}
null <- vector("list", opt$n_random)
t0 <- Sys.time()
for (i in seq_len(opt$n_random)) {
  r <- shuffle_once(opt$seed * 100000L + i)
  null[[i]] <- rbind(data.table(rep = i, annotation = CATS, n = count_categories(r)),
                     if (length(marks)) data.table(rep = i, annotation = names(marks), n = count_marks(r)))
  null[[i]][, n_regions := length(r)]
  if (i == 1) cat(sprintf("[info] ~%.1f min for all sets\n", as.numeric(difftime(Sys.time(), t0, units = "mins")) * opt$n_random))
}
null <- rbindlist(null)
fwrite(null, o("random_set_counts.tsv"), sep = "\t")

fig <- fread(opt$figvals)
s10 <- summarise_null(obs_cat, null[annotation %in% CATS])
c10 <- compare_to_published(s10, fig, "S10")
fwrite(c10, o("S10_reproduction.tsv"), sep = "\t")
plot_paper_style(s10, "FDR5 caQTL Peak Genome Annotations (reproduced)", o("S10_reproduced"))
plot_vs_published(c10, "Fig S10: reproduced vs published", o("S10_vs_published"))
if (length(marks)) {
  s11 <- summarise_null(obs_mark, null[annotation %in% names(marks)])
  c11 <- compare_to_published(s11, fig, "S11")
  fwrite(c11, o("S11_reproduction.tsv"), sep = "\t")
  plot_paper_style(s11, "Liver FDR5 caQTL Peaks ENCODE Roadmap Epigenomic Marks (reproduced)", o("S11_reproduced"))
  plot_vs_published(c11, "Fig S11: reproduced vs published", o("S11_vs_published"))
}
show <- function(x) print(x[, .(annotation, observed, obs_approx, null_median = round(null_median), null_median_approx,
                                fold = round(fold, 2), p_two_sided, direction_matches_figure, call_matches_caption)])
cat("\n== Fig S10 ==\n"); show(c10)
if (length(marks)) { cat("\n== Fig S11 ==\n"); show(c11) }

wr <- function(gr, f) fwrite(as.data.table(reduce(gr, ignore.strand = TRUE))[, .(seqnames, start - 1L, end)],
                             f, sep = "\t", col.names = FALSE)
for (ct in CATS) wr(annots[annots$type == ct], o("annotations", "categories", paste0(ct, ".bed")))
for (m in names(marks)) wr(marks[[m]], o("annotations", "marks", paste0(m, ".bed")))
man <- rbind(
  data.table(figure = "S10", annotation = CATS, annotation_bed = o("annotations", "categories", paste0(CATS, ".bed"))),
  if (length(marks)) data.table(figure = "S11", annotation = names(marks), annotation_bed = o("annotations", "marks", paste0(names(marks), ".bed"))))
man[, `:=`(pair_id = paste(figure, annotation, sep = ":"),
           query_bed = normalizePath(file.path(opt$regions, "caqtl_peaks.bed")),
           universe_bed = normalizePath(file.path(opt$regions, "all_peaks.bed")),
           mask_bed = file.path(normalizePath(opt$regions), "mask_gaps_blacklist.bed"),
           chrom_sizes = file.path(normalizePath(opt$regions), "chrom_sizes.autosomes.txt"),
           authors_null = switch(opt$null,
             bedtools = sprintf("bedtools shuffle -chrom -excl %s (gaps not excluded); %d iterations", basename(opt$blacklist), opt$n_random),
             annotatr = "annotatr::randomize_regions -> regioneR::randomizeRegions; per.chromosome=TRUE; allow.overlaps=TRUE; no mask"),
           overlap_definition = sprintf("distinct query regions with >= %d bp overlap", opt$minoverlap))]
res_all <- rbind(c10, if (length(marks)) c11, fill = TRUE)
man <- merge(man, res_all[, .(annotation, reproduced_observed = observed, reproduced_null_median = null_median,
                              reproduced_fold = fold, reproduced_p = p_two_sided, published_obs_approx = obs_approx,
                              published_null_approx = null_median_approx, published_direction = direction_in_figure,
                              published_caption_claim = caption_claim)], by = "annotation", all.x = TRUE)
setcolorder(man, c("pair_id", "figure", "annotation"))
fwrite(man, o("region_pairs_manifest.tsv"), sep = "\t")

writeLines(c(sprintf("null=%s minoverlap=%d n_random=%d seed=%d blacklist=%s", opt$null, opt$minoverlap,
                     opt$n_random, opt$seed, ifelse(is.null(opt$blacklist), "NA", opt$blacklist)),
             capture.output(sessionInfo()),
             sprintf("annotatr randomize_regions deprecated in this version: %s",
                     any(grepl("Deprecated", deparse(body(annotatr::randomize_regions)))))),
           o("sessionInfo.txt"))
cat("\n[done]", opt$out, ": S10/S11 tables and figures, region_pairs_manifest.tsv, annotations/\n")
