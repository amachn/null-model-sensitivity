#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
opt <- list(peaks = NULL, caqtl = NULL, chrom_sizes = NULL, gaps = NULL, blacklist = NULL, out = "regions",
            cis_window = 10000, id_col = NULL, chr_col = NULL, pos_col = NULL)
a <- commandArgs(trailingOnly = TRUE)
for (i in seq(1, length(a), by = 2)) opt[[sub("^--", "", a[i])]] <- a[i + 1]
opt$cis_window <- as.numeric(opt$cis_window)
dir.create(opt$out, showWarnings = FALSE, recursive = TRUE)
o <- function(f) file.path(opt$out, f)
AUTO <- paste0("chr", 1:22)
ucsc <- function(x) { x <- as.character(x); ifelse(grepl("^chr", x), x, paste0("chr", x)) }
wbed <- function(d, f) fwrite(d, o(f), sep = "\t", col.names = FALSE)

pk <- fread(opt$peaks, header = FALSE, sep = "\t", select = 1:4, colClasses = list(character = c(1, 4)))
setnames(pk, c("chr", "start", "end", "name"))
pk[, `:=`(chr = ucsc(chr), row1 = .I, row0 = .I - 1L)]
cat(sprintf("[info] %s peaks (paper: 2,518,633)\n", format(nrow(pk), big.mark = ",")))

RQ <- c("feature", "rsid", "chr", "pos", "ref", "alt", "af", "hwe_chisq", "ia", "log10_q", "chisq", "pi", "delta",
        "phi", "overdisp", "snp_idx", "n_fsnp", "n_tested", "iter_null", "iter_alt", "tie", "loglik_null",
        "converge", "r2_fsnp", "r2_rsnp")
cq <- fread(opt$caqtl)
if (!is.null(opt$id_col)) {
  setnames(cq, c(opt$id_col, opt$chr_col, opt$pos_col), c("feature", "chr", "pos"))
} else if (!all(c("feature", "chr", "pos") %in% names(cq))) {
  if (ncol(cq) < 25) stop("Unrecognised caQTL columns: ", paste(names(cq), collapse = " "),
                          "\nSee Files_column_description.rtf and pass --id_col/--chr_col/--pos_col.")
  if (!grepl("^V[0-9]+$", names(cq)[1])) cq <- fread(opt$caqtl, header = FALSE, skip = 1)
  setnames(cq, 1:25, RQ)
}
cq[, `:=`(feature = as.character(feature), chr = ucsc(chr), pos = as.numeric(pos))]
if ("log10_q" %in% names(cq)) setorder(cq, feature, log10_q, -chisq)
cq <- unique(cq, by = "feature")
key <- function(x) sub("^peak_?", "", as.character(x))
cq[, k := key(feature)]
cat(sprintf("[info] %s caQTL peaks (paper: 14,076)\n", format(nrow(cq), big.mark = ",")))

sc <- rbindlist(lapply(c("name", "row1", "row0"), function(s) {
  m <- merge(cq[, .(k, lchr = chr, pos)], pk[, .(k = key(get(s)), chr, start, end)], by = "k")
  data.table(scheme = s, matched = nrow(m),
             frac_in_window = m[, sum(lchr == chr & pos > start - opt$cis_window & pos <= end + opt$cis_window)] / nrow(cq))
}))
print(sc)
best <- sc[which.max(frac_in_window)]
if (best$frac_in_window < 0.95) stop("No ID scheme puts >=95% of lead SNPs within +/-", opt$cis_window, " bp of their peak.")
cat(sprintf("[info] ID scheme '%s' (%.2f%% of leads in window)\n", best$scheme, 100 * best$frac_in_window))
pk[, id := key(get(best$scheme))]
pk[, is_caqtl := id %in% cq$k]

pk <- pk[chr %in% AUTO]
setorder(pk, chr, start)
wbed(pk[, .(chr, start, end, id)], "all_peaks.bed")
wbed(pk[is_caqtl == TRUE, .(chr, start, end, id)], "caqtl_peaks.bed")
fwrite(merge(cq[, intersect(c("k", "chr", "pos", "rsid", "log10_q", "chisq", "pi"), names(cq)), with = FALSE],
             pk[, .(k = id, peak_chr = chr, start, end)], by = "k"), o("caqtl_leads.tsv"), sep = "\t")
cat(sprintf("[info] caQTL peaks written: %d (median width %d bp)\n", pk[is_caqtl == TRUE, .N],
            pk[is_caqtl == TRUE, as.integer(median(end - start))]))

if (!is.null(opt$chrom_sizes)) {
  cs <- fread(opt$chrom_sizes, header = FALSE)[V1 %in% AUTO]
  fwrite(cs, o("chrom_sizes.autosomes.txt"), sep = "\t", col.names = FALSE)
}
if (!is.null(opt$gaps) && !is.null(opt$blacklist)) {
  g <- fread(opt$gaps, header = FALSE)[, .(chr = ucsc(V2), start = V3, end = V4)]
  b <- fread(opt$blacklist, header = FALSE, select = 1:3)[, .(chr = ucsc(V1), start = V2, end = V3)]
  m <- rbind(g, b)[chr %in% AUTO][order(chr, start)]
  wbed(m, "mask.unmerged.bed")
  if (system(sprintf("bedtools merge -i %s > %s && rm %s", o("mask.unmerged.bed"), o("mask_gaps_blacklist.bed"),
                     o("mask.unmerged.bed"))) != 0) cat("[warn] bedtools not found; mask left unmerged\n")
}
cat("[done]", opt$out, "\n")
