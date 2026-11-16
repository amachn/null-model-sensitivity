suppressPackageStartupMessages({
  library(rtracklayer); library(GenomicRanges); library(GenomeInfoDb)
  library(IRanges); library(S4Vectors)
})

# shared helper: genome -> GRanges that always has seqlengths
as_genome_gr <- function(genome) {
  g <- if (methods::is(genome, "GRanges")) genome else regioneR::getGenome(genome)
  GenomicRanges::strand(g) <- "*"
  if (anyNA(GenomeInfoDb::seqlengths(g))) {
    # fall back to the range ends as chromosome lengths
    GenomeInfoDb::seqlengths(g) <- as.integer(
      tapply(end(g), as.character(GenomicRanges::seqnames(g)), max)[GenomeInfoDb::seqlevels(g)])
  }
  g
}

conform <- function(gr, genome_gr) {
  if (!any(seqlevels(gr) %in% seqlevels(genome_gr)))
    tryCatch(seqlevelsStyle(gr) <- "UCSC", error = function(e) NULL)
  keep <- intersect(seqlevels(gr), seqlevels(genome_gr))
  if (length(keep) == 0) stop("No chromosome names match genome_gr")
  n0 <- length(gr)
  gr <- keepSeqlevels(gr, keep, pruning.mode = "coarse")
  if (length(gr) < n0) message("Dropped ", n0 - length(gr), " ranges not in genome_gr")
  suppressWarnings({
    seqlevels(gr) <- seqlevels(genome_gr)
    seqinfo(gr)   <- seqinfo(genome_gr)
    gr <- trim(gr)
  })
  strand(gr) <- "*"
  gr
}
read_ccre <- function(path, class_col = 10) {
  df <- read.delim(path, header = FALSE, quote = "", comment.char = "",
                   stringsAsFactors = FALSE)
  gr <- GRanges(df[[1]], IRanges(df[[2]] + 1L, df[[3]]))   # BED is 0-based
  gr$class <- df[[class_col]]
  gr
}



read_any <- function(path, chr = "Chr", start = "Start", end = "End",
                     header = TRUE) {
  ext <- tolower(tools::file_ext(sub("\\.gz$", "", path)))
  if (ext %in% c("bed", "narrowpeak", "broadpeak", "bigbed", "bb", "gff", "gtf")) {
    return(rtracklayer::import(path))
  }
  df <- read.delim(path, header = header, quote = "", comment.char = "",
                   check.names = FALSE, stringsAsFactors = FALSE)
  if (!header) {                       # positional: chrom, start, end, name, class
    chr <- 1; start <- 2; end <- 3
  } else stopifnot(all(c(chr, start, end) %in% names(df)))
  chrs <- as.character(df[[chr]])
  if (!any(grepl("^chr", chrs))) chrs <- paste0("chr", chrs)
  gr <- GRanges(chrs, IRanges(as.integer(df[[start]]), as.integer(df[[end]])))
  if (!header && ncol(df) >= 5) {
    gr$name  <- df[[4]]
    gr$class <- df[[5]]
  }
  gr
}
