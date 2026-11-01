# ~~~ src/nullsens/r/nulls/mappability.R ~~~
# - prepares forbidden genomic regions for mappability-constrained randomization
suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(IRanges)
  library(S4Vectors)
})

prepare_mappability <- function(mappability_bed, genome, threshold = NULL, invert = FALSE) {
  if (is.null(mappability_bed) || !file.exists(mappability_bed)) {
    stop("Mappability reference not found: ", mappability_bed)
  }

  track <- rtracklayer::import(mappability_bed)

  has_score <- "score" %in% names(S4Vectors::mcols(track))
  if (!is.null(threshold) && has_score) {
    track <- track[S4Vectors::mcols(track)$score < threshold]
  }

  GenomicRanges::strand(track) <- "*"
  track <- GenomicRanges::reduce(track)

  genome_gr <- if (methods::is(genome, "GRanges")) genome else regioneR::getGenome(genome)
  GenomicRanges::strand(genome_gr) <- "*"

  genome_seqs <- GenomeInfoDb::seqlevels(genome_gr)
  track_seqs  <- GenomeInfoDb::seqlevels(track)
  common      <- intersect(track_seqs, genome_seqs)

  if (length(common) == 0) {
    stop(
      "Mappability reference shares no chromosome names with the genome.\n",
      "  reference chromosomes: ", paste(utils::head(track_seqs, 5), collapse = ", "), " ...\n",
      "  genome chromosomes:    ", paste(utils::head(genome_seqs, 5), collapse = ", "), " ...\n",
      "Check genome build and chr-prefix consistency for '", mappability_bed, "'."
    )
  }

  dropped <- setdiff(track_seqs, genome_seqs)
  
  if (length(dropped) > 0) {
    warning(
      length(dropped), " chromosome(s) in '", mappability_bed,
      "' are not in the genome and will be dropped: ",
      paste(utils::head(dropped, 5), collapse = ", "),
      if (length(dropped) > 5) " ..." else ""
    )
  }

  track <- track[GenomicRanges::seqnames(track) %in% genome_seqs]
  GenomeInfoDb::seqlevels(track) <- genome_seqs

  mask <- if (invert) GenomicRanges::setdiff(genome_gr, track) else track

  list(
    mask               = mask,
    source             = mappability_bed,
    threshold          = threshold,
    invert             = invert,
    n_intervals        = length(mask),
    frac_genome_masked = sum(as.numeric(GenomicRanges::width(mask))) /
                           sum(as.numeric(GenomicRanges::width(genome_gr)))
  )
}

validate_mappability <- function(context, original, candidate) {
  !IRanges::overlapsAny(candidate, context$mask, ignore.strand = TRUE)
}
