suppressPackageStartupMessages({
  library(regioneR)
  library(rtracklayer)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(S4Vectors)
})

build_mappability_mask <- function(mappability_bed, genome, threshold = NULL, invert = FALSE) {
  track <- rtracklayer::import(mappability_bed)
  
  has_score <- "score" %in% names(S4Vectors::mcols(track))
  if (!is.null(threshold) && has_score) {
    track <- track[S4Vectors::mcols(track)$score < threshold]
  }
  
  track <- GenomicRanges::reduce(track)
  
  genome_gr   <- regioneR::getGenome(genome)
  genome_seqs <- GenomeInfoDb::seqlevels(genome_gr)
  mask_seqs   <- GenomeInfoDb::seqlevels(track)
  
  common <- intersect(mask_seqs, genome_seqs)
  if (length(common) == 0) {
    stop(
      "Mappability mask shares no chromosome names with genome '", genome, "'.\n",
      "  mask chromosomes:   ", paste(utils::head(mask_seqs, 5), collapse = ", "), " ...\n",
      "  genome chromosomes: ", paste(utils::head(genome_seqs, 5), collapse = ", "), " ...\n",
      "Check genome build and chr-prefix consistency between '", mappability_bed,
      "' and the requested genome."
    )
  }
  
  dropped <- setdiff(mask_seqs, genome_seqs)
  if (length(dropped) > 0) {
    warning(
      length(dropped), " chromosome(s) in '", mappability_bed,
      "' are not present in genome '", genome, "' and will be dropped: ",
      paste(utils::head(dropped, 5), collapse = ", "),
      if (length(dropped) > 5) " ..." else ""
    )
  }
  
  track <- track[GenomicRanges::seqnames(track) %in% common]
  GenomeInfoDb::seqlevels(track) <- common
  
  genome_gr <- genome_gr[GenomicRanges::seqnames(genome_gr) %in% common]
  GenomeInfoDb::seqlevels(genome_gr) <- common
  
  if (invert) {
    track <- GenomicRanges::setdiff(genome_gr, track)
  }
  
  track
}


get_mappability_randomizer <- function(
    masked, mappability_bed = NULL, genome = NULL, threshold = NULL, invert = FALSE
) {
  if (!masked || is.null(mappability_bed)) {
    return(function(A, genome, mask = NULL, ...) {
      regioneR::randomizeRegions(A = A, genome = genome, mask = mask, ...)
    })
  }
  
  mappability_mask <- build_mappability_mask(mappability_bed, genome, threshold, invert)
  
  function(A, genome, mask = NULL, ...) {
    combined_mask <- if (is.null(mask)) {
      mappability_mask
    } else {
      GenomicRanges::reduce(c(mask, mappability_mask))
    }
    regioneR::randomizeRegions(A = A, genome = genome, mask = combined_mask, ...)
  }
}