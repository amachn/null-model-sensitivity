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

prepare_chromatin <- function(chromatin_bed, genome, pad = 2000,
                              label_col = NULL, priority = NULL) {
  if (is.null(chromatin_bed) || !file.exists(chromatin_bed))
    stop("Chromatin reference not found: ", chromatin_bed)
  genome_gr <- as_genome_gr(genome)
  context   <- conform(rtracklayer::import(chromatin_bed), genome_gr)
  prepare_chromatin_null(context, genome_gr, pad, label_col, priority)
  strand(context) <- "*"
  if (pad > 0)
    context <- suppressWarnings(trim(resize(context, width(context) + 2 * pad, fix = "center")))
  groups <- if (is.null(label_col)) list(`in` = context)
  else as.list(split(context, mcols(context)[[label_col]]))
  if (!is.null(priority)) groups <- groups[intersect(priority, names(groups))]
  
  states <- list(); covered <- GRanges()
  for (nm in names(groups)) {
    g <- GenomicRanges::setdiff(
      GenomicRanges::intersect(reduce(groups[[nm]]), genome_gr), covered)
    states[[nm]] <- g
    covered <- reduce(c(covered, g))
  }
  states[["out"]] <- GenomicRanges::setdiff(genome_gr, covered)
  list(states = states, genome_gr = genome_gr)
}

# ---- validate: per-candidate accept/reject. Candidate i replaces original i.
validate_chromatin_null <- function(prep, original, candidate, strict = FALSE) {
  lens   <- seqlengths(prep$genome_gr)
  in_bnd <- as.character(seqnames(candidate)) %in% names(lens) &
    start(candidate) >= 1 &
    end(candidate) <= lens[as.character(seqnames(candidate))]
  
  orig_state <- rep("out", length(original))
  omid <- resize(original, 1, fix = "center")
  for (nm in setdiff(names(prep$states), "out"))
    orig_state[overlapsAny(omid, prep$states[[nm]])] <- nm
  
  ok <- logical(length(candidate))
  for (s in unique(orig_state)) {
    i <- which(orig_state == s & in_bnd)
    if (!length(i)) next
    ok[i] <- if (strict) overlapsAny(candidate[i], prep$states[[s]], type = "within")
    else        overlapsAny(resize(candidate[i], 1, fix = "center"), prep$states[[s]])
  }
  ok
}

