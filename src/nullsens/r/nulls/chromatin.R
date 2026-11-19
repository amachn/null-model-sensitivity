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
  context   <- conform(read_ccre(chromatin_bed), genome_gr) #will match chrom names and do additional processing
  if (pad > 0)
    context <- suppressWarnings(trim(resize(context, width(context) + 2 * pad, fix = "center")))
  groups <- if (is.null(label_col)) list(`in` = context)
  else as.list(split(context, mcols(context)[[label_col]]))
  if (!is.null(priority)) groups <- groups[intersect(priority, names(groups))]
  
  states <- list(); covered <- GRanges(seqinfo = seqinfo(genome_gr))
  for (nm in names(groups)) {
    g <- GenomicRanges::setdiff(
      GenomicRanges::intersect(reduce(groups[[nm]]), genome_gr), covered)
    states[[nm]] <- g
    covered <- reduce(c(covered, g))
  }
  states[["out"]] <- GenomicRanges::setdiff(genome_gr, covered)
  list(states = states, genome_gr = genome_gr)
}

# which state each region's midpoint falls in
state_of <- function(ctx, gr) {
  mid <- resize(gr, 1, fix = "center")
  st <- rep("out", length(gr))
  for (nm in setdiff(names(ctx$states), "out"))
    st[overlapsAny(mid, ctx$states[[nm]], ignore.strand = TRUE)] <- nm
  st
}

# PROPOSE: generate candidates already inside each original's own state
propose_chromatin <- function(ctx, A) {
  st <- state_of(ctx, A); w <- width(A); lens <- seqlengths(ctx$genome_gr)
  chr <- character(length(A)); s0 <- integer(length(A))
  for (s in unique(st)) {
    i <- which(st == s); pool <- ctx$states[[s]]
    if (!length(pool)) stop("State '", s, "' has regions in A but no sequence")
    k   <- sample(length(pool), length(i), replace = TRUE, prob = as.numeric(width(pool)))
    mid <- start(pool)[k] + floor(runif(length(i)) * width(pool)[k])
    chr[i] <- as.character(seqnames(pool))[k]
    s0[i]  <- as.integer(mid - floor(w[i] / 2))
  }
  s0 <- pmax(1L, pmin(s0, as.integer(lens[chr] - w + 1)))   # clamp so widths are kept
  GRanges(chr, IRanges(s0, width = w), seqinfo = seqinfo(ctx$genome_gr))
}


# ---- validate: per-candidate accept/reject. Candidate i replaces original i.
validate_chromatin <- function(prep, original, candidate, strict = FALSE) {
  lens   <- seqlengths(prep$genome_gr)
  chr <- as.character(seqnames(candidate))
  in_bnd <- chr %in% names(lens) & start(candidate) >= 1 & end(candidate) <= lens[chr]
  st <- state_of(prep, original)
  ok <- logical(length(candidate))
  for (s in unique(st)) {
    i <- which(st == s & in_bnd)
    if (!length(i)) next
    ok[i] <- if (strict)
      overlapsAny(candidate[i], prep$states[[s]], type = "within", ignore.strand = TRUE)
    else
      overlapsAny(resize(candidate[i], 1, fix = "center"), prep$states[[s]], ignore.strand = TRUE)
  }
  ok
}


