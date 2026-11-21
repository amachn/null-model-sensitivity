
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(IRanges)
  library(Biostrings)
})

prepare_gc <- function(genome, tolerance = 0.02, flank = 0L, max_n = 0.1) {
  if (!is.numeric(tolerance) || length(tolerance) != 1L || tolerance <= 0) {
    stop("GC tolerance must be a single number > 0, got: ", paste(tolerance, collapse = ", "))
  }
  if (!is.numeric(flank) || length(flank) != 1L || flank < 0) {
    stop("GC flank must be a single number >= 0, got: ", paste(flank, collapse = ", "))
  }
  if (!is.numeric(max_n) || length(max_n) != 1L || max_n < 0 || max_n > 1) {
    stop("GC max_n must be between 0 and 1, got: ", paste(max_n, collapse = ", "))
  }

  # GC is computed from the genome sequence, so a BSgenome is required
  sequences <- if (inherits(genome, "BSgenome") || methods::is(genome, "DNAStringSet")) {
    genome
  } else {
    if (!requireNamespace("BSgenome", quietly = TRUE)) {
      stop("GC matching requires the BSgenome package.")
    }
    tryCatch(
      suppressWarnings(
        suppressMessages(
          BSgenome::getBSgenome(genome, masked = FALSE)
        )
      ),
      error = function(err) {
        stop(
          "GC matching could not load a BSgenome for '", genome, "'.\n",
          conditionMessage(err), "\n",
          "Run: nullsens genomes install ", genome
        )
      }
    )
  }

  genome_seqs <- if (inherits(sequences, "BSgenome")) {
    GenomeInfoDb::seqlengths(sequences)
  } else {
    stats::setNames(as.numeric(Biostrings::width(sequences)), names(sequences))
  }

  # chromosome sequences are loaded on first use and cached for later permutations
  cache <- new.env(parent = emptyenv())
  fetch <- function(chr) {
    if (is.null(cache[[chr]])) {
      seq <- sequences[[chr]]
      if (methods::is(seq, "MaskedDNAString")) seq <- Biostrings::unmasked(seq)
      cache[[chr]] <- seq
    }
    cache[[chr]]
  }

  list(
    fetch       = fetch,
    seqnames    = names(genome_seqs),
    source      = if (is.character(genome)) genome else class(genome)[1],
    tolerance   = tolerance,
    flank       = as.integer(flank),
    max_n       = max_n,
    n_sequences = length(genome_seqs)
  )
}

# GC fraction and N fraction of each interval (+/- flank), order preserved
gc_content <- function(context, regions) {
  n <- length(regions)
  gc <- rep(NA_real_, n)
  n_frac <- rep(1, n)

  chroms <- as.character(GenomicRanges::seqnames(regions))
  starts <- GenomicRanges::start(regions) - context$flank
  ends <- GenomicRanges::end(regions) + context$flank

  for (chr in intersect(unique(chroms), context$seqnames)) {
    idx <- which(chroms == chr)
    seq <- context$fetch(chr)
    views <- IRanges::Views(
      seq,
      start = as.integer(pmax(1, starts[idx])),
      end = as.integer(pmin(length(seq), ends[idx]))
    )
    widths <- IRanges::width(views)
    n_count <- Biostrings::letterFrequency(views, "N")[, 1]
    gc_count <- Biostrings::letterFrequency(views, "GC")[, 1]
    called <- widths - n_count

    gc[idx] <- ifelse(called > 0, gc_count / called, NA_real_)
    n_frac[idx] <- n_count / widths
  }

  list(gc = gc, n_frac = n_frac)
}

validate_gc <- function(context, original, candidate, original_gc = NULL) {
  if (length(original) != length(candidate)) {
    stop("validate_gc: original and candidate must have the same length and order.")
  }

  # calculate original GC only when it wasn't precomputed
  if (is.null(original_gc)){
    original_gc <- gc_content(context, original)$gc
  }

  if (length(original_gc) != length(original)) {
    stop("validate_gc: original_gc must contain one value per original region.")
  }
  
  # candidate sequence changes every proposal, so this still needs calculation
  candidate_content <- gc_content(context, candidate)
  
  valid <- (
    abs(candidate_content$gc - original_gc) <= context$tolerance &
    candidate_content$n_frac <= context$max_n
  )

  # an original with no called bases has no GC to match; any low-N candidate is fine
  no_gc <- is.na(original_gc)
  valid[no_gc] <- candidate_content$n_frac[no_gc] <= context$max_n

  valid[is.na(valid)] <- FALSE
  valid
}
