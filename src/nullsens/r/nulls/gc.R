# ~~~ src/nullsens/r/nulls/gc.R ~~~
# - prepares genome sequence indexes for GC matching
# - computes GC and N fractions for genomic intervals
# - validates candidate regions against the GC content of their originals

# resolve native GC backend relative to this file
nulls_dir <- dirname(
  normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = TRUE))

# compiled GC-index functions are cached for the lifetime of the R process
gc_lib <- new.env(parent = emptyenv())

gc_get_lib <- function() {
  if (is.null(gc_lib$fns)) {
    if (!requireNamespace("Rcpp", quietly = TRUE)) {
      stop("GC matching requires the Rcpp package.", call. = FALSE)
    }

    env <- new.env(parent = globalenv())

    tryCatch(
      Rcpp::sourceCpp(file = file.path(nulls_dir, "gc_index.cpp"), env = env),
      error = function(err) {
        stop(
          "GC matching requires a working C++ compiler for Rcpp.\n",
          conditionMessage(err), call. = FALSE
        )
      }
    )

    gc_lib$fns <- env
  }

  gc_lib$fns
}

# window indexes, shared by all randomizers in one R session
gc_index_store <- new.env(parent = emptyenv())

prepare_gc <- function(genome, tolerance = 0.02, flank = 0L, max_n = 0.1) {
  if (!is.numeric(tolerance) || length(tolerance) != 1L || tolerance <= 0) {
    stop(
      "GC tolerance must be a single number > 0, got: ",
      paste(tolerance, collapse = ", ")
    )
  }

  if (!is.numeric(flank) || length(flank) != 1L || flank < 0) {
    stop(
      "GC flank must be a single number >= 0, got: ",
      paste(flank, collapse = ", ")
    )
  }

  if (!is.numeric(max_n) || length(max_n) != 1L || max_n < 0 || max_n > 1) {
    stop(
      "GC max_n must be between 0 and 1, got: ",
      paste(max_n, collapse = ", ")
    )
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

  # compile the C++ backend during GC setup
  lib <- gc_get_lib()

  # chromosome indexes are built lazily and cached
  index_cache <- new.env(parent = emptyenv())
  source_name <- if (is.character(genome)) genome else NULL

  fetch_index <- function(chr) {
    key <- if (is.null(source_name)) NULL else paste(source_name, chr, sep = "|")
    hit <- if (is.null(key)) index_cache[[chr]] else gc_index_store[[key]]

    if (is.null(hit)) {
      seq <- sequences[[chr]]
      if (methods::is(seq, "MaskedDNAString")) seq <- Biostrings::unmasked(seq)
      hit <- lib$gc_build_index(charToRaw(as.character(seq)))
      if (is.null(key)) index_cache[[chr]] <- hit else gc_index_store[[key]] <- hit
    }

    hit
  }

  list(
    fetch_index  = fetch_index,
    window_stats = lib$gc_window_stats,
    seqnames     = names(genome_seqs),
    seq_lengths  = genome_seqs,
    source       = if (is.character(genome)) genome else class(genome)[1],
    tolerance    = tolerance,
    flank        = as.integer(flank),
    max_n        = max_n,
    n_sequences  = length(genome_seqs)
  )
}

# compute GC and N fractions for each interval (+/- flank), preserving order
gc_content <- function(context, regions) {
  n <- length(regions)
  gc <- rep(NA_real_, n)
  n_frac <- rep(1, n)
  if (n == 0L) return(list(gc = gc, n_frac = n_frac))

  chroms <- as.character(GenomicRanges::seqnames(regions))
  chrom_lengths <- unname(context$seq_lengths[chroms])

  left <- pmax(1L, GenomicRanges::start(regions) - context$flank)
  right <- pmin(chrom_lengths, GenomicRanges::end(regions) + context$flank)

  # skip windows that are off the chromosome
  usable <- !is.na(chrom_lengths) & left <= right
  if (!any(usable)) return(list(gc = gc, n_frac = n_frac))

  used_chroms <- unique(chroms[usable])
  indexes <- lapply(used_chroms, context$fetch_index)
  chr_ids <- match(chroms[usable], used_chroms)
  counts <- context$window_stats(
    indexes, chr_ids, as.integer(left[usable]), as.integer(right[usable])
  )

  widths <- right[usable] - left[usable] + 1L
  called <- widths - counts$n

  gc[usable] <- ifelse(called > 0, counts$gc / called, NA_real_)
  n_frac[usable] <- counts$n / widths

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
