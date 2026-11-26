# ~~~ src/nullsens/r/nulls/chromatin.R ~~~
# - prepares chromatin-state references
# - proposes candidates within the original region's chromatin state
# - validates chromatin-state preservation

# import a chromatin-state BED and normalize its selected state label
# - standard BED files (e.g. BED4 `name`) are imported through rtracklayer
# - `class` preserves compatibility with the ENCODE cCRE-style reference
read_chromatin_reference <- function(path, label_col = NULL) {
  if (identical(label_col, "class")) {
    class_col <- 10L
    data <- utils::read.delim(
      path, header = FALSE, quote = "",
      comment.char = "", stringsAsFactors = FALSE
    )
    
    if (ncol(data) < class_col) {
      stop(
        "Chromatin reference doesn't contain cCRE class column",
        class_col, ": ", path
      )
    }
    
    regions <- GenomicRanges::GRanges(
      seqnames = data[[1]],
      ranges = IRanges::IRanges(
        start = data[[2]] + 1L, # BED is 0-based
        end = data[[3]]
      )
    )
    
    labels <- as.character(data[[class_col]])
  } else {
    regions <- rtracklayer::import(path)
    
    if (is.null(label_col)) {
      S4Vectors::mcols(regions) <- NULL
      return(regions)
    }
    
    available <- names(S4Vectors::mcols(regions))
    
    if (!label_col %in% available) {
      stop(
        "Chromatin label column '", label_col, "' wasn't found in ", path,
        ". Available metadata columns: ",
        if (length(available) > 0L) {
          paste(available, collapse = ", ")
        } else { "none" }
      )
    }
    
    labels <- as.character(S4Vectors::mcols(regions)[[label_col]])
  }
  
  if (anyNA(labels) || any(!nzchar(labels))) {
    stop("Chromatin reference contains missing/empty state labels: ", path)
  }
  
  S4Vectors::mcols(regions) <- S4Vectors::DataFrame(chromatin_label = labels)
  
  regions
}

# prepare mutually exclusive chromatin-state regions
prepare_chromatin <- function(
  chromatin_bed, genome, padding = 0L, label_col = NULL, priority = NULL
) {
  if (is.null(chromatin_bed) || !file.exists(chromatin_bed)) {
    stop("Chromatin reference not found: ", chromatin_bed)
  }

  genome_gr <- as_genome_gr(genome)
  context <- conform_regions(
    read_chromatin_reference(chromatin_bed, label_col),
    genome_gr
  )
  
  # optionally enlarge chromatin-state intervals before constructing the map
  if (padding > 0L) {
    context <- suppressWarnings(
      GenomicRanges::trim(
        GenomicRanges::resize(
          context, GenomicRanges::width(context) + 2L * padding, fix = "center"
        )
      )
    )
  }
  
  # construct either one generic "in" state or one state per supplied class
  if (is.null(label_col)) {
    groups <- list(`in` = context)
  } else {
    groups <- as.list(split(context, context$chromatin_label))
  }
  
  # priority determines which state wins if state annotations overlap
  if (!is.null(priority)) groups <- groups[intersect(priority, names(groups))]
  
  # build mutually exclusive state intervals
  states <- list()
  covered <- GenomicRanges::GRanges(seqinfo = GenomeInfoDb::seqinfo(genome_gr))
  
  for (name in names(groups)) {
    state <- GenomicRanges::setdiff(
      GenomicRanges::intersect(
        GenomicRanges::reduce(groups[[name]]), genome_gr), 
      covered
    )
    
    states[[name]] <- state
    covered <- GenomicRanges::reduce(c(covered, state))
  }
  
  # everything not assigned to another state becomes "out"
  states[["out"]] <- GenomicRanges::setdiff(genome_gr, covered)
  
  list(states = states, genome_gr = genome_gr)
}

# determine the chromatin state containing each region's midpoint
state_of <- function(context, regions) {
  midpoint <- GenomicRanges::resize(regions, 1L, fix = "center")
  state <- rep("out", length(regions))
  
  for (name in setdiff(names(context$states), "out")) {
    state[IRanges::overlapsAny(
      midpoint, context$states[[name]], ignore.strand = TRUE
    )] <- name
  }
  
  state
}

# propose candidates directly in each original region's chromatin state
propose_chromatin <- function(context, regions) {
  states <- state_of(context, regions)
  widths <- GenomicRanges::width(regions)
  chromosome_lengths <- GenomeInfoDb::seqlengths(context$genome_gr)
  chromosomes <- character(length(regions))
  starts <- integer(length(regions))
  
  for (state in unique(states)) {
    indices <- which(states == state)
    pool <- context$states[[state]]
    
    if (length(pool) == 0L) {
      stop(
        "Chromatin state '", state, "' has original regions but no sequence."
      )
    }
    
    pool_chromosomes <- as.character(GenomicRanges::seqnames(pool))
    pool_chromosome_lengths <- as.numeric(chromosome_lengths[pool_chromosomes])
    
    # sample state intervals proportional to their widths
    selected <- sample(
      length(pool), length(indices), 
      replace = TRUE, prob = as.numeric(GenomicRanges::width(pool))
    )
    
    # a region can't be placed on a chromosome shorter than itself
    # - resample those rare proposals from only eligible chromosomes
    too_short <- (pool_chromosome_lengths[selected] < widths[indices])
    
    if (any(too_short)) {
      for (position in which(too_short)) {
        region_width <- widths[indices[position]]
        
        eligible <- which(
          !is.na(pool_chromosome_lengths) &
          pool_chromosome_lengths >= region_width
        )
        
        if (length(eligible) == 0L) {
          stop(
            "Chromatin state '",
            state,
            "' has no chromosome long enough for a ",
            region_width,
            " bp region."
          )
        }
        
        selected[position] <- sample(
          eligible, 1L,
          prob = as.numeric(GenomicRanges::width(pool)[eligible])
        )
      }
    }
    
    midpoint <- (
      GenomicRanges::start(pool)[selected] + 
      floor(
        stats::runif(length(indices)) * 
        GenomicRanges::width(pool)[selected]
      )
    )
    
    chromosomes[indices] <- as.character(
      GenomicRanges::seqnames(pool))[selected]
    
    starts[indices] <- as.integer(midpoint - floor(widths[indices] / 2))
  }
  
  # keep each candidate within chromosome bounds while preserving width
  starts <- pmax(
    1L, pmin(
      starts, as.integer(chromosome_lengths[chromosomes] - widths + 1L)
    )
  )
  
  GenomicRanges::GRanges(
    chromosomes,
    IRanges::IRanges(starts, width = widths),
    seqinfo = GenomeInfoDb::seqinfo(context$genome_gr)
  )
}

# validate chromatin-state preservation
# - strict = FALSE -> candidate midpoint must be in original state
# - strict = TRUE  -> entire candidate interval must be in original state
validate_chromatin <- function(context, original, candidate, strict = FALSE) {
  chromosome_lengths <- GenomeInfoDb::seqlengths(context$genome_gr)
  chromosomes <- as.character(GenomicRanges::seqnames(candidate))
  in_bounds <- chromosomes %in% names(chromosome_lengths) &
    GenomicRanges::start(candidate) >= 1L &
    GenomicRanges::end(candidate) <= chromosome_lengths[chromosomes]
  
  original_states <- state_of(context, original)
  valid <- logical(length(candidate))
  
  for (state in unique(original_states)) {
    indices <- which(original_states == state & in_bounds)
    if (length(indices) == 0L) { next }
    
    if (strict) {
      valid[indices] <- IRanges::overlapsAny(
        candidate[indices], context$states[[state]],
        type = "within", ignore.strand = TRUE
      )
    } else {
      valid[indices] <- IRanges::overlapsAny(
        GenomicRanges::resize(candidate[indices], 1L, fix = "center"),
        context$states[[state]],
        ignore.strand = TRUE
      )
    }
  }
  
  valid
}
