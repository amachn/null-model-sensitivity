# ~~~ src/nullsens/r/nulls/build_randomizer.R ~~~
# - builds one regioneR-compatible randomization function from the
#   requested null constraints
#
# - at build time:
#   - prepare mappability mask, if enabled
#   - prepare GC context, if enabled
#   - prepare chromatin context, if enabled
#   - construct validator functions
# - at permutation time:
#   - generate exact-width randomized regions
#   - enforce static workspace masks during randomization
#   - validate GC/chromatin constraints
#   - resample only failed regions until every region passes

# runs every enabled per-region validator
# - each validator must return exactly one T/F value for each candidate region:
#   - TRUE = candidate satisfies the constraint
#   - FALSE = candidate must be resampled
# - a candidate is accepted only when it passes every enabled validator
validate_candidates <- function(original, candidate, validators) {
  valid <- rep(TRUE, length(candidate))
  
  for (name in names(validators)) {
    result <- validators[[name]](original, candidate)
    
    if (!is.logical(result) || length(result) != length(candidate)) {
      stop(
        "Null validator '", name, "' must return one TRUE/FALSE value ",
        "per candidate region."
      )
    }
    
    # an indeterminate result can't be accepted
    result[is.na(result)] <- FALSE
    valid <- valid & result
  }
  
  valid
}

# verify invariants required of every nullsens randomization
# - the randomized set must:
#   - contain exactly one candidate per original region
#   - preserve the exact width of each corresponding original region
# - maintaining index order is important because GC/chromatin matching
#   compare candidate[i] against original[i]
validate_randomization <- function(original, randomized) {
  if (length(original) != length(randomized)) {
    stop("Randomization changed the number of regions.")
  }
  
  if (
    !all(GenomicRanges::width(original) == GenomicRanges::width(randomized))
  ) {
    stop("Randomization failed to preserve exact region widths.")
  }
  
  invisible(TRUE)
}

# build one randomization function for one null-model config
# - `regions` is kept because it is the original region set
#   - future impls may use it to precompute original-region properties once
#     rather than repeatedly
build_randomizer <- function(
  config, refs, regions, genome, max_attempts = 100L
) {
  if (max_attempts < 1L) { stop("max_attempts must be an integer >= 1.") }
  
  # ----
  # 1 -> prepare context, static mask and candidate validators
  # ----
  # these context objects need to exist outside of their if statements
  mappability_context <- NULL
  chromatin_context <- NULL
  
  # NA means nullsens explicitly has no baseline forbidden mask
  proposal_mask <- NA
  
  # GC and chrom validators will be stored here
  validators <- list()
  
  # ----
  # 2 -> prepare mappability
  # ----
  # mappability is primarily a static constraint rather than a validator.
  # its forbidden regions can thus be handed directly to regioneR.
  if (isTRUE(config$mappability$enabled)) {
    if (is.null(refs$mappability)) {
      stop("Mappability masking is enabled but no reference was provided.")
    }
    
    mappability_context <- prepare_mappability(
      mappability_bed     = refs$mappability,
      genome              = genome,
      threshold           = config$mappability$threshold,
      inclusion_workspace = config$mappability$inclusion_workspace
    )
    
    # normally, mappability is enforced directly by randomizeRegions()
    proposal_mask <- mappability_context$mask
  }
  
  # ----
  # 3 -> prepare GC matching
  # ----
  # gc is region-specific: 
  # candidate[i] must have approx. the same GC fraction as original[i].
  # because of this, it becomes a validator used during rejection/resampling.
  if (isTRUE(config$gc$enabled)) {
    gc_context <- prepare_gc(
      genome,
      tolerance = config$gc$tolerance,
      flank     = config$gc$flank,
      max_n     = config$gc$max_n
    )
    
    # original regions never change between permutations, so calculate
    # their GC content only one time
    gc_original_content <- gc_content(gc_context, regions)
    
    # validators don't always receive the full original region set
    # - may receive only A[failed], and batched retries can contain multiple
    #   copies of the same original region
    # - this creates a stable key for each prepared region so those 
    #   subsets/copies can be mapped back
    gc_original_keys <- paste(
      as.character(GenomicRanges::seqnames(regions)),
      GenomicRanges::start(regions),
      GenomicRanges::end(regions),
      sep = ":"
    )
    
    validators$gc <- function(original, candidate) {
      # build the same coordinate-based keys for the original regions in this
      # particular validation call
      original_keys <- paste(
        as.character(GenomicRanges::seqnames(original)),
        GenomicRanges::start(original),
        GenomicRanges::end(original),
        sep = ":"
      )
      
      # map each region back to its position in the complete original set
      original_indices <- match(original_keys, gc_original_keys)
      
      validate_gc(
        context     = gc_context, 
        original    = original,
        candidate   = candidate,
        original_gc = gc_original_content$gc[original_indices]
      )
    }
  }
  
  # ----
  # 4 -> prepare chromatin-context matching
  # ----
  # like gc, chromatin depends on the corresponding original interval, so it
  # is also represented by a validator rather than a static mask.
  if (isTRUE(config$chromatin$enabled)) {
    if (is.null(refs$chromatin)) {
      stop("Chromatin matching is enabled but no reference was provided.")
    }
    
    chromatin_context <- prepare_chromatin(
      refs$chromatin, genome,
      padding   = config$chromatin$padding,
      label_col = config$chromatin$label_col,
      priority  = config$chromatin$priority
    )
    
    # keep validation even though the proposer is already chromatin-aware.
    validators$chromatin <- function(original, candidate) {
      validate_chromatin(
        context   = chromatin_context,
        original  = original,
        candidate = candidate,
        strict    = isTRUE(config$chromatin$strict)
      )
    }
    
    # propose_chromatin() doesn't understand regioneR masks - thus, when
    # mappability + chromatin are both enabled, mappability must move from
    # proposal-time enforcement to post-proposal validation
    if (!is.null(mappability_context)) {
      validators$mappability <- function(original, candidate) {
        validate_mappability(
          context   = mappability_context,
          original  = original,
          candidate = candidate
        )
      }
      
      # the chromatin proposer will be used instead of randomizeRegions(),
      # so this static mask is no longer used
      proposal_mask <- NA
    }
  }
  
  # ----
  # 5 -> define the randomization proposer
  # ----
  propose <- function(A, genome, ...) {
    # chromatin has its own state-aware proposal mechanism
    if (!is.null(chromatin_context)) {
      return(propose_chromatin(chromatin_context, A))
    }
    
    # otherwise, use regioneR's normal exact-width randomization
    args <- list(...)
    args$mask <- NULL # ignore any externally supplied mask
    
    do.call(
      regioneR::randomizeRegions,
      c(list(A = A, genome = genome, mask = proposal_mask), args)
    )
  }
  
  # everything above happens ONCE when the analysis is constructed
  # ----
  # 6 -> return the regioneR-compatible randomization closure
  # ----
  # everything below happens repeatedly, once per permutation
  function(A, genome, ...) {
    # -----
    # 6A -> generate one candidate for every original region
    # -----
    randomized <- propose(A = A, genome = genome, ...)
    validate_randomization(A, randomized)
    
    # -----
    # 6B -> if no region-specific constraints are enabled, finish here
    # -----
    if (length(validators) == 0L) {
      return(randomized)
    }
    
    # -----
    # 6C -> check the first set of proposals
    # -----
    valid <- validate_candidates(
      original = A, candidate = randomized, validators = validators
    )
    
    # initial proposal counts as attempt 1
    attempts <- rep(1L, length(A))
    
    # -----
    # 6D -> resample only failed candidates
    # -----
    while(any(!valid)) {
      failed <- which(!valid)
      
      # number of attempts each failed region still has available
      remaining <- max_attempts - attempts[failed]
      
      # stop before making another proposal for regions that have already 
      # reached the allowed attempt count
      exceeded <- failed[remaining <= 0L]
      if (length(exceeded) > 0L) {
        failed_regions <- A[exceeded]
        
        failed_labels <- paste0(
          as.character(GenomicRanges::seqnames(failed_regions)), ":",
          GenomicRanges::start(failed_regions), "-",
          GenomicRanges::end(failed_regions)
        )
        
        stop(
          "Unable to place ", length(exceeded), " region(s) after ",
          max_attempts, " attempts under the requested null-model constraints.",
          "\nFailed region(s): ",
          paste(utils::head(failed_labels, 10L), collapse = ", "),
          if (length(exceeded) > 10L) " ..." else ""
        )
      }
      
      # number of proposals generated for each currently failed region
      batch_sizes <- ifelse(
        attempts[failed] < 25L, # if <25 attempts, only 1 proposal
        1L,
        ifelse(
          attempts[failed] < 100L, # if <100 attempts, 8 proposals
          8L,
          ifelse(
            attempts[failed] < 500L, # if <500 attempts, 32 proposals
            32L,
            64L # otherwise, 64 proposals at a time
          )
        )
      )
      batch_sizes <- pmin(batch_sizes, remaining)
      
      # repeat each failed original region batch_sizes[i] times
      source_indices <- rep(failed, times = batch_sizes)
      batch_original <- A[source_indices]
      
      # generate all retry candidates in one proposal call
      replacements <- propose(A = batch_original, genome = genome, ...)
      validate_randomization(batch_original, replacements)
      
      # test the replacement candidates against every enabled constraint
      replacement_valid <- validate_candidates(
        original = batch_original,
        candidate = replacements,
        validators = validators
      )
      
      # collect successful mappings first instead of modifying the randomized
      # GRanges once per accepted region
      accepted_original <- integer(length(failed))
      accepted_replacement <- integer(length(failed))
      attempt_increments <- integer(length(failed))
      
      n_accepted <- 0L
      offset <- 0L
      
      for (group_idx in seq_along(failed)) {
        batch_size <- batch_sizes[group_idx]
        
        # candidates belonging to this one original region
        candidate_indices <- seq.int(offset + 1L, offset + batch_size)
        offset <- offset + batch_size
        
        valid_positions <- which(replacement_valid[candidate_indices])
        
        if (length(valid_positions) > 0L) {
          # rejection sampling stops at the first successful proposal
          first_valid <- valid_positions[1]
          n_accepted <- n_accepted + 1L
          
          accepted_original[n_accepted] <- failed[group_idx]
          accepted_replacement[n_accepted] <- candidate_indices[first_valid]
          attempt_increments[group_idx] <- first_valid
        } else {
          # every proposal in this batch failed
          attempt_increments[group_idx] <- batch_size
        }
      }
      
      # update attempt counts in one vectorized operation
      attempts[failed] <- attempts[failed] + attempt_increments
      
      # perform one GRanges replacement for all successes in this retry cycle
      if (n_accepted > 0L) {
        accepted_indices <- seq_len(n_accepted)
        original_indices <- accepted_original[accepted_indices]
        replacement_indices <- accepted_replacement[accepted_indices]
        
        randomized[original_indices] <- replacements[replacement_indices]
        valid[original_indices] <- TRUE
      }
    }
    
    # -----
    # 6E -> final invariant check and return
    # -----
    validate_randomization(A, randomized)
    
    randomized
  }
}
